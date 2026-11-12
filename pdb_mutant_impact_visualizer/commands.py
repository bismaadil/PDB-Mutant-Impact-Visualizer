"""The ``run`` command (single mutation, single structure)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import click

from .config import RunConfig
from .models import Mutation
from .pipeline import run_pipeline
from .report import render_text, write_json




@click.command(name="run", context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--config", "config_file", type=click.Path(exists=True, dir_okay=False, path_type=Path), help="YAML/JSON run configuration; command-line options override it.")
@click.option("--pdb", "pdb_id", help="Four-character RCSB PDB identifier; mutually exclusive with --structure.")
@click.option("--structure", "structure_path", type=click.Path(exists=True, dir_okay=False, path_type=Path), help="Local PDB coordinate file.")
@click.option("--chain", help="PDB author chain identifier.")
@click.option("--mutation", help="Single substitution, e.g. S50V.")
@click.option("--insertion-code", default=None, help="Residue insertion code, if present.")
@click.option("--neighbor-cutoff", default=None, type=click.FloatRange(min=0.01), help="C-alpha to heavy-atom cutoff in Angstroms [default 10].")
@click.option("--contact-cutoff", default=None, type=click.FloatRange(min=0.01), help="Side-chain contact cutoff in Angstroms [default 4.5].")
@click.option("--stability", "stability_backend", type=click.Choice(["none", "dynamut2", "foldx"], case_sensitive=False), default=None, help="Remote DynaMut2 uploads your coordinates; opt-in. [default none]")
@click.option("--foldx-exe", type=click.Path(exists=True, dir_okay=False, path_type=Path), help="Full path to your licensed FoldX executable.")
@click.option("--conservation", "conservation_backend", type=click.Choice(["none", "ebi"], case_sensitive=False), default=None, help="[default none]")
@click.option("--email", help="Contact email required for EMBL-EBI programmatic jobs (never written to reports).")
@click.option("--uniprot", "uniprot_accession", help="Curated UniProt accession (for uploaded structures).")
@click.option("--uniprot-position", type=click.IntRange(min=1), help="Explicit 1-based UniProt position.")
@click.option("--esm-model", help="Enable the ESM-2 language-model score, e.g. facebook/esm2_t12_35M_UR50D (needs the 'ml' extra).")
@click.option("--ml-model", "ml_model_path", type=click.Path(exists=True, dir_okay=False, path_type=Path), help="Model JSON from `pdbmi train` (research use).")
@click.option("--timeout", default=None, type=click.IntRange(min=1), help="Maximum wait per remote job or FoldX command, seconds [default 900].")
@click.option("--poll-interval", default=None, type=click.IntRange(min=1), help="Seconds between remote status checks [default 10].")
@click.option("--output", "output_path", type=click.Path(dir_okay=False, path_type=Path), help="Write the complete JSON report.")
@click.option("--refresh-structure", is_flag=True, help="Re-download the RCSB structure instead of using the cache.")
@click.option("--strict", is_flag=True, help="Exit with an error if any requested optional stage fails (default: record and continue).")
@click.option("--data-dir", default=None, help="Directory for cached structures and FoldX runs [default ./data].")
def run_command(
    config_file: Path | None, pdb_id: str | None, structure_path: Path | None, chain: str | None,
    mutation: str | None, insertion_code: str | None, neighbor_cutoff: float | None,
    contact_cutoff: float | None, stability_backend: str | None, foldx_exe: Path | None,
    conservation_backend: str | None, email: str | None, uniprot_accession: str | None,
    uniprot_position: int | None, esm_model: str | None, ml_model_path: Path | None,
    timeout: int | None, poll_interval: int | None, output_path: Path | None,
    refresh_structure: bool, strict: bool, data_dir: str | None,
) -> None:
    """Analyze one real PDB structure and one single amino-acid substitution (executes the pipeline DAG)."""
    overrides: dict[str, Any] = {
        "pdb_id": pdb_id, "structure_path": str(structure_path) if structure_path else None, "chain": chain,
        "mutation": mutation, "insertion_code": insertion_code, "neighbor_cutoff": neighbor_cutoff,
        "contact_cutoff": contact_cutoff, "stability_backend": stability_backend,
        "foldx_executable": str(foldx_exe) if foldx_exe else None, "conservation_backend": conservation_backend,
        "email": email, "uniprot_accession": uniprot_accession, "uniprot_position": uniprot_position,
        "esm_model": esm_model, "ml_model_path": str(ml_model_path) if ml_model_path else None,
        "poll_timeout_s": timeout, "poll_interval_s": poll_interval, "output_root": data_dir,
        "refresh_structure": True if refresh_structure else None,
    }
    try:
        if config_file:
            config = RunConfig.from_file(config_file, **overrides)
        else:
            if not chain or not mutation:
                raise click.UsageError("--chain and --mutation are required (or supply --config).")
            config = RunConfig(**{k: v for k, v in overrides.items() if v is not None})
        Mutation.parse(config.mutation)
        result = run_pipeline(config, strict=strict)
    except ImportError as exc:
        raise click.ClickException(
            f"A required package is missing ({exc.name}). Install the project first: python -m pip install -e '.[app]'"
        ) from exc
    except (ValueError, RuntimeError) as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(render_text(result))
    if output_path:
        click.echo(f"\nFull JSON result: {write_json(result, output_path)}")
