"""The pipeline's stage registry: the single source of truth for the DAG.

Each stage reads upstream outputs from ``ctx`` (keyed by stage name) and the
:class:`RunConfig` from ``ctx["config"]``. Heavy / Biopython-dependent imports
happen inside the stage functions so that the graph can be built, drawn and
tested without a structure parser installed.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Mapping

from .config import RunConfig
from .dag import Graph, Stage
from .models import Mutation

Ctx = Mapping[str, Any]

METHOD_URLS = {
    "rcsb": "https://www.rcsb.org/docs/programmatic-access/file-download-services",
    "dynamut2": "https://doi.org/10.1002/pro.3942",
    "foldx": "https://doi.org/10.1093/bioinformatics/btz184",
}

EBI_REFERENCES: tuple[str, str] = (
    "https://doi.org/10.1093/nar/gkae241",
    "https://www.ebi.ac.uk/pdbe/coordinates/",
)


def _cfg(ctx: Ctx) -> RunConfig:
    config: RunConfig = ctx["config"]
    return config


def _mut(ctx: Ctx) -> Mutation:
    return Mutation.parse(_cfg(ctx).mutation)


# ------------------------------------------------------------------ stages
def fetch_structure(ctx: Ctx) -> Path:
    cfg = _cfg(ctx)
    if cfg.pdb_id:
        from .fetch_structure import download_pdb

        return download_pdb(cfg.pdb_id, Path(cfg.output_root) / "structures", cfg.refresh_structure)
    assert cfg.structure_path is not None
    return Path(cfg.structure_path)


def parse_structure(ctx: Ctx) -> Any:
    from .fetch_structure import load_structure

    cfg = _cfg(ctx)
    path = ctx["fetch_structure"]
    return load_structure(path, cfg.pdb_id or Path(path).stem)


def validate_mutation(ctx: Ctx) -> dict[str, Any]:
    from .local_structure import validate_mutation as _validate

    cfg = _cfg(ctx)
    residue = _validate(ctx["parse_structure"], cfg.chain, _mut(ctx), cfg.insertion_code)
    return {"residue": f"{cfg.chain}:{residue.id[1]}", "resname": residue.resname}


def local_geometry(ctx: Ctx) -> dict[str, Any]:
    from .local_structure import analyze_local_structure

    cfg = _cfg(ctx)
    return analyze_local_structure(
        ctx["parse_structure"], cfg.chain, _mut(ctx).position, cfg.neighbor_cutoff, cfg.insertion_code
    )


def solvent_accessibility(ctx: Ctx) -> dict[str, Any]:
    from .descriptors import solvent_accessibility as _sasa

    cfg = _cfg(ctx)
    return _sasa(ctx["parse_structure"], cfg.chain, _mut(ctx).position, cfg.insertion_code)


def substitution(ctx: Ctx) -> dict[str, Any]:
    from .descriptors import substitution_descriptors

    return substitution_descriptors(_mut(ctx))


def structure_warnings(ctx: Ctx) -> list[str]:
    from .descriptors import structure_warnings as _warn

    cfg = _cfg(ctx)
    return _warn(ctx["parse_structure"], cfg.chain, _mut(ctx).position, cfg.insertion_code, cfg.neighbor_cutoff)


def contacts(ctx: Ctx) -> dict[str, Any]:
    from .contacts import atoms_from_structure, contact_profile

    cfg = _cfg(ctx)
    table = atoms_from_structure(ctx["parse_structure"])
    return contact_profile(
        table, cfg.chain, _mut(ctx).position, cfg.insertion_code.strip(), side_chain_cutoff=cfg.contact_cutoff
    )


def language_model(ctx: Ctx) -> dict[str, Any]:
    """ESM-2 masked-marginal score using deposited polymer sequence when available."""
    from .esm_backend import (
        ESMError,
        HuggingFaceESM,
        masked_marginal_llr,
        score_chain_mutation,
    )
    from .local_structure import THREE_TO_ONE

    cfg = _cfg(ctx)
    assert cfg.esm_model
    mutation = _mut(ctx)
    structure = ctx["parse_structure"]

    residues = [
        (int(r.id[1]), str(r.id[2]).strip(), THREE_TO_ONE.get(r.resname.upper(), "X"))
        for r in structure[0][cfg.chain] if r.id[0] == " "
    ]

    factory = ctx.get("esm_factory") or HuggingFaceESM
    scorer = factory(cfg.esm_model)

    # A PDB ID gives us access to RCSB polymer metadata, including
    # unresolved residues. A user-supplied local PDB does not.
    if cfg.pdb_id:
        from .fetch_structure import download_mmcif
        from .mmcif_sequence import reconstruct_chain_sequence

        structure_path = Path(ctx["fetch_structure"])
        cif_path = download_mmcif(
            cfg.pdb_id,
            out_dir=structure_path.parent,
            refresh=cfg.refresh_structure,
        )
        sequence, residue_to_index, metadata = reconstruct_chain_sequence(
            cif_path, cfg.chain
        )

        key = (mutation.position, cfg.insertion_code.strip())
        if key not in residue_to_index:
            raise ESMError(
                f"Residue {mutation.position}{cfg.insertion_code.strip()} "
                f"has no author-number mapping in the deposited polymer sequence."
            )

        index = residue_to_index[key]
        result = masked_marginal_llr(
            scorer, sequence, index, mutation.wild_type, mutation.mutant
        )
        result.update({
            "backend": "ESM-2 masked marginal",
            "model": cfg.esm_model,
            "chain_sequence_length": len(sequence),
            "sequence_metadata": metadata,
            "coordinate_residue_count": len(residues),
            "references": [
                "https://doi.org/10.1101/2021.07.09.450648",
                "https://doi.org/10.1126/science.ade2574",
            ],
        })
        return result

    # Local PDB fallback: only coordinate-observed residues are available.
    result = score_chain_mutation(
        scorer, residues, mutation.position, cfg.insertion_code.strip(),
        mutation.wild_type, mutation.mutant, cfg.esm_model,
    )
    result["sequence_source"] = "coordinate-observed residues from local PDB"
    result["warning"] = (
        (result.get("warning", "") + " ")
        + "Local PDB input does not provide polymer metadata; unresolved residues "
        "cannot be restored and are absent from the ESM-2 sequence."
    ).strip()
    return result


def stability(ctx: Ctx) -> dict[str, Any]:
    from .stability import predict_dynamut2, predict_foldx

    cfg = _cfg(ctx)
    path, mutation = ctx["fetch_structure"], _mut(ctx)
    if cfg.stability_backend == "dynamut2":
        return predict_dynamut2(path, cfg.chain, mutation, cfg.poll_timeout_s, cfg.poll_interval_s)
    assert cfg.foldx_executable is not None
    return predict_foldx(
        cfg.foldx_executable, path, cfg.chain, mutation,
        work_root=Path(cfg.output_root) / "foldx-runs", timeout_s=cfg.poll_timeout_s,
    )


def conservation(ctx: Ctx) -> dict[str, Any]:
    from .conservation import analyze_conservation

    cfg = _cfg(ctx)
    assert cfg.email is not None
    return analyze_conservation(
        cfg.pdb_id or "", cfg.chain, _mut(ctx), cfg.email, insertion_code=cfg.insertion_code.strip(),
        uniprot_accession=cfg.uniprot_accession, uniprot_position=cfg.uniprot_position,
        timeout_s=cfg.poll_timeout_s, poll_interval_s=cfg.poll_interval_s,
    )


def local_block(ctx: Ctx) -> dict[str, Any]:
    """Merge geometry, accessibility and contacts into the report's local_structure block."""
    block = dict(ctx["local_geometry"])
    block["solvent_accessibility"] = ctx["solvent_accessibility"]
    if ctx.get("contacts") is not None:
        block["contacts"] = ctx["contacts"]
    return block


def ml_predict(ctx: Ctx) -> dict[str, Any]:
    from .features import extract_features
    from .ml import load_model, predict_from_features

    cfg = _cfg(ctx)
    assert cfg.ml_model_path
    model, card = load_model(cfg.ml_model_path)
    report_like = {
        "mutation": cfg.mutation, "substitution": ctx["substitution"], "local_structure": local_block(ctx),
        "language_model": ctx.get("language_model"),
    }
    return predict_from_features(model, card, extract_features(report_like))


def assemble(ctx: Ctx) -> Any:
    from .models import AnalysisResult

    cfg = _cfg(ctx)
    path = Path(ctx["fetch_structure"])
    from datetime import datetime, timezone

    from .descriptors import REFERENCES as DESCRIPTOR_REFERENCES
    from .fetch_structure import file_sha256
    from .provenance import dependency_versions

    refs: list[str] = [METHOD_URLS["rcsb"]] if cfg.pdb_id else []
    if cfg.stability_backend == "dynamut2":
        refs.append(METHOD_URLS["dynamut2"])
    elif cfg.stability_backend == "foldx":
        refs.append(METHOD_URLS["foldx"])
    if cfg.conservation_backend == "ebi":
        refs.extend(EBI_REFERENCES)
    if cfg.esm_model:
        from .esm_backend import ESM2_REFERENCE, REFERENCE

        refs.extend([REFERENCE, ESM2_REFERENCE])
    refs.extend(DESCRIPTOR_REFERENCES)
    policy = (
        "refreshed from RCSB for this run" if cfg.pdb_id and cfg.refresh_structure
        else "existing local cache reused when available" if cfg.pdb_id else "user-supplied file"
    )
    return AnalysisResult.create(
        structure={
            "pdb_id": cfg.pdb_id,
            "file_name": cfg.source_name or path.name,
            "sha256": file_sha256(path),
            "source": f"https://files.rcsb.org/download/{cfg.pdb_id}.pdb" if cfg.pdb_id else "user-supplied local PDB file",
            "coordinate_file_mtime_utc": datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat(),
            "cache_policy": policy,
        },
        chain=cfg.chain,
        mutation=str(_mut(ctx)),
        stability=ctx.get("stability"),
        conservation=ctx.get("conservation"),
        local_structure=local_block(ctx),
        substitution=ctx["substitution"],
        method_references=refs,
        dependency_versions=dependency_versions(),
        warnings=list(ctx["structure_warnings"]),
        language_model=ctx.get("language_model"),
        ml_prediction=ctx.get("ml_predict"),
    )


# ------------------------------------------------------------------- graph
def build_graph() -> Graph:
    def on(field: str) -> Callable[[Ctx], bool]:
        return lambda ctx: bool(getattr(ctx["config"], field))

    def backend(field: str) -> Callable[[Ctx], bool]:
        return lambda ctx: bool(getattr(ctx["config"], field) != "none")

    return Graph([
        Stage("fetch_structure", fetch_structure, kind="input", remote=True,
              label="Fetch coordinates<br/>(RCSB cache or local file)"),
        Stage("parse_structure", parse_structure, deps=("fetch_structure",), label="Parse PDB"),
        Stage("validate_mutation", validate_mutation, deps=("parse_structure",),
              label="Validate chain,<br/>residue, wild type"),
        Stage("local_geometry", local_geometry, deps=("validate_mutation",), label="Neighbors and B-factors"),
        Stage("solvent_accessibility", solvent_accessibility, deps=("validate_mutation",), label="SASA (Shrake-Rupley)"),
        Stage("substitution", substitution, deps=("validate_mutation",), label="BLOSUM62 / hydropathy /<br/>volume / charge"),
        Stage("structure_warnings", structure_warnings, deps=("validate_mutation",), label="Structure warnings"),
        Stage("contacts", contacts, deps=("validate_mutation",), optional=True, label="Contacts: H-bond, salt bridge,<br/>disulfide, interface"),
        Stage("language_model", language_model, deps=("validate_mutation",), optional=True, enabled=on("esm_model"),
              label="ESM-2 masked-marginal LLR"),
        Stage("stability", stability, deps=("validate_mutation",), optional=True, remote=True,
              enabled=backend("stability_backend"), label="Stability: DynaMut2 / FoldX"),
        Stage("conservation", conservation, deps=("validate_mutation",), optional=True, remote=True,
              enabled=backend("conservation_backend"), label="Conservation: SIFTS, BLAST,<br/>Clustal Omega, entropy"),
        Stage("ml_predict", ml_predict, deps=("local_geometry", "solvent_accessibility", "substitution", "contacts"),
              after=("language_model",), optional=True, enabled=on("ml_model_path"),
              label="Research ML model<br/>(user-trained ridge)"),
        Stage("assemble", assemble,
              deps=("fetch_structure", "local_geometry", "solvent_accessibility", "substitution",
                    "structure_warnings"),
              after=("contacts", "language_model", "stability", "conservation", "ml_predict"),
              kind="output", label="AnalysisResult<br/>schema 2.0 + manifest"),
    ])
