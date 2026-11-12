
"""Run one mutation analysis by executing the stage DAG (see ``stages.py``)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .config import RunConfig
from .dag import RunOutcome, StageFailure
from .models import AnalysisResult, Mutation
from .provenance import environment, reproduce_command
from .stages import build_graph

OPTIONAL_STAGE_LABELS = {
    "contacts": "Contact analysis",
    "stability": "Stability prediction",
    "conservation": "Conservation analysis",
    "language_model": "ESM-2 language-model score",
    "ml_predict": "Research ML model",
}


def run_pipeline(
    config: RunConfig,
    *,
    strict: bool = False,
    extra_context: dict[str, Any] | None = None,
) -> AnalysisResult:
    """Validate ``config``, execute the DAG, and attach the run manifest.

    ``strict=True`` turns any failed *requested* optional stage into an error.
    Otherwise the failure is recorded in ``warnings`` and in ``run.stages`` and
    the corresponding report field stays ``null`` (never a made-up value).

    If a required stage fails, the original exception type is preserved.
    The exception also receives the completed DAG records and failed-stage
    name so callers can inspect what happened before the failure.
    """
    config.validate()
    Mutation.parse(config.mutation)  # fail fast, before any network or disk work
    graph = build_graph()

    try:
        outcome: RunOutcome = graph.run(
            {"config": config, **(extra_context or {})}
        )
    except StageFailure as exc:
        cause = exc.cause

        # Preserve diagnostic information without wrapping the original error.
        setattr(
            cause,
            "stage_records",
            [record.to_dict() for record in exc.records],
        )
        setattr(cause, "failed_stage", exc.stage)

        raise cause from exc

    result: AnalysisResult = outcome.outputs["assemble"]
    failures = [record for record in outcome.records if record.status == "failed"]

    for record in failures:
        label = OPTIONAL_STAGE_LABELS.get(record.name, record.name)
        result.warnings.append(
            f"{label} FAILED and is reported as missing: {record.error}"
        )

    if strict and failures:
        raise RuntimeError(
            "; ".join(f"{record.name}: {record.error}" for record in failures)
        )

    public = config.public_dict()
    result.run = {
        "config_sha256": config.sha256(),
        "config": public,
        "stages": [record.to_dict() for record in outcome.records],
        "execution_order": graph.order,
        "environment": environment(),
        "reproduce_command": reproduce_command(public),
        "seed": config.seed,
    }
    return result


def analyze_mutation(
    *,
    mutation: Mutation,
    chain: str,
    pdb_id: str | None = None,
    structure_path: str | Path | None = None,
    insertion_code: str = " ",
    neighbor_cutoff: float = 10.0,
    stability_backend: str = "none",
    foldx_executable: str | Path | None = None,
    conservation_backend: str = "none",
    email: str | None = None,
    uniprot_accession: str | None = None,
    uniprot_position: int | None = None,
    output_root: str | Path = "data",
    poll_timeout_s: int = 900,
    poll_interval_s: int = 10,
    refresh_structure: bool = False,
    source_name: str | None = None,
    contact_cutoff: float = 4.5,
    esm_model: str | None = None,
    ml_model_path: str | Path | None = None,
    strict: bool = False,
    seed: int = 0,
) -> AnalysisResult:
    """Backward-compatible function API; builds a :class:`RunConfig` and runs the DAG."""
    config = RunConfig(
        chain=chain,
        mutation=str(mutation),
        pdb_id=pdb_id,
        structure_path=str(structure_path) if structure_path is not None else None,
        insertion_code=insertion_code,
        neighbor_cutoff=neighbor_cutoff,
        contact_cutoff=contact_cutoff,
        stability_backend=stability_backend,
        foldx_executable=str(foldx_executable) if foldx_executable else None,
        conservation_backend=conservation_backend,
        email=email,
        uniprot_accession=uniprot_accession,
        uniprot_position=uniprot_position,
        esm_model=esm_model,
        ml_model_path=str(ml_model_path) if ml_model_path else None,
        output_root=str(output_root),
        poll_timeout_s=poll_timeout_s,
        poll_interval_s=poll_interval_s,
        refresh_structure=refresh_structure,
        source_name=source_name,
        seed=seed,
    )
    return run_pipeline(config, strict=strict)
