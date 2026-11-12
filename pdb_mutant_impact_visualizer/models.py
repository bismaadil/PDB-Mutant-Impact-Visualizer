"""Validated, provenance-bearing data models for mutation analyses."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from typing import Any

AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
SCHEMA_VERSION = "2.0"


@dataclass(frozen=True)
class Mutation:
    wild_type: str
    position: int
    mutant: str

    @classmethod
    def parse(cls, value: str) -> "Mutation":
        import re

        match = re.fullmatch(r"([A-Z])(\d+)([A-Z])", value.strip().upper())
        if not match:
            raise ValueError("Mutation must look like A123V (one-letter amino-acid codes).")
        wild_type, position, mutant = match.groups()
        if wild_type not in AA or mutant not in AA:
            raise ValueError("Use standard amino-acid codes (ACDEFGHIKLMNPQRSTVWY).")
        if int(position) < 1:
            raise ValueError("Residue positions start at 1.")
        if wild_type == mutant:
            raise ValueError("The mutant residue must differ from the wild-type residue.")
        return cls(wild_type, int(position), mutant)

    def __str__(self) -> str:
        return f"{self.wild_type}{self.position}{self.mutant}"


@dataclass
class AnalysisResult:
    schema_version: str
    software: dict[str, str]
    dependency_versions: dict[str, str]
    analyzed_at_utc: str
    structure: dict[str, Any]
    chain: str
    mutation: str
    stability: dict[str, Any] | None
    conservation: dict[str, Any] | None
    local_structure: dict[str, Any]
    substitution: dict[str, Any] | None = None
    method_references: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    language_model: dict[str, Any] | None = None
    ml_prediction: dict[str, Any] | None = None
    run: dict[str, Any] | None = None  # DAG manifest: config hash, per-stage status/timing, environment

    @classmethod
    def create(
        cls, *, structure: dict[str, Any], chain: str, mutation: str,
        stability: dict[str, Any] | None, conservation: dict[str, Any] | None,
        local_structure: dict[str, Any], dependency_versions: dict[str, str],
        substitution: dict[str, Any] | None = None,
        method_references: list[str] | None = None, warnings: list[str] | None = None,
        language_model: dict[str, Any] | None = None, ml_prediction: dict[str, Any] | None = None,
        run: dict[str, Any] | None = None,
    ) -> "AnalysisResult":
        from . import __version__

        return cls(
            schema_version=SCHEMA_VERSION,
            software={"name": "PDB Mutant Impact Visualizer", "version": __version__},
            analyzed_at_utc=datetime.now(timezone.utc).isoformat(),
            structure=structure,
            chain=chain,
            mutation=mutation,
            stability=stability,
            conservation=conservation,
            local_structure=local_structure,
            substitution=substitution,
            dependency_versions=dependency_versions,
            method_references=method_references or [],
            warnings=warnings or [],
            language_model=language_model,
            ml_prediction=ml_prediction,
            run=run,
        )

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)
