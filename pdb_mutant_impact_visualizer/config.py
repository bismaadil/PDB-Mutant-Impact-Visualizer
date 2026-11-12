"""Run configuration: one validated, serialisable object describes a whole run.

The configuration hash (``config_sha256``) is recorded in every report. Two
runs with the same coordinate SHA-256, configuration hash and software/image
version are intended to produce the same *offline* measurements. The contact
email is deliberately excluded from the hash and from the stored config
(personal data); only ``email_provided`` is kept.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Any

STABILITY_BACKENDS = ("none", "dynamut2", "foldx")
CONSERVATION_BACKENDS = ("none", "ebi")
SECRET_FIELDS = ("email",)


@dataclass
class RunConfig:
    chain: str
    mutation: str
    pdb_id: str | None = None
    structure_path: str | None = None
    insertion_code: str = " "
    neighbor_cutoff: float = 10.0
    contact_cutoff: float = 4.5
    stability_backend: str = "none"
    foldx_executable: str | None = None
    conservation_backend: str = "none"
    email: str | None = None
    uniprot_accession: str | None = None
    uniprot_position: int | None = None
    esm_model: str | None = None          # e.g. "facebook/esm2_t12_35M_UR50D"; None = disabled
    ml_model_path: str | None = None      # JSON model produced by `pdbmi train`
    output_root: str = "data"
    poll_timeout_s: int = 900
    poll_interval_s: int = 10
    refresh_structure: bool = False
    source_name: str | None = None
    seed: int = 0
    extras: dict[str, Any] = field(default_factory=dict)

    def validate(self) -> "RunConfig":
        self.chain = self.chain.strip()
        if self.pdb_id is not None:
            self.pdb_id = self.pdb_id.strip().upper() or None
        self.insertion_code = self.insertion_code or " "
        if (self.pdb_id is None) == (self.structure_path is None):
            raise ValueError("Provide exactly one of a PDB ID or a local PDB structure file.")
        if not self.chain:
            raise ValueError("Chain identifier cannot be empty.")
        if self.stability_backend not in STABILITY_BACKENDS:
            raise ValueError("Stability backend must be none, dynamut2, or foldx.")
        if self.conservation_backend not in CONSERVATION_BACKENDS:
            raise ValueError("Conservation backend must be none or ebi.")
        if self.stability_backend == "foldx" and not self.foldx_executable:
            raise ValueError("FoldX was selected; provide its executable with --foldx-exe.")
        if self.conservation_backend == "ebi" and (not self.email or "@" not in self.email):
            raise ValueError("EBI conservation requires a valid contact email address (--email).")
        if self.conservation_backend == "ebi" and not self.pdb_id and not (
            self.uniprot_accession and self.uniprot_position
        ):
            raise ValueError("EBI conservation for an uploaded file requires both a UniProt accession and position.")
        if self.neighbor_cutoff <= 0 or self.contact_cutoff <= 0:
            raise ValueError("Distance cutoffs must be greater than zero.")
        if self.poll_timeout_s <= 0 or self.poll_interval_s <= 0:
            raise ValueError("Polling timeout and interval must be positive.")
        return self

    # -- serialisation -------------------------------------------------
    def public_dict(self) -> dict[str, Any]:
        """Config as stored in reports: no secrets, no machine-specific paths."""
        data = asdict(self)
        for name in SECRET_FIELDS:
            data["email_provided"] = bool(data.pop(name))
        data["structure_path"] = Path(data["structure_path"]).name if data["structure_path"] else None
        data["foldx_executable"] = Path(data["foldx_executable"]).name if data["foldx_executable"] else None
        data["output_root"] = None  # where files land does not change the science
        data["ml_model_path"] = Path(data["ml_model_path"]).name if data["ml_model_path"] else None
        for volatile in ("source_name", "refresh_structure", "poll_timeout_s", "poll_interval_s"):
            data.pop(volatile, None)
        return data

    def sha256(self) -> str:
        payload = json.dumps(self.public_dict(), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()

    @classmethod
    def from_mapping(cls, data: dict[str, Any]) -> "RunConfig":
        known = {f.name for f in fields(cls)}
        unknown = set(data) - known
        if unknown:
            raise ValueError(f"Unknown configuration field(s): {', '.join(sorted(unknown))}")
        return cls(**data)

    @classmethod
    def from_file(cls, path: str | Path, **overrides: Any) -> "RunConfig":
        """Load YAML or JSON. Command-line overrides win over the file."""
        text = Path(path).read_text(encoding="utf-8")
        if str(path).lower().endswith((".yml", ".yaml")):
            import yaml  # PyYAML is a declared dependency

            data = yaml.safe_load(text) or {}
        else:
            data = json.loads(text)
        data.update({k: v for k, v in overrides.items() if v is not None})
        return cls.from_mapping(data)
