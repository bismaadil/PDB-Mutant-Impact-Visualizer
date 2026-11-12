"""Run manifest: everything needed to re-run or audit an analysis."""

from __future__ import annotations

import os
import platform
import subprocess
from importlib.metadata import PackageNotFoundError, version
from typing import Any, Mapping

TRACKED_DISTRIBUTIONS = (
    "biopython", "numpy", "scipy", "requests", "click", "PyYAML", "streamlit", "plotly",
    "torch", "transformers",
)

# Report fields that must be identical when an *offline* analysis is repeated.
DETERMINISTIC_PATHS = (
    ("structure", "sha256"), ("chain",), ("mutation",),
    ("local_structure", "neighbor_residue_count"),
    ("local_structure", "solvent_accessibility", "sasa_angstrom2"),
    ("local_structure", "contacts", "side_chain_contact_residue_count"),
    ("local_structure", "contacts", "hydrogen_bond_candidate_count"),
    ("substitution", "blosum62_score"),
    ("run", "config_sha256"),
)


def _git_commit() -> str | None:
    env = os.environ.get("GIT_COMMIT")
    if env:
        return env
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, timeout=3, check=False)
        return out.stdout.strip() or None
    except (OSError, subprocess.SubprocessError):
        return None


def dependency_versions() -> dict[str, str]:
    versions = {"python": platform.python_version()}
    for name in TRACKED_DISTRIBUTIONS:
        try:
            versions[name] = version(name)
        except PackageNotFoundError:
            continue
    return versions


def environment() -> dict[str, Any]:
    return {
        "platform": platform.platform(),
        "git_commit": _git_commit(),
        "container_image_digest": os.environ.get("IMAGE_DIGEST"),
        "in_container": os.path.exists("/.dockerenv") or bool(os.environ.get("PDBMI_IN_CONTAINER")),
    }


def reproduce_command(config_public: Mapping[str, Any]) -> str:
    parts = ["pdbmi", "run"]
    if config_public.get("pdb_id"):
        parts += ["--pdb", str(config_public["pdb_id"])]
    elif config_public.get("structure_path"):
        parts += ["--structure", f"<{config_public['structure_path']}>"]
    parts += ["--chain", str(config_public["chain"]), "--mutation", str(config_public["mutation"])]
    parts += ["--stability", config_public["stability_backend"], "--conservation", config_public["conservation_backend"]]
    if config_public.get("esm_model"):
        parts += ["--esm-model", config_public["esm_model"]]
    if config_public.get("ml_model_path"):
        parts += ["--ml-model", f"<{config_public['ml_model_path']}>"]
    parts += ["--neighbor-cutoff", str(config_public["neighbor_cutoff"])]
    return " ".join(parts)


def _dig(data: Mapping[str, Any], path: tuple[str, ...]) -> Any:
    cur: Any = data
    for key in path:
        if not isinstance(cur, Mapping) or key not in cur:
            return None
        cur = cur[key]
    return cur


def compare_reports(a: Mapping[str, Any], b: Mapping[str, Any]) -> list[str]:
    """Return human-readable differences in deterministic fields (empty = reproduced)."""
    problems = []
    for path in DETERMINISTIC_PATHS:
        va, vb = _dig(a, path), _dig(b, path)
        if va != vb:
            problems.append(f"{'.'.join(path)}: {va!r} != {vb!r}")
    return problems
