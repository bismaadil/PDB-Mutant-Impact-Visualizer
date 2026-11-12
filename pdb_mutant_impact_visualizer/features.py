"""Numeric feature vectors from an analysis report (for the research/ML layer).

Only measurements already present in the report are used. ``None`` is kept as
``None`` here; callers decide whether to drop or reject rows. Feature groups
let experiments run clean ablations (e.g. substitution-only vs +structure
vs +language model).
"""

from __future__ import annotations

from typing import Any, Mapping

from .tables import HYDROPATHY, VOLUME
from .models import Mutation

FEATURE_GROUPS: dict[str, list[str]] = {
    "substitution": [
        "blosum62", "d_hydropathy", "d_volume", "d_charge", "proline_introduced", "glycine_removed",
        "wt_hydropathy", "wt_volume",
    ],
    "structure": [
        "relative_sasa", "neighbor_count_10A", "cb_neighbors_10A", "side_chain_contacts",
        "hbond_candidates", "salt_bridge_candidates", "is_interface", "has_disulfide", "site_b_ratio",
    ],
    "language_model": ["esm_llr"],
}
ALL_FEATURES: list[str] = [f for names in FEATURE_GROUPS.values() for f in names]


def _flag(value: Any) -> float | None:
    return None if value is None else float(bool(value))


def extract_features(report: Mapping[str, Any]) -> dict[str, float | None]:
    """Map a report dict (``AnalysisResult.to_dict()``) to the feature dictionary."""
    mutation = Mutation.parse(report["mutation"])
    sub = report.get("substitution") or {}
    local = report.get("local_structure") or {}
    sasa = local.get("solvent_accessibility") or {}
    contacts = local.get("contacts") or {}
    esm = report.get("language_model") or {}
    site_b, protein_b = local.get("site_mean_b_factor_angstrom2"), local.get("protein_mean_b_factor_angstrom2")
    b_ratio = round(site_b / protein_b, 4) if site_b is not None and protein_b else None
    wt = mutation.wild_type
    return {
        "blosum62": sub.get("blosum62_score"),
        "d_hydropathy": sub.get("hydropathy_change_kyte_doolittle"),
        "d_volume": sub.get("volume_change_angstrom3"),
        "d_charge": sub.get("charge_change_ph7_approx"),
        "proline_introduced": _flag(sub.get("proline_introduced")),
        "glycine_removed": _flag(sub.get("glycine_removed")),
        "wt_hydropathy": HYDROPATHY[wt],
        "wt_volume": VOLUME[wt],
        "relative_sasa": sasa.get("relative_sasa"),
        "neighbor_count_10A": local.get("neighbor_residue_count") if local.get("neighbor_cutoff_angstrom") == 10.0 else None,
        "cb_neighbors_10A": contacts.get("cb_neighbors_within_10A"),
        "side_chain_contacts": contacts.get("side_chain_contact_residue_count"),
        "hbond_candidates": contacts.get("hydrogen_bond_candidate_count"),
        "salt_bridge_candidates": contacts.get("salt_bridge_candidate_count"),
        "is_interface": _flag(contacts.get("is_interface_residue")),
        "has_disulfide": None if not contacts else float(bool(contacts.get("disulfide_bonds"))),
        "site_b_ratio": b_ratio,
        "esm_llr": esm.get("esm_llr"),
    }


def select(features: Mapping[str, float | None], groups: list[str]) -> list[str]:
    unknown = [g for g in groups if g not in FEATURE_GROUPS]
    if unknown:
        raise ValueError(f"Unknown feature group(s): {unknown}. Choose from {list(FEATURE_GROUPS)}.")
    return [name for g in groups for name in FEATURE_GROUPS[g]]
