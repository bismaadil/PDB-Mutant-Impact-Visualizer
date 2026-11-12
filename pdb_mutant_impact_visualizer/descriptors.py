"""Deterministic, offline descriptors for a substitution and its structural site.

Every value comes from a published table or algorithm and stays a separate
measurement; nothing here is combined into an impact score.
"""

from __future__ import annotations

import copy
from typing import Any

from Bio.Align import substitution_matrices
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Structure import Structure

from .local_structure import THREE_TO_ONE, get_residue
from .models import Mutation

from .tables import CHARGE_PH7, HYDROPATHY, MAX_ASA, VOLUME  # noqa: F401  (re-exported)

REFERENCES = [
    "https://doi.org/10.1016/0022-2836(82)90515-0",   # Kyte & Doolittle
    "https://doi.org/10.1073/pnas.89.22.10915",       # BLOSUM62
    "https://doi.org/10.1371/journal.pone.0080635",   # Tien max ASA
    "https://doi.org/10.1016/0022-2836(73)90011-9",   # Shrake-Rupley
]


def substitution_descriptors(mutation: Mutation) -> dict[str, Any]:
    wt, mt = mutation.wild_type, mutation.mutant
    blosum = substitution_matrices.load("BLOSUM62")
    return {
        "blosum62_score": int(blosum[wt][mt]),
        "hydropathy_change_kyte_doolittle": round(HYDROPATHY[mt] - HYDROPATHY[wt], 3),
        "volume_change_angstrom3": round(VOLUME[mt] - VOLUME[wt], 1),
        "charge_change_ph7_approx": CHARGE_PH7.get(mt, 0) - CHARGE_PH7.get(wt, 0),
        "proline_introduced": mt == "P",
        "glycine_removed": wt == "G",
        "note": "Table-based properties of the two amino acids; independent of structural context.",
    }


def solvent_accessibility(
    structure: Structure, chain_id: str, position: int, insertion_code: str = " "
) -> dict[str, Any]:
    """Shrake-Rupley SASA of the target residue in the first model (protein atoms only)."""
    target = get_residue(structure, chain_id, position, insertion_code)
    aa = THREE_TO_ONE[target.resname.upper()]
    work = copy.deepcopy(structure[0])
    for chain in work:
        for res in [r for r in chain if r.id[0] != " "]:
            chain.detach_child(res.id)
    ShrakeRupley(probe_radius=1.40, n_points=960).compute(work, level="R")
    residue = work[chain_id][target.id]
    sasa = float(residue.sasa)
    return {
        "sasa_angstrom2": round(sasa, 2),
        "relative_sasa": round(min(sasa / MAX_ASA[aa], 1.0), 3),
        "probe_radius_angstrom": 1.40,
        "context": "protein chains in the first model only; ligands, waters and ions removed",
        "note": "Relative SASA = SASA / theoretical maximum for the residue type (capped at 1).",
    }


def structure_warnings(structure: Structure, chain_id: str, position: int, insertion_code: str = " ",
                       cutoff: float = 10.0) -> list[str]:
    warnings: list[str] = []
    if len(structure) > 1:
        warnings.append(f"Structure has {len(structure)} models (e.g. NMR); only the first model was analysed.")
    target = get_residue(structure, chain_id, position, insertion_code)
    if target.is_disordered() or any(a.is_disordered() for a in target):
        warnings.append("Target residue has alternate conformations; the highest-occupancy conformer was used.")
    ca = target["CA"] if "CA" in target else None
    if ca is not None:
        ligands = {
            res.resname for res in structure[0].get_residues()
            if res.id[0].startswith("H_") and any((a - ca) <= cutoff for a in res)
        }
        if ligands:
            warnings.append("Hetero groups within the cutoff are not part of neighbor counts: " + ", ".join(sorted(ligands)))
    if all(float(a.bfactor) == 0.0 for a in target):
        warnings.append("All B-factors at the site are 0; the file may be a model without real B-factors.")
    return warnings
