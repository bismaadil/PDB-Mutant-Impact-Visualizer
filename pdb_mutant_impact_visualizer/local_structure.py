"""Coordinate-derived measurements; no biological interpretation is inferred."""

from __future__ import annotations

from typing import Any

from Bio.PDB import NeighborSearch
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from .models import AA, Mutation

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def get_residue(structure: Structure, chain_id: str, position: int, insertion_code: str = " ") -> Residue:
    if chain_id not in structure[0]:
        raise ValueError(f"Chain {chain_id!r} was not found in the structure.")
    # Only standard polymer residues (hetero flag " "); waters/ligands may share a number.
    polymer = [r for r in structure[0][chain_id] if r.id[0] == " "]
    code = (insertion_code or " ").strip() or " "
    matches = [r for r in polymer if r.id[1] == position and r.id[2] == code]
    if not matches and code == " ":
        matches = [r for r in polymer if r.id[1] == position]
    if not matches:
        raise ValueError(f"Residue {chain_id}:{position}{insertion_code.strip()} was not found.")
    if len(matches) > 1:
        raise ValueError("Residue location is ambiguous; supply its insertion code.")
    return matches[0]


def validate_mutation(structure: Structure, chain_id: str, mutation: Mutation, insertion_code: str = " ") -> Residue:
    residue = get_residue(structure, chain_id, mutation.position, insertion_code)
    actual = THREE_TO_ONE.get(residue.resname.upper())
    if actual is None:
        raise ValueError(f"Residue is non-standard ({residue.resname}); mutation modelling is unsupported.")
    if actual != mutation.wild_type:
        raise ValueError(
            f"{mutation} does not match the structure: {chain_id}:{mutation.position} is {actual}, "
            f"not {mutation.wild_type}. Check chain, numbering, and the selected structure."
        )
    if mutation.mutant not in AA:
        raise ValueError("Mutant must be a standard amino acid.")
    return residue


def analyze_local_structure(
    structure: Structure, chain_id: str, position: int, cutoff: float = 10.0,
    insertion_code: str = " ",
) -> dict[str, Any]:
    if cutoff <= 0:
        raise ValueError("Neighbor cutoff must be greater than zero.")
    residue = get_residue(structure, chain_id, position, insertion_code)
    if "CA" not in residue:
        raise ValueError(f"Residue {chain_id}:{position} has no alpha-carbon atom.")
    target = residue["CA"]
    model = structure[0]  # NMR/multi-model files: analyse the first model only
    is_protein = lambda res: res.id[0] == " " and THREE_TO_ONE.get(res.resname.upper()) is not None  # noqa: E731
    atoms = [a for a in model.get_atoms() if a.element not in {"H", "D"} and is_protein(a.get_parent())]
    neighbors = set()
    for neighbor in NeighborSearch(atoms).search(target.coord, cutoff, level="R"):
        if neighbor is residue:
            continue
        parent_chain = neighbor.get_parent().id
        neighbors.add((str(parent_chain), int(neighbor.id[1]), str(neighbor.id[2]).strip()))
    site_b = [float(atom.bfactor) for atom in residue if atom.element not in {"H", "D"}]
    protein_b = [float(a.bfactor) for a in atoms]
    return {
        "neighbor_cutoff_angstrom": float(cutoff),
        "neighbor_residue_count": len(neighbors),
        "neighbors": [
            {"chain": c, "author_residue_number": n, "insertion_code": ins}
            for c, n, ins in sorted(neighbors)
        ],
        "site_mean_b_factor_angstrom2": round(sum(site_b) / len(site_b), 3) if site_b else None,
        "protein_mean_b_factor_angstrom2": round(sum(protein_b) / len(protein_b), 3) if protein_b else None,
        "scope": "first model, standard amino-acid residues only (waters, ligands and ions excluded)",
        "b_factor_note": "Experimental/model B-factors are reported as supplied; they are not a direct measure of mutation-induced dynamics.",
    }
