"""Non-covalent contact environment of one residue, from raw atom coordinates.

All functions work on plain arrays (``AtomTable``) so they are fast, unit
testable without a structure parser, and independent of Biopython internals.
Distances are geometric *candidates*, not energetic claims: a donor/acceptor
pair within 3.5 A is reported as a hydrogen-bond *candidate* because hydrogen
positions and angles are not evaluated.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

import numpy as np
from scipy.spatial import cKDTree

BACKBONE = frozenset({"N", "CA", "C", "O", "OXT"})
POLAR_ELEMENTS = frozenset({"N", "O"})
ACID_ATOMS = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2")}
BASE_ATOMS = {("LYS", "NZ"), ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2")}
HBOND_MAX = 3.5
SALT_BRIDGE_MAX = 4.0
DISULFIDE_MAX = 2.5
INTERFACE_CUTOFF = 5.0
CB_SPHERE = 10.0


@dataclass(frozen=True)
class AtomTable:
    """Columnar atom data. Index ``i`` describes one heavy atom."""

    chain: np.ndarray      # str
    resnum: np.ndarray     # int
    icode: np.ndarray      # str ("" when absent)
    resname: np.ndarray    # str, upper case 3-letter
    atom: np.ndarray       # str, e.g. "CA"
    element: np.ndarray    # str, upper case
    coord: np.ndarray      # float (n, 3)

    def __len__(self) -> int:
        return int(self.coord.shape[0])

    @classmethod
    def from_rows(cls, rows: Iterable[tuple[str, int, str, str, str, str, tuple[float, float, float]]]) -> "AtomTable":
        rows = list(rows)
        if not rows:
            raise ValueError("No atoms supplied.")
        c, n, i, r, a, e, xyz = zip(*rows, strict=True)
        return cls(
            np.array(c, dtype=str), np.array(n, dtype=int), np.array(i, dtype=str),
            np.array([x.upper() for x in r], dtype=str), np.array(a, dtype=str),
            np.array([x.upper() for x in e], dtype=str), np.array(xyz, dtype=float),
        )

    def residue_mask(self, chain: str, resnum: int, icode: str) -> np.ndarray:
        return (self.chain == chain) & (self.resnum == resnum) & (self.icode == icode)


def atoms_from_structure(structure: Any) -> AtomTable:
    """Adapter: first model, standard amino-acid residues, heavy atoms only."""
    from .local_structure import THREE_TO_ONE  # local import keeps this module parser-free

    rows: list[tuple[str, int, str, str, str, str, tuple[float, float, float]]] = []
    for residue in structure[0].get_residues():
        if residue.id[0] != " " or residue.resname.upper() not in THREE_TO_ONE:
            continue
        chain = str(residue.get_parent().id)
        for atom in residue:
            element = str(atom.element or atom.get_id()[0]).strip().upper()
            if element in {"H", "D"}:
                continue
            coord = (float(atom.coord[0]), float(atom.coord[1]), float(atom.coord[2]))
            rows.append((
                chain,
                int(residue.id[1]),
                str(residue.id[2]).strip(),
                str(residue.resname),
                str(atom.get_id()).strip(),
                element,
                coord,
            ))
    return AtomTable.from_rows(rows)


def _residue_keys(table: AtomTable, indices: np.ndarray) -> set[tuple[str, int, str]]:
    return {(str(table.chain[i]), int(table.resnum[i]), str(table.icode[i])) for i in indices}


def _are_sequence_adjacent(a: tuple[str, int, str], b: tuple[str, int, str]) -> bool:
    return a[0] == b[0] and abs(a[1] - b[1]) <= 1


def contact_profile(
    table: AtomTable, chain: str, resnum: int, icode: str = "",
    side_chain_cutoff: float = 4.5,
) -> dict[str, Any]:
    """Describe the environment of the residue ``chain:resnum[icode]``."""
    mask = table.residue_mask(chain, resnum, icode)
    if not mask.any():
        raise ValueError(f"Residue {chain}:{resnum}{icode} has no atoms in the table.")
    key = (chain, resnum, icode)
    target_idx = np.flatnonzero(mask)
    others = np.flatnonzero(~mask)
    resname = str(table.resname[target_idx[0]])
    tree = cKDTree(table.coord[others])

    def partners(atom_i: int, radius: float) -> list[int]:
        return [int(others[j]) for j in tree.query_ball_point(table.coord[atom_i], radius)]

    # side-chain heavy-atom contacts (target side chain vs any heavy atom of other residues)
    side = [i for i in target_idx if table.atom[i] not in BACKBONE] or [i for i in target_idx if table.atom[i] == "CA"]
    contact_res: set[tuple[str, int, str]] = set()
    for i in side:
        contact_res |= _residue_keys(table, np.array(partners(int(i), side_chain_cutoff), dtype=int))
    contact_res.discard(key)

    # hydrogen-bond candidates: polar heavy atoms of target vs polar atoms elsewhere
    hbonds: set[tuple[str, str, int, str]] = set()
    for i in target_idx:
        if table.element[i] not in POLAR_ELEMENTS:
            continue
        for j in partners(int(i), HBOND_MAX):
            if table.element[j] not in POLAR_ELEMENTS:
                continue
            other = (str(table.chain[j]), int(table.resnum[j]), str(table.icode[j]))
            if _are_sequence_adjacent(key, other) and table.atom[i] in BACKBONE and table.atom[j] in BACKBONE:
                continue  # trivial peptide-bond geometry, not a contact of interest
            hbonds.add((str(table.atom[i]), str(table.atom[j]), other[1], other[0]))

    # salt-bridge candidates
    def charge_class(idx: int) -> str | None:
        pair = (str(table.resname[idx]), str(table.atom[idx]))
        return "acid" if pair in ACID_ATOMS else "base" if pair in BASE_ATOMS else None

    salt: set[tuple[str, int, str]] = set()
    for i in target_idx:
        mine = charge_class(int(i))
        if mine is None:
            continue
        for j in partners(int(i), SALT_BRIDGE_MAX):
            theirs = charge_class(j)
            if theirs is not None and theirs != mine:
                salt.add((str(table.chain[j]), int(table.resnum[j]), str(table.icode[j])))

    # disulfide
    disulfide: list[dict[str, Any]] = []
    if resname == "CYS":
        for i in target_idx:
            if table.atom[i] != "SG":
                continue
            for j in partners(int(i), DISULFIDE_MAX):
                if table.atom[j] == "SG" and table.resname[j] == "CYS":
                    d = float(np.linalg.norm(table.coord[i] - table.coord[j]))
                    disulfide.append({"chain": str(table.chain[j]), "residue": int(table.resnum[j]), "distance_angstrom": round(d, 2)})

    # chain-chain interface
    other_chain = [j for j in others if table.chain[j] != chain]
    interface_partners: set[tuple[str, int, str]] = set()
    for i in target_idx:
        for j in partners(int(i), INTERFACE_CUTOFF):
            if table.chain[j] != chain:
                interface_partners |= _residue_keys(table, np.array([j]))

    # CB-sphere neighbour count (CA stands in for Gly); standard burial proxy
    def representative(res_mask: np.ndarray) -> np.ndarray | None:
        idx = np.flatnonzero(res_mask)
        for name in ("CB", "CA"):
            hit = [i for i in idx if table.atom[i] == name]
            if hit:
                return table.coord[hit[0]]
        return None

    target_rep = representative(mask)
    cb_neighbors: int | None = None
    if target_rep is not None:
        reps = set()
        count = 0
        for res_key in sorted({(str(table.chain[i]), int(table.resnum[i]), str(table.icode[i])) for i in others}):
            if res_key in reps:
                continue
            reps.add(res_key)
            rep = representative(table.residue_mask(*res_key))
            if rep is not None and float(np.linalg.norm(rep - target_rep)) <= CB_SPHERE:
                count += 1
        cb_neighbors = count

    return {
        "side_chain_contact_cutoff_angstrom": side_chain_cutoff,
        "side_chain_contact_residue_count": len(contact_res),
        "hydrogen_bond_candidate_count": len(hbonds),
        "hydrogen_bond_candidate_max_distance_angstrom": HBOND_MAX,
        "salt_bridge_candidate_count": len(salt),
        "salt_bridge_candidate_max_distance_angstrom": SALT_BRIDGE_MAX,
        "disulfide_bonds": disulfide,
        "is_interface_residue": bool(interface_partners) and bool(other_chain),
        "interface_partner_residue_count": len(interface_partners),
        "interface_cutoff_angstrom": INTERFACE_CUTOFF,
        "cb_neighbors_within_10A": cb_neighbors,
        "scope": "first model, standard amino-acid heavy atoms only (waters, ligands, ions excluded)",
        "note": "Geometric candidates only; hydrogen positions, angles and energies are not evaluated.",
    }
