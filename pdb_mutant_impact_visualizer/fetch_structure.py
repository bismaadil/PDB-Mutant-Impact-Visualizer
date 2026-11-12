"""Download and parse experimental PDB coordinates with explicit provenance."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

import requests
from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure


class StructureError(RuntimeError):
    pass


def download_pdb(pdb_id: str, out_dir: str | Path = "data/structures", refresh: bool = False) -> Path:
    pdb_id = pdb_id.strip().upper()
    if not re.fullmatch(r"[A-Z0-9]{4}", pdb_id):
        raise StructureError("PDB identifier must contain exactly four letters or digits.")
    destination = Path(out_dir) / f"{pdb_id.lower()}.pdb"
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not refresh and destination.exists() and destination.stat().st_size > 0:
        return destination
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        response = requests.get(url, timeout=(10, 60))
        response.raise_for_status()
    except requests.RequestException as exc:
        raise StructureError(f"Could not download {pdb_id} from RCSB PDB: {exc}") from exc
    if not response.text.startswith(("HEADER", "TITLE", "REMARK", "ATOM", "MODEL", "COMPND")):
        raise StructureError(f"RCSB response for {pdb_id} did not look like a PDB coordinate file.")
    destination.write_bytes(response.content)
    return destination


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_structure(path: str | Path, structure_id: str = "input") -> Structure:
    path = Path(path)
    if not path.is_file():
        raise StructureError(f"Structure file does not exist: {path}")
    try:
        return PDBParser(QUIET=True).get_structure(structure_id, str(path))
    except Exception as exc:
        raise StructureError(f"Could not parse PDB coordinates in {path}: {exc}") from exc


def download_mmcif(
    pdb_id: str,
    out_dir: str | Path = "data/structures",
    refresh: bool = False,
) -> Path:
    """Download the RCSB mmCIF file used for polymer-sequence reconstruction."""
    pdb_id = pdb_id.strip().upper()
    if not re.fullmatch(r"[A-Z0-9]{4}", pdb_id):
        raise StructureError("PDB identifier must contain exactly four letters or digits.")

    destination = Path(out_dir) / f"{pdb_id.lower()}.cif"
    destination.parent.mkdir(parents=True, exist_ok=True)

    if not refresh and destination.exists() and destination.stat().st_size > 0:
        return destination

    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    try:
        response = requests.get(url, timeout=(10, 60))
        response.raise_for_status()
    except requests.RequestException as exc:
        raise StructureError(f"Could not download mmCIF for {pdb_id} from RCSB: {exc}") from exc

    if not response.text.lstrip().lower().startswith("data_"):
        raise StructureError(
            f"RCSB response for {pdb_id} did not look like an mmCIF file."
        )

    destination.write_bytes(response.content)
    return destination
