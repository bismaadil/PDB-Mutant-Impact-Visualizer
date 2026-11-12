"""Reconstruct reference-aligned protein sequences from RCSB mmCIF metadata."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from .local_structure import THREE_TO_ONE


class MMCIFSequenceError(RuntimeError):
    """Raised when a complete, unambiguous chain sequence cannot be reconstructed."""


def _values(data: dict[str, Any], category: str, field: str) -> list[str]:
    key = f"{category}.{field}"
    if key not in data:
        raise MMCIFSequenceError(f"Required mmCIF field is missing: {key}")
    value = data[key]
    return value if isinstance(value, list) else [value]


def _clean(value: str) -> str | None:
    value = str(value).strip()
    return None if value in {"", ".", "?"} else value


def reconstruct_chain_sequence(
    cif_path: str | Path,
    chain_id: str,
) -> tuple[str, dict[tuple[int, str], int], dict[str, Any]]:
    """Return aligned construct sequence, author-residue mapping and provenance.

    The mapping values are zero-based sequence indices. Unresolved polymer
    residues remain in the sequence but have no author-residue mapping.
    """
    try:
        data = MMCIF2Dict(str(cif_path))
    except Exception as exc:
        raise MMCIFSequenceError(f"Could not parse mmCIF file {cif_path}: {exc}") from exc

    # Find reference alignments explicitly associated with the requested chain.
    ref_chains = _values(data, "_struct_ref_seq", "pdbx_strand_id")
    ref_begins = _values(data, "_struct_ref_seq", "seq_align_beg")
    ref_ends = _values(data, "_struct_ref_seq", "seq_align_end")
    ref_accessions = _values(data, "_struct_ref_seq", "pdbx_db_accession")

    alignments: list[tuple[int, int, str | None]] = []
    for strand, begin, end, accession in zip(
        ref_chains, ref_begins, ref_ends, ref_accessions, strict=True
    ):
        strands = {part.strip() for part in strand.split(",")}
        if chain_id in strands:
            try:
                start, stop = int(begin), int(end)
            except ValueError as exc:
                raise MMCIFSequenceError("Invalid reference-alignment boundaries.") from exc
            if start < 1 or stop < start:
                raise MMCIFSequenceError("Invalid reference-alignment boundaries.")
            alignments.append((start, stop, _clean(accession)))

    if len(alignments) != 1:
        raise MMCIFSequenceError(
            f"Expected exactly one reference alignment for chain {chain_id!r}; "
            f"found {len(alignments)}. Construct boundaries are ambiguous."
        )

    start, stop, reference_accession = alignments[0]

    # Find the polymer-scheme rows belonging to the requested author chain.
    fields = (
        "asym_id", "entity_id", "seq_id", "mon_id",
        "auth_seq_num", "pdb_ins_code", "pdb_strand_id",
    )
    columns = {
        field: _values(data, "_pdbx_poly_seq_scheme", field)
        for field in fields
    }

    groups: dict[tuple[str, str], list[dict[str, str]]] = {}
    for scheme_values in zip(*(columns[field] for field in fields), strict=True):
        record = dict(zip(fields, scheme_values, strict=True))
        if record["pdb_strand_id"] == chain_id:
            key = (record["asym_id"], record["entity_id"])
            groups.setdefault(key, []).append(record)

    candidates = []
    for (asym_id, entity_id), rows in groups.items():
        rows.sort(key=lambda row: int(row["seq_id"]))
        positions = [int(row["seq_id"]) for row in rows]

        if len(positions) != len(set(positions)):
            continue
        aligned = [
            row for row in rows
            if start <= int(row["seq_id"]) <= stop
        ]
        if (
            not aligned
            or int(aligned[0]["seq_id"]) != start
            or int(aligned[-1]["seq_id"]) != stop
            or len(aligned) != stop - start + 1
        ):
            continue

        try:
            full_sequence = "".join(THREE_TO_ONE[row["mon_id"].upper()] for row in rows)
        except KeyError:
            continue

        entity_ids = _values(data, "_entity_poly", "entity_id")
        entity_sequences = _values(data, "_entity_poly", "pdbx_seq_one_letter_code")
        matching_entity_sequences = [
            re.sub(r"\s+", "", sequence).upper()
            for entity, sequence in zip(entity_ids, entity_sequences, strict=True)
            if entity == entity_id
        ]
        if matching_entity_sequences != [full_sequence]:
            continue

        candidates.append((asym_id, entity_id, rows, aligned))

    if len(candidates) != 1:
        raise MMCIFSequenceError(
            f"Could not identify exactly one polymer sequence for chain {chain_id!r} "
            f"and alignment {start}–{stop}; found {len(candidates)} valid candidates."
        )

    asym_id, entity_id, all_rows, aligned_rows = candidates[0]
    sequence = "".join(
        THREE_TO_ONE[scheme_row["mon_id"].upper()] for scheme_row in aligned_rows
    )

    # Map only real author residue identifiers. Never invent IDs for missing residues.
    residue_to_index: dict[tuple[int, str], int] = {}
    for index, aligned_row in enumerate(aligned_rows):
        author_number = _clean(aligned_row["auth_seq_num"])
        if author_number is None:
            continue

        try:
            number = int(author_number)
        except ValueError as exc:
            raise MMCIFSequenceError(
                f"Invalid author residue number {author_number!r}."
            ) from exc

        insertion_code = _clean(aligned_row["pdb_ins_code"]) or ""
        key = (number, insertion_code)
        if key in residue_to_index:
            raise MMCIFSequenceError(
                f"Duplicate author residue identifier {key!r} in chain {chain_id!r}."
            )
        residue_to_index[key] = index

    metadata = {
        "sequence_source": "RCSB mmCIF polymer sequence",
        "chain": chain_id,
        "label_asym_id": asym_id,
        "entity_id": entity_id,
        "reference_accession": reference_accession,
        "polymer_alignment_start": start,
        "polymer_alignment_end": stop,
        "sequence_length": len(sequence),
        "unresolved_residues_included": sum(
            1 for row in aligned_rows if _clean(row["auth_seq_num"]) is None
        ),
    }
    return sequence, residue_to_index, metadata
