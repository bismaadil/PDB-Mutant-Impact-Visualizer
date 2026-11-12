from pathlib import Path

import pytest

from pdb_mutant_impact_visualizer.mmcif_sequence import (
    MMCIFSequenceError,
    reconstruct_chain_sequence,
)


def write_cif(path: Path, *, duplicate_alignment: bool = False) -> None:
    alignments = """\
1 A 2 6 P12345
"""
    if duplicate_alignment:
        alignments += "2 A 2 6 Q67890\n"

    path.write_text(
        f"""\
data_test
#
loop_
_entity_poly.entity_id
_entity_poly.pdbx_seq_one_letter_code
1
;MAGSTV
;
#
loop_
_struct_ref_seq.align_id
_struct_ref_seq.pdbx_strand_id
_struct_ref_seq.seq_align_beg
_struct_ref_seq.seq_align_end
_struct_ref_seq.pdbx_db_accession
{alignments}#
loop_
_pdbx_poly_seq_scheme.asym_id
_pdbx_poly_seq_scheme.entity_id
_pdbx_poly_seq_scheme.seq_id
_pdbx_poly_seq_scheme.mon_id
_pdbx_poly_seq_scheme.auth_seq_num
_pdbx_poly_seq_scheme.pdb_ins_code
_pdbx_poly_seq_scheme.pdb_strand_id
A 1 1 MET -1 . A
A 1 2 ALA 10 . A
A 1 3 GLY 11 . A
A 1 4 SER ? . A
A 1 5 THR 12 B A
A 1 6 VAL 13 . A
#
""",
        encoding="utf-8",
    )


def test_reconstructs_aligned_sequence_and_keeps_unresolved_residue(tmp_path):
    cif = tmp_path / "test.cif"
    write_cif(cif)

    sequence, mapping, metadata = reconstruct_chain_sequence(cif, "A")

    assert sequence == "AGSTV"
    assert (10, "") in mapping
    assert (11, "") in mapping
    assert (12, "B") in mapping
    assert (13, "") in mapping
    assert len(mapping) == 4
    assert metadata["unresolved_residues_included"] == 1
    assert metadata["polymer_alignment_start"] == 2
    assert metadata["polymer_alignment_end"] == 6


def test_excludes_tag_and_maps_author_residue_numbers(tmp_path):
    cif = tmp_path / "test.cif"
    write_cif(cif)

    sequence, mapping, _ = reconstruct_chain_sequence(cif, "A")

    assert sequence[0] == "A"
    assert mapping[(10, "")] == 0
    assert mapping[(12, "B")] == 3


def test_rejects_ambiguous_reference_alignments(tmp_path):
    cif = tmp_path / "ambiguous.cif"
    write_cif(cif, duplicate_alignment=True)

    with pytest.raises(MMCIFSequenceError, match="exactly one reference alignment"):
        reconstruct_chain_sequence(cif, "A")
