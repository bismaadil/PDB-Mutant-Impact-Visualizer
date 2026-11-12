import pytest

from pdb_mutant_impact_visualizer.descriptors import solvent_accessibility, structure_warnings, substitution_descriptors
from pdb_mutant_impact_visualizer.fetch_structure import load_structure
from pdb_mutant_impact_visualizer.local_structure import analyze_local_structure, get_residue
from pdb_mutant_impact_visualizer.models import Mutation
from pdb_mutant_impact_visualizer.pipeline import analyze_mutation

ATOMS = """ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 10.00           N
ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00 10.00           C
ATOM      3  C   ALA A   1       2.000   1.400   0.000  1.00 10.00           C
ATOM      4  N   GLY A   2       3.300   1.600   0.000  1.00 12.00           N
ATOM      5  CA  GLY A   2       3.900   2.900   0.000  1.00 12.00           C
ATOM      6  C   GLY A   2       5.400   2.900   0.000  1.00 12.00           C
"""


def write(tmp_path, body, name="s.pdb"):
    path = tmp_path / name
    path.write_text(body)
    return path


def test_water_with_same_number_does_not_make_residue_ambiguous(tmp_path):
    water = "HETATM  900  O   HOH A   1       2.000   2.000   2.000  1.00 30.00           O\n"
    structure = load_structure(write(tmp_path, ATOMS + water + "END\n"), "t")
    assert get_residue(structure, "A", 1).resname == "ALA"


def test_waters_and_ligands_are_not_counted_as_neighbors(tmp_path):
    water = "HETATM  900  O   HOH A 101       2.000   2.000   2.000  1.00 30.00           O\n"
    plain = load_structure(write(tmp_path, ATOMS + "END\n", "a.pdb"), "a")
    wet = load_structure(write(tmp_path, ATOMS + water + "END\n", "b.pdb"), "b")
    assert (analyze_local_structure(plain, "A", 1)["neighbor_residue_count"]
            == analyze_local_structure(wet, "A", 1)["neighbor_residue_count"])


def test_only_first_model_is_used_and_warned(tmp_path):
    body = "MODEL        1\n" + ATOMS + "ENDMDL\nMODEL        2\n" + ATOMS + "ENDMDL\nEND\n"
    structure = load_structure(write(tmp_path, body), "nmr")
    assert analyze_local_structure(structure, "A", 1)["neighbor_residue_count"] == 1
    assert any("models" in w for w in structure_warnings(structure, "A", 1))


def test_substitution_descriptors_are_deterministic():
    d = substitution_descriptors(Mutation.parse("A50V"))
    assert d["blosum62_score"] == 0
    assert d["volume_change_angstrom3"] == pytest.approx(51.4)
    assert d["hydropathy_change_kyte_doolittle"] == pytest.approx(2.4)
    assert substitution_descriptors(Mutation.parse("D5K"))["charge_change_ph7_approx"] == 2
    assert substitution_descriptors(Mutation.parse("G5P"))["proline_introduced"] is True


def test_sasa_is_reported_and_bounded(structure):
    s = solvent_accessibility(structure, "A", 50)
    assert s["sasa_angstrom2"] >= 0 and 0 <= s["relative_sasa"] <= 1


def test_ebi_email_checked_before_any_remote_work(fixture_path):
    with pytest.raises(ValueError, match="email"):
        analyze_mutation(mutation=Mutation.parse("A50V"), chain="A", structure_path=fixture_path,
                         conservation_backend="ebi", uniprot_accession="P00698", uniprot_position=50,
                         stability_backend="dynamut2")
