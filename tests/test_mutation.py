import pytest

from pdb_mutant_impact_visualizer.local_structure import validate_mutation
from pdb_mutant_impact_visualizer.models import Mutation


@pytest.mark.parametrize("text,expected", [("a50v", "A50V"), ("S7G", "S7G")])
def test_parse_mutation(text, expected):
    assert str(Mutation.parse(text)) == expected


@pytest.mark.parametrize("text", ["A0V", "A50A", "B4V", "A-1V", "A50"])
def test_reject_invalid_mutation(text):
    with pytest.raises(ValueError):
        Mutation.parse(text)


def test_validate_reference_amino_acid(structure):
    residue = validate_mutation(structure, "A", Mutation.parse("A50V"))
    assert residue.resname == "ALA"
    with pytest.raises(ValueError, match="does not match"):
        validate_mutation(structure, "A", Mutation.parse("S50V"))


def test_missing_chain_and_position_are_actionable(structure):
    with pytest.raises(ValueError, match="Chain"):
        validate_mutation(structure, "Z", Mutation.parse("A50V"))
    with pytest.raises(ValueError, match="not found"):
        validate_mutation(structure, "A", Mutation.parse("A999V"))
