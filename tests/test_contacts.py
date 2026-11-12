import pytest

from pdb_mutant_impact_visualizer.contacts import AtomTable, contact_profile


def A(chain, num, res, atom, xyz, el=None):
    return (chain, num, "", res, atom, el or atom[0], xyz)


def test_hbond_salt_bridge_and_isolation():
    table = AtomTable.from_rows([
        A("A", 10, "ASP", "CA", (0, 0, 0)), A("A", 10, "ASP", "CB", (1, 0, 0)),
        A("A", 10, "ASP", "OD1", (2, 0, 0)),
        A("A", 20, "LYS", "CA", (5, 0, 0)), A("A", 20, "LYS", "NZ", (4.0 - 0.1, 0, 0)),  # 1.9 A from OD1 ... salt + hbond
        A("A", 30, "ALA", "CA", (80, 80, 80)), A("A", 30, "ALA", "CB", (81, 80, 80)),
    ])
    p = contact_profile(table, "A", 10)
    assert p["salt_bridge_candidate_count"] == 1
    assert p["hydrogen_bond_candidate_count"] == 1
    assert p["side_chain_contact_residue_count"] == 1
    assert p["cb_neighbors_within_10A"] == 1  # LYS only; ALA far away
    far = contact_profile(table, "A", 30)
    assert far["side_chain_contact_residue_count"] == 0 and far["hydrogen_bond_candidate_count"] == 0


def test_adjacent_backbone_pairs_are_not_reported_as_hbonds():
    table = AtomTable.from_rows([
        A("A", 1, "ALA", "N", (0, 0, 0)), A("A", 1, "ALA", "CA", (1, 0, 0)), A("A", 1, "ALA", "O", (0, 3, 0)),
        A("A", 2, "GLY", "N", (0, 3.0, 0.5)), A("A", 2, "GLY", "CA", (3, 3, 0)),
    ])
    assert contact_profile(table, "A", 1)["hydrogen_bond_candidate_count"] == 0


def test_disulfide_and_interface_detection():
    table = AtomTable.from_rows([
        A("A", 5, "CYS", "CA", (0, 0, 0)), A("A", 5, "CYS", "CB", (1, 0, 0)), A("A", 5, "CYS", "SG", (2, 0, 0)),
        A("A", 40, "CYS", "CA", (5, 0, 0)), A("A", 40, "CYS", "SG", (4.05, 0, 0)),
        A("B", 7, "ASN", "CA", (2, 4.5, 0)), A("B", 7, "ASN", "ND2", (2, 4.0, 0)),
    ])
    p = contact_profile(table, "A", 5)
    assert p["disulfide_bonds"] == [{"chain": "A", "residue": 40, "distance_angstrom": 2.05}]
    assert p["is_interface_residue"] is True
    assert p["interface_partner_residue_count"] == 1


def test_missing_residue_is_an_error():
    table = AtomTable.from_rows([A("A", 1, "ALA", "CA", (0, 0, 0))])
    with pytest.raises(ValueError, match="no atoms"):
        contact_profile(table, "A", 99)
