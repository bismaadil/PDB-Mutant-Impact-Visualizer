from pdb_mutant_impact_visualizer.features import ALL_FEATURES, FEATURE_GROUPS, extract_features, select


def report(**over):
    base = {
        "mutation": "A50V",
        "substitution": {"blosum62_score": 0, "hydropathy_change_kyte_doolittle": 2.4, "volume_change_angstrom3": 51.4,
                         "charge_change_ph7_approx": 0, "proline_introduced": False, "glycine_removed": False},
        "local_structure": {
            "neighbor_residue_count": 12, "neighbor_cutoff_angstrom": 10.0,
            "site_mean_b_factor_angstrom2": 20.0, "protein_mean_b_factor_angstrom2": 25.0,
            "solvent_accessibility": {"relative_sasa": 0.12},
            "contacts": {"cb_neighbors_within_10A": 14, "side_chain_contact_residue_count": 7,
                         "hydrogen_bond_candidate_count": 2, "salt_bridge_candidate_count": 0,
                         "is_interface_residue": False, "disulfide_bonds": []},
        },
    }
    base.update(over)
    return base


def test_full_report_yields_every_feature():
    f = extract_features(report(language_model={"esm_llr": -1.5}))
    assert set(f) == set(ALL_FEATURES)
    assert all(v is not None for v in f.values())
    assert f["site_b_ratio"] == 0.8 and f["esm_llr"] == -1.5 and f["wt_volume"] == 88.6


def test_missing_measurements_stay_none_not_zero():
    f = extract_features(report(language_model=None))
    assert f["esm_llr"] is None
    r = report(); r["local_structure"].pop("contacts")
    f = extract_features(r)
    assert f["side_chain_contacts"] is None and f["has_disulfide"] is None


def test_nonstandard_cutoff_does_not_masquerade_as_10A_count():
    r = report(); r["local_structure"]["neighbor_cutoff_angstrom"] = 8.0
    assert extract_features(r)["neighbor_count_10A"] is None


def test_group_selection():
    assert select({}, ["substitution", "language_model"]) == FEATURE_GROUPS["substitution"] + ["esm_llr"]
