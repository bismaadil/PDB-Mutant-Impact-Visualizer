import types
import zlib

import numpy as np
import pytest

from pdb_mutant_impact_visualizer import research
from pdb_mutant_impact_visualizer.features import FEATURE_GROUPS


def fake_result(cfg, extra_context=None):
    """Fake pipeline: report features depend deterministically on the mutation string."""
    rng = np.random.default_rng(zlib.crc32(f"{cfg.pdb_id}{cfg.mutation}".encode()))
    pos = int(cfg.mutation[1:-1])
    rsasa = (pos % 10) / 10
    report = dict(
        schema_version="2.0", software={}, dependency_versions={}, analyzed_at_utc="t",
        structure={"sha256": "ab" * 32, "pdb_id": cfg.pdb_id}, chain=cfg.chain, mutation=cfg.mutation,
        stability=None, conservation=None,
        local_structure={
            "neighbor_residue_count": 10 + pos % 7, "neighbor_cutoff_angstrom": 10.0,
            "site_mean_b_factor_angstrom2": 20.0, "protein_mean_b_factor_angstrom2": 25.0,
            "solvent_accessibility": {"relative_sasa": rsasa},
            "contacts": {"cb_neighbors_within_10A": 10, "side_chain_contact_residue_count": 5 + pos % 4,
                         "hydrogen_bond_candidate_count": 1, "salt_bridge_candidate_count": 0,
                         "is_interface_residue": False, "disulfide_bonds": []}},
        substitution={"blosum62_score": -1, "hydropathy_change_kyte_doolittle": float(rng.normal()),
                      "volume_change_angstrom3": float(rng.normal() * 30), "charge_change_ph7_approx": 0,
                      "proline_introduced": False, "glycine_removed": False},
        language_model={"esm_llr": float(rng.normal())},
        run={"config_sha256": cfg.sha256(), "config": cfg.public_dict()},
    )
    return types.SimpleNamespace(to_dict=lambda: report)


def test_feature_building_records_failures_without_aborting(tmp_path):
    rows = [{"pdb_id": "1ABC", "chain": "A", "mutation": "A10V", "ddg": "1.2"},
            {"pdb_id": "1ABC", "chain": "A", "mutation": "BAD", "ddg": "0.1"}]
    good, bad = research.build_feature_rows(rows, run_fn=fake_result)
    assert len(good) == 1 and len(bad) == 1 and "BAD" in bad[0]["row"]
    path = research.write_csv(good, tmp_path / "f.csv", research.FEATURE_TABLE_COLUMNS)
    back = research.read_feature_table(path)
    assert back[0]["ddg"] == "1.2" and back[0]["structure_sha256"] == "ab" * 32


def test_training_pipeline_learns_planted_signal_with_grouped_cv_and_counts_drops():
    rows = []
    for p in range(14):
        for j in range(12):
            r = {"pdb_id": "1ABC", "chain": "A", "mutation": f"A{10 + p * 12 + j}V", "group": f"P{p}"}
            good, _ = research.build_feature_rows([r], run_fn=fake_result)
            f = good[0]
            f["group"] = f"P{p}"
            f["ddg"] = str(3.0 * (1 - float(f["relative_sasa"])) - 0.5 * float(f["esm_llr"]))  # planted signal
            rows.append(f)
    rows[0]["ddg"] = ""          # dropped: no target
    rows[1]["relative_sasa"] = ""  # dropped: missing feature
    model, card = research.train_from_rows(rows, ["substitution", "structure", "language_model"], k=4, n_boot=100,
                                           target_convention="positive = destabilising", dataset_sha256="abc")
    assert card["dropped_rows"] == {"no_target": 1, "missing_feature": 1}
    assert card["cross_validation"]["model"]["spearman"] > 0.8
    assert card["cross_validation"]["model"]["rmse"] < card["cross_validation"]["mean_baseline"]["rmse"]
    assert card["n_proteins"] == 14 and card["training_data_sha256"] == "abc"
    table = research.ablation_markdown(research.run_ablation(rows, list(FEATURE_GROUPS), k=4, n_boot=50,
                                                             target_convention="x"))
    assert "substitution+structure+language_model" in table and "language_model_only" in table


def test_training_refuses_tiny_data():
    with pytest.raises(ValueError, match="at least 10"):
        research.train_from_rows([{"ddg": "1", "group": "a"}] * 3, ["substitution"], target_convention="x")


def test_positions_parsing():
    assert research.parse_positions("50") == [50]
    assert research.parse_positions("45-47,50,46") == [45, 46, 47, 50]
    for bad in ("", "a", "9-3", "0-2"):
        with pytest.raises(ValueError):
            research.parse_positions(bad)


def test_saturation_scan_covers_19_substitutions_per_position(tmp_path):
    from pdb_mutant_impact_visualizer.config import RunConfig
    base = RunConfig(chain="A", mutation="A1V", pdb_id="1ABC")
    rows = research.saturation_scan(base, [12, 13], lambda p: "A" if p == 12 else "G", run_fn=fake_result)
    assert len(rows) == 38 and not any(r["error"] for r in rows)
    assert {r["mutant"] for r in rows if r["position"] == 12} == set("CDEFGHIKLMNPQRSTVWY")
    png = research.scan_heatmap(rows, "esm_llr", tmp_path / "h.png")
    assert png.stat().st_size > 1000


def test_reproduce_detects_drift_and_matches_when_identical(tmp_path):
    from pdb_mutant_impact_visualizer.config import RunConfig
    cfg = RunConfig(chain="A", mutation="A10V", pdb_id="1ABC").validate()
    original = fake_result(cfg).to_dict()
    out = research.reproduce(original, run_fn=lambda c: fake_result(c))
    assert out["reproduced"] is True
    original["structure"]["sha256"] = "00" * 32  # coordinates changed upstream
    out = research.reproduce(original, run_fn=lambda c: fake_result(c))
    assert out["reproduced"] is False and "structure.sha256" in out["differences"][0]
    with pytest.raises(ValueError, match="no run manifest"):
        research.config_from_report({"mutation": "A1V"})
