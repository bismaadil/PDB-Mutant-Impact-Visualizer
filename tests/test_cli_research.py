import csv
import json
from pathlib import Path

import numpy as np
from click.testing import CliRunner

from pdb_mutant_impact_visualizer.features import ALL_FEATURES
from pdb_mutant_impact_visualizer.research_cli import cli, render_dag_markdown
from pdb_mutant_impact_visualizer.stages import build_graph

ROOT = Path(__file__).resolve().parents[1]


def test_dag_doc_is_in_sync_with_code():
    committed = (ROOT / "docs" / "pipeline-dag.md").read_text(encoding="utf-8")
    assert committed == render_dag_markdown(build_graph()), "Run: pdbmi dag --write docs/pipeline-dag.md"


def test_cli_help_lists_every_command():
    out = CliRunner().invoke(cli, ["--help"]).output
    for name in ("run", "dag", "scan", "features", "train", "reproduce", "validate"):
        assert name in out


def test_run_requires_chain_and_rejects_bad_mutation_without_traceback():
    res = CliRunner().invoke(cli, ["run", "--pdb", "1AKI", "--chain", "A", "--mutation", "bad"])
    assert res.exit_code != 0 and "Traceback" not in res.output and "Mutation must" in res.output
    res = CliRunner().invoke(cli, ["run", "--pdb", "1AKI"])
    assert res.exit_code != 0 and "--chain and --mutation" in res.output


def synthetic_feature_csv(path, n_proteins=12, per=10, seed=0):
    rng = np.random.default_rng(seed)
    cols = ["pdb_id", "chain", "mutation", "group", "ddg", "structure_sha256", *ALL_FEATURES]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
        for p in range(n_proteins):
            for j in range(per):
                row = {c: round(float(rng.normal()), 4) for c in ALL_FEATURES}
                for flag in ("proline_introduced", "glycine_removed", "is_interface", "has_disulfide"):
                    row[flag] = 0
                row.update(pdb_id=f"P{p:03d}", chain="A", mutation=f"A{j + 1}V", group=f"P{p}", structure_sha256="0" * 64,
                           ddg=round(1.2 * row["relative_sasa"] - 0.6 * row["esm_llr"] + float(rng.normal(scale=0.2)), 4))
                w.writerow(row)


def test_train_command_end_to_end_on_synthetic_data(tmp_path):
    feats = tmp_path / "features.csv"
    synthetic_feature_csv(feats)
    out = tmp_path / "model"
    res = CliRunner().invoke(cli, ["train", "--features", str(feats), "--feature-groups", "substitution,structure,language_model",
                                   "--target-convention", "positive = destabilising", "--k", "4", "--bootstrap", "50", "--out", str(out)])
    assert res.exit_code == 0, res.output
    card = json.loads((out / "model.json").read_text())["model_card"]
    assert card["n_proteins"] == 12 and card["cross_validation"]["model"]["spearman"] > 0.7
    assert len(card["training_data_sha256"]) == 64
    assert "language_model_only" in (out / "ablation.md").read_text()


def test_train_requires_explicit_sign_convention(tmp_path):
    feats = tmp_path / "f.csv"; synthetic_feature_csv(feats)
    res = CliRunner().invoke(cli, ["train", "--features", str(feats)])
    assert res.exit_code != 0 and "target-convention" in res.output


def test_validate_command(tmp_path):
    good = {"schema_version": "2.0", "software": {}, "dependency_versions": {}, "analyzed_at_utc": "t",
            "structure": {"sha256": "a" * 64}, "chain": "A", "mutation": "A1V", "stability": None, "conservation": None,
            "local_structure": {}, "substitution": {}, "warnings": [], "run": {"config_sha256": "b" * 64}}
    p = tmp_path / "r.json"; p.write_text(json.dumps(good))
    assert CliRunner().invoke(cli, ["validate", str(p)]).exit_code == 0
    good["structure"]["sha256"] = "short"; p.write_text(json.dumps(good))
    res = CliRunner().invoke(cli, ["validate", str(p)])
    assert res.exit_code != 0 and "sha256" in res.output
