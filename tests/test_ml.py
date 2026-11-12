import numpy as np
import pytest

from pdb_mutant_impact_visualizer import ml


def synthetic(n_groups=12, per_group=15, seed=1, noise=0.3):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n_groups * per_group, 3))
    y = 1.5 * X[:, 0] - 0.8 * X[:, 1] + rng.normal(scale=noise, size=len(X))
    groups = [f"P{i}" for i in range(n_groups) for _ in range(per_group)]
    return X, y, groups, ["f0", "f1", "f2"]


def test_ridge_recovers_signal():
    X, y, _, names = synthetic()
    model = ml.RidgeModel.fit(X, y, names, alpha=0.1)
    assert ml.pearson(y, model.predict(X)) > 0.95
    assert abs(model.coefficients()["f2"]) < 0.2


def test_nan_training_data_is_rejected():
    X, y, _, names = synthetic()
    X[0, 0] = np.nan
    with pytest.raises(ValueError, match="NaN"):
        ml.RidgeModel.fit(X, y, names)


def test_grouped_folds_never_split_a_protein():
    _, _, groups, _ = synthetic()
    folds = ml.grouped_folds(groups, 4, seed=3)
    arr = np.array(groups)
    seen = []
    for f in folds:
        seen.append(set(arr[f]))
    assert sum(len(s) for s in seen) == len(set(groups))
    assert all(a.isdisjoint(b) for i, a in enumerate(seen) for b in seen[i + 1:])
    assert sorted(np.concatenate(folds).tolist()) == list(range(len(groups)))
    with pytest.raises(ValueError, match="distinct protein groups"):
        ml.grouped_folds(["a", "a", "b"], 5)


def test_cv_beats_mean_baseline_and_is_deterministic():
    X, y, groups, names = synthetic()
    a = ml.cross_validate(X, y, groups, names, k=4, seed=7)
    b = ml.cross_validate(X, y, groups, names, k=4, seed=7)
    assert a["model"] == b["model"]
    assert a["model"]["rmse"] < a["mean_baseline"]["rmse"]
    assert a["model"]["spearman"] > 0.8


def test_noise_labels_show_no_skill():
    X, _, groups, names = synthetic()
    y = np.random.default_rng(99).normal(size=len(X))
    cv = ml.cross_validate(X, y, groups, names, k=4, seed=0)
    assert abs(cv["model"]["spearman"]) < 0.3  # guards against leakage / overfitting claims


def test_metrics_and_bootstrap():
    y = np.array([0.0, 1.0, 2.0, 3.0, 4.0]); p = np.array([0.1, 0.9, 2.2, 2.8, 4.1])
    assert ml.spearman(y, p) == pytest.approx(1.0)
    assert ml.auroc(np.array([0, 0, 1, 1]), np.array([0.1, 0.4, 0.35, 0.8])) == pytest.approx(0.75)
    X, yy, groups, names = synthetic()
    cv = ml.cross_validate(X, yy, groups, names, k=4)
    lo, hi = ml.cluster_bootstrap_ci(yy, cv["out_of_fold_predictions"], groups, n_boot=200)
    assert lo <= cv["model"]["spearman"] <= hi


def test_save_load_roundtrip_and_prediction(tmp_path):
    X, y, _, names = synthetic()
    model = ml.RidgeModel.fit(X, y, names)
    path = ml.save_model(model, tmp_path / "m.json", {"target_sign_convention": "positive = destabilising"})
    loaded, card = ml.load_model(path)
    assert np.allclose(loaded.predict(X), model.predict(X))
    out = ml.predict_from_features(loaded, card, {"f0": 1.0, "f1": 0.0, "f2": 0.0})
    assert out["target_sign_convention"] == "positive = destabilising"
    with pytest.raises(ValueError, match="needs feature"):
        ml.predict_from_features(loaded, card, {"f0": 1.0, "f1": None, "f2": 0.0})
    path.write_text('{"schema": "other"}')
    with pytest.raises(ValueError, match="Unsupported"):
        ml.load_model(path)
