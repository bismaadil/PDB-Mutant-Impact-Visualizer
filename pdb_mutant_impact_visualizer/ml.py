"""Small, auditable supervised-learning toolkit for the research layer.

Design choices (deliberate, defensible in a methods section):
* NumPy-only ridge regression: closed form, deterministic, no hidden state.
* **Protein-grouped** cross-validation: all variants of one protein stay in one
  fold, which prevents the standard leakage that inflates stability-predictor
  scores.
* **Cluster bootstrap** confidence intervals (resampling proteins, not rows).
* Baselines are first-class citizens: a model must beat a mean predictor and
  the best single feature to be called useful.
* A saved model records training-data SHA-256, feature order, CV metrics and
  the sign convention of the target. No pretrained weights are bundled.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Sequence

import numpy as np
from scipy import stats

MODEL_SCHEMA = "pdbmi-ridge-1"


# ---------------------------------------------------------------- metrics
def pearson(y: np.ndarray, p: np.ndarray) -> float:
    if len(y) < 3 or np.allclose(y, y[0]) or np.allclose(p, p[0]):
        return float("nan")
    return float(np.corrcoef(y, p)[0, 1])


def spearman(y: np.ndarray, p: np.ndarray) -> float:
    if len(y) < 3 or np.allclose(y, y[0]) or np.allclose(p, p[0]):
        return float("nan")
    return float(stats.spearmanr(y, p).statistic)


def rmse(y: np.ndarray, p: np.ndarray) -> float:
    return float(np.sqrt(np.mean((y - p) ** 2)))


def mae(y: np.ndarray, p: np.ndarray) -> float:
    return float(np.mean(np.abs(y - p)))


def auroc(labels: np.ndarray, scores: np.ndarray) -> float:
    """Rank-based AUROC (Mann-Whitney U); ties get average rank."""
    labels = np.asarray(labels).astype(bool)
    pos, neg = int(labels.sum()), int((~labels).sum())
    if pos == 0 or neg == 0:
        return float("nan")
    ranks = stats.rankdata(scores)
    return float((ranks[labels].sum() - pos * (pos + 1) / 2) / (pos * neg))


METRICS: dict[str, Callable[[np.ndarray, np.ndarray], float]] = {
    "pearson": pearson, "spearman": spearman, "rmse": rmse, "mae": mae,
}


def summarize(y: np.ndarray, p: np.ndarray, destabilising_threshold: float = 1.0) -> dict[str, float]:
    out = {name: round(fn(y, p), 4) for name, fn in METRICS.items()}
    out["auroc_destabilising"] = round(auroc(y >= destabilising_threshold, p), 4)
    out["n"] = int(len(y))
    return out


# ------------------------------------------------------------------ model
@dataclass
class RidgeModel:
    feature_names: list[str]
    mean: np.ndarray
    scale: np.ndarray
    weights: np.ndarray
    intercept: float
    alpha: float

    @classmethod
    def fit(cls, X: np.ndarray, y: np.ndarray, feature_names: Sequence[str], alpha: float = 1.0) -> "RidgeModel":
        if X.ndim != 2 or X.shape[0] != len(y):
            raise ValueError("X must be 2-D with one row per target value.")
        if not np.isfinite(X).all() or not np.isfinite(y).all():
            raise ValueError("Training data contain NaN/inf; drop or impute rows explicitly before fitting.")
        mean = X.mean(axis=0)
        scale = X.std(axis=0)
        scale[scale == 0] = 1.0
        Z = (X - mean) / scale
        y_mean = float(y.mean())
        gram = Z.T @ Z + alpha * np.eye(Z.shape[1])
        w = np.linalg.solve(gram, Z.T @ (y - y_mean))
        return cls(list(feature_names), mean, scale, w, y_mean, alpha)

    def predict(self, X: np.ndarray) -> np.ndarray:
        return ((np.atleast_2d(X) - self.mean) / self.scale) @ self.weights + self.intercept

    def coefficients(self) -> dict[str, float]:
        """Standardised coefficients (comparable magnitudes across features)."""
        return {n: round(float(w), 5) for n, w in zip(self.feature_names, self.weights, strict=True)}


class MeanBaseline:
    """Predicts the training mean; the floor every model must clear."""

    def __init__(self, value: float):
        self.value = value

    @classmethod
    def fit(cls, X: np.ndarray, y: np.ndarray, *_: Any, **__: Any) -> "MeanBaseline":
        return cls(float(np.mean(y)))

    def predict(self, X: np.ndarray) -> np.ndarray:
        return np.full(len(np.atleast_2d(X)), self.value)


# ------------------------------------------------------------- evaluation
def grouped_folds(groups: Sequence[str], k: int, seed: int = 0) -> list[np.ndarray]:
    """Assign whole groups to ``k`` folds, balancing sizes greedily (deterministic)."""
    unique = sorted(set(groups))
    if k < 2:
        raise ValueError("k must be at least 2.")
    if len(unique) < k:
        raise ValueError(f"Need at least {k} distinct protein groups for {k}-fold grouped CV; found {len(unique)}.")
    rng = np.random.default_rng(seed)
    order = list(rng.permutation(len(unique)))
    sizes = {g: sum(1 for x in groups if x == g) for g in unique}
    order.sort(key=lambda i: -sizes[unique[i]])  # largest groups first (stable after shuffle)
    fold_members: list[list[str]] = [[] for _ in range(k)]
    fold_size = [0] * k
    for i in order:
        g = unique[i]
        target = int(np.argmin(fold_size))
        fold_members[target].append(g)
        fold_size[target] += sizes[g]
    arr = np.array(groups)
    return [np.flatnonzero(np.isin(arr, members)) for members in fold_members]


def cross_validate(
    X: np.ndarray, y: np.ndarray, groups: Sequence[str], feature_names: Sequence[str],
    k: int = 5, seed: int = 0, alpha: float = 1.0,
) -> dict[str, Any]:
    folds = grouped_folds(groups, k, seed)
    oof = np.full(len(y), np.nan)
    base = np.full(len(y), np.nan)
    for test_idx in folds:
        train_idx = np.setdiff1d(np.arange(len(y)), test_idx)
        model = RidgeModel.fit(X[train_idx], y[train_idx], feature_names, alpha)
        oof[test_idx] = model.predict(X[test_idx])
        base[test_idx] = MeanBaseline.fit(X[train_idx], y[train_idx]).predict(X[test_idx])
    single = {
        name: round(spearman(y, X[:, j]), 4) for j, name in enumerate(feature_names)
    }  # sign-agnostic ranking: |rho| matters; no fitting, so no leakage
    return {
        "k": k, "seed": seed, "alpha": alpha,
        "model": summarize(y, oof), "mean_baseline": summarize(y, base),
        "single_feature_spearman_with_target": single,
        "out_of_fold_predictions": oof,
        "fold_sizes": [int(len(f)) for f in folds],
    }


def cluster_bootstrap_ci(
    y: np.ndarray, p: np.ndarray, groups: Sequence[str], metric: str = "spearman",
    n_boot: int = 1000, seed: int = 0, level: float = 0.95,
) -> tuple[float, float]:
    """Resample proteins (clusters) with replacement; percentile interval."""
    fn = METRICS[metric]
    arr = np.array(groups)
    unique = sorted(set(groups))
    index = {g: np.flatnonzero(arr == g) for g in unique}
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(n_boot):
        pick = rng.choice(len(unique), size=len(unique), replace=True)
        rows = np.concatenate([index[unique[i]] for i in pick])
        v = fn(y[rows], p[rows])
        if not math.isnan(v):
            values.append(v)
    if not values:
        return float("nan"), float("nan")
    lo, hi = np.percentile(values, [(1 - level) / 2 * 100, (1 + level) / 2 * 100])
    return round(float(lo), 4), round(float(hi), 4)


# ----------------------------------------------------------- persistence
def file_sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def save_model(model: RidgeModel, path: str | Path, card: dict[str, Any]) -> Path:
    payload = {
        "schema": MODEL_SCHEMA,
        "feature_names": model.feature_names,
        "mean": model.mean.tolist(), "scale": model.scale.tolist(),
        "weights": model.weights.tolist(), "intercept": model.intercept, "alpha": model.alpha,
        "model_card": card,
    }
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return out


def load_model(path: str | Path) -> tuple[RidgeModel, dict[str, Any]]:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if data.get("schema") != MODEL_SCHEMA:
        raise ValueError(f"Unsupported model schema {data.get('schema')!r}; expected {MODEL_SCHEMA!r}.")
    model = RidgeModel(
        data["feature_names"], np.array(data["mean"]), np.array(data["scale"]),
        np.array(data["weights"]), float(data["intercept"]), float(data["alpha"]),
    )
    return model, data["model_card"]


def predict_from_features(
    model: RidgeModel, card: dict[str, Any], features: dict[str, float | None]
) -> dict[str, Any]:
    values: list[float] = []
    missing: list[str] = []
    for name in model.feature_names:
        value = features.get(name)
        if value is None:
            missing.append(name)
        else:
            values.append(value)
    if missing:
        raise ValueError(f"Model needs feature(s) that this run did not produce: {missing}")
    x = np.array([values])
    return {
        "backend": "pdbmi ridge (research model)",
        "ddg_kcal_per_mol": round(float(model.predict(x)[0]), 4),
        "target_sign_convention": card.get("target_sign_convention", "unspecified"),
        "trained_on_sha256": card.get("training_data_sha256"),
        "cv_spearman_grouped": card.get("cross_validation", {}).get("model", {}).get("spearman"),
        "feature_names": model.feature_names,
        "note": "Research model trained by the user on a stated dataset; see its model card. Not validated for clinical use.",
    }
