"""Research workflows built on the pipeline: features, scans, training, reproduction.

Every function takes an injectable ``run_fn`` (default: the real pipeline) so the
logic can be tested offline and so that alternative runners can be swapped in.
"""

from __future__ import annotations

import csv
import dataclasses
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence

import numpy as np

from . import ml
from .config import RunConfig
from .features import ALL_FEATURES, FEATURE_GROUPS, extract_features
from .models import AA
from .provenance import compare_reports

RunFn = Callable[..., Any]
REQUIRED_COLUMNS = ("pdb_id", "chain", "mutation")


def _default_run() -> RunFn:
    from .pipeline import run_pipeline

    return run_pipeline


# ---------------------------------------------------------------- datasets
def read_mutation_table(path: str | Path) -> list[dict[str, str]]:
    """CSV with columns pdb_id, chain, mutation and optionally ddg, group."""
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = [c for c in REQUIRED_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"Mutation table lacks required column(s): {missing}")
        rows = [{k.strip(): (v or "").strip() for k, v in row.items() if k} for row in reader]
    if not rows:
        raise ValueError("Mutation table is empty.")
    return rows


def build_feature_rows(
    rows: Iterable[dict[str, str]], base: dict[str, Any] | None = None,
    run_fn: RunFn | None = None, extra_context: dict[str, Any] | None = None,
    progress: Callable[[str], None] | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    """Run the pipeline for each row; return (feature rows, failures). Never aborts the batch."""
    run = run_fn or _default_run()
    good: list[dict[str, Any]] = []
    bad: list[dict[str, str]] = []
    for row in rows:
        label = f"{row['pdb_id']}:{row['chain']}:{row['mutation']}"
        try:
            cfg = RunConfig(chain=row["chain"], mutation=row["mutation"], pdb_id=row["pdb_id"], **(base or {}))
            result = run(cfg, extra_context=extra_context) if extra_context else run(cfg)
            report = result.to_dict()
            feats = extract_features(report)
            good.append({
                "pdb_id": row["pdb_id"].upper(), "chain": row["chain"], "mutation": row["mutation"],
                "group": row.get("group") or row["pdb_id"].upper(),
                "ddg": row.get("ddg", ""), "structure_sha256": report["structure"]["sha256"],
                **feats,
            })
            if progress:
                progress(f"ok     {label}")
        except Exception as exc:  # noqa: BLE001 - recorded per row
            bad.append({"row": label, "error": f"{type(exc).__name__}: {exc}"})
            if progress:
                progress(f"FAILED {label}: {exc}")
    return good, bad


FEATURE_TABLE_COLUMNS = ["pdb_id", "chain", "mutation", "group", "ddg", "structure_sha256", *ALL_FEATURES]


def write_csv(rows: Sequence[dict[str, Any]], path: str | Path, columns: Sequence[str]) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if row.get(k) is None else row.get(k)) for k in columns})
    return out


def read_feature_table(path: str | Path) -> list[dict[str, Any]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return [dict(r) for r in csv.DictReader(handle)]


def _num(value: Any) -> float | None:
    if value in ("", None):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------- training
def ablation_sets(available_groups: Sequence[str]) -> dict[str, list[str]]:
    sets: dict[str, list[str]] = {}
    order = [g for g in ("substitution", "structure", "language_model") if g in available_groups]
    for i in range(1, len(order) + 1):
        sets["+".join(order[:i])] = order[:i]
    if "language_model" in order and len(order) > 1:
        sets["language_model_only"] = ["language_model"]
    return sets


def train_from_rows(
    table: list[dict[str, Any]], groups_to_use: Sequence[str], *, k: int = 5, seed: int = 0, alpha: float = 1.0,
    n_boot: int = 1000, target_convention: str, dataset_sha256: str = "",
) -> tuple[ml.RidgeModel, dict[str, Any]]:
    """Fit + evaluate. Rows missing the target or any selected feature are dropped *and counted*."""
    names = [n for g in groups_to_use for n in FEATURE_GROUPS[g]]
    kept, X_rows, y, grp = [], [], [], []
    dropped_no_target = dropped_missing_feature = 0
    for row in table:
        target = _num(row.get("ddg"))
        if target is None:
            dropped_no_target += 1
            continue
        values = [_num(row.get(n)) for n in names]
        if any(v is None for v in values):
            dropped_missing_feature += 1
            continue
        kept.append(row); X_rows.append(values); y.append(target); grp.append(str(row.get("group") or row.get("pdb_id")))
    if len(y) < 10:
        raise ValueError(f"Only {len(y)} complete rows; at least 10 are required (hundreds recommended).")
    X, yv = np.array(X_rows, dtype=float), np.array(y, dtype=float)
    cv = ml.cross_validate(X, yv, grp, names, k=k, seed=seed, alpha=alpha)
    oof = cv.pop("out_of_fold_predictions")
    ci = {m: ml.cluster_bootstrap_ci(yv, oof, grp, m, n_boot=n_boot, seed=seed) for m in ("spearman", "pearson")}
    model = ml.RidgeModel.fit(X, yv, names, alpha=alpha)
    card = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "feature_groups": list(groups_to_use),
        "n_rows": len(yv), "n_proteins": len(set(grp)),
        "dropped_rows": {"no_target": dropped_no_target, "missing_feature": dropped_missing_feature},
        "target_sign_convention": target_convention,
        "training_data_sha256": dataset_sha256,
        "cross_validation": cv,
        "cluster_bootstrap_95ci": {m: list(v) for m, v in ci.items()},
        "standardised_coefficients": model.coefficients(),
        "intended_use": "Research ranking of substitutions within the training distribution. Not for clinical use.",
        "limitations": [
            "Trained only on the supplied dataset; performance outside its protein families is unknown.",
            "Linear model on interpretable descriptors; it will miss non-additive effects.",
            "Features use wild-type structures only; no mutant structure is modelled.",
        ],
    }
    return model, card


def run_ablation(
    table: list[dict[str, Any]], available_groups: Sequence[str], **kw: Any
) -> dict[str, Any]:
    results = {}
    for label, groups in ablation_sets(available_groups).items():
        try:
            _, card = train_from_rows(table, groups, **kw)
            results[label] = {
                "n_rows": card["n_rows"], "cv_model": card["cross_validation"]["model"],
                "cv_mean_baseline": card["cross_validation"]["mean_baseline"],
                "spearman_95ci": card["cluster_bootstrap_95ci"]["spearman"],
            }
        except ValueError as exc:
            results[label] = {"error": str(exc)}
    return results


def ablation_markdown(results: dict[str, Any]) -> str:
    lines = ["| Feature set | n | Spearman (95% CI) | Pearson | RMSE | Mean-baseline RMSE | AUROC (destab.) |", "|---|---|---|---|---|---|---|"]
    for label, r in results.items():
        if "error" in r:
            lines.append(f"| {label} | – | {r['error']} | | | | |")
            continue
        m, b, ci = r["cv_model"], r["cv_mean_baseline"], r["spearman_95ci"]
        lines.append(f"| {label} | {r['n_rows']} | {m['spearman']} ({ci[0]}, {ci[1]}) | {m['pearson']} | {m['rmse']} | {b['rmse']} | {m['auroc_destabilising']} |")
    return "\n".join(lines)


# -------------------------------------------------------------------- scan
def parse_positions(text: str) -> list[int]:
    out: list[int] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        m = re.fullmatch(r"(\d+)(?:-(\d+))?", part)
        if not m:
            raise ValueError(f"Bad position specification {part!r}; use e.g. 50 or 45-60,70.")
        lo, hi = int(m.group(1)), int(m.group(2) or m.group(1))
        if hi < lo or lo < 1:
            raise ValueError(f"Bad position range {part!r}.")
        out.extend(range(lo, hi + 1))
    if not out:
        raise ValueError("No positions given.")
    return sorted(set(out))


def saturation_scan(
    base: RunConfig, positions: Sequence[int], wild_type_at: Callable[[int], str],
    run_fn: RunFn | None = None, extra_context: dict[str, Any] | None = None,
    progress: Callable[[str], None] | None = None,
) -> list[dict[str, Any]]:
    """All 19 substitutions at each position, offline by default; one row per variant."""
    run = run_fn or _default_run()
    rows: list[dict[str, Any]] = []
    for pos in positions:
        wt = wild_type_at(pos)
        for mt in sorted(AA - {wt}):
            cfg = dataclasses.replace(base, mutation=f"{wt}{pos}{mt}")
            try:
                result = run(cfg, extra_context=extra_context) if extra_context else run(cfg)
                report = result.to_dict()
                feats = extract_features(report)
                ml_pred = (report.get("ml_prediction") or {}).get("ddg_kcal_per_mol")
                rows.append({"position": pos, "wild_type": wt, "mutant": mt, "mutation": f"{wt}{pos}{mt}",
                             "ml_ddg": ml_pred, "error": "", **feats})
            except Exception as exc:  # noqa: BLE001
                rows.append({"position": pos, "wild_type": wt, "mutant": mt, "mutation": f"{wt}{pos}{mt}",
                             "ml_ddg": None, "error": f"{type(exc).__name__}: {exc}"})
            if progress:
                progress(f"{wt}{pos}{mt}")
    return rows


SCAN_COLUMNS = ["mutation", "position", "wild_type", "mutant", "ml_ddg", "error", *ALL_FEATURES]


def scan_heatmap(rows: Sequence[dict[str, Any]], metric: str, path: str | Path) -> Path:
    """Positions x amino acids heatmap of one column (matplotlib is an optional dependency)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    aas = "ACDEFGHIKLMNPQRSTVWY"
    positions = sorted({r["position"] for r in rows})
    grid = np.full((len(positions), len(aas)), np.nan)
    wts = {}
    for r in rows:
        v = r.get(metric)
        if v is not None and v != "":
            grid[positions.index(r["position"]), aas.index(r["mutant"])] = float(v)
        wts[r["position"]] = r["wild_type"]
    fig, ax = plt.subplots(figsize=(8, max(2.5, 0.28 * len(positions) + 1.5)))
    lim = np.nanmax(np.abs(grid)) if np.isfinite(grid).any() else 1.0
    im = ax.imshow(grid, aspect="auto", cmap="RdBu", vmin=-lim, vmax=lim)
    ax.set_xticks(range(len(aas)), list(aas)); ax.set_yticks(range(len(positions)), [f"{wts[p]}{p}" for p in positions])
    ax.set_xlabel("Mutant amino acid"); ax.set_title(f"Saturation scan: {metric} (blank = wild type / unavailable)")
    fig.colorbar(im, ax=ax, label=metric)
    out = Path(path); out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    return out


# --------------------------------------------------------------- reproduce
def config_from_report(report: dict[str, Any], structure_path: str | None = None, email: str | None = None) -> RunConfig:
    run = report.get("run")
    if not run or "config" not in run:
        raise ValueError("Report has no run manifest (schema < 2.0); it cannot be re-run automatically.")
    data = dict(run["config"])
    data.pop("email_provided", None)
    data["email"] = email
    if data.get("structure_path"):
        if not structure_path:
            raise ValueError("This report used a local structure file; pass --structure with the same file.")
        data["structure_path"] = structure_path
    data["foldx_executable"] = None  # machine-specific path; FoldX runs are never auto-repeated
    data.pop("extras", None)
    data["output_root"] = "data"
    return RunConfig.from_mapping({k: v for k, v in data.items() if v is not None or k in {"email"}})


def reproduce(report: dict[str, Any], structure_path: str | None = None, run_fn: RunFn | None = None,
              email: str | None = None) -> dict[str, Any]:
    cfg = config_from_report(report, structure_path, email)
    if cfg.stability_backend != "none" or cfg.conservation_backend != "none":
        cfg.stability_backend = "none"; cfg.conservation_backend = "none"
        remote_skipped = True
    else:
        remote_skipped = False
    fresh = (run_fn or _default_run())(cfg).to_dict()
    # The config hash legitimately differs when remote backends are skipped; compare the science only.
    a, b = json.loads(json.dumps(report)), json.loads(json.dumps(fresh))
    if remote_skipped:
        a["run"]["config_sha256"] = b["run"]["config_sha256"] = ""
    differences = compare_reports(a, b)
    return {"reproduced": not differences, "differences": differences,
            "remote_backends_skipped": remote_skipped, "fresh_report": fresh}
