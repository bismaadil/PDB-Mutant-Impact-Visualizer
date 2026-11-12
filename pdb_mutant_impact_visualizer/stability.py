"""External stability backends with explicit errors and raw provenance."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

import requests

from .models import Mutation

DYNAMUT2_URL = "https://biosig.unimelb.edu.au/dynamut2/api/prediction_single"


class StabilityError(RuntimeError):
    pass


def predict_dynamut2(
    pdb_path: str | Path, chain: str, mutation: Mutation, timeout_s: int = 900,
    poll_interval_s: int = 10, session: requests.Session | None = None,
) -> dict[str, Any]:
    if timeout_s <= 0 or poll_interval_s <= 0:
        raise ValueError("Polling timeout and interval must be positive.")
    client = session or requests.Session()
    try:
        with Path(pdb_path).open("rb") as handle:
            response = client.post(
                DYNAMUT2_URL,
                data={"mutation": str(mutation), "chain": chain},
                files={"pdb_file": (Path(pdb_path).name, handle, "chemical/x-pdb")},
                timeout=(15, 120),
            )
        response.raise_for_status()
        job_id = response.json().get("job_id")
    except (OSError, requests.RequestException, ValueError) as exc:
        raise StabilityError(f"Could not submit the DynaMut2 job: {exc}") from exc
    if not job_id:
        raise StabilityError("DynaMut2 did not return a job ID.")

    elapsed = 0
    while elapsed <= timeout_s:
        try:
            response = client.get(DYNAMUT2_URL, params={"job_id": job_id}, timeout=(15, 60))
            response.raise_for_status()
            payload = response.json()
        except (requests.RequestException, ValueError) as exc:
            raise StabilityError(f"Could not retrieve DynaMut2 job {job_id}: {exc}") from exc
        prediction = payload.get("prediction")
        if prediction is not None:
            try:
                value = float(prediction)
            except (TypeError, ValueError) as exc:
                raise StabilityError(f"DynaMut2 returned a non-numeric prediction: {prediction!r}") from exc
            return {
                "backend": "DynaMut2",
                "job_id": str(job_id),
                "ddg_kcal_per_mol": value,
                "sign_convention": "DynaMut2: positive ddG = stabilising, negative = destabilising (opposite of FoldX). Service prediction; see the method paper.",
                "results_url": payload.get("results_page"),
                "raw_response": payload,
            }
        state = str(payload.get("status", payload.get("message", ""))).upper()
        if state not in {"RUNNING", "QUEUED", "PENDING", "PROCESSING"}:
            raise StabilityError(f"DynaMut2 job {job_id} did not return a prediction: {payload}")
        if elapsed == timeout_s:
            break
        import time

        time.sleep(min(poll_interval_s, timeout_s - elapsed))
        elapsed += poll_interval_s
    raise StabilityError(f"DynaMut2 job {job_id} timed out after {timeout_s} seconds.")


def _write_mutation_file(work_dir: Path, chain: str, mutation: Mutation) -> Path:
    path = work_dir / "individual_list.txt"
    path.write_text(f"{mutation.wild_type}{chain}{mutation.position}{mutation.mutant};\n", encoding="ascii")
    return path


def _parse_foldx_difference(path: str | Path) -> float:
    path = Path(path)
    if not path.is_file():
        raise StabilityError(f"FoldX did not create its expected difference file: {path.name}")
    lines = [line for line in path.read_text(errors="replace").splitlines() if line.strip()]
    header_i = next((i for i, line in enumerate(lines) if line.split("\t")[0].strip() == "Pdb"), None)
    if header_i is None:
        raise StabilityError(f"Could not find a FoldX data header in {path.name}.")
    header = [cell.strip() for cell in lines[header_i].split("\t")]
    col = next((i for i, name in enumerate(header) if name.lower() == "total energy"), None)
    if col is None:
        raise StabilityError(f"No 'total energy' field in FoldX output header: {header}")
    rows = [line.split("\t") for line in lines[header_i + 1:] if line.strip()]
    if not rows or len(rows[0]) <= col:
        raise StabilityError(f"FoldX difference file has no usable energy row: {path.name}")
    try:
        return float(rows[0][col])
    except ValueError as exc:
        raise StabilityError(f"Invalid FoldX energy value {rows[0][col]!r}.") from exc


def predict_foldx(
    executable: str | Path, pdb_path: str | Path, chain: str, mutation: Mutation,
    work_root: str | Path = "data/foldx-runs", timeout_s: int = 600,
) -> dict[str, Any]:
    exe = Path(executable).expanduser().resolve()
    structure_path = Path(pdb_path).resolve()
    if not exe.is_file():
        raise StabilityError(
            f"FoldX executable was not found at {exe}. Install FoldX under its licence terms "
            "and pass its full path with --foldx-exe."
        )
    if not structure_path.is_file():
        raise StabilityError(f"PDB file does not exist: {structure_path}")
    root = Path(work_root)
    root.mkdir(parents=True, exist_ok=True)
    run_dir = Path(tempfile.mkdtemp(prefix="foldx-", dir=root)).resolve()
    pdb_name = "input_structure.pdb"
    shutil.copy2(structure_path, run_dir / pdb_name)

    def run(args: list[str]) -> None:
        try:
            completed = subprocess.run(
                [str(exe), *args], cwd=run_dir, text=True, capture_output=True,
                timeout=timeout_s, check=False,
            )
        except subprocess.TimeoutExpired as exc:
            raise StabilityError(f"FoldX exceeded the {timeout_s}s command timeout.") from exc
        if completed.returncode:
            tail = (completed.stdout + "\n" + completed.stderr)[-4000:]
            raise StabilityError(f"FoldX exited with status {completed.returncode}:\n{tail}")

    try:
        run(["--command=RepairPDB", f"--pdb={pdb_name}"])
        repaired = run_dir / "input_structure_Repair.pdb"
        if not repaired.is_file():
            raise StabilityError("FoldX RepairPDB completed without creating input_structure_Repair.pdb.")
        _write_mutation_file(run_dir, chain, mutation)
        run(["--command=BuildModel", f"--pdb={repaired.name}", "--mutant-file=individual_list.txt"])
        ddg = _parse_foldx_difference(run_dir / "Dif_input_structure_Repair.fxout")
    except StabilityError:
        raise
    except OSError as exc:
        raise StabilityError(f"Could not start FoldX: {exc}") from exc
    return {
        "backend": "FoldX",
        "ddg_kcal_per_mol": ddg,
        "output_file": str(run_dir / "Dif_input_structure_Repair.fxout"),
        "work_directory": str(run_dir),
        "sign_convention": "FoldX: positive ddG = destabilising. DynaMut2 reports the opposite sign (positive = stabilising); do not compare raw values across backends.",
        "note": "Energy difference reported by FoldX BuildModel; inspect the preserved raw output before interpretation.",
    }
