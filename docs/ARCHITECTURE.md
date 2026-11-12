# Architecture

## Layers

```
cli / pdbmi / app.py            presentation (no science here)
        │
pipeline.py  run_pipeline()     validate RunConfig → execute DAG → attach run manifest
        │
dag.py + stages.py              the DAG: stage registry = executable definition = diagram
        │
 ┌──────┴───────────────────────────────────────────────────────────────┐
 fetch_structure  local_structure  descriptors  contacts  esm_backend    measurement modules
 stability        conservation     tables       features  ml             (pure functions, one job each)
 └──────────────────────────────────────────────────────────────────────┘
research.py                     batch features, saturation scan, training, ablation, reproduce
provenance.py  config.py        manifest, hashes, environment, validated run configuration
```

## The DAG contract (`dag.py`)

- A `Stage` has a name, a function `f(ctx) -> output`, hard dependencies (`deps`), soft inputs (`after`), and flags (`optional`, `remote`, `enabled`).
- Stages run in deterministic topological order. `ctx` holds the `RunConfig` plus every upstream output.
- **Required** stage fails → the run aborts with the original exception type.
- **Optional** stage fails → status `failed`, its dependents are `skipped`, everything independent continues, and the report field stays `null`.
- `enabled` false → status `disabled` (e.g. no stability backend chosen).
- `after` inputs order the schedule but never block (the ML stage works with or without the language-model score).
- `Graph.to_mermaid()` renders the registry; `tests/test_cli_research.py::test_dag_doc_is_in_sync_with_code` fails if `docs/pipeline-dag.md` is stale.

Adding a measurement = write one function, add one `Stage(...)` line in `stages.py`, add a field in the assembler, add tests, run `make dag`.

## Reproducibility design

| Concern | Mechanism |
|---|---|
| Same inputs? | `structure.sha256` of the exact coordinate file; `cache_policy` recorded. |
| Same settings? | `run.config_sha256` over the canonical public config (no email, no local paths). |
| Same code/environment? | `dependency_versions`, `run.environment.git_commit`, `container_image_digest` (set `IMAGE_DIGEST`), `constraints.txt` via `make lock`. |
| Same result? | `pdbmi reproduce report.json` re-runs offline stages and diffs deterministic fields (SASA, counts, descriptors, hashes). |
| Randomness | None in the offline pipeline. Cross-validation and bootstrap take an explicit `--seed`. |
| Remote services | Job IDs and raw responses are stored; remote results are **not** re-run by `reproduce` (they can change over time). |

A container standardises the environment but does not by itself guarantee bit-identical floats across CPU architectures. Archive the image digest with any published result.

## Containers and HPC

One Dockerfile, build arg `EXTRAS`: `app` (default) or `app,ml,plots` (CPU PyTorch). Runs as UID 10001; writable state lives in the mounted `data/` volume (structures, results, Hugging Face cache).

```bash
docker build -t pdbmi:2.0.0 .
docker save pdbmi:2.0.0 -o pdbmi.tar
apptainer build pdbmi.sif docker-archive:pdbmi.tar        # on a machine with Apptainer
apptainer exec --bind "$PWD/data:/opt/pdb-mutant-impact-visualizer/data" pdbmi.sif pdbmi run --pdb 1AKI --chain A --mutation S50V
```
Check your cluster's container, network and data policies first.

## Report schema 2.0

Top-level keys: `schema_version, software, dependency_versions, analyzed_at_utc, structure, chain, mutation, stability, conservation, local_structure (with solvent_accessibility, contacts), substitution, language_model, ml_prediction, method_references, warnings, run`.
`run` = `{config_sha256, config, stages[{name,status,seconds,error,remote,optional,depends_on}], execution_order, environment, reproduce_command, seed}`.
`pdbmi validate` checks structure; the `email` is never stored.

## Security notes

- Web app: FoldX path only from `FOLDX_EXE`; uploads capped (`.streamlit/config.toml`, 50 MB); temp files deleted.
- Remote submission (DynaMut2, EBI) is opt-in and described in `SECURITY.md`.
- Do not expose the app publicly without authentication.
