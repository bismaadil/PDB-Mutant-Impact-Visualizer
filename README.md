# PDB Mutant Impact Visualizer (v2.0)

A **reproducible, DAG-orchestrated workflow** for studying one amino-acid substitution in one real PDB structure, plus a **research layer** (saturation scans, feature tables, protein-grouped cross-validation, ablations, an optional ESM-2 language-model score).

It is an integration, measurement and evaluation framework. It is **not** a new state-of-the-art predictor, does not build a mutant structure, and deliberately produces **no single composite "impact score"** (see [limitations](docs/LIMITATIONS_AND_VALIDATION.md)).

```
PDB ID / file ─► fetch ─► parse ─► validate mutation ─┬─► neighbors, B-factors
                                                      ├─► SASA (Shrake–Rupley)
                                                      ├─► BLOSUM62 / hydropathy / volume / charge
                                                      ├─► contacts: H-bond, salt bridge, disulfide, interface
                                                      ├─► [opt] ESM-2 masked-marginal LLR
                                                      ├─► [opt, remote] DynaMut2 / FoldX
                                                      ├─► [opt, remote] SIFTS + BLAST + Clustal Ω entropy
                                                      └─► [opt] user-trained ridge model ─► AnalysisResult 2.0 + run manifest
```
The full, auto-generated diagram is in [`docs/pipeline-dag.md`](docs/pipeline-dag.md). It is produced from the same registry that executes the pipeline, and a test fails if the two drift apart.

## Why it is built this way

| Goal | How |
|---|---|
| **DAG** | `dag.py` executes `stages.py`: topological order, failure isolation (a failed remote service never corrupts offline results), per-stage status/timing in every report. |
| **Modularity** | Each measurement is one stage and one module; backends are swappable; heavy imports are lazy. |
| **Reproducibility** | Coordinate SHA-256, config SHA-256 (email excluded), dependency versions, git commit, container digest, seed, and a `pdbmi reproduce` command that re-runs a report and diffs its deterministic fields. |
| **Containers** | One Dockerfile, non-root, `app` and `research` (CPU PyTorch) targets, `compose.yaml`, Codespaces devcontainer, optional pinned `constraints.txt` (`make lock`). |
| **ML / AI** | Zero-shot ESM-2 score and a transparent ridge model with **protein-grouped CV**, **cluster-bootstrap CIs**, **ablations** and mean/single-feature baselines. No pretrained weights are bundled; you train on a dataset you cite. |
| **Honesty** | Missing = `null`, never a default. Sign conventions are recorded. Contacts are geometric *candidates*. |

## Quick start

```bash
python -m venv .venv && source .venv/bin/activate        # Windows: .venv\Scripts\Activate.ps1
python -m pip install -e ".[app,dev,plots]"
make verify            # tests, mypy, ruff, DAG-doc check, offline smoke test
pdbmi run --pdb 1AKI --chain A --mutation S50V --output data/results/1aki-S50V.json
streamlit run app.py
```

### Commands

| Command | Purpose |
|---|---|
| `pdbmi run` | One structure, one mutation (also `python cli.py …`, kept for compatibility). YAML/JSON `--config` supported. |
| `pdbmi dag [--write FILE]` | Print/generate the pipeline diagram from code. |
| `pdbmi scan` | Saturation scan (19 substitutions × positions) → CSV + heatmap. |
| `pdbmi features` | Mutation table CSV → feature table CSV (failures logged, never silently dropped). |
| `pdbmi train` | Grouped CV, bootstrap CIs, ablations → `model.json` + model card + `ablation.md`. |
| `pdbmi reproduce REPORT` | Re-run and compare deterministic fields. |
| `pdbmi validate REPORT` | Structural check of a report (schema 2.0). |

Remote services are opt-in: `--stability dynamut2|foldx`, `--conservation ebi --email you@uni.edu`, `--esm-model facebook/esm2_t12_35M_UR50D` (needs `pip install -e ".[ml]"`). A failed optional stage is recorded under `warnings` and `run.stages`; add `--strict` to make it fatal.

## Containers

```bash
make docker-build && make docker-smoke          # CLI inside the image on the bundled fixture
docker compose up app                           # web app on http://localhost:8501
docker compose --profile research run --rm research scan --pdb 1AKI --chain A --positions 45-55
```
FoldX is licensed separately and is never baked into an image; in the web app it is read only from the server-side `FOLDX_EXE` variable. For HPC (Apptainer/Singularity) see [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md#containers-and-hpc).

## Documentation

- [Architecture and reproducibility design](docs/ARCHITECTURE.md)
- [Research protocol (question, datasets, evaluation, threats to validity)](docs/RESEARCH_PROTOCOL.md)
- [Scientific methods and citations](docs/scientific-methods.md)
- [Limitations and validation status](docs/LIMITATIONS_AND_VALIDATION.md)
- [Codespaces testing and GitHub publishing guide](docs/GITHUB_CODESPACES_GUIDE.md)
- [Presenting this project to a selection committee](docs/PRESENTING_THIS_PROJECT.md)

## Cite

See [`CITATION.cff`](CITATION.cff). Please also cite the underlying methods listed in [`docs/scientific-methods.md`](docs/scientific-methods.md).

## License

MIT for this project's code. FoldX, DynaMut2, ESM-2 weights and all datasets have their own terms.
