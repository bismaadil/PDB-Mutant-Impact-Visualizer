# Changelog

## 2.0.0
### Added
- **DAG executor** (`dag.py`) and stage registry (`stages.py`); failure isolation, per-stage status/timing; Mermaid diagram generated from code (`pdbmi dag`) with a drift test.
- **RunConfig** (YAML/JSON, validation, SHA-256 excluding secrets) and **run manifest** in every report (schema 2.0).
- **Contacts**: side-chain contacts, H-bond/salt-bridge candidates, disulfides, interface, Cβ burial.
- **ESM-2** masked-marginal LLR backend (optional `ml` extra; lazy imports; testable maths).
- **Research layer**: `pdbmi scan | features | train | reproduce | validate`; protein-grouped CV, cluster-bootstrap CIs, ablations, model cards.
- Containers: multi-target Dockerfile, `compose.yaml`, devcontainer for Codespaces, optional pinned `constraints.txt`.
- CI: tests (3.10/3.12), smoke test, Docker build, weekly live-service workflow. Makefile, scripts, examples, documentation set.
### Changed
- `analyze_mutation()` keeps its signature but now executes the DAG; failed *optional* stages are recorded as warnings and `null` fields (`--strict` restores fail-fast).
- Substitution tables moved to `tables.py` (no third-party imports).
- Web app: FoldX path read only from `FOLDX_EXE`; upload size capped; shared report renderer; DAG status table.
### Security
- Removed the ability to execute an arbitrary path typed into the web UI.

## 1.2.0
### Fixed
- Waters/ligands sharing a residue number no longer make a residue "ambiguous".
- Neighbor counts and protein B-factor mean exclude waters/ligands and use the first model only (NMR files were mixing models).
- Streamlit: dependent fields (email, UniProt, FoldX path) were hidden inside `st.form`; widgets now react immediately. Uploaded file names are recorded; temp files always cleaned up.
- EBI email is validated before long DynaMut2/FoldX jobs start.
- PDB ID is normalised once; insertion codes reach the SIFTS mapper; query accession matched exactly in the MSA.
- mypy errors in `app.py` (CI was failing).
### Changed
- CLI/app default stability backend is now `none` (remote upload is opt-in).
- Report schema 1.1; FoldX/DynaMut2 sign conventions recorded.
### Added
- SASA, BLOSUM62/hydropathy/volume/charge descriptors, structure warnings, backbone trace in 3D view, pipeline DAG, regression tests.
