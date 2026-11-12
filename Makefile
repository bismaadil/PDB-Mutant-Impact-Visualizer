.PHONY: help install test typecheck lint dag dag-check smoke live verify lock docker-build docker-run docker-smoke clean
PY ?= python
help:            ## list targets
	@grep -E '^[a-z-]+:.*##' $(MAKEFILE_LIST) | sed 's/:.*##/\t/'
install:         ## editable install with app + dev + plots extras
	$(PY) -m pip install -e ".[app,dev,plots]"
test:            ## offline unit tests (no network)
	$(PY) -m pytest
typecheck:       ## mypy
	$(PY) -m mypy
lint:            ## ruff
	$(PY) -m ruff check .
dag:             ## regenerate docs/pipeline-dag.md from the code
	pdbmi dag --write docs/pipeline-dag.md
dag-check:       ## fail if the committed DAG diagram is stale
	$(PY) -m pytest tests/test_cli_research.py -k dag_doc
smoke:           ## offline end-to-end run on the bundled fixture + reproduce + validate
	bash scripts/smoke_test.sh
live:            ## network tests against RCSB (1AKI)
	$(PY) -m pytest -m live
verify: test typecheck lint dag-check smoke   ## everything a reviewer should be able to run
lock:            ## write constraints.txt (fully pinned versions) for reproducible Docker builds
	$(PY) -m pip install uv && uv pip compile pyproject.toml --extra app --extra plots --extra ml --python-version 3.12 --output-file constraints.txt
docker-build:    ## build app image
	docker build --build-arg GIT_COMMIT=$$(git rev-parse HEAD) -t pdbmi:2.0.0 .
docker-run:      ## run the web app on :8501
	mkdir -p data && docker run --rm -p 8501:8501 --user "$$(id -u):$$(id -g)" -v "$$(pwd)/data:/opt/pdb-mutant-impact-visualizer/data" pdbmi:2.0.0
docker-smoke:    ## run the CLI inside the container on the bundled fixture
	mkdir -p data && docker run --rm --user "$$(id -u):$$(id -g)" -v "$$(pwd)/tests/fixtures:/fx:ro" -v "$$(pwd)/data:/opt/pdb-mutant-impact-visualizer/data" pdbmi:2.0.0 \
	  pdbmi run --structure /fx/tiny_test.pdb --chain A --mutation A50V --output data/results/docker-smoke.json
clean:
	rm -rf .pytest_cache .mypy_cache .ruff_cache build dist *.egg-info data/results
