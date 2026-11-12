# Contributing

Thanks for helping improve the workflow.

- Keep remote-service calls out of the default test suite. Mock HTTP responses and job polling so tests are deterministic and do not burden public services.
- Do not add static biological predictions or silently replace failed analyses with defaults.
- Include source, parameter, coordinate hash, and method details for new measurements.
- Add tests for valid input, invalid input, service failures, and missing results.
- Document scientific assumptions, limitations, and relevant citations in `docs/scientific-methods.md`.
- Before opening a pull request, run `make verify` (tests, mypy, ruff, DAG-doc check, smoke test). If you change a stage, run `make dag`.
- New measurements must be a `Stage` in `stages.py` with tests; new ML claims need grouped-CV evidence.

Do not include user structures, remote-service credentials, FoldX binaries, or licensed FoldX data in contributions.
