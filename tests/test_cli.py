import json

from click.testing import CliRunner

from cli import main


def test_cli_runs_local_pipeline_and_writes_complete_json(fixture_path, tmp_path):
    output = tmp_path / "result.json"
    result = CliRunner().invoke(main, [
        "--structure", str(fixture_path), "--chain", "A", "--mutation", "A50V",
        "--stability", "none", "--conservation", "none", "--output", str(output),
    ])
    assert result.exit_code == 0, result.output
    assert "Structure SHA-256:" in result.output
    payload = json.loads(output.read_text())
    assert payload["mutation"] == "A50V"
    assert payload["dependency_versions"]["python"]
    assert payload["dependency_versions"]["biopython"]
    assert payload["stability"] is None
    assert payload["conservation"] is None
    assert payload["software"]["version"]


def test_cli_rejects_mutation_format_without_traceback(fixture_path):
    result = CliRunner().invoke(main, ["--structure", str(fixture_path), "--chain", "A", "--mutation", "bad"])
    assert result.exit_code != 0
    assert "Traceback" not in result.output
