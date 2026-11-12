
import pytest

from pdb_mutant_impact_visualizer.models import Mutation
from pdb_mutant_impact_visualizer.stability import StabilityError, _parse_foldx_difference, _write_mutation_file, predict_dynamut2, predict_foldx


def test_foldx_mutation_syntax(tmp_path):
    path = _write_mutation_file(tmp_path, "A", Mutation.parse("S50V"))
    assert path.read_text() == "SA50V;\n"


def test_foldx_difference_parser_reads_named_energy_column(tmp_path):
    path = tmp_path / "Dif_test.fxout"
    path.write_text("banner\nPdb\ttotal energy\tvdw\nmutant\t1.25\t.4\n")
    assert _parse_foldx_difference(path) == 1.25


def test_foldx_difference_parser_rejects_malformed_file(tmp_path):
    path = tmp_path / "bad.fxout"
    path.write_text("not a FoldX output\n")
    with pytest.raises(StabilityError, match="header"):
        _parse_foldx_difference(path)


def test_missing_foldx_executable_is_clear(fixture_path, tmp_path):
    with pytest.raises(StabilityError, match="not found"):
        predict_foldx(tmp_path / "missing-foldx", fixture_path, "A", Mutation.parse("A50V"))


class FakeResponse:
    def __init__(self, payload):
        self.payload = payload

    def raise_for_status(self):
        pass

    def json(self):
        return self.payload


class FakeDynaMutSession:
    def post(self, *args, **kwargs):
        assert kwargs["data"]["mutation"] == "A50V"
        return FakeResponse({"job_id": "test-job-42"})

    def get(self, *args, **kwargs):
        return FakeResponse({"prediction": "0.82", "results_page": "https://example.invalid/job/42"})


def test_dynamut_returns_service_result_and_provenance(fixture_path):
    result = predict_dynamut2(fixture_path, "A", Mutation.parse("A50V"), session=FakeDynaMutSession())
    assert result["ddg_kcal_per_mol"] == 0.82
    assert result["job_id"] == "test-job-42"
    assert result["results_url"].endswith("42")
