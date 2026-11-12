import json

from pdb_mutant_impact_visualizer.local_structure import analyze_local_structure
from pdb_mutant_impact_visualizer.models import Mutation
from pdb_mutant_impact_visualizer.pipeline import analyze_mutation
from pdb_mutant_impact_visualizer.report import render_text, write_json


def test_local_metrics_are_raw_coordinate_observations(structure):
    result = analyze_local_structure(structure, "A", 50, cutoff=10)
    assert result["neighbor_residue_count"] == 3
    assert result["site_mean_b_factor_angstrom2"] == 12.4
    assert "not a direct measure" in result["b_factor_note"]


def test_pipeline_keeps_unrequested_backends_missing(fixture_path, tmp_path):
    result = analyze_mutation(
        mutation=Mutation.parse("A50V"), chain="A", structure_path=fixture_path,
        stability_backend="none", conservation_backend="none", output_root=tmp_path,
    )
    assert result.stability is None
    assert result.conservation is None
    assert result.structure["sha256"]
    assert "no composite impact score" in render_text(result).lower()
    output = write_json(result, tmp_path / "analysis.json")
    saved = json.loads(output.read_text())
    assert saved["stability"] is None
    assert saved["conservation"] is None


def test_pipeline_requires_exactly_one_structure_source(fixture_path):
    import pytest
    with pytest.raises(ValueError, match="exactly one"):
        analyze_mutation(mutation=Mutation.parse("A50V"), chain="A", structure_path=fixture_path, pdb_id="1ABC")
