"""Network tests. Run with: pytest -m live  (skipped by default)."""

import pytest

from pdb_mutant_impact_visualizer.pipeline import analyze_mutation
from pdb_mutant_impact_visualizer.models import Mutation

pytestmark = pytest.mark.live


def test_real_structure_1aki_offline_measurements_are_sane(tmp_path):
    result = analyze_mutation(mutation=Mutation.parse("S50V"), chain="A", pdb_id="1AKI", output_root=tmp_path, strict=True)
    data = result.to_dict()
    assert len(data["structure"]["sha256"]) == 64
    assert 0.0 <= data["local_structure"]["solvent_accessibility"]["relative_sasa"] <= 1.0
    assert data["local_structure"]["neighbor_residue_count"] > 5
    assert data["local_structure"]["contacts"]["side_chain_contact_residue_count"] >= 1
    assert all(s["status"] in {"ok", "disabled"} for s in data["run"]["stages"])
