from pdb_mutant_impact_visualizer.config import RunConfig
from pdb_mutant_impact_visualizer.provenance import compare_reports, dependency_versions, reproduce_command


def test_reproduce_command_has_no_secrets():
    cfg = RunConfig(chain="A", mutation="S50V", pdb_id="1AKI", conservation_backend="ebi", email="me@uni.edu").validate()
    cmd = reproduce_command(cfg.public_dict())
    assert "me@uni.edu" not in cmd and "--pdb 1AKI" in cmd and "--conservation ebi" in cmd


def test_compare_reports_flags_deterministic_differences():
    a = {"structure": {"sha256": "x"}, "chain": "A", "mutation": "S50V", "run": {"config_sha256": "c"}}
    assert compare_reports(a, dict(a)) == []
    b = {**a, "structure": {"sha256": "y"}}
    assert any("structure.sha256" in p for p in compare_reports(a, b))


def test_dependency_versions_include_python():
    assert dependency_versions()["python"]
