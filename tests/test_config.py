import json

import pytest

from pdb_mutant_impact_visualizer.config import RunConfig


def cfg(**kw):
    base = dict(chain="A", mutation="S50V", pdb_id="1aki")
    base.update(kw)
    return RunConfig(**base).validate()


def test_validation_normalises_and_rejects():
    assert cfg().pdb_id == "1AKI"
    with pytest.raises(ValueError, match="exactly one"):
        cfg(structure_path="x.pdb")
    with pytest.raises(ValueError, match="valid contact email"):
        cfg(conservation_backend="ebi")
    with pytest.raises(ValueError, match="FoldX"):
        cfg(stability_backend="foldx")
    with pytest.raises(ValueError, match="greater than zero"):
        cfg(neighbor_cutoff=0)


def test_hash_is_stable_and_ignores_email_and_paths():
    a = cfg(conservation_backend="ebi", email="a@b.org", output_root="/tmp/x")
    b = cfg(conservation_backend="ebi", email="other@c.org", output_root="/elsewhere")
    assert a.sha256() == b.sha256()
    assert "email" not in a.public_dict() and a.public_dict()["email_provided"] is True
    assert cfg(neighbor_cutoff=8.0).sha256() != cfg().sha256()


def test_load_json_with_override_and_unknown_field(tmp_path):
    path = tmp_path / "run.json"
    path.write_text(json.dumps({"chain": "A", "mutation": "S50V", "pdb_id": "1AKI"}))
    assert RunConfig.from_file(path, neighbor_cutoff=7.5).neighbor_cutoff == 7.5
    path.write_text(json.dumps({"chain": "A", "mutation": "S50V", "typo_field": 1}))
    with pytest.raises(ValueError, match="Unknown"):
        RunConfig.from_file(path)
