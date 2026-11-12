
"""Checks the DAG wiring end to end with fake parser/service modules.

These tests need neither Biopython nor network access: they prove that stages
are ordered correctly, outputs reach the report, failures are isolated and the
manifest is attached. Real-structure behaviour is covered by the other tests.
"""

import contextlib
import json
import sys
import types

import pytest

import pdb_mutant_impact_visualizer.conservation as conservation_module
from pdb_mutant_impact_visualizer.config import RunConfig
from pdb_mutant_impact_visualizer.contacts import AtomTable, contact_profile
from pdb_mutant_impact_visualizer.models import Mutation
from pdb_mutant_impact_visualizer.pipeline import run_pipeline
from pdb_mutant_impact_visualizer.stages import build_graph

PKG = "pdb_mutant_impact_visualizer"


@contextlib.contextmanager
def fake_modules(structure_file, stability_error=None):
    table = AtomTable.from_rows([
        ("A", 50, "", "ALA", "CA", "C", (0.0, 0.0, 0.0)),
        ("A", 50, "", "ALA", "CB", "C", (1.0, 0.0, 0.0)),
        ("A", 51, "", "VAL", "CA", "C", (4.0, 0.0, 0.0)),
    ])

    def mod(name, **attrs):
        m = types.ModuleType(f"{PKG}.{name}")
        m.__dict__.update(attrs)
        return m

    class Residue:
        id = (" ", 50, " ")
        resname = "ALA"

    fakes = {
        "fetch_structure": mod(
            "fetch_structure",
            download_pdb=lambda *a, **k: structure_file,
            load_structure=lambda p, i: object(),
            file_sha256=lambda p: "f" * 64,
        ),
        "local_structure": mod(
            "local_structure",
            THREE_TO_ONE={"ALA": "A"},
            validate_mutation=lambda s, c, m, i: Residue(),
            analyze_local_structure=lambda s, c, p, cut, i: {
                "neighbor_cutoff_angstrom": cut,
                "neighbor_residue_count": 3,
                "site_mean_b_factor_angstrom2": 12.4,
                "protein_mean_b_factor_angstrom2": 20.0,
            },
        ),
        "descriptors": mod(
            "descriptors",
            REFERENCES=["https://doi.org/x"],
            solvent_accessibility=lambda s, c, p, i: {
                "sasa_angstrom2": 10.0,
                "relative_sasa": 0.1,
            },
            substitution_descriptors=lambda m: {"blosum62_score": 0},
            structure_warnings=lambda *a: ["fake warning"],
        ),
        "contacts": mod(
            "contacts",
            atoms_from_structure=lambda s: table,
            contact_profile=contact_profile,
        ),
    }

    def predict_dynamut2(*a, **k):
        if stability_error:
            raise RuntimeError(stability_error)
        return {"backend": "DynaMut2", "ddg_kcal_per_mol": -0.5}

    fakes["stability"] = mod(
        "stability",
        predict_dynamut2=predict_dynamut2,
        predict_foldx=None,
    )

    saved = {
        f"{PKG}.{k}": sys.modules.get(f"{PKG}.{k}")
        for k in fakes
    }

    try:
        for k, v in fakes.items():
            sys.modules[f"{PKG}.{k}"] = v
        yield
    finally:
        for name, original in saved.items():
            if original is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = original


def make_config(tmp_path, **kw):
    f = tmp_path / "s.pdb"
    f.write_text("ATOM\n")
    base = dict(
        chain="A",
        mutation="A50V",
        structure_path=str(f),
        output_root=str(tmp_path),
    )
    base.update(kw)
    return RunConfig(**base), f


def test_graph_shape_matches_design():
    g = build_graph()
    assert g.order[0] == "fetch_structure" and g.order[-1] == "assemble"
    assert g.order.index("language_model") < g.order.index("ml_predict")
    assert "-.-> assemble" in g.to_mermaid()
    first_levels = g.levels()
    assert first_levels[0] == ["fetch_structure"]
    assert {
        "local_geometry",
        "solvent_accessibility",
        "substitution",
        "contacts",
    } <= set(first_levels[3])


def test_offline_run_populates_report_and_manifest(tmp_path):
    cfg, f = make_config(tmp_path)

    with fake_modules(f):
        result = run_pipeline(cfg)

    data = result.to_dict()

    assert data["schema_version"] == "2.0"
    assert data["stability"] is None
    assert data["conservation"] is None
    assert data["language_model"] is None
    assert (
        data["local_structure"]["contacts"]["side_chain_contact_residue_count"]
        == 1
    )
    assert "fake warning" in data["warnings"]

    statuses = {s["name"]: s["status"] for s in data["run"]["stages"]}
    assert statuses["stability"] == "disabled"
    assert statuses["contacts"] == "ok"
    assert statuses["assemble"] == "ok"

    assert data["run"]["config_sha256"] == cfg.sha256()
    assert len(data["run"]["config_sha256"]) == 64
    assert "pdbmi run" in data["run"]["reproduce_command"]

    json.dumps(data)  # fully serialisable


def test_remote_failure_is_isolated_and_reported_as_missing(tmp_path):
    cfg, f = make_config(tmp_path, stability_backend="dynamut2")

    with fake_modules(f, stability_error="service down"):
        result = run_pipeline(cfg)

    assert result.stability is None
    assert any(
        "Stability prediction FAILED" in w and "service down" in w
        for w in result.warnings
    )
    assert (
        result.local_structure["solvent_accessibility"]["relative_sasa"]
        == 0.1
    )  # rest unaffected

    with fake_modules(f, stability_error="service down"):
        with pytest.raises(RuntimeError, match="service down"):
            run_pipeline(cfg, strict=True)


def test_successful_remote_stage_and_language_model_injection(tmp_path):
    cfg, f = make_config(
        tmp_path,
        stability_backend="dynamut2",
        esm_model="fake/esm",
    )

    # structure[0][chain] is needed by the ESM stage: give the fake structure
    # that shape.
    class R:
        def __init__(s, n, name):
            s.id = (" ", n, " ")
            s.resname = name

    class S:
        def __getitem__(s, i):
            return {"A": [R(49, "ALA"), R(50, "ALA"), R(51, "ALA")]}

    def fake_esm(name):
        return lambda masked, idx: {
            a: (-0.1 if a == "A" else -2.0)
            for a in "ACDEFGHIKLMNPQRSTVWY"
        }

    with fake_modules(f):
        sys.modules[f"{PKG}.fetch_structure"].load_structure = (
            lambda p, i: S()
        )
        result = run_pipeline(
            cfg,
            extra_context={"esm_factory": fake_esm},
        )

    assert result.stability["ddg_kcal_per_mol"] == -0.5
    assert result.language_model["esm_llr"] == pytest.approx(-1.9)


def test_bad_inputs_fail_before_any_work(tmp_path):
    cfg, f = make_config(tmp_path, mutation="A50A")

    with pytest.raises(ValueError, match="differ"):
        run_pipeline(cfg)

    with pytest.raises(ValueError, match="Mutation must"):
        Mutation.parse("junk")


def test_conservation_pipeline_preserves_submitted_ebi_job_ids(
    tmp_path,
    monkeypatch,
):
    """Exercise real conservation logic and preserve both IDs on poll failure."""
    cfg, f = make_config(
        tmp_path,
        conservation_backend="ebi",
        email="researcher@example.org",
        uniprot_accession="P12345",
        uniprot_position=50,
    )

    # A valid-format UniProt sequence with alanine at position 50.
    monkeypatch.setattr(
        conservation_module,
        "fetch_uniprot_sequence",
        lambda accession, session: (
            f"{accession} reviewed protein",
            "A" * 100,
        ),
    )

    def fake_submit(tool, email, session, **params):
        if tool == "ncbiblast":
            return "blast-test-id"
        if tool == "clustalo":
            return "clustalo-test-id"
        raise AssertionError(f"Unexpected EBI tool: {tool}")

    monkeypatch.setattr(conservation_module, "_submit", fake_submit)

    def fake_wait(tool, job_id, timeout_s, interval_s, session):
        if tool == "ncbiblast":
            assert job_id == "blast-test-id"
            return

        if tool == "clustalo":
            assert job_id == "clustalo-test-id"
            raise conservation_module.ConservationError(
                "Simulated Clustal Omega polling failure."
            )

        raise AssertionError(f"Unexpected EBI polling tool: {tool}")

    monkeypatch.setattr(conservation_module, "_wait", fake_wait)

    monkeypatch.setattr(
        conservation_module,
        "_result",
        lambda tool, job_id, kind, session: (
            "sp:P23456\n"
            if tool == "ncbiblast"
            else pytest.fail("Clustal result must not be requested after failure")
        ),
    )

    monkeypatch.setattr(
        conservation_module,
        "_fetch_homologs",
        lambda accessions, session: (
            ">sp|P12345|QUERY\n"
            + "A" * 100
            + "\n>sp|P23456|HOMOLOG\n"
            + "A" * 100
            + "\n"
        ),
    )

    # The pipeline, DAG, stage wrapper, and actual analyze_conservation()
    # function execute. Only external-service interactions are simulated.
    with fake_modules(f):
        result = run_pipeline(cfg)

    data = result.to_dict()

    assert result.conservation is None

    conservation_record = next(
        stage
        for stage in data["run"]["stages"]
        if stage["name"] == "conservation"
    )

    assert conservation_record["status"] == "failed"
    assert conservation_record["metadata"] == {
        "blast_job_id": "blast-test-id",
        "clustalo_job_id": "clustalo-test-id",
    }

    assert "blast-test-id" in conservation_record["error"]
    assert "clustalo-test-id" in conservation_record["error"]
    assert "Simulated Clustal Omega polling failure" in (
        conservation_record["error"]
    )

    # A failed optional conservation stage must not prevent assembly.
    assert data["local_structure"] is not None

    statuses = {
        stage["name"]: stage["status"]
        for stage in data["run"]["stages"]
    }
    assert statuses["assemble"] == "ok"
    assert statuses["contacts"] == "ok"

    # The final report, including failure metadata, must remain serialisable.
    json.dumps(data)
