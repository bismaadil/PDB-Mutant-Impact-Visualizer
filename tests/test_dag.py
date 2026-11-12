import pytest

from pdb_mutant_impact_visualizer.dag import DAGError, Graph, Stage, StageFailure


def _ok(value):
    return lambda ctx: value


def _boom(ctx):
    raise RuntimeError("service down")


def test_topological_order_respects_dependencies():
    g = Graph([
        Stage("c", _ok(3), deps=("a", "b")),
        Stage("a", _ok(1)),
        Stage("b", _ok(2), deps=("a",)),
    ])
    assert g.order == ["a", "b", "c"]
    assert g.levels() == [["a"], ["b"], ["c"]]


def test_independent_stages_share_a_level():
    g = Graph([Stage("root", _ok(0)), Stage("x", _ok(1), deps=("root",)), Stage("y", _ok(2), deps=("root",))])
    assert g.levels() == [["root"], ["x", "y"]]


def test_cycle_and_unknown_dependency_are_rejected():
    with pytest.raises(DAGError, match="Cycle"):
        Graph([Stage("a", _ok(1), deps=("b",)), Stage("b", _ok(2), deps=("a",))])
    with pytest.raises(DAGError, match="unknown"):
        Graph([Stage("a", _ok(1), deps=("ghost",))])
    with pytest.raises(DAGError, match="Duplicate"):
        Graph([Stage("a", _ok(1)), Stage("a", _ok(2))])


def test_outputs_flow_downstream():
    g = Graph([Stage("a", _ok(2)), Stage("b", lambda ctx: ctx["a"] * 10, deps=("a",))])
    out = g.run({})
    assert out.outputs == {"a": 2, "b": 20}
    assert all(r.status == "ok" for r in out.records)


def test_optional_failure_is_isolated_and_dependents_skipped():
    g = Graph([
        Stage("base", _ok(1)),
        Stage("remote", _boom, deps=("base",), optional=True, remote=True),
        Stage("after_remote", _ok(9), deps=("remote",), optional=True),
        Stage("local", _ok(5), deps=("base",)),
    ])
    out = g.run({})
    assert out.record("remote").status == "failed"
    assert "service down" in out.record("remote").error
    assert out.record("after_remote").status == "skipped"
    assert out.record("local").status == "ok"
    assert "remote" not in out.outputs  # never fabricate a missing result


def test_required_failure_aborts_with_partial_records():
    g = Graph([Stage("a", _ok(1)), Stage("b", _boom, deps=("a",)), Stage("c", _ok(3), deps=("b",))])
    with pytest.raises(StageFailure, match="Stage 'b' failed") as info:
        g.run({})
    assert [r.name for r in info.value.records] == ["a", "b"]


def test_disabled_stage_and_its_dependents_are_not_run():
    ran = []
    g = Graph([
        Stage("a", _ok(1)),
        Stage("opt", lambda c: ran.append("opt"), deps=("a",), optional=True, enabled=lambda c: c["flag"]),
        Stage("down", lambda c: ran.append("down"), deps=("opt",), optional=True),
    ])
    out = g.run({"flag": False})
    assert ran == []
    assert out.record("opt").status == "disabled"
    assert out.record("down").status == "disabled"


def test_mermaid_is_generated_from_the_registry():
    g = Graph([Stage("a", _ok(1), kind="input", label="Input"), Stage("b", _ok(2), deps=("a",), optional=True, remote=True, label="Remote")])
    text = g.to_mermaid()
    assert text.startswith("flowchart TD")
    assert "a -->" not in text and "a -.-> b" in text
    assert "class b remote;" in text


def test_soft_inputs_order_but_never_block():
    seen = {}
    g = Graph([
        Stage("base", _ok(1)),
        Stage("esm", _boom, deps=("base",), optional=True),
        Stage("off", _ok(0), deps=("base",), optional=True, enabled=lambda c: False),
        Stage("model", lambda c: seen.update(esm=c.get("esm"), off=c.get("off")) or "ran",
              deps=("base",), after=("esm", "off"), optional=True),
    ])
    assert g.order.index("esm") < g.order.index("model")
    out = g.run({})
    assert out.record("esm").status == "failed" and out.record("off").status == "disabled"
    assert out.record("model").status == "ok" and seen == {"esm": None, "off": None}
