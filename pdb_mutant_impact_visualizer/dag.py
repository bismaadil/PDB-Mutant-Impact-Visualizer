
"""A small, dependency-free DAG executor that *is* the pipeline definition.

Every analysis step is a :class:`Stage` with declared dependencies. The same
registry drives (a) execution order, (b) failure isolation, (c) per-stage
provenance (status, timing, error text and metadata) and (d) the Mermaid
diagram in the documentation, so the diagram cannot drift away from the code.

Semantics
---------
* Stages run in a deterministic topological order (ties broken by name).
* A *required* stage that fails aborts the run with :class:`StageFailure`.
* An *optional* stage that fails is recorded as ``failed``; stages that depend
  on it are recorded as ``skipped`` and the rest of the DAG keeps going.
* A stage whose ``enabled`` predicate is false is recorded as ``disabled``;
  stages that hard-depend (``deps``) on it are ``disabled`` too.
* ``after`` lists soft inputs: they are scheduled first when present, but a
  disabled/failed soft input never blocks the stage (it reads
  ``ctx.get(name)`` and copes with absence).
* Outputs are never invented: a stage that did not run has no output.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any, Callable, Iterable, Mapping

Status = str  # "ok" | "failed" | "skipped" | "disabled"


class DAGError(ValueError):
    """The graph definition itself is invalid (cycle, unknown dependency...)."""


class StageFailure(RuntimeError):
    """A required stage failed; carries the partial run record."""

    def __init__(self, stage: str, cause: BaseException, records: list["StageRecord"]):
        super().__init__(f"Stage '{stage}' failed: {cause}")
        self.stage = stage
        self.cause = cause
        self.records = records


@dataclass(frozen=True)
class Stage:
    name: str
    func: Callable[[Mapping[str, Any]], Any]
    deps: tuple[str, ...] = ()
    after: tuple[str, ...] = ()  # soft inputs: run earlier if enabled, never block this stage
    optional: bool = False
    enabled: Callable[[Mapping[str, Any]], bool] | None = None
    label: str = ""
    kind: str = "compute"  # input | compute | remote | output  (diagram styling)
    remote: bool = False   # contacts a network service / external program

    def display(self) -> str:
        return self.label or self.name


@dataclass
class StageRecord:
    name: str
    status: Status
    seconds: float = 0.0
    error: str | None = None
    remote: bool = False
    optional: bool = False
    depends_on: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "status": self.status,
            "seconds": round(self.seconds, 4),
            "error": self.error,
            "remote": self.remote,
            "optional": self.optional,
            "depends_on": list(self.depends_on),
            "metadata": dict(self.metadata),
        }


@dataclass
class RunOutcome:
    outputs: dict[str, Any]
    records: list[StageRecord]

    def record(self, name: str) -> StageRecord:
        return next(r for r in self.records if r.name == name)

    def ok(self, name: str) -> bool:
        return self.record(name).status == "ok"


class Graph:
    def __init__(self, stages: Iterable[Stage]):
        self.stages: dict[str, Stage] = {}
        for stage in stages:
            if stage.name in self.stages:
                raise DAGError(f"Duplicate stage name: {stage.name}")
            self.stages[stage.name] = stage
        for stage in self.stages.values():
            for dep in (*stage.deps, *stage.after):
                if dep not in self.stages:
                    raise DAGError(f"Stage '{stage.name}' depends on unknown stage '{dep}'")
        self._order = self._toposort()

    def _toposort(self) -> list[str]:
        indegree = {n: len(s.deps) + len(s.after) for n, s in self.stages.items()}
        children: dict[str, list[str]] = {n: [] for n in self.stages}
        for s in self.stages.values():
            for d in (*s.deps, *s.after):
                children[d].append(s.name)
        ready = sorted(n for n, k in indegree.items() if k == 0)
        order: list[str] = []
        while ready:
            node = ready.pop(0)
            order.append(node)
            for child in sorted(children[node]):
                indegree[child] -= 1
                if indegree[child] == 0:
                    ready.append(child)
            ready.sort()
        if len(order) != len(self.stages):
            cyclic = sorted(set(self.stages) - set(order))
            raise DAGError(f"Cycle detected involving: {', '.join(cyclic)}")
        return order

    @property
    def order(self) -> list[str]:
        return list(self._order)

    def levels(self) -> list[list[str]]:
        """Stages grouped by longest-path depth (stages in a level are independent)."""
        depth: dict[str, int] = {}
        for name in self._order:
            deps = (*self.stages[name].deps, *self.stages[name].after)
            depth[name] = 1 + max((depth[d] for d in deps), default=-1)
        grouped: dict[int, list[str]] = {}
        for name, d in depth.items():
            grouped.setdefault(d, []).append(name)
        return [sorted(grouped[k]) for k in sorted(grouped)]

    def run(self, context: Mapping[str, Any]) -> RunOutcome:
        """Execute the graph. ``context`` is read-only input; stage outputs are
        passed to downstream stages as ``ctx[<stage name>]`` merged over it."""
        outputs: dict[str, Any] = {}
        records: list[StageRecord] = []
        status: dict[str, Status] = {}
        for name in self._order:
            stage = self.stages[name]
            base = StageRecord(
                name=name,
                status="ok",
                remote=stage.remote,
                optional=stage.optional,
                depends_on=list(stage.deps) + list(stage.after),
            )
            ctx: dict[str, Any] = {**context, **outputs}
            if stage.enabled is not None and not stage.enabled(ctx):
                base.status = "disabled"
            elif any(status[d] in {"failed", "skipped"} for d in stage.deps):
                bad = [d for d in stage.deps if status[d] in {"failed", "skipped"}]
                base.status = "skipped"
                base.error = f"upstream stage(s) did not complete: {', '.join(bad)}"
            elif any(status[d] == "disabled" for d in stage.deps):
                base.status = "disabled"
            else:
                started = time.perf_counter()
                try:
                    outputs[name] = stage.func(ctx)
                except Exception as exc:  # noqa: BLE001 - recorded, then re-raised if required
                    base.status = "failed"
                    base.error = f"{type(exc).__name__}: {exc}"
                    base.seconds = time.perf_counter() - started

                    # Remote stages can attach structured diagnostic metadata
                    # to their exceptions. Keep it in the stage record rather
                    # than relying on parsing the human-readable error string.
                    exception_metadata = getattr(exc, "metadata", None)
                    if isinstance(exception_metadata, Mapping):
                        base.metadata = dict(exception_metadata)

                    records.append(base)
                    status[name] = "failed"
                    if not stage.optional:
                        raise StageFailure(name, exc, records) from exc
                    continue
                base.seconds = time.perf_counter() - started
            records.append(base)
            status[name] = base.status
        return RunOutcome(outputs=outputs, records=records)

    def to_mermaid(self) -> str:
        lines = ["flowchart TD"]
        shapes = {"input": ('[/"', '"/]'), "output": ('(["', '"])'), "compute": ('["', '"]'), "remote": ('["', '"]')}
        for name in self._order:
            s = self.stages[name]
            left, right = shapes.get(s.kind, shapes["compute"])
            suffix = "<br/>(optional, remote)" if s.remote and s.optional else "<br/>(optional)" if s.optional else ""
            lines.append(f"  {name}{left}{s.display()}{suffix}{right}")
        for name in self._order:
            stage = self.stages[name]
            for dep in stage.deps:
                lines.append(f"  {dep} {'-.->' if stage.optional else '-->'} {name}")
            for dep in stage.after:
                lines.append(f"  {dep} -.-> {name}")
        remote = [n for n in self._order if self.stages[n].remote]
        optional = [n for n in self._order if self.stages[n].optional and not self.stages[n].remote]
        lines.append("  classDef remote fill:#fde8e0,stroke:#c2522f,stroke-dasharray: 4 3;")
        lines.append("  classDef optional fill:#eef2fb,stroke:#5b7fc7,stroke-dasharray: 2 2;")
        if remote:
            lines.append(f"  class {','.join(remote)} remote;")
        if optional:
            lines.append(f"  class {','.join(optional)} optional;")
        return "\n".join(lines)
