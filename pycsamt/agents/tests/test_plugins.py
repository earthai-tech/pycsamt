"""Tests for :mod:`pycsamt.agents.plugins` and the MTpy example."""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

from pycsamt.agents import _workflow_plan, _workflows, orchestrator
from pycsamt.agents import router
from pycsamt.agents._workflow_plan import WorkflowPlan
from pycsamt.agents.plugins import (
    DEFAULT_PLUGIN,
    DomainPlugin,
    PlanValidator,
    RecordWriter,
    WorkflowRegistry,
    available_plugins,
    get_plugin,
    register_plugin,
)

REPO = Path(__file__).resolve().parents[3]
EXAMPLES = REPO / "examples" / "agent_plugins"
E2_REQUESTS = Path(os.environ.get(
    "PYCSAMT_E2_REQUESTS",
    REPO / "PAPERS" / "CAGEO" / "R2" / "E2" / "frozen"
    / "requests.jsonl"))

ALL_IDS = sorted(set(_workflows.WORKFLOW_DESCRIPTIONS)
                 | set(orchestrator._WORKFLOW_STEPS)
                 | set(_workflow_plan.VALID_WORKFLOWS))
SAMPLE_TEXTS = [
    "run QC on L26PLT", "estimate geoelectric strike for L26PLT",
    "denoise the survey data", "plot the phase tensor map",
    "invert line L22 with Occam 2D", "what is static shift?",
    "write code to load EDI files", "", "xyzzy",
]


def _e2_texts():
    if not E2_REQUESTS.is_file():
        return []
    return [json.loads(line)["query"] for line in
            E2_REQUESTS.read_text(encoding="utf-8").splitlines()
            if line.strip()]


@pytest.fixture()
def builtin():
    return get_plugin()


# ── built-in plugin ───────────────────────────────────────────────


def test_builtin_satisfies_protocols(builtin, tmp_path):
    assert builtin.name == DEFAULT_PLUGIN
    assert isinstance(builtin.registry, WorkflowRegistry)
    assert isinstance(builtin.validator, PlanValidator)
    rec = builtin.record_factory({"request": "qc"})
    assert isinstance(rec, RecordWriter)
    rec.emit("probe", value=1)
    rec.persist(tmp_path / "rec.json")
    data = json.loads((tmp_path / "rec.json").read_text("utf-8"))
    assert any(e["event"] == "probe" for e in data["events"])
    assert DEFAULT_PLUGIN in available_plugins()


def test_builtin_tables_are_the_legacy_objects(builtin):
    reg = builtin.registry
    assert reg.descriptions is _workflows.WORKFLOW_DESCRIPTIONS
    assert reg.valid_workflows is _workflow_plan.VALID_WORKFLOWS
    assert reg.dispatchable == frozenset(orchestrator._WORKFLOW_STEPS)


@pytest.mark.parametrize("wf", ALL_IDS)
def test_builtin_registry_matches_legacy_per_id(builtin, wf):
    reg = builtin.registry
    assert reg.normalise(wf) == _workflows.normalise_workflow(wf)
    assert reg.family(wf) == _workflows.workflow_family(wf)
    assert reg.steps(wf) is orchestrator._WORKFLOW_STEPS.get(wf)


@pytest.mark.parametrize("wf", sorted(_workflow_plan.VALID_WORKFLOWS))
def test_builtin_validator_matches_legacy_per_id(builtin, wf,
                                                 tmp_path):
    for cfg in ({"workflow": wf},
                {"workflow": wf, "data_path": str(tmp_path)}):
        plan = WorkflowPlan.from_config(cfg, request="run " + wf)
        assert (builtin.validator.errors(plan)
                == _workflow_plan.validate_workflow_plan(plan))
        assert (builtin.validator.needs_review(plan)
                is plan.requires_human_review)


def _texts():
    return SAMPLE_TEXTS + _e2_texts()


def test_builtin_matches_legacy_on_requests(builtin):
    # The 120 E2 requests are added when the (untracked) file exists.
    texts = _texts()
    reg, val = builtin.registry, builtin.validator
    for text in texts:
        assert reg.classify(text) == _workflows.classify_workflow(text)
        assert (reg.classify(text, default="qc")
                == _workflows.classify_workflow(text, default="qc"))
        for wf in ALL_IDS:
            for intent in ("workflow", "plot", "question"):
                assert (val.cross_check(text, intent, wf)
                        == router.route_disagreement(text, intent, wf))
            assert (val.route_flags(text, wf)
                    == orchestrator._route_cross_check(text, wf))


def test_e2_requests_present_when_required():
    if os.environ.get("PYCSAMT_REQUIRE_E2"):
        assert len(_e2_texts()) == 120


def test_steps_follow_patched_table(builtin, monkeypatch):
    fake = {"qc": [("first", "Dummy", None, "dummy")]}
    monkeypatch.setattr(orchestrator, "_WORKFLOW_STEPS", fake)
    assert builtin.registry.steps("qc") is fake["qc"]
    assert builtin.registry.dispatchable == frozenset({"qc"})


# ── registry of plugins ───────────────────────────────────────────


def test_register_rejects_builtin_and_duplicates(builtin):
    with pytest.raises(ValueError):
        register_plugin(builtin)
    clone = DomainPlugin("dup_test", builtin.registry,
                         builtin.validator, builtin.record_factory)
    register_plugin(clone, replace=True)
    with pytest.raises(ValueError):
        register_plugin(clone)
    with pytest.raises(TypeError):
        register_plugin("not a plugin")


def test_unknown_plugin_raises():
    with pytest.raises(KeyError):
        get_plugin("no_such_plugin_xyz")


def test_domain_plugin_checks_protocols(builtin):
    with pytest.raises(TypeError):
        DomainPlugin("bad", object(), builtin.validator,
                     builtin.record_factory)


# ── MTpy example plugin ───────────────────────────────────────────


@pytest.fixture()
def mtpy_plugin(monkeypatch):
    monkeypatch.syspath_prepend(str(EXAMPLES))
    sys.modules.pop("mtpy_example", None)
    sys.modules.pop("mtpy_example.plugin", None)
    import mtpy_example

    return register_plugin(mtpy_example.make_plugin(), replace=True)


def test_example_imports_without_mtpy(mtpy_plugin):
    import mtpy_example

    assert get_plugin("mtpy_example") is mtpy_plugin
    assert isinstance(mtpy_plugin.registry, WorkflowRegistry)
    assert isinstance(mtpy_plugin.validator, PlanValidator)
    if not mtpy_example.mtpy_available():
        assert "mtpy" not in sys.modules


@pytest.mark.parametrize("text, wf", [
    ("plot apparent resistivity and phase for all sites",
     "read_edi_plot_rhophi"),
    ("map the phase tensor ellipses", "phase_tensor"),
    ("estimate geoelectric strike", "strike"),
    ("what is the weather", None),
])
def test_example_classifies(mtpy_plugin, text, wf):
    assert mtpy_plugin.registry.classify(text) == wf


def _plan(wf, path, **params):
    return WorkflowPlan(request="r", workflow_type=wf,
                        data_path=str(path), parameters=params)


def test_example_validation(mtpy_plugin, tmp_path):
    val = mtpy_plugin.validator
    empty = tmp_path / "empty"
    empty.mkdir()
    (tmp_path / "S01.edi").write_text(">HEAD\n>END\n")
    assert val.errors(_plan("strike", tmp_path)) == []
    assert any("does not exist" in e for e in
               val.errors(_plan("strike", tmp_path / "nope")))
    assert any(".edi" in e for e in val.errors(_plan("strike", empty)))
    assert val.errors(_plan("bogus", tmp_path))
    good = _plan("phase_tensor", tmp_path, freq_band=(0.1, 1e4),
                 strike_tolerance_deg=10, n_periods=30)
    assert val.errors(good) == []
    for bad in ({"freq_band": (100.0, 1.0)},
                {"freq_band": (1e-9, 1.0)},
                {"freq_band": (float("nan"), 1.0)},
                {"freq_band": "abc"},
                {"strike_tolerance_deg": 0},
                {"n_periods": 0}, {"n_periods": True}):
        assert val.errors(_plan("strike", tmp_path, **bad)), bad


def test_example_cross_check_and_review(mtpy_plugin):
    val = mtpy_plugin.validator
    text = "estimate geoelectric strike"
    assert val.cross_check(text, "workflow", "phase_tensor") is None
    assert val.cross_check(text, "workflow", "read_edi_plot_rhophi")
    assert val.route_flags(text, "strike") == []
    plan = WorkflowPlan(request=text, workflow_type="strike",
                        risk_flags=["x"])
    assert val.needs_review(plan) is True


def test_example_end_to_end_with_fake_executor(mtpy_plugin, tmp_path):
    """Classify, validate, gate, dispatch and record, no MTpy."""
    (tmp_path / "S01.edi").write_text(">HEAD\n>END\n")
    text = "estimate geoelectric strike for this survey"
    reg, val = mtpy_plugin.registry, mtpy_plugin.validator
    wf = reg.classify(text)
    plan = _plan(wf, tmp_path, freq_band=(1.0, 1e3))
    plan.risk_flags.extend(val.route_flags(text, wf))
    assert val.errors(plan) == []
    assert val.needs_review(plan) is False

    record = mtpy_plugin.record_factory({"request": text})
    assert isinstance(record, RecordWriter)
    assert record.corpus_manifest() is None
    results = {}
    for name, agent, upstream, _desc in reg.steps(wf):
        record.emit("dispatch_started", step=name, agent=agent)
        results[name] = {"agent": agent,
                         "input": results.get(upstream)}
        record.emit("dispatch_finished", step=name)
    out = tmp_path / "run_record.json"
    record.persist(out)
    events = json.loads(out.read_text("utf-8"))["events"]
    started = [e["step"] for e in events
               if e["event"] == "dispatch_started"]
    assert started == ["read", "phase_tensor", "strike"]
    assert results["strike"]["input"]["agent"] == "MTpyPhaseTensor"


def test_example_real_mtpy_step(mtpy_plugin):
    import mtpy_example
    from mtpy_example.plugin import STEP_FUNCTIONS

    if not mtpy_example.mtpy_available():
        with pytest.raises(ImportError, match="MTpy"):
            STEP_FUNCTIONS["read"]({})
        pytest.skip("MTpy not installed")
    assert "collection" in STEP_FUNCTIONS["read"]({})
