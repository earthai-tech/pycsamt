"""Fast branch tests for agent modules that are otherwise integration-heavy."""

from __future__ import annotations

import sys
from types import SimpleNamespace

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from pycsamt.agents._base import AgentResult, BaseAgent


class _Agent(BaseAgent):
    def __init__(self, **kwargs):
        super().__init__("test", **kwargs)

    def execute(self, input_data):
        return AgentResult("success", "ok", data=input_data)


def test_base_result_helpers_and_validation(tmp_path, monkeypatch):
    result = AgentResult(
        "success", "ok", data={"answer": 42}, warnings=["warn"],
        cost_estimate_usd=0.25,
    )
    assert result["answer"] == 42
    assert "answer" in result and result.get("missing", 7) == 7
    assert result and "warnings" in repr(result) and "$0.2500" in repr(result)
    failed = AgentResult.failed("boom", hint="fix it", elapsed=1.5)
    assert not failed and failed.error_fix_hint == "fix it"

    with pytest.raises(ValueError):
        _Agent(llm_provider="unknown")
    agent = _Agent()
    assert agent.require_keys({"a": 1}, "a", "b") == ["b"]
    assert agent.extract_json("") is None
    assert agent.extract_json("prefix [1, 2] suffix") == [1, 2]
    assert agent.extract_json("not {valid] json") is None
    assert agent._citation_paths(None) == []
    context = SimpleNamespace(citations=[
        {"source_path": "a.py"}, {"source_path": "a.py"},
        {"source_path": "b.py"}, {"other": "ignored"},
    ])
    assert agent._citation_paths(context, limit=2) == ["a.py", "b.py"]
    assert agent._save_figure(object(), None, "x") is None
    warnings = []
    monkeypatch.setattr(agent._plot_cfg, "save", lambda *a, **k: (_ for _ in ()).throw(OSError("disk")))
    assert agent._save_figure(object(), str(tmp_path), "x", warnings_list=warnings) is None
    assert "disk" in warnings[0]
    assert "_Agent" in repr(agent)


@pytest.mark.parametrize("provider", ["openai", "deepseek", "minimax"])
def test_base_openai_compatible_backends(provider, monkeypatch):
    captured = {}

    class Client:
        def __init__(self, **kwargs):
            captured.update(kwargs)
            self.chat = SimpleNamespace(completions=self)

        def create(self, **kwargs):
            captured.update(kwargs)
            usage = SimpleNamespace(prompt_tokens=10, completion_tokens=5)
            choice = SimpleNamespace(message=SimpleNamespace(content="answer"))
            return SimpleNamespace(choices=[choice], usage=usage)

    monkeypatch.setitem(sys.modules, "openai", SimpleNamespace(OpenAI=Client))
    agent = _Agent(api_key="key", llm_provider=provider)
    text, cost = getattr(agent, f"_query_{provider}")("p", "s", 0.1, 12)
    assert text == "answer" and cost >= 0
    if provider != "openai":
        assert provider in captured["base_url"]


def test_base_claude_and_gemini_backends(monkeypatch):
    class Messages:
        def create(self, **kwargs):
            return SimpleNamespace(
                content=[SimpleNamespace(text="claude")],
                usage=SimpleNamespace(input_tokens=4, output_tokens=2),
            )

    anthropic = SimpleNamespace(
        Anthropic=lambda **kwargs: SimpleNamespace(messages=Messages())
    )
    monkeypatch.setitem(sys.modules, "anthropic", anthropic)
    agent = _Agent(api_key="key")
    assert agent._query_claude("p", "s", 0.2, 8)[0] == "claude"

    class Model:
        def __init__(self, *args, **kwargs):
            pass

        def generate_content(self, prompt):
            return SimpleNamespace(text="gemini response")

    genai = SimpleNamespace(
        configure=lambda **kwargs: None,
        GenerationConfig=lambda **kwargs: kwargs,
        GenerativeModel=Model,
    )
    google = SimpleNamespace(generativeai=genai)
    monkeypatch.setitem(sys.modules, "google", google)
    monkeypatch.setitem(sys.modules, "google.generativeai", genai)
    gemini = _Agent(api_key="key", llm_provider="gemini")
    assert gemini._query_gemini("two words", "s", 0.2, 8)[0] == "gemini response"


def test_query_llm_retry_and_non_rate_failure(monkeypatch):
    agent = _Agent(api_key="key", llm_provider="openai")
    attempts = {"count": 0}

    def rate_then_ok(*args):
        attempts["count"] += 1
        if attempts["count"] == 1:
            raise RuntimeError("429 rate")
        return "ok", 0.5

    monkeypatch.setattr(agent, "_query_openai", rate_then_ok)
    monkeypatch.setattr("pycsamt.agents._base.time.sleep", lambda delay: None)
    assert agent.query_llm("prompt") == "ok"
    assert agent._last_cost == 0.5
    monkeypatch.setattr(agent, "_query_openai", lambda *a: (_ for _ in ()).throw(ValueError("bad")))
    assert agent.query_llm("prompt") is None


def test_topography_validation_and_explicit_arrays():
    from pycsamt.agents._topography import resolve_agent_topography

    common = dict(sites=None, station_names=["A", "B", "C"])
    warnings = []
    assert resolve_agent_topography("yes", coords_m=None, warnings_list=warnings, **common) is None
    assert warnings
    for cfg in ({"exaggeration": 0}, {"interp_method": "bad"}, {"elevation_m": [1]}):
        assert resolve_agent_topography(cfg, coords_m=None, warnings_list=[], **common) is None
    assert resolve_agent_topography(
        {"elevation_m": [1, 2, 3]}, coords_m=None, warnings_list=[], **common
    ) is None
    meta = resolve_agent_topography(
        {"elevation_m": [10, 20, 30], "exaggeration": 2},
        coords_m=np.array([[0, 0], [300, 400], [900, 400]]),
        warnings_list=[], **common,
    )
    assert meta["applied"] and meta["source"] == "array"
    assert np.allclose(meta["chainage_km"], [0, 0.5, 1.1])
    zeros = resolve_agent_topography(
        {"elevation_m": [0, 0, 0], "chainage_km": [0, 1, 2]},
        coords_m=None, warnings_list=[], **common,
    )
    assert zeros["applied"] is False
    assert resolve_agent_topography(
        {"elevation_m": [1, 2, 3], "chainage_km": [0, 2, 1]},
        coords_m=None, warnings_list=[], **common,
    ) is None


def test_topography_sites_alignment(monkeypatch):
    from pycsamt.agents._topography import resolve_agent_topography
    import pycsamt.topo.extract as extract

    monkeypatch.setattr(extract, "extract_station_names", lambda sites: [" A ", "b"])
    monkeypatch.setattr(extract, "extract_elevation", lambda sites: [100, 120])
    monkeypatch.setattr(extract, "extract_chainage", lambda sites: [5, 6])
    meta = resolve_agent_topography(True, sites=object(), station_names=["a", "B"],
                                    coords_m=None, warnings_list=[])
    assert meta["applied"] and np.allclose(meta["chainage_km"], [0, 1])
    warnings = []
    meta = resolve_agent_topography(True, sites=object(), station_names=["missing"],
                                    coords_m=None, warnings_list=warnings)
    assert meta["applied"] is False and warnings


def test_ai_inversion_from_pretrained(monkeypatch, tmp_path):
    from pycsamt.agents.ai_inversion import AIInversionAgent
    import pycsamt.ai._zoo as zoo

    checkpoint = tmp_path / "model.pt"
    monkeypatch.setattr(zoo, "get_pretrained_info", lambda name: {"arch": "cnn1d", "n_layers": 4})
    monkeypatch.setattr(zoo, "download_checkpoint", lambda *args, **kwargs: checkpoint)
    agent = AIInversionAgent.from_pretrained("tiny", cache_dir=str(tmp_path), force_download=True)
    assert agent.arch == "cnn1d" and agent.n_layers == 4
    assert agent.pretrained == str(checkpoint)


def test_context_regex_normalisation_validation_and_llm(monkeypatch, tmp_path):
    from pycsamt.agents.context import (
        ContextInputAgent, _normalise_config, _regex_extract, _validate_config,
    )

    cfg = _regex_extract(
        'load "/data/input" save to /tmp/out frequency 10 to 100 hz '
        "off-diagonal occam2d depth_max=500 m station=S01"
    )
    assert cfg["output_dir"] and cfg["period_range"] == [0.01, 0.1]
    assert cfg["component"] == "off_diagonal"
    assert cfg["inversion_code"] == "occam2d"
    assert cfg["depth_max_km"] == 0.5 and cfg["station"] == "S01"

    normal = _normalise_config(
        {"workflow": "", "data_path": tmp_path, "output_dir": tmp_path,
         "period_range": [10, 1], "component": "invalid"}, "request",
    )
    assert normal["workflow"] == "qc" and normal["period_range"] == [1, 10]
    assert normal["component"] == "xy"
    warnings = _validate_config({"period_range": [0, 0]})
    assert len(warnings) == 3

    rag = SimpleNamespace(context_text="known symbols", citations=[{"source_path": "api.py"}])
    agent = ContextInputAgent(api_key="key")
    monkeypatch.setattr(agent, "_retrieve_context", lambda *a, **k: rag)
    answers = iter(['{"workflow":"forward","data_path":"/data/x"}', "summary"])
    monkeypatch.setattr(agent, "query_llm", lambda *a, **k: next(answers))
    result = agent.execute({"request": "run frequency decimation on /data/x"})
    assert result.status == "success"
    assert result.data["config"]["workflow"] == "freq_decimation"
    assert result.data["rag_citations"] == ["api.py"]
    assert result.llm_interpretation == "summary"


def test_code_generation_all_workflow_blocks(tmp_path, monkeypatch):
    from pycsamt.agents.code_gen import CodeGenerationAgent

    layered = SimpleNamespace(resistivity=[10, 20], thickness=[100])
    results = {name: AgentResult("success", "ok") for name in (
        "static_shift", "phase_analysis", "pre_inversion", "tipper",
        "sensitivity", "ensemble",
    )}
    results["forward"] = AgentResult("success", "ok", {"layered_model": layered})
    results["ai_inv"] = AgentResult("success", "ok", {"n_layers": 3})
    agent = CodeGenerationAgent(api_key="key")
    monkeypatch.setattr(agent, "query_llm", lambda *a, **k: "```python\n# refined\n```" )
    result = agent.execute({
        "workflow_config": {"workflow": "full_ai_workflow", "n_members": 2},
        "results": results, "output_dir": str(tmp_path), "rag_context": "real API",
    })
    assert result.status == "success" and result.data["script_path"]
    assert result.data["code"] == "# refined"


def test_batch_survey_edge_paths(monkeypatch):
    import pycsamt.agents.batch_survey as mod

    agent = mod.BatchSurveyAgent()
    assert agent.execute({"profiles": {}}).status == "failed"
    assert mod._run_pipeline([], {}, api_key=None, model=None, llm_provider="claude").status == "failed"
    good = AgentResult("success", "ok", {"n_stations": 2}, warnings=["w"])
    monkeypatch.setattr(mod, "_run_pipeline", lambda *a, **k: good)
    monkeypatch.setattr(mod, "_plot_batch_summary", lambda results: None)
    result = agent.execute({"profiles": ["a", "b"], "workflow": "unknown", "extra": 1})
    assert result.status == "success" and result.data["n_success"] == 2
    assert result.data["workflow"] == "qc" and result.warnings
    assert mod._agent_module("SomethingAgent") == "somethingagent"


def test_anomaly_agent_failure_branches(monkeypatch):
    import pycsamt.agents.anomaly_agent as mod
    import pycsamt.backends as backends
    import pycsamt.emtools._core as core

    monkeypatch.setattr(backends, "get_backend_instance", lambda: None)
    assert mod.AnomalyDetectionAgent().execute({"sites": object()}).status == "failed"
    monkeypatch.setattr(backends, "get_backend_instance", lambda: object())
    assert mod.AnomalyDetectionAgent().execute({}).status == "failed"
    monkeypatch.setattr(core, "ensure_sites", lambda value, verbose=0: (_ for _ in ()).throw(ValueError("bad sites")))
    assert "bad sites" in mod.AnomalyDetectionAgent().execute({"sites": object()}).error

    class Detector:
        def __init__(self, **kwargs):
            pass

        def fit(self, *args, **kwargs):
            raise RuntimeError("training broke")

    import pycsamt.ai.processing.anomaly as anomaly_impl
    import pycsamt.agents.ai_inversion as ai_mod

    monkeypatch.setattr(anomaly_impl, "AnomalyDetector", Detector)
    monkeypatch.setattr(core, "ensure_sites", lambda value, verbose=0: [1, 2, 3])
    monkeypatch.setattr(core, "_iter_items", iter)
    monkeypatch.setattr(core, "_name", lambda item, i: f"S{i}")
    monkeypatch.setattr(core, "_get_z_block", lambda item: (None, np.ones((2, 2, 2)), np.array([1, 2])))
    monkeypatch.setattr(ai_mod, "_z_to_features", lambda *args: np.ones(80))
    result = mod.AnomalyDetectionAgent().execute({"sites": object()})
    assert result.status == "failed" and "training broke" in result.error


def test_coordinator_required_optional_resume_and_cleanup(tmp_path):
    from pycsamt.agents.coordinator import AgentCoordinator

    coord = AgentCoordinator("edges", checkpoint_dir=str(tmp_path), verbose=True)
    agent = _Agent()
    coord.add_step("optional", agent, required=False,
                   input_fn=lambda results: (_ for _ in ()).throw(ValueError("skip")))
    coord.add_step("run", agent)
    with pytest.raises(ValueError):
        coord.add_step("run", agent)
    first = coord.execute({"x": 1})
    assert first.status == "success" and first.warnings
    resumed = coord.execute({"x": 2}, resume=True)
    assert resumed.status == "success" and "run" in resumed.data
    assert "steps=2" in repr(coord)
    coord.reset_checkpoints()
    assert not tmp_path.exists()

    broken = AgentCoordinator("broken", checkpoint_dir=str(tmp_path / "b"))
    broken.add_step("boom", _Agent(), input_fn=lambda r: 1 / 0)
    assert broken.execute({}).status == "failed"

    class Raising(_Agent):
        def execute(self, input_data):
            raise RuntimeError("explode")

    crashed = AgentCoordinator("crashed", checkpoint_dir=str(tmp_path / "c"))
    crashed.add_step("boom", Raising())
    assert crashed.execute({}).status == "failed"
