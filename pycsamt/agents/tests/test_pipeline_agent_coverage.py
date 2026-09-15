# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.agents.pipeline_agent`.

Runs the real ``pycsamt.pipeline`` machinery against the bundled 3-EDI
dataset (fast: basic_qc completes in well under a second) so the config
resolution, guided-mode LLM recommendation parsing, pipeline-build
branches, param-override application, and the run-exception guard are
all exercised for real. ``query_llm`` is stubbed on the agent instance
to avoid any network call.
"""

from __future__ import annotations

import pytest

from pycsamt.agents.pipeline_agent import PipelineAgent, _count

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


def test_ensure_sites_exception_fails(monkeypatch):
    import pycsamt.emtools._core as core

    monkeypatch.setattr(
        core,
        "ensure_sites",
        lambda *a, **k: (_ for _ in ()).throw(ValueError("bad sites")),
    )
    agent = PipelineAgent()
    result = agent.execute({"path": "/nope"})
    assert result.status == "failed"
    assert "bad sites" in result.error


def test_no_sites_or_path_fails():
    agent = PipelineAgent()
    result = agent.execute({})
    assert result.status == "failed"
    assert "No 'sites' or 'path'" in result.error


def test_explicit_steps_build_custom_pipeline(loaded_sites):
    agent = PipelineAgent()
    result = agent.execute(
        {"sites": loaded_sites, "steps": ["NR001", "QC001"]}
    )
    assert result.status in ("success", "needs_review")
    assert result["preset_used"] == "custom"
    assert result["steps_run"] == ["NR001", "QC001"]


def test_constructor_default_preset_used(loaded_sites):
    agent = PipelineAgent(preset="basic_qc")
    result = agent.execute({"sites": loaded_sites})
    assert result["preset_used"] == "basic_qc"


def test_no_preset_or_steps_falls_back_to_basic_qc(loaded_sites):
    agent = PipelineAgent(preset=None)
    result = agent.execute({"sites": loaded_sites, "preset": None})
    assert result["preset_used"] == "basic_qc"
    assert any(
        "defaulting to 'basic_qc'" in w for w in result.warnings
    )


def test_invalid_preset_name_fails(loaded_sites):
    agent = PipelineAgent()
    result = agent.execute(
        {"sites": loaded_sites, "preset": "not_a_real_preset_xyz"}
    )
    assert result.status == "failed"
    assert "Could not build pipeline" in result.error


def test_param_overrides_are_applied(loaded_sites):
    agent = PipelineAgent()
    result = agent.execute(
        {
            "sites": loaded_sites,
            "preset": "basic_qc",
            "param_overrides": {"NR001": {"mains_hz": 60}},
        }
    )
    assert result.status in ("success", "needs_review")
    # sanity: the override didn't crash the run
    assert result["n_sites_in"] > 0


def test_param_override_exception_is_recorded(monkeypatch, loaded_sites):
    import pycsamt.pipeline as pipe_pkg

    def _boom_replace(self, label, step):
        raise RuntimeError("replace boom")

    monkeypatch.setattr(pipe_pkg.Pipeline, "replace", _boom_replace)
    agent = PipelineAgent()
    result = agent.execute(
        {
            "sites": loaded_sites,
            "preset": "basic_qc",
            "param_overrides": {"NR001": {"mains_hz": 60}},
        }
    )
    assert any(
        "Could not apply override for step" in w for w in result.warnings
    )


def test_pipeline_run_exception_fails(monkeypatch, loaded_sites):
    import pycsamt.pipeline as pipe_pkg

    monkeypatch.setattr(
        pipe_pkg.Pipeline,
        "run",
        lambda self, *a, **k: (_ for _ in ()).throw(
            RuntimeError("run boom")
        ),
    )
    agent = PipelineAgent()
    result = agent.execute({"sites": loaded_sites, "preset": "basic_qc"})
    assert result.status == "failed"
    assert "Pipeline run failed" in result.error


def test_guided_mode_recommendation_parsed_and_applied(loaded_sites):
    calls = {"n": 0}

    def _fake_query(prompt, **k):
        calls["n"] += 1
        if calls["n"] == 1:
            return (
                '{"preset": "basic_qc", "steps": [], '
                '"param_overrides": {"NR001": {"mains_hz": 60}}, '
                '"rationale": "50/60Hz grid noise suspected"}'
            )
        return "mocked interpretation"

    agent = PipelineAgent()
    agent.query_llm = _fake_query
    result = agent.execute(
        {
            "sites": loaded_sites,
            "request": "50 Hz grid noise, possible static shift",
        }
    )
    assert result["recommendation"]["preset"] == "basic_qc"
    assert result.llm_interpretation == "mocked interpretation"


def test_guided_mode_unparseable_recommendation_warns(loaded_sites):
    calls = {"n": 0}

    def _fake_query(prompt, **k):
        calls["n"] += 1
        if calls["n"] == 1:
            return "not valid json at all"
        return "mocked interpretation"

    agent = PipelineAgent()
    agent.query_llm = _fake_query
    result = agent.execute(
        {"sites": loaded_sites, "request": "clean this up please"}
    )
    assert any(
        "could not be parsed as JSON" in w for w in result.warnings
    )
    assert result["preset_used"] == "basic_qc"


# ── _count direct tests ──────────────────────────────────────────────────────


def test_count_uses_n_sites_attribute():
    class _Obj:
        n_sites = 5

    assert _count(_Obj()) == 5


def test_count_falls_back_to_len():
    assert _count([1, 2, 3]) == 3


def test_count_returns_zero_when_everything_fails():
    class _Useless:
        pass

    assert _count(_Useless()) == 0
