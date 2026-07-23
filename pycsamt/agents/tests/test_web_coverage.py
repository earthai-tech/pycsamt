"""Fast unit coverage for the optional Gradio front end.

The web callbacks are deliberately tested with result doubles: their job is
formatting agent output and maintaining session state, not re-testing the
numerical agents or starting a web server.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from pycsamt.agents import web
from pycsamt.agents._base import AgentResult


class _Figure:
    def savefig(self, path, **_kwargs):
        self.saved = path


class _Agent:
    def __init__(self, result):
        self.result = result
        self.calls = []

    def execute(self, payload):
        self.calls.append(payload)
        return self.result


def _result(**data):
    return AgentResult(
        status="success",
        summary="completed",
        data=data,
        warnings=["check"],
        llm_interpretation="interpreted",
    )


@pytest.fixture(autouse=True)
def _clean_session(monkeypatch, tmp_path):
    web._SESSION.clear()
    monkeypatch.setattr(web.tempfile, "mkdtemp", lambda **_: str(tmp_path))
    monkeypatch.setattr("matplotlib.pyplot.close", lambda *_: None)
    yield
    web._SESSION.clear()


def test_agent_cache(monkeypatch):
    import pycsamt.agents as package
    import pycsamt.agents.orchestrator as orch

    class Stub:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

    names = [
        "AIInversionAgent", "DataQCAgent", "EnsembleAgent",
        "ForwardModelAgent", "Inv2DAgent", "Inv3DAgent",
        "JointInversionAgent", "ModelZooAgent", "MTLoaderAgent",
        "PhaseAnalysisAgent", "ReportAgent", "StaticShiftAgent",
    ]
    for name in names:
        monkeypatch.setattr(package, name, Stub)
    monkeypatch.setattr(orch, "WorkflowOrchestratorAgent", Stub)

    first = web._get_agent("qc", "secret", "openai")
    assert first.kwargs["api_key"] == "secret"
    assert web._get_agent("qc", "secret", "openai") is first
    assert web._get_agent("qc", "new", "openai") is not first


def test_chat_and_load_qc_callbacks(monkeypatch, tmp_path):
    chat = _Agent(_result(workflow_type="qc"))
    monkeypatch.setattr(web, "_get_agent", lambda name, *_: chat)
    assert "Workflow" in web._chat_run("qc", "input", "", "", "claude")

    assert web._load_qc_run("missing", "", "claude")[1:] == (None, None)
    data = tmp_path / "edis"
    data.mkdir()
    loader = _Agent(_result(
        sites=["S1"], n_stations=1, summary_stats={"mean_qc_score": 91.0}
    ))
    qc = _Agent(_result(
        n_flagged=1,
        flagged_stations=["S1"],
        figures={"confidence_section": _Figure(), "confidence_profile": _Figure()},
    ))
    monkeypatch.setattr(
        web, "_get_agent", lambda name, *_: loader if name == "loader" else qc
    )
    text, fig1, fig2 = web._load_qc_run(str(data), "", "claude")
    assert "1 stations" in text and fig1 and fig2
    assert web._SESSION["_last_sites"] == ["S1"]

    loader.result = AgentResult.failed("bad input")
    assert "Loader failed" in web._load_qc_run(str(data), "", "claude")[0]


def test_processing_callbacks(monkeypatch):
    assert "Load data" in web._ss_run("ama", "", "claude")[0]
    assert "Load data" in web._pt_run("", "claude")[0]
    assert "Invalid" in web._fwd_run("x", "1", "", "claude")[0]
    assert "Load data" in web._ai_inv_run("mlp", 1, "", "claude")[0]
    assert "Load data" in web._inv2d_run(1, "", "claude")[0]
    assert "Load data" in web._ensemble_run(2, 1, "", "claude")[0]
    assert "Load data" in web._joint_run("tem", 1, "", "claude")[0]
    assert "Load data" in web._inv3d_run(1, 2.0, 3, "", "claude")[0]

    web._SESSION["_last_sites"] = ["S1"]
    figures = {
        "ss_comparison": _Figure(),
        "pt_psection": _Figure(),
        "strike_analysis": _Figure(),
        "survey_fingerprint": _Figure(),
        "response_and_model": _Figure(),
        "ai_section": _Figure(),
        "inv2d_section": _Figure(),
        "uncertainty_section": _Figure(),
        "uncertainty_profile": _Figure(),
        "joint_section": _Figure(),
        "depth_slices": _Figure(),
        "resistivity_section": _Figure(),
    }
    result = _result(
        corrected_sites=["S2"], delta_stats={"mean": 1, "max": 2, "n_shifted": 1},
        strike_consensus=30, strike_iqr=4, n_1d=1, n_2d=2, n_3d=3,
        rms_global=0.2, coverage=0.9, n_edges=4, figures=figures,
    )
    agent = _Agent(result)
    monkeypatch.setattr(web, "_get_agent", lambda *_: agent)
    monkeypatch.setattr(web, "AIInversionAgent_lazy", lambda **_: agent)

    assert web._ss_run("ama", "", "claude")[1]
    assert all(web._pt_run("", "claude")[1:])
    assert web._fwd_run("100, 10", "200", "", "claude")[1]
    assert web._ai_inv_run("mlp", 1, "", "claude")[1]
    assert web._inv2d_run(1, "", "claude")[1]
    assert all(web._ensemble_run(2, 1, "", "claude")[1:])
    assert web._joint_run("tem", 1, "", "claude")[1]
    assert all(web._inv3d_run(1, 2.0, 3, "", "claude")[1:])


def test_zoo_and_report_callbacks(monkeypatch):
    listing = _Agent(_result(details=[{
        "name": "tiny", "arch": "mlp", "n_layers": 4,
        "solver": "adam", "description": "A tiny model",
    }]))
    monkeypatch.setattr(web, "_get_agent", lambda *_: listing)
    assert "`tiny`" in web._zoo_list_run("", "claude")
    assert "Enter a model" in web._zoo_predict_run(" ", "", "claude")[0]
    assert "Load data" in web._zoo_predict_run("tiny", "", "claude")[0]
    assert "Run at least" in web._report_run("", "", "claude")[0]

    web._SESSION["_last_sites"] = ["S1"]
    listing.result = _result(figures={"ai_section": _Figure()})
    assert web._zoo_predict_run("tiny", "", "claude")[1]

    web._SESSION["_last_qc"] = _result()
    listing.result = _result(report_path_md="report.md")
    assert web._report_run("Title", "", "claude")[1] == "report.md"


class _Widget:
    def __init__(self, *args, **kwargs):
        self.args, self.kwargs = args, kwargs

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return False

    def click(self, *args, **kwargs):
        return self


class _Gradio:
    Blocks = Row = Column = Tabs = Tab = TabItem = Accordion = _Widget
    Markdown = Textbox = Dropdown = Button = Plot = Image = Slider = Number = _Widget
    Radio = File = _Widget
    themes = SimpleNamespace(Soft=lambda: object())


def test_build_and_launch_with_gradio_double(monkeypatch):
    monkeypatch.setattr(web, "_require_gradio", lambda: _Gradio)
    app = web.build_app()
    assert isinstance(app, _Widget)
    launched = {}
    app.launch = lambda **kwargs: launched.update(kwargs)
    monkeypatch.setattr(web, "build_app", lambda: app)
    web.launch(server_name="127.0.0.1", server_port=9999)
    assert launched["server_port"] == 9999


def test_require_gradio_import_error(monkeypatch):
    import builtins

    real_import = builtins.__import__

    def deny(name, *args, **kwargs):
        if name == "gradio":
            raise ImportError
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", deny)
    with pytest.raises(ImportError, match="pip install gradio"):
        web._require_gradio()
