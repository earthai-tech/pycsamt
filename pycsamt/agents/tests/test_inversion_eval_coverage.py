# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.agents.inversion_eval`.

Mocks ``ensure_sites``/``_iter_items``/``_name``/``_get_z_block``
(imported locally inside ``execute()``) and
``build_phase_tensor_table``/``plot_station_response`` with fast fakes
so the RMS loop, residual phase-tensor merge, and figure branches can
all be exercised without real EDI data: the ``ensure_sites`` exception
for both observed and model sources, the "no model response" branch,
the per-station RMS skip/continue/exception paths, the residual-PT
merge and exception branches, the figure success/None/exception
branches, and the LLM interpretation path.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from pycsamt.agents.inversion_eval import InversionEvaluationAgent

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


class _FakeED:
    def __init__(self, name, z=None, fr=None):
        self.name = name
        self.z = z
        self.fr = fr


def _z_block(rho_xy, freqs):
    """Build a fake (n, 2, 2) complex Z array whose |Zxy|^2 * 0.2/f == rho_xy."""
    z = np.zeros((len(freqs), 2, 2), dtype=complex)
    amp = np.sqrt(np.asarray(rho_xy) / (0.2 / freqs))
    z[:, 0, 1] = amp
    return z


def _mock_core(monkeypatch, *, identity=True):
    import pycsamt.emtools._core as core

    if identity:
        monkeypatch.setattr(core, "ensure_sites", lambda x, **k: x)
    monkeypatch.setattr(core, "_iter_items", lambda sites: list(sites))
    monkeypatch.setattr(core, "_name", lambda ed, i: ed.name)
    monkeypatch.setattr(
        core, "_get_z_block", lambda ed, **k: (None, ed.z, ed.fr)
    )


def _mock_pt(monkeypatch, *, raises=False, empty=False, frames=None):
    import pycsamt.emtools.tensor as tensor_mod

    calls = {"n": 0}

    def _fake(sites, verbose=0):
        if raises:
            raise RuntimeError("pt boom")
        if frames is not None:
            df = frames[calls["n"]]
            calls["n"] += 1
            return df
        if empty:
            return pd.DataFrame(columns=["station", "period"])
        return pd.DataFrame(
            {
                "station": ["S1", "S1"],
                "period": [1.0, 10.0],
                "skew": [1.0, 2.0],
                "ellipt": [0.1, 0.2],
                "theta": [10.0, 20.0],
            }
        )

    monkeypatch.setattr(tensor_mod, "build_phase_tensor_table", _fake)


# ── ensure_sites branches ───────────────────────────────────────────────────


def test_no_obs_source_fails():
    agent = InversionEvaluationAgent()
    result = agent.execute({})
    assert result.status == "failed"


def test_ensure_sites_obs_exception_fails(monkeypatch):
    import pycsamt.emtools._core as core

    monkeypatch.setattr(
        core,
        "ensure_sites",
        lambda *a, **k: (_ for _ in ()).throw(ValueError("bad obs")),
    )
    agent = InversionEvaluationAgent()
    result = agent.execute({"sites_obs": object()})
    assert result.status == "failed"
    assert "bad obs" in result.error


def test_no_model_response_warns_and_needs_review(monkeypatch):
    _mock_core(monkeypatch)
    agent = InversionEvaluationAgent()
    sites_obs = [_FakeED("S1")]
    result = agent.execute({"sites_obs": sites_obs})
    assert result.status == "needs_review"
    assert any("No model response provided" in w for w in result.warnings)


def test_ensure_sites_mod_exception_recorded_as_warning(monkeypatch):
    import pycsamt.emtools._core as core

    _mock_core(monkeypatch)
    real_ensure = core.ensure_sites

    def _flaky(x, **k):
        if isinstance(x, str) and x == "bad-model-path":
            raise RuntimeError("bad model")
        return real_ensure(x, **k)

    monkeypatch.setattr(core, "ensure_sites", _flaky)
    agent = InversionEvaluationAgent()
    result = agent.execute(
        {"sites_obs": [_FakeED("S1")], "sites_mod": "bad-model-path"}
    )
    assert any("Could not load model response" in w for w in result.warnings)
    assert result.status == "needs_review"


# ── per-station RMS loop branches ───────────────────────────────────────────


def test_rms_loop_skip_missing_and_bad_mask_and_exception(monkeypatch):
    _mock_core(monkeypatch)
    _mock_pt(monkeypatch, empty=True)
    freqs = np.array([1.0, 10.0, 100.0])

    # S1: matches, good data -> contributes to rms
    z1 = _z_block([100.0, 100.0, 100.0], freqs)
    ed_o1 = _FakeED("S1", z=z1, fr=freqs)
    ed_m1 = _FakeED("S1", z=z1, fr=freqs)

    # S2: only in model, not in obs -> "ed_o is None: continue"
    ed_m2 = _FakeED("S2", z=z1, fr=freqs)

    # S3: matches but mostly non-finite/negative rho -> mask.sum() < 2
    z3 = _z_block([-1.0, -1.0, -1.0], freqs)
    ed_o3 = _FakeED("S3", z=z3, fr=freqs)
    ed_m3 = _FakeED("S3", z=z3, fr=freqs)

    # S4: raises inside the try block (z with wrong shape)
    ed_o4 = _FakeED("S4", z=np.zeros((3,)), fr=freqs)
    ed_m4 = _FakeED("S4", z=np.zeros((3,)), fr=freqs)

    # S5: matches but one side has no Z block at all -> "z is None: continue"
    ed_o5 = _FakeED("S5", z=None, fr=freqs)
    ed_m5 = _FakeED("S5", z=z1, fr=freqs)

    sites_obs = [ed_o1, ed_o3, ed_o4, ed_o5]
    sites_mod = [ed_m1, ed_m2, ed_m3, ed_m4, ed_m5]

    agent = InversionEvaluationAgent()
    result = agent.execute({"sites_obs": sites_obs, "sites_mod": sites_mod})
    assert result.status == "success"
    assert "S1" in result["rms_per_station"]
    assert "S3" not in result["rms_per_station"]
    assert "S4" not in result["rms_per_station"]
    assert "S5" not in result["rms_per_station"]
    assert any("RMS for S4" in w for w in result.warnings)


# ── residual PT branches ────────────────────────────────────────────────────


def test_residual_pt_exception_recorded(monkeypatch):
    _mock_core(monkeypatch)
    _mock_pt(monkeypatch, raises=True)
    freqs = np.array([1.0, 10.0])
    z = _z_block([100.0, 100.0], freqs)
    ed = _FakeED("S1", z=z, fr=freqs)
    agent = InversionEvaluationAgent()
    result = agent.execute({"sites_obs": [ed], "sites_mod": [ed]})
    assert any("Residual PT computation" in w for w in result.warnings)


def test_residual_pt_empty_frames_skip_merge(monkeypatch):
    _mock_core(monkeypatch)
    _mock_pt(monkeypatch, empty=True)
    freqs = np.array([1.0, 10.0])
    z = _z_block([100.0, 100.0], freqs)
    ed = _FakeED("S1", z=z, fr=freqs)
    agent = InversionEvaluationAgent()
    result = agent.execute({"sites_obs": [ed], "sites_mod": [ed]})
    assert result["residual_pt_table"] is None


def test_residual_pt_merge_with_partial_columns(monkeypatch):
    _mock_core(monkeypatch)
    df_obs = pd.DataFrame(
        {"station": ["S1"], "period": [1.0], "skew": [1.0]}
    )
    df_mod = pd.DataFrame(
        {"station": ["S1"], "period": [1.0], "skew": [0.5]}
    )
    _mock_pt(monkeypatch, frames=[df_obs, df_mod])
    freqs = np.array([1.0, 10.0])
    z = _z_block([100.0, 100.0], freqs)
    ed = _FakeED("S1", z=z, fr=freqs)
    agent = InversionEvaluationAgent()
    result = agent.execute({"sites_obs": [ed], "sites_mod": [ed]})
    merged = result["residual_pt_table"]
    assert merged is not None
    assert "d_skew" in merged.columns
    assert "d_ellipt" not in merged.columns


# ── figure branches ──────────────────────────────────────────────────────


def test_plot_station_response_exception_recorded(monkeypatch):
    _mock_core(monkeypatch)
    _mock_pt(monkeypatch, empty=True)
    import pycsamt.emtools.inspect as inspect_mod

    def _raise(*a, **k):
        raise RuntimeError("plot boom")

    monkeypatch.setattr(inspect_mod, "plot_station_response", _raise)
    freqs = np.array([1.0, 10.0])
    z = _z_block([100.0, 100.0], freqs)
    ed = _FakeED("S1", z=z, fr=freqs)
    agent = InversionEvaluationAgent()
    result = agent.execute({"sites_obs": [ed], "sites_mod": [ed]})
    assert any("plot_station_response" in w for w in result.warnings)


def test_plot_station_response_none_skips_figure(monkeypatch):
    _mock_core(monkeypatch)
    _mock_pt(monkeypatch, empty=True)
    import pycsamt.emtools.inspect as inspect_mod

    monkeypatch.setattr(inspect_mod, "plot_station_response", lambda *a, **k: None)
    freqs = np.array([1.0, 10.0])
    z = _z_block([100.0, 100.0], freqs)
    ed = _FakeED("S1", z=z, fr=freqs)
    agent = InversionEvaluationAgent()
    result = agent.execute({"sites_obs": [ed], "sites_mod": [ed]})
    assert result["figures"] == {}


def test_plot_station_response_neither_savefig_nor_get_figure(monkeypatch):
    _mock_core(monkeypatch)
    _mock_pt(monkeypatch, empty=True)
    import pycsamt.emtools.inspect as inspect_mod

    monkeypatch.setattr(
        inspect_mod, "plot_station_response", lambda *a, **k: object()
    )
    freqs = np.array([1.0, 10.0])
    z = _z_block([100.0, 100.0], freqs)
    ed = _FakeED("S1", z=z, fr=freqs)
    agent = InversionEvaluationAgent()
    result = agent.execute({"sites_obs": [ed], "sites_mod": [ed]})
    assert result["figures"] == {}


def test_plot_station_response_savefig_without_output_dir(monkeypatch):
    _mock_core(monkeypatch)
    _mock_pt(monkeypatch, empty=True)
    import pycsamt.emtools.inspect as inspect_mod

    fig, ax = plt.subplots()
    monkeypatch.setattr(inspect_mod, "plot_station_response", lambda *a, **k: fig)
    freqs = np.array([1.0, 10.0])
    z = _z_block([100.0, 100.0], freqs)
    ed = _FakeED("S1", z=z, fr=freqs)
    agent = InversionEvaluationAgent()
    result = agent.execute({"sites_obs": [ed], "sites_mod": [ed]})
    assert "station_response" in result["figures"]
    assert "station_response" not in result["figure_paths"]


def test_plot_station_response_savefig_object_saved(monkeypatch, tmp_output):
    _mock_core(monkeypatch)
    _mock_pt(monkeypatch, empty=True)
    import pycsamt.emtools.inspect as inspect_mod

    fig, ax = plt.subplots()
    monkeypatch.setattr(inspect_mod, "plot_station_response", lambda *a, **k: fig)
    freqs = np.array([1.0, 10.0])
    z = _z_block([100.0, 100.0], freqs)
    ed = _FakeED("S1", z=z, fr=freqs)
    agent = InversionEvaluationAgent()
    result = agent.execute(
        {
            "sites_obs": [ed],
            "sites_mod": [ed],
            "output_dir": str(tmp_output),
        }
    )
    assert result["figure_paths"].get("station_response")


# ── LLM interpretation branch ───────────────────────────────────────────────


def test_llm_interpretation_branch(monkeypatch):
    _mock_core(monkeypatch)
    _mock_pt(monkeypatch, empty=True)
    freqs = np.array([1.0, 10.0])
    z = _z_block([100.0, 100.0], freqs)
    ed = _FakeED("S1", z=z, fr=freqs)
    agent = InversionEvaluationAgent(api_key="fake-key")
    agent.query_llm = lambda *a, **k: "mocked interpretation"
    result = agent.execute({"sites_obs": [ed], "sites_mod": [ed]})
    assert result.llm_interpretation == "mocked interpretation"
