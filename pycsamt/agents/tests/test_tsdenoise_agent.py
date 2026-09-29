# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for :mod:`pycsamt.agents.tsdenoise_agent`.

Real data: ``data/MT/TS/kap103as.ts`` (station kap103, 5 channels, 5 s
sampling), the same file used throughout
``docs/scripts/generate_user_guide_ai_processing_tsdenoise_figures.py``.
A short real slice around a known EY interference burst keeps the test
fast without resorting to synthetic data.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")

from pycsamt.agents.tsdenoise_agent import TimeSeriesDenoisingAgent

_ROOT = Path(__file__).resolve().parents[3]
_TS_PATH = _ROOT / "data" / "MT" / "TS" / "kap103as.ts" / "kap103as.ts"
_HAS_TS = _TS_PATH.exists()


def _load_ey_slice():
    from pycsamt.ts import TSData, read_ts

    ts_full = read_ts(str(_TS_PATH))
    s, e = 341100, 341900  # ~66-minute window around a real EY burst
    ey = ts_full.get("EY")[s:e]
    ts_sub = TSData(dt=ts_full.dt, station=ts_full.station)
    ts_sub.add_channel("EY", ey)
    return ts_sub


def test_no_ts_or_ts_path(no_llm_kw):
    agent = TimeSeriesDenoisingAgent(**no_llm_kw)
    result = agent.execute({})
    assert result.status == "failed"
    assert "No 'ts' or 'ts_path'" in result.error


def test_sites_only_input_is_not_mistaken_for_a_time_series(no_llm_kw):
    """AgentWorker always sets input_data['sites'] -- confirm this agent
    genuinely needs its own 'ts'/'ts_path' key and does not silently
    accept the desktop app's default EDI-Sites payload."""
    agent = TimeSeriesDenoisingAgent(**no_llm_kw)
    result = agent.execute({"sites": object()})
    assert result.status == "failed"
    assert "No 'ts' or 'ts_path'" in result.error


@pytest.mark.skipif(not _HAS_TS, reason="kap103as.ts dataset not bundled")
def test_ts_path_reads_and_denoises_a_real_file(no_llm_kw, tmp_output, monkeypatch):
    """Exercises the ts_path -> read_ts() wiring against a real (small,
    pre-sliced) record -- read_ts() is monkeypatched to skip re-parsing
    the full ~27-day/461747-sample file a second time (the 'ts=' test
    below already covers that same real slice end-to-end); the data
    itself is still real, not synthetic."""
    ts_sub = _load_ey_slice()
    # read_ts is imported function-locally inside execute() (from ..ts
    # import read_ts), so patching the source module is what actually
    # takes effect at call time.
    monkeypatch.setattr("pycsamt.ts.read_ts", lambda path: ts_sub)

    agent = TimeSeriesDenoisingAgent(
        **no_llm_kw, mmf_size=241, win_seconds=300.0
    )
    result = agent.execute(
        {"ts_path": str(_TS_PATH), "output_dir": str(tmp_output)}
    )

    assert result.status == "success"
    assert result["n_windows"] > 0
    assert 0 <= result["n_flagged"] <= result["n_windows"]
    assert "EY" in result["denoised_ts"].channels()
    assert "ts_denoise_summary" in result["figures"]


@pytest.mark.skipif(not _HAS_TS, reason="kap103as.ts dataset not bundled")
def test_ts_object_input_denoises_a_real_burst(no_llm_kw):
    ts_sub = _load_ey_slice()
    agent = TimeSeriesDenoisingAgent(
        **no_llm_kw, mmf_size=241, win_seconds=300.0
    )
    result = agent.execute({"ts": ts_sub})

    assert result.status == "success"
    assert result["diagnostics"] is not None
    assert len(result["diagnostics"]) == result["n_windows"]
    raw = ts_sub.get("EY")
    denoised = result["denoised_ts"].get("EY")
    assert raw.shape == denoised.shape


def test_missing_dt_fails_cleanly(no_llm_kw):
    from pycsamt.ts import TSData

    ts = TSData(dt=None, station="no_dt")
    ts.add_channel("EX", [1.0, 2.0, 3.0])
    agent = TimeSeriesDenoisingAgent(**no_llm_kw)
    result = agent.execute({"ts": ts})
    assert result.status == "failed"
    assert "dt" in result.error.lower()


def test_no_channels_fails_cleanly(no_llm_kw):
    from pycsamt.ts import TSData

    ts = TSData(dt=5.0, station="empty")
    agent = TimeSeriesDenoisingAgent(**no_llm_kw)
    result = agent.execute({"ts": ts})
    assert result.status == "failed"
    assert "no channels" in result.error.lower()
