# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.agents.tensor_rotation`.

Targets branches left uncovered by the broader agent battery: the
``seg.ops`` import guard, ``ensure_sites`` failure, missing
``output_dir``, per-station Z-missing / rotation-failure / tipper-failure
/ existing-file / write-failure branches, the rotation-summary figure
exception, the LLM interpretation path, and the private
``_write_rotated_edi`` / ``_plot_rotation_summary`` helpers directly.
"""

from __future__ import annotations

import copy
import sys
import types

import numpy as np
import pytest

from pycsamt.agents.tensor_rotation import (
    TensorRotationAgent,
    _plot_rotation_summary,
    _write_rotated_edi,
)

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


# ── fake site helpers ────────────────────────────────────────────────────────


class _FakeZ:
    def __init__(self, z, freq):
        self.z = z
        self.freq = freq


class _FakeSite:
    def __init__(self, z=None, freq=None, tipper=None, station="st"):
        self.station = station
        if z is not None:
            self.Z = _FakeZ(z, freq)
        if tipper is not None:
            self.Tip = types.SimpleNamespace(tipper=tipper)


def _zblock(n=4, seed=0):
    rng = np.random.default_rng(seed)
    z = rng.normal(size=(n, 2, 2)) + 1j * rng.normal(size=(n, 2, 2))
    freq = np.linspace(1.0, 100.0, n)
    return z, freq


def _passthrough_ensure_sites(monkeypatch, sites):
    """Make ``ensure_sites`` return *sites* unchanged.

    The real ``ensure_sites`` converts its input into a :class:`Sites`
    container that silently drops objects it does not recognise (see
    :func:`pycsamt.emtools._core.ensure_sites`), which would discard the
    lightweight fake site doubles used below. Bypassing it keeps the
    fakes flowing through ``execute`` untouched.
    """
    import pycsamt.emtools._core as core

    monkeypatch.setattr(core, "ensure_sites", lambda x, **k: sites)


# ── execute() guard-clause tests ─────────────────────────────────────────────


def test_seg_ops_import_error(no_llm_kw, monkeypatch, tmp_output):
    fake_mod = types.ModuleType("pycsamt.seg.ops")
    monkeypatch.setitem(sys.modules, "pycsamt.seg.ops", fake_mod)
    agent = TensorRotationAgent(**no_llm_kw)
    result = agent.execute(
        {"sites": [_FakeSite()], "output_dir": str(tmp_output)}
    )
    assert result.status == "failed"
    assert "pycsamt.seg.ops not available" in result.error


def test_no_sites_or_path(no_llm_kw, tmp_output):
    agent = TensorRotationAgent(**no_llm_kw)
    result = agent.execute({"output_dir": str(tmp_output)})
    assert result.status == "failed"
    assert "No 'sites' or 'path'" in result.error


def test_ensure_sites_raises(no_llm_kw, monkeypatch, tmp_output):
    import pycsamt.emtools._core as core

    monkeypatch.setattr(
        core, "ensure_sites", lambda *a, **k: (_ for _ in ()).throw(
            ValueError("bad")
        )
    )
    agent = TensorRotationAgent(**no_llm_kw)
    result = agent.execute(
        {"sites": [_FakeSite()], "output_dir": str(tmp_output)}
    )
    assert result.status == "failed"
    assert "bad" in result.error


def test_no_output_dir(no_llm_kw):
    agent = TensorRotationAgent(**no_llm_kw)
    result = agent.execute({"sites": [_FakeSite()]})
    assert result.status == "failed"
    assert "output_dir" in result.error


def test_station_with_no_z_is_skipped(no_llm_kw, monkeypatch, tmp_output):
    site = _FakeSite(z=None, station="no_z")
    _passthrough_ensure_sites(monkeypatch, [site])
    agent = TensorRotationAgent(**no_llm_kw, strike_deg=20.0)
    result = agent.execute(
        {"sites": [site], "output_dir": str(tmp_output)}
    )
    assert result.status == "needs_review"
    assert "no_z" in result["failed_stations"]
    assert any("no Z data" in w for w in result.warnings)


def test_rotate_impedance_failure_is_recorded(
    no_llm_kw, monkeypatch, tmp_output
):
    import pycsamt.seg.ops as ops

    z, fr = _zblock()
    site = _FakeSite(z=z, freq=fr, station="s1")
    _passthrough_ensure_sites(monkeypatch, [site])
    monkeypatch.setattr(
        ops,
        "rotate_impedance",
        lambda z, theta: (_ for _ in ()).throw(RuntimeError("rot boom")),
    )
    agent = TensorRotationAgent(**no_llm_kw, strike_deg=10.0)
    result = agent.execute(
        {"sites": [site], "output_dir": str(tmp_output)}
    )
    assert "s1" in result["failed_stations"]
    assert any("rotate_impedance failed" in w for w in result.warnings)


def test_tipper_rotation_failure_is_recorded_but_write_continues(
    no_llm_kw, monkeypatch, tmp_output
):
    import pycsamt.seg.ops as ops

    z, fr = _zblock()
    tipper = np.ones((len(fr), 2), dtype=complex)
    site = _FakeSite(z=z, freq=fr, tipper=tipper, station="s1")
    _passthrough_ensure_sites(monkeypatch, [site])
    monkeypatch.setattr(
        ops,
        "rotate_tipper",
        lambda t, theta: (_ for _ in ()).throw(RuntimeError("tip boom")),
    )
    agent = TensorRotationAgent(**no_llm_kw, strike_deg=10.0)
    result = agent.execute(
        {"sites": [site], "output_dir": str(tmp_output)}
    )
    assert any("tipper rotation failed" in w for w in result.warnings)
    # tipper failure does not abort the station's own EDI write attempt
    assert "s1" not in result["failed_stations"] or result["written_paths"] == []


def test_existing_output_file_skipped_without_overwrite(
    no_llm_kw, monkeypatch, tmp_output
):
    z, fr = _zblock()
    site = _FakeSite(z=z, freq=fr, station="dup")
    _passthrough_ensure_sites(monkeypatch, [site])
    (tmp_output / "dup_rot.edi").write_text("placeholder")
    agent = TensorRotationAgent(**no_llm_kw, strike_deg=5.0)
    result = agent.execute(
        {"sites": [site], "output_dir": str(tmp_output)}
    )
    assert "dup" in result["failed_stations"]
    assert any("skipped (set overwrite=True)" in w for w in result.warnings)


def test_write_returns_none_marks_failed(no_llm_kw, monkeypatch, tmp_output):
    import pycsamt.agents.tensor_rotation as tr

    z, fr = _zblock()
    site = _FakeSite(z=z, freq=fr, station="s2")
    _passthrough_ensure_sites(monkeypatch, [site])
    monkeypatch.setattr(tr, "_write_rotated_edi", lambda *a, **k: None)
    agent = TensorRotationAgent(**no_llm_kw, strike_deg=5.0)
    result = agent.execute(
        {"sites": [site], "output_dir": str(tmp_output)}
    )
    assert "s2" in result["failed_stations"]
    assert result["n_written"] == 0


def test_write_rotated_edi_raising_is_caught(
    no_llm_kw, monkeypatch, tmp_output
):
    import pycsamt.agents.tensor_rotation as tr

    def _boom(*a, **k):
        raise RuntimeError("write crashed")

    z, fr = _zblock()
    site = _FakeSite(z=z, freq=fr, station="s3")
    _passthrough_ensure_sites(monkeypatch, [site])
    monkeypatch.setattr(tr, "_write_rotated_edi", _boom)
    agent = TensorRotationAgent(**no_llm_kw, strike_deg=5.0)
    result = agent.execute(
        {"sites": [site], "output_dir": str(tmp_output)}
    )
    assert "s3" in result["failed_stations"]
    assert any("EDI write failed: write crashed" in w for w in result.warnings)


def test_rotation_summary_figure_exception_is_captured(
    no_llm_kw, monkeypatch, tmp_output
):
    import pycsamt.agents.tensor_rotation as tr

    z, fr = _zblock()
    site = _FakeSite(z=z, freq=fr, station="s4")
    _passthrough_ensure_sites(monkeypatch, [site])
    monkeypatch.setattr(
        tr,
        "_plot_rotation_summary",
        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("plot boom")),
    )
    monkeypatch.setattr(tr, "_write_rotated_edi", lambda *a, **k: "ok.edi")
    agent = TensorRotationAgent(**no_llm_kw, strike_deg=5.0)
    result = agent.execute(
        {"sites": [site], "output_dir": str(tmp_output)}
    )
    assert any("Rotation summary figure" in w for w in result.warnings)


def test_llm_interpretation_called_when_api_key_set(
    monkeypatch, tmp_output
):
    import pycsamt.agents.tensor_rotation as tr

    z, fr = _zblock()
    site = _FakeSite(z=z, freq=fr, station="s5")
    _passthrough_ensure_sites(monkeypatch, [site])
    monkeypatch.setattr(tr, "_write_rotated_edi", lambda *a, **k: "ok.edi")
    monkeypatch.setattr(
        tr, "_plot_rotation_summary", lambda *a, **k: None
    )
    agent = TensorRotationAgent(api_key="fake-key", strike_deg=5.0)
    agent.query_llm = lambda *a, **k: "mocked interpretation"
    result = agent.execute(
        {"sites": [site], "output_dir": str(tmp_output)}
    )
    assert result.llm_interpretation == "mocked interpretation"


# ── _write_rotated_edi direct tests ──────────────────────────────────────────


class _FakeZWrapper:
    def __init__(self, z):
        self.z = z


class _FakeTipWrapper:
    def __init__(self, tipper):
        self.tipper = tipper


class _EdWriteNewEdi:
    def __init__(self, z, tip=None):
        self.Z = _FakeZWrapper(z)
        self.Tip = tip
        self.calls = []

    def write_new_edi(self, *, edi_fn, savepath, Z, Tipper):
        self.calls.append((edi_fn, savepath, Z, Tipper))
        return str((__import__("pathlib").Path(savepath) / edi_fn))


class _EdWriteOnly:
    def __init__(self, z, tip=None):
        self.Z = _FakeZWrapper(z)
        self.Tip = tip

    def write(self, *, savepath, new_edifn):
        return str((__import__("pathlib").Path(savepath) / new_edifn))


class _EdNoWriter:
    pass


class _EdWriteNewEdiRaises(_EdWriteNewEdi):
    def write_new_edi(self, **k):
        raise IOError("disk full")


def test_write_rotated_edi_prefers_write_new_edi(tmp_output):
    z, _ = _zblock()
    ed = _EdWriteNewEdi(z)
    out_path = str(tmp_output / "s1_rot.edi")
    warnings: list[str] = []
    result = _write_rotated_edi(ed, z, None, out_path, "s1", 10.0, warnings)
    assert result is not None
    assert ed.calls and ed.calls[0][2].z is z


def test_write_rotated_edi_includes_tipper(tmp_output):
    z, fr = _zblock()
    tip_rot = np.ones((len(fr), 2), dtype=complex)
    ed = _EdWriteNewEdi(z, tip=_FakeTipWrapper(np.zeros((len(fr), 1, 2))))
    out_path = str(tmp_output / "s2_rot.edi")
    warnings: list[str] = []
    result = _write_rotated_edi(
        ed, z, tip_rot, out_path, "s2", 10.0, warnings
    )
    assert result is not None
    assert ed.calls[0][3] is not None


def test_write_rotated_edi_falls_back_to_write(tmp_output):
    z, _ = _zblock()
    ed = _EdWriteOnly(z)
    out_path = str(tmp_output / "s3_rot.edi")
    warnings: list[str] = []
    result = _write_rotated_edi(ed, z, None, out_path, "s3", 10.0, warnings)
    assert result is not None


def test_write_rotated_edi_no_writer_warns(tmp_output):
    ed = _EdNoWriter()
    out_path = str(tmp_output / "s4_rot.edi")
    warnings: list[str] = []
    result = _write_rotated_edi(ed, None, None, out_path, "s4", 10.0, warnings)
    assert result is None
    assert any("no write method" in w for w in warnings)


def test_write_rotated_edi_exception_is_caught(tmp_output):
    z, _ = _zblock()
    ed = _EdWriteNewEdiRaises(z)
    out_path = str(tmp_output / "s5_rot.edi")
    warnings: list[str] = []
    result = _write_rotated_edi(ed, z, None, out_path, "s5", 10.0, warnings)
    assert result is None
    assert any("disk full" in w for w in warnings)


# ── _plot_rotation_summary direct tests ──────────────────────────────────────


def test_plot_rotation_summary_no_stations_returns_none():
    fig = _plot_rotation_summary([_FakeSite(z=None)], 10.0, [])
    assert fig is None


def test_plot_rotation_summary_builds_bar_chart():
    z, fr = _zblock(n=6, seed=3)
    site = _FakeSite(z=z, freq=fr, station="s1")
    fig = _plot_rotation_summary([site], 15.0, [])
    assert fig is not None
    assert fig.axes


def test_plot_rotation_summary_skips_all_zero_off_diag():
    z = np.zeros((4, 2, 2), dtype=complex)
    z[:, 0, 0] = 1.0
    fr = np.linspace(1, 4, 4)
    site = _FakeSite(z=z, freq=fr, station="flat")
    fig = _plot_rotation_summary([site], 15.0, [])
    assert fig is None
