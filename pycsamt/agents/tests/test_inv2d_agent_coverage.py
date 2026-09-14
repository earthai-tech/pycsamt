# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.agents.inv2d_agent`.

Targets the cheap, deterministic guard clauses at the top of
``Inv2DAgent.execute()``: the no-DL-backend guard, the ``ensure_sites``
exception, the invalid-``freqs`` / invalid-``depth_max`` guards, the
per-station "no Z" / "bad data" skip branches, and the
fewer-than-3-usable-stations guard. The heavier ``mt2d`` / ``mt2d_tri``
training paths are left to the existing inversion battery tests, which
already require and exercise a real DL backend.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.agents.inv2d_agent import Inv2DAgent

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


class _FakeZ:
    def __init__(self, z, freq):
        self.z = z
        self.freq = freq


class _FakeSite:
    def __init__(self, z=None, freq=None, station="st"):
        self.station = station
        if z is not None:
            self.Z = _FakeZ(z, freq)


def _zblock(n=8, seed=0):
    rng = np.random.default_rng(seed)
    z = rng.normal(size=(n, 2, 2)) + 1j * rng.normal(size=(n, 2, 2))
    freq = np.linspace(1.0, 100.0, n)
    return z, freq


def _passthrough_ensure_sites(monkeypatch, sites):
    import pycsamt.emtools._core as core

    monkeypatch.setattr(core, "ensure_sites", lambda x, **k: sites)


def _mock_backend(monkeypatch, available=True):
    import pycsamt.backends as backends

    monkeypatch.setattr(
        backends,
        "get_backend_instance",
        lambda: object() if available else None,
    )


def test_no_dl_backend_fails(monkeypatch):
    _mock_backend(monkeypatch, available=False)
    agent = Inv2DAgent()
    result = agent.execute({"sites": [object()]})
    assert result.status == "failed"
    assert "requires PyTorch or TensorFlow" in result.error


def test_ensure_sites_exception_fails(monkeypatch):
    import pycsamt.emtools._core as core

    _mock_backend(monkeypatch, available=True)
    monkeypatch.setattr(
        core,
        "ensure_sites",
        lambda *a, **k: (_ for _ in ()).throw(ValueError("bad sites")),
    )
    agent = Inv2DAgent()
    result = agent.execute({"path": "/nope"})
    assert result.status == "failed"
    assert "bad sites" in result.error


def test_invalid_freqs_fails(monkeypatch):
    z, fr = _zblock()
    sites = [_FakeSite(z, fr, f"s{i}") for i in range(4)]
    _mock_backend(monkeypatch, available=True)
    _passthrough_ensure_sites(monkeypatch, sites)
    agent = Inv2DAgent()
    result = agent.execute({"sites": sites, "freqs": [1.0]})
    assert result.status == "failed"
    assert "freqs" in result.error


def test_invalid_depth_max_fails(monkeypatch):
    z, fr = _zblock()
    sites = [_FakeSite(z, fr, f"s{i}") for i in range(4)]
    _mock_backend(monkeypatch, available=True)
    _passthrough_ensure_sites(monkeypatch, sites)
    agent = Inv2DAgent()
    result = agent.execute({"sites": sites, "depth_max": -10.0})
    assert result.status == "failed"
    assert "depth_max" in result.error


def test_stations_with_no_z_or_bad_data_are_skipped(monkeypatch):
    z_ok, fr_ok = _zblock(n=8, seed=1)
    no_z = _FakeSite(z=None, station="no_z")
    # a single-sample Z makes `_z_to_features` return None (needs >= 2
    # finite samples) -> the "skipped (bad data)" branch.
    z_bad = np.ones((1, 2, 2), dtype=complex)
    fr_bad = np.array([1.0])
    bad_data = _FakeSite(z_bad, fr_bad, station="bad_data")
    good = [
        _FakeSite(z_ok, fr_ok, station=f"s{i}") for i in range(3)
    ]
    sites = [no_z, bad_data, *good]
    _mock_backend(monkeypatch, available=True)
    _passthrough_ensure_sites(monkeypatch, sites)
    agent = Inv2DAgent()
    result = agent.execute({"sites": sites, "freqs": np.logspace(-2, 2, 8)})
    # A real backend is not installed/mocked here, so training itself
    # fails downstream -- but the per-station skip warning collected
    # earlier must still surface on that later failure (see the
    # AgentResult.failed(..., warnings=warnings) fix applied across this
    # file's guard clauses).
    assert any("skipped (bad data)" in w for w in result.warnings)


def test_fewer_than_three_usable_stations_fails(monkeypatch):
    z_ok, fr_ok = _zblock(n=8, seed=2)
    sites = [_FakeSite(z_ok, fr_ok, station="only_one")]
    _mock_backend(monkeypatch, available=True)
    _passthrough_ensure_sites(monkeypatch, sites)
    agent = Inv2DAgent()
    result = agent.execute({"sites": sites, "freqs": np.logspace(-2, 2, 8)})
    assert result.status == "failed"
    assert "Fewer than 3 usable stations" in result.error
