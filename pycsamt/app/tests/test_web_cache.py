# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.web.cache."""

from __future__ import annotations

import pytest

from pycsamt.app.web import cache as cache_mod


@pytest.fixture(autouse=True)
def _clean_mem():
    """Isolate the in-memory fallback dict across tests."""
    saved = dict(cache_mod._MEM)
    cache_mod._MEM.clear()
    yield
    cache_mod._MEM.clear()
    cache_mod._MEM.update(saved)


def test_cache_set_get_roundtrip_disk_backed():
    cache_mod.cache_set("sess-1", {"n_stations": 5})
    assert cache_mod.cache_get("sess-1") == {"n_stations": 5}


def test_cache_set_noop_without_session_id():
    cache_mod.cache_set("", {"x": 1})
    cache_mod.cache_set(None, {"x": 1})
    assert cache_mod._MEM == {}


def test_cache_get_returns_none_without_session_id():
    assert cache_mod.cache_get("") is None
    assert cache_mod.cache_get(None) is None


def test_cache_get_falls_back_to_memory_when_disk_cache_absent(monkeypatch):
    monkeypatch.setattr(cache_mod, "cache", None)
    cache_mod.cache_set("sess-mem", {"n_stations": 3})
    assert cache_mod.cache_get("sess-mem") == {"n_stations": 3}


def test_cache_set_get_swallow_disk_errors(monkeypatch):
    class _BoomCache:
        def set(self, *a, **k):
            raise RuntimeError("disk full")

        def get(self, *a, **k):
            raise RuntimeError("disk error")

    monkeypatch.setattr(cache_mod, "cache", _BoomCache())
    cache_mod.cache_set("sess-err", {"n_stations": 1})
    # falls back to the in-memory copy despite the disk error
    assert cache_mod.cache_get("sess-err") == {"n_stations": 1}


def test_cache_get_disk_miss_falls_back_to_memory(monkeypatch):
    class _MissCache:
        def set(self, *a, **k):
            pass

        def get(self, *a, **k):
            return None

    monkeypatch.setattr(cache_mod, "cache", _MissCache())
    cache_mod._MEM["sess-miss"] = {"n_stations": 9}
    assert cache_mod.cache_get("sess-miss") == {"n_stations": 9}


def test_has_diskcache_reflects_bg_manager_flag(monkeypatch):
    monkeypatch.setattr(cache_mod, "HAS_BG_MANAGER", True)
    assert cache_mod.has_diskcache() is True
    monkeypatch.setattr(cache_mod, "HAS_BG_MANAGER", False)
    assert cache_mod.has_diskcache() is False


# ── cache_merge_sites ────────────────────────────────────────────────────────


class _FakeEdi:
    def __init__(self, station):
        self.station = station


class _FakeSites:
    def __init__(self, edic):
        self.edic = list(edic)


def test_cache_merge_sites_stores_new_when_nothing_cached(monkeypatch):
    monkeypatch.setattr(cache_mod, "cache", None)
    new_sites = _FakeSites([_FakeEdi("S1")])
    result = cache_mod.cache_merge_sites("sess-new", new_sites)
    assert result is new_sites
    assert cache_mod.cache_get("sess-new") is new_sites


def _make_edi(tmp_path_factory, station: str):
    """Build a minimal synthetic EDI file for *station* and parse it."""
    import numpy as np

    from pycsamt.seg.edi import EDIFile

    tmp = tmp_path_factory.mktemp("edi_" + station)
    edi_path = tmp / f"{station}.edi"

    nfreq = 8
    freqs = np.logspace(2, -1, nfreq)
    z_real = np.ones(nfreq) * 10.0

    lines = [
        ">HEAD",
        f" DATAID={station}",
        " LAT=48:30:0.0",
        " LONG=7:45:0.0",
        " ELEV=200.0",
        ">INFO",
        " MAXINFO=999",
        ">=DEFINEMEAS",
        " MAXCHAN=7",
        " MAXRUN=999",
        " MAXMEAS=9999",
        " UNITS=M",
        " REFTYPE=CART",
        " REFLAT=48:30:0.0",
        " REFLONG=7:45:0.0",
        " REFELEV=200.0",
        ">=MTSECT",
        f" SECTID={station}",
        f" NFREQ={nfreq}",
        " HX=1001.001",
        " HY=1002.001",
        " HZ=1003.001",
        " EX=1004.001",
        " EY=1005.001",
        f">FREQ // {nfreq}",
    ]
    lines.append("  " + "  ".join(f"{f:.6E}" for f in freqs))
    for comp in (
        "ZXXR",
        "ZXXI",
        "ZXYR",
        "ZXYI",
        "ZYXR",
        "ZYXI",
        "ZYYR",
        "ZYYI",
    ):
        lines.append(f">{comp} // {nfreq}")
        lines.append("  " + "  ".join(f"{v:.6E}" for v in z_real))
    lines.append(">END")

    edi_path.write_text("\n".join(lines), encoding="utf-8")
    return EDIFile(str(edi_path))


def test_cache_merge_sites_deduplicates_by_station_name(tmp_path_factory):
    from pycsamt.site.base import Sites

    existing = Sites(
        [_make_edi(tmp_path_factory, "S1"), _make_edi(tmp_path_factory, "S2")]
    )
    cache_mod.cache_set("sess-merge", existing)

    new_sites = Sites(
        [_make_edi(tmp_path_factory, "S2"), _make_edi(tmp_path_factory, "S3")]
    )

    result = cache_mod.cache_merge_sites("sess-merge", new_sites)

    assert isinstance(result, Sites)
    names = sorted(getattr(e, "station", None) for e in result.edic)
    assert names == ["S1", "S2", "S3"]


def test_cache_merge_sites_falls_back_to_new_on_merge_failure(monkeypatch):
    monkeypatch.setattr(cache_mod, "cache", None)
    existing = _FakeSites([_FakeEdi("S1")])
    cache_mod.cache_set("sess-broken", existing)

    class _BrokenNewSites:
        # no .edic attribute -> AttributeError inside the merge try-block
        pass

    broken = _BrokenNewSites()
    result = cache_mod.cache_merge_sites("sess-broken", broken)
    assert result is broken
    assert cache_mod.cache_get("sess-broken") is broken


# ── forward-response cache ────────────────────────────────────────────────────


def test_cache_fwd_roundtrip_and_empty_session():
    cache_mod.cache_set_fwd("sess-1", "2d", {"resp": 1})
    assert cache_mod.cache_get_fwd("sess-1", "2d") == {"resp": 1}
    assert cache_mod.cache_get_fwd("sess-1", "3d") is None

    cache_mod.cache_set_fwd("", "2d", {"resp": 1})
    assert cache_mod.cache_get_fwd("", "2d") is None


# ── inversion-result cache ────────────────────────────────────────────────────


def test_cache_inversion_result_roundtrip_and_empty_session():
    cache_mod.cache_set_inversion_result("sess-1", {"rms": 1.2})
    assert cache_mod.cache_get_inversion_result("sess-1") == {"rms": 1.2}

    cache_mod.cache_set_inversion_result("", {"rms": 1.2})
    assert cache_mod.cache_get_inversion_result("") is None


# ── frequency-edit cache ─────────────────────────────────────────────────────


def test_cache_freq_edit_roundtrip_and_empty_session():
    cache_mod.cache_set_freq_edit("sess-1", {"edited": True})
    assert cache_mod.cache_get_freq_edit("sess-1") == {"edited": True}

    cache_mod.cache_set_freq_edit("", {"edited": True})
    assert cache_mod.cache_get_freq_edit("") is None
