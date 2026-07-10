"""Tests for pycsamt.map.inversion — importing ModEM 3-D results."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.map._core import StationRecord
from pycsamt.map.inversion import (
    group_modem_stations,
    load_modem_lines,
)
from pycsamt.map.view import MapView


class _Site:
    """Minimal duck-typed EDI-like site for ModEmData.from_edi."""

    def __init__(self, name, lat, lon, elev, freq, z):
        self.name = name
        self.coords = (lat, lon, elev)
        self.freq = freq
        self.z = z
        self.z_err = None


def _synthetic_sites(n_per_line=5, lat0=32.12, lon0=119.12):
    """Two parallel lines, named like real survey stations (survey-line-id)."""
    freq = np.array([100.0, 33.0, 10.0, 3.3, 1.0])
    n_freq = freq.size
    sites = []
    for line_idx, line in enumerate(("18", "22")):
        lon = lon0 + line_idx * 0.002
        for j in range(n_per_line):
            lat = lat0 + j * 0.0009
            name = f"23-{line}-{j:03d}A"
            z = np.zeros((n_freq, 2, 2), dtype=complex)
            z[:, 0, 1] = 100.0 + 50j
            z[:, 1, 0] = -(100.0 + 50j)
            sites.append(_Site(name, lat, lon, 300.0 + j, freq, z))
    return sites


def _write_modem_pair(tmp_path):
    """Build+write a small synthetic ModEM data + 3-D model pair."""
    from pycsamt.models.modem.config import ModEmConfig
    from pycsamt.models.modem.data import ModEmData
    from pycsamt.models.modem.model3d import ModEmModel3D

    cfg = ModEmConfig(
        component_type="Full_Impedance",
        nz=6,
        n_airlayers=1,
        n_padding_xy=2,
        cell_size_h=200.0,
        cell_size_v_top=20.0,
        depth_scale=1.5,
        initial_rho=100.0,
    )
    data = ModEmData.from_edi(_synthetic_sites(), config=cfg)
    model = ModEmModel3D.halfspace(data, config=cfg)
    # Non-uniform resistivity so the curtain slicer can be verified
    # against a known, position-dependent value.
    model.rho_loge = np.log(
        100.0 + 10.0 * np.arange(model.rho_loge.size).reshape(model.shape)
    )

    # InversionResult._scan() recognises the old m0/d0 naming convention.
    dat_path = tmp_path / "d0.dat"
    rho_path = tmp_path / "m0.rho"
    data.write(dat_path)
    model.write(rho_path)
    return tmp_path, data, model


# ---------------------------------------------------------------------------
# group_modem_stations
# ---------------------------------------------------------------------------

def test_group_modem_stations_name_token_fallback():
    names = ["23-18-001A", "23-18-002A", "23-22-001A"]
    groups = group_modem_stations(names)
    assert set(groups) == {"18", "22"}
    assert groups["18"] == ["23-18-001A", "23-18-002A"]
    assert groups["22"] == ["23-22-001A"]


def test_group_modem_stations_prefers_known_stations():
    names = ["23-18-001A", "23-22-001A"]
    known = [
        StationRecord(id="23-18-001A", line="LINE-A"),
        StationRecord(id="23-22-001A", line="LINE-B"),
    ]
    groups = group_modem_stations(names, known_stations=known)
    assert set(groups) == {"LINE-A", "LINE-B"}


# ---------------------------------------------------------------------------
# station_curtain
# ---------------------------------------------------------------------------

def test_station_curtain_samples_real_positions(tmp_path):
    from pycsamt.models.modem.section import station_curtain

    _, data, model = _write_modem_pair(tmp_path)
    line18 = [n for n in data.site_names if n.split("-")[1] == "18"]
    curtain = station_curtain(model, data, line18)

    assert curtain.rho.shape == (curtain.z.size, len(line18))
    assert len(curtain.station_names) == len(line18)
    assert np.isfinite(curtain.rho).all()
    # Depth axis is earth-layers only (air layer excluded) and increasing.
    assert curtain.z.size == model.nz - model.n_air
    assert np.all(np.diff(curtain.z) > 0)


def test_station_curtain_skips_unknown_stations(tmp_path):
    from pycsamt.models.modem.section import station_curtain

    _, data, model = _write_modem_pair(tmp_path)
    curtain = station_curtain(model, data, ["not-a-real-station"])
    assert curtain.rho.shape[1] == 0
    assert curtain.station_names == []


# ---------------------------------------------------------------------------
# load_modem_lines / MapView.from_inversion_results (end-to-end)
# ---------------------------------------------------------------------------

def test_load_modem_lines_builds_multiline_mapdata(tmp_path):
    folder, data, _ = _write_modem_pair(tmp_path)
    map_data = load_modem_lines(folder, fetch_elevation=False)

    assert set(map_data.lines) == {"18", "22"}
    assert len(map_data.stations) == len(data.site_names)
    assert map_data.has_geo  # native ModEM lat/lon resolved
    sections = map_data.metadata["sections"]
    assert set(sections) == {"18", "22"}
    for _line, section in sections.items():
        assert section["rho"].shape[1] == section["stations"].size


def test_mapview_from_inversion_results_renders_fence(tmp_path):
    folder, _, _ = _write_modem_pair(tmp_path)
    mv = MapView.from_inversion_results(folder, fetch_elevation=False)

    assert set(mv.lines) == {"18", "22"}
    fig = mv.map3d(mode="fence")
    surfaces = [t for t in fig.data if t.type == "surface"]
    assert len(surfaces) == 2
    # Real-geometry offsets: the two parallel lines must not collapse
    # onto the same cross-line (y) position.
    y_values = {float(np.nanmean(np.asarray(t.y))) for t in surfaces}
    assert len(y_values) == 2


def test_load_modem_lines_missing_folder_raises(tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(ValueError):
        load_modem_lines(empty)


# ---------------------------------------------------------------------------
# fetch_elevation (mocked — no live network call in tests)
# ---------------------------------------------------------------------------

def test_fetch_elevation_fills_only_missing_stations(tmp_path, monkeypatch):
    folder, data, _ = _write_modem_pair(tmp_path)

    def _fake_fetch(probe_data, **_kw):
        return {s.id: 999.0 for s in probe_data.stations}

    monkeypatch.setattr(
        "pycsamt.map.topo.fetch_elevations", _fake_fetch
    )

    line18 = [n for n in data.site_names if n.split("-")[1] == "18"]
    known = [StationRecord(id=line18[0], elevation=42.0)]

    map_data = load_modem_lines(
        folder, known_stations=known, fetch_elevation=True,
    )
    by_id = {s.id: s for s in map_data.stations}
    # Already-known elevation is preserved, not overwritten by the fetch.
    assert by_id[line18[0]].elevation == 42.0
    # Every other station got the (mocked) fetched value.
    others = [s for s in map_data.stations if s.id != line18[0]]
    assert others and all(s.elevation == 999.0 for s in others)


def test_fetch_elevation_false_skips_fetch(tmp_path, monkeypatch):
    folder, _, _ = _write_modem_pair(tmp_path)
    calls = []
    monkeypatch.setattr(
        "pycsamt.map.topo.fetch_elevations",
        lambda *a, **k: calls.append(1) or {},
    )
    map_data = load_modem_lines(folder, fetch_elevation=False)
    assert not calls
    assert all(s.elevation is None for s in map_data.stations)
