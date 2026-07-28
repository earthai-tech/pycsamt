# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
# ruff: noqa: E501
"""Coverage for map topography export/fetch, geometry helpers, and
figure export writers (no plotly/kaleido required)."""

from __future__ import annotations

import io as _io
import json
from pathlib import Path

import numpy as np
import pytest

from pycsamt.map import MapView, apply_elevations, parse_elevation_file
from pycsamt.map._core import MapData, StationRecord
from pycsamt.map.export import (
    ExportOptions,
    export_figure,
    save_png,
    write_dict,
    write_html,
    write_image,
    write_json,
)
from pycsamt.map.geometry import (
    equirect_xy,
    fit_strike,
    normalize_offsets,
    resolve_offset,
    survey_uv,
)
from pycsamt.map.topo import export_elevations, fetch_elevations

# ── fake survey (no file IO) ───────────────────────────


class _Z:
    freq = [10.0, 1.0]

    def __init__(self) -> None:
        self.resistivity = np.ones((2, 2, 2)) * 100.0
        self.phase = np.ones((2, 2, 2)) * 45.0


class _Edi:
    def __init__(self, station: str) -> None:
        self.station = station
        self.Z = _Z()


def _data(with_elev: bool = True) -> MapData:
    e0 = 120.0 if with_elev else None
    e1 = 140.0 if with_elev else None
    return MapData(
        sites=[_Edi("S00"), _Edi("S01")],
        stations=(
            StationRecord("S00", 28.00, 102.00, e0, "L1", 0),
            StationRecord("S01", 28.01, 102.02, e1, "L1", 1),
        ),
    )


# ── export_elevations / parse_elevation_file roundtrips ─────────────────


def test_export_elevations_csv_roundtrip(tmp_path):
    out = export_elevations(_data(), tmp_path / "topo.csv")
    assert out.exists()
    parsed = parse_elevation_file(out.read_bytes(), out.name)
    assert parsed == {"S00": 120.0, "S01": 140.0}


def test_export_elevations_h5_roundtrip(tmp_path):
    pytest.importorskip("h5py")
    out = export_elevations(_data(), tmp_path / "topo.h5")
    assert out.exists()
    parsed = parse_elevation_file(out.read_bytes(), out.name)
    assert parsed == {"S00": 120.0, "S01": 140.0}


def test_export_elevations_fmt_override(tmp_path):
    out = export_elevations(_data(), tmp_path / "topo.txt", fmt="csv")
    text = out.read_text(encoding="utf-8")
    assert "station" in text and "S00" in text


def test_export_elevations_errors(tmp_path):
    with pytest.raises(ValueError, match="elevation"):
        export_elevations(_data(with_elev=False), tmp_path / "t.csv")
    with pytest.raises(ValueError, match="Unsupported"):
        export_elevations(_data(), tmp_path / "t.xyz", fmt="xyz")


def test_parse_npz_elevation_file():
    buf = _io.BytesIO()
    np.savez(
        buf,
        station=np.asarray(["S00", "S01"]),
        elevation=np.asarray([12.5, 13.5]),
    )
    parsed = parse_elevation_file(buf.getvalue(), "topo.npz")
    assert parsed == {"S00": 12.5, "S01": 13.5}


def test_parse_unknown_extension_returns_empty():
    assert parse_elevation_file(b"junk", "topo.weird") == {}
    assert parse_elevation_file(1234, "topo.csv") == {}


# ── fetch_elevations (API stubbed) ───────────────────────────────────────


def test_fetch_elevations_maps_ids(monkeypatch):
    import pycsamt.gis.utils as gis_utils

    def fake_api(lats, lons, api_name="open_meteo"):
        return [100.0 + i for i in range(len(list(lats)))]

    monkeypatch.setattr(gis_utils, "get_elevation_from_api", fake_api)
    out = fetch_elevations(_data(), api_name="open_meteo")
    assert out == {"S00": 100.0, "S01": 101.0}


def test_fetch_elevations_skips_nonfinite(monkeypatch):
    import pycsamt.gis.utils as gis_utils

    monkeypatch.setattr(
        gis_utils,
        "get_elevation_from_api",
        lambda lats, lons, api_name="x": [float("nan"), 55.0],
    )
    out = fetch_elevations(_data())
    assert out == {"S01": 55.0}


def test_fetch_elevations_no_coords():
    data = MapData(
        sites=[_Edi("S00")],
        stations=(StationRecord("S00", None, None, None, "L1", 0),),
    )
    assert fetch_elevations(data) == {}


# ── MapView delegation ───────────────────────────────────────────────────


def test_view_with_elevations_returns_new_view():
    view = MapView(_data(with_elev=False))
    new = view.with_elevations({"S00": 5.0, "S01": 6.0})
    assert new is not view
    assert [s.elevation for s in new.data.stations] == [5.0, 6.0]
    # original untouched
    assert [s.elevation for s in view.data.stations] == [None, None]


def test_view_fetch_and_export_topography(tmp_path, monkeypatch):
    import pycsamt.gis.utils as gis_utils

    monkeypatch.setattr(
        gis_utils,
        "get_elevation_from_api",
        lambda lats, lons, api_name="x": [1.0, 2.0],
    )
    view = MapView(_data())
    fetched = view.fetch_elevations()
    assert set(fetched) == {"S00", "S01"}

    out = view.export_topography(tmp_path / "sta_topo.csv")
    assert Path(out).exists()


# ── geometry helpers ─────────────────────────────────────────────────────


def test_equirect_xy_scales_with_latitude():
    x, y = equirect_xy([0.0, 1.0], [0.0, 1.0])
    assert y[1] - y[0] == pytest.approx(111_320, rel=0.05)
    assert x[1] > 0


def test_fit_strike_recovers_line_direction():
    rng = np.random.default_rng(0)
    t = np.linspace(0.0, 1.0, 12)
    pts = np.column_stack([t, 0.5 * t]) + rng.normal(0, 1e-4, (12, 2))
    groups = np.asarray(["L1"] * 12, dtype=object)
    d = fit_strike(pts - pts.mean(axis=0), groups)
    expected = np.array([1.0, 0.5]) / np.hypot(1.0, 0.5)
    assert d is not None
    assert abs(float(np.dot(d, expected))) == pytest.approx(1.0, abs=1e-3)


def test_fit_strike_degenerate_returns_none():
    pts = np.zeros((1, 2))
    assert fit_strike(pts, np.asarray(["L1"], dtype=object)) is None


def test_fit_strike_singleton_groups_fall_back():
    # every group has < 2 points -> global SVD fallback
    pts = np.array([[0.0, 0.0], [1.0, 1.0]])
    groups = np.asarray(["A", "B"], dtype=object)
    d = fit_strike(pts - pts.mean(axis=0), groups)
    assert d is not None
    assert np.hypot(*d) == pytest.approx(1.0)


def test_survey_uv_projection_and_empty_cases():
    uv = survey_uv(
        ["a", "b", "c"],
        [28.0, 28.0, 28.0],
        [102.00, 102.01, 102.02],
        ["L1", "L1", "L1"],
    )
    assert set(uv) == {"a", "b", "c"}
    u_vals = sorted(v[0] for v in uv.values())
    assert u_vals[0] < u_vals[1] < u_vals[2]
    # v (cross-strike) ~ 0 for a straight EW line
    assert all(abs(v[1]) < 1.0 for v in uv.values())

    assert survey_uv(["a"], [28.0], [102.0], ["L1"]) == {}
    assert survey_uv(["a", "b"], [28.0, 28.0], [102.0, 102.0], ["L1", "L1"]) == {}


def test_normalize_offsets_and_resolve():
    assert normalize_offsets({}) is None
    assert normalize_offsets({"L1": 10.0, "L2": None}) is None
    assert normalize_offsets({"L1": 10.0, "L2": float("nan")}) is None
    out = normalize_offsets({"L1": 10.0, "L2": 25.0})
    assert out == {"L1": 0.0, "L2": 15.0}

    assert resolve_offset("L2", 3, out, 100.0, 2.0) == pytest.approx(30.0)
    assert resolve_offset("L2", 3, None, 100.0, 2.0) == pytest.approx(600.0)


# ── figure export writers (fake figures, no kaleido) ─────────────────────


class _FakeFig:
    """Duck-typed Plotly-like figure without any plotly dependency."""

    def __init__(self):
        self.saved = {}

    def to_json(self):
        return json.dumps({"data": [], "layout": {"title": "fake"}})

    def to_dict(self):
        return {"data": [], "layout": {"title": "fake"}}

    def to_html(self, **kwargs):
        return "<html><body>fake</body></html>"

    def write_html(self, path, **kwargs):
        Path(path).write_text(self.to_html(), encoding="utf-8")


def test_write_json_variants(tmp_path):
    fig = _FakeFig()
    out = write_json(fig, tmp_path / "fig.json")
    assert json.loads(out.read_text(encoding="utf-8"))["layout"]["title"] == "fake"


def test_write_dict_and_export_figure_dispatch(tmp_path):
    fig = _FakeFig()
    out = write_dict(fig, tmp_path / "fig.dict")
    assert "layout" in out.read_text(encoding="utf-8")

    via_json = export_figure(fig, ExportOptions(path=tmp_path / "a.json"))
    assert via_json.suffix == ".json"

    via_html = export_figure(fig, ExportOptions(path=tmp_path / "a.html"))
    assert via_html.exists()

    with pytest.raises(ValueError, match="extension"):
        export_figure(fig, ExportOptions(path=tmp_path / "noext"))


def test_write_html_requires_write_html_method(tmp_path):
    class _HtmlOnly:
        def to_html(self, **kwargs):  # pragma: no cover - never called
            return "<html>x</html>"

    with pytest.raises(TypeError, match="write_html"):
        write_html(_HtmlOnly(), tmp_path / "fig.html")


def test_save_png_matplotlib(tmp_path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    out = save_png(fig, tmp_path / "fig.png")
    assert out.exists() and out.stat().st_size > 0
    plt.close(fig)


def test_write_image_requires_plotly_like(tmp_path):
    with pytest.raises(TypeError, match="write_image"):
        write_image(object(), tmp_path / "x.png")


# ── station map option branches ──────────────────────────────────────────


class _Sites:
    def as_list(self):
        return [_Edi("S00"), _Edi("S01"), _Edi("S02")]


def _geo_data() -> MapData:
    return MapData(
        sites=_Sites(),
        stations=(
            StationRecord("S00", 28.00, 102.00, 10.0, "L1", 0),
            StationRecord("S01", 28.01, 102.01, 20.0, "L1", 1),
            StationRecord("S02", 28.00, 102.02, 30.0, "L2", 2),
        ),
    )


def test_station_map_skin_depth_overlay():
    from pycsamt.map.config import StationMapOptions
    from pycsamt.map.station import build_station_map

    opts = StationMapOptions(overlay="skin_depth", frequency=10.0, log_color=True)
    fig = build_station_map(_geo_data(), opts)
    assert fig.data


def test_station_map_contour_image_and_filter():
    from pycsamt.map.config import StationMapOptions
    from pycsamt.map.station import build_station_map

    opts = StationMapOptions(
        overlay="elevation",
        contour_image=True,
        show_contours=True,
        contour_grid_res=40,
        line_filter=("L1",),
        show_labels=False,
    )
    fig = build_station_map(_geo_data(), opts)
    assert fig is not None


def test_station_map_unknown_overlay_falls_back_to_index():
    from pycsamt.map.config import StationMapOptions
    from pycsamt.map.station import build_station_map

    fig = build_station_map(_geo_data(), StationMapOptions(overlay="not-a-column"))
    assert fig.data


# ── apply_elevations extra paths ─────────────────────────────────────────


def test_apply_elevations_normalized_id_fallback():
    data = MapData(
        sites=[_Edi("S00")],
        stations=(StationRecord("S00", 1.0, 2.0, None, "L1", 0),),
    )
    # keyed by a case/format variant of the id
    out = apply_elevations(data, {"s00": 77.0})
    assert out.stations[0].elevation == 77.0
