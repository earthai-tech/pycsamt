# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Unit tests for :mod:`pycsamt.gis.coord_correction`.

Covers the pure DataFrame-based ``df_*`` correction functions, the
low-level numeric helpers (``_pca_azimuth``, ``_loess_smooth``,
``_to_utm_arrays``/``_to_ll_arrays``), the dispatcher, and the
EDI/Sites-backed ``correct_*`` wrappers using minimal synthetic EDI
files (mirrors the pattern used in ``site/tests/test_site_base.py``).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pycsamt.gis.coord_correction import (
    _get_coords_df,
    _loess_smooth,
    _pca_azimuth,
    _to_ll_arrays,
    _to_utm_arrays,
    apply_coord_correction_to_df,
    apply_coords_df_to_sites,
    correct_coordinate_shift,
    correct_elevation_smooth,
    correct_interpolate_missing,
    correct_outlier_snap,
    correct_profile_projection,
    correct_spacing_regularize,
    df_elevation_smooth,
    df_interpolate_missing,
    df_outlier_snap,
    df_profile_projection,
    df_shift,
    df_spacing_regularize,
)
from pycsamt.seg.edi import EDIFile
from pycsamt.site.base import Sites

# ─────────────────────────────────────────────────────────────────────────────
# Synthetic EDI fixtures: a 5-station roughly N-S profile at lon=-118.34,
# with station S02 offset 0.005 deg east (a cross-profile "outlier") and an
# elevation spike, so profile-projection / outlier-snap / elevation-smooth
# all have something real to correct.
# ─────────────────────────────────────────────────────────────────────────────

_EDI_TEMPLATE = "\n".join(
    [
        ">HEAD",
        "  DATAID={dataid}",
        "  LAT={lat}",
        "  LONG={lon}",
        "  ELEV={elev}",
        "  STDVERS=SEG 1.0",
        "",
        ">INFO",
        "  PROJECT=SIM",
        "  PROCESSEDBY=pyCSAMT",
        "  PROCESSINGSOFTWARE=pyCSAMT",
        "",
        ">=MTSECT",
        "  SECTID={dataid}",
        "  NFREQ=2",
        "",
        ">!****FREQUENCIES****!",
        ">FREQ  //2",
        "  1.000000E+02  2.000000E+02",
        "",
        ">ZXXR ROT=ZROT  //2",
        "  1.000000E+00  1.000000E+00",
        "",
        ">END",
    ]
)

_STATIONS = [
    ("S00", 34.00, -118.340, 490.0),
    ("S01", 34.01, -118.340, 500.0),
    ("S02", 34.02, -118.335, 700.0),  # cross-profile outlier + elevation spike
    ("S03", 34.03, -118.340, 510.0),
    ("S04", 34.04, -118.340, 520.0),
]


def _make_edi(tmp_path: Path, stem: str, lat: float, lon: float, elev: float):
    p = tmp_path / f"{stem}.edi"
    p.write_text(
        _EDI_TEMPLATE.format(dataid=stem, lat=lat, lon=lon, elev=elev),
        encoding="utf-8",
    )
    return EDIFile(p)


@pytest.fixture()
def profile_sites(tmp_path: Path) -> Sites:
    edis = [_make_edi(tmp_path, *row) for row in _STATIONS]
    return Sites(edis)


@pytest.fixture()
def profile_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "station": [row[0] for row in _STATIONS],
            "lat": [row[1] for row in _STATIONS],
            "lon": [row[2] for row in _STATIONS],
            "elev": [row[3] for row in _STATIONS],
        }
    )


# ------------------------------ pure helpers -----------------------------


def test_pca_azimuth_north_south_line():
    east = np.zeros(5)
    north = np.array([0.0, 100.0, 200.0, 300.0, 400.0])
    az = _pca_azimuth(east, north)
    assert az == pytest.approx(0.0, abs=1e-6)


def test_pca_azimuth_east_west_line():
    east = np.array([0.0, 100.0, 200.0, 300.0, 400.0])
    north = np.zeros(5)
    az = _pca_azimuth(east, north)
    assert az == pytest.approx(90.0, abs=1e-6)


def test_loess_smooth_reduces_noise_variance():
    rng = np.random.default_rng(0)
    x = np.arange(50, dtype=float)
    trend = 0.5 * x
    noisy = trend + rng.normal(0, 3.0, size=x.size)
    smoothed = _loess_smooth(x, noisy, span=5)
    # Smoothed series should track the trend more closely than raw noise.
    assert np.std(smoothed - trend) < np.std(noisy - trend)


def test_loess_smooth_constant_series_is_unchanged():
    x = np.arange(10, dtype=float)
    y = np.full(10, 5.0)
    out = _loess_smooth(x, y, span=3)
    assert out == pytest.approx(y)


def test_to_utm_and_to_ll_arrays_roundtrip():
    lat = np.array([34.00, 34.01, 34.02])
    lon = np.array([-118.34, -118.34, -118.335])
    east, north, zone = _to_utm_arrays(lat, lon)
    assert zone == "11S"
    lat2, lon2 = _to_ll_arrays(east, north, zone)
    assert lat2 == pytest.approx(lat, abs=1e-6)
    assert lon2 == pytest.approx(lon, abs=1e-6)


def test_to_utm_arrays_fallback_on_exception(monkeypatch):
    import pycsamt.gis.coord_correction as cc

    def _boom(*a, **k):
        raise RuntimeError("no gis backend")

    monkeypatch.setattr(cc, "_to_utm_arrays", cc._to_utm_arrays)
    import pycsamt.gis.utils as gu

    monkeypatch.setattr(gu, "to_utm", _boom)
    lat = np.array([34.00, 34.01, 34.02])
    lon = np.array([-118.34, -118.34, -118.335])
    east, north, zone = cc._to_utm_arrays(lat, lon)
    assert zone is None
    # Flat-earth approximation: near-constant longitude -> small easting
    # spread relative to the northing spread.
    assert np.ptp(north) > np.ptp(east)


def test_to_ll_arrays_fallback_on_exception(monkeypatch):
    import pycsamt.gis.utils as gu

    def _boom(*a, **k):
        raise RuntimeError("no gis backend")

    monkeypatch.setattr(gu, "to_ll", _boom)
    east = np.array([1.0, 2.0, 3.0])
    north = np.array([4.0, 5.0, 6.0])
    lat, lon = _to_ll_arrays(east, north, "11S")
    # Fallback returns the inputs unchanged (caller's responsibility).
    assert lat is east
    assert lon is north


# ------------------------------- df_* functions ---------------------------


def test_df_profile_projection_snaps_outlier_onto_line(profile_df):
    out = df_profile_projection(profile_df)
    # All stations converge onto (nearly) the same longitude.
    assert out["lon"].std() < profile_df["lon"].std()
    assert out["lon"].nunique() <= 2  # essentially collinear now


def test_df_profile_projection_too_few_points_returns_copy():
    df = pd.DataFrame({"station": ["A"], "lat": [1.0], "lon": [2.0], "elev": [0.0]})
    out = df_profile_projection(df)
    pd.testing.assert_frame_equal(out, df)
    assert out is not df


def test_df_profile_projection_explicit_azimuth(profile_df):
    out = df_profile_projection(profile_df, azimuth=0.0)
    assert out["lon"].std() < profile_df["lon"].std()


def test_df_spacing_regularize_preserve_extent(profile_df):
    out = df_spacing_regularize(profile_df, spacing_m=200.0, preserve_extent=True)
    assert len(out) == len(profile_df)
    # First/last station chainage endpoints preserved (lat extent roughly
    # matches original), only interior spacing changes.
    assert out["lat"].min() == pytest.approx(profile_df["lat"].min(), abs=1e-3)
    assert out["lat"].max() == pytest.approx(profile_df["lat"].max(), abs=1e-3)


def test_df_spacing_regularize_strict_spacing(profile_df):
    out = df_spacing_regularize(profile_df, spacing_m=50.0, preserve_extent=False)
    assert len(out) == len(profile_df)


def test_df_spacing_regularize_too_few_points_returns_copy():
    df = pd.DataFrame({"station": ["A"], "lat": [1.0], "lon": [2.0], "elev": [0.0]})
    out = df_spacing_regularize(df)
    pd.testing.assert_frame_equal(out, df)


def test_df_outlier_snap_moves_only_outliers(profile_df):
    out = df_outlier_snap(profile_df, threshold_m=200.0)
    # S02 (index 2) was the deliberate outlier; its lon should change.
    assert out.loc[2, "lon"] != pytest.approx(profile_df.loc[2, "lon"])
    # In-line stations remain effectively untouched.
    assert out.loc[0, "lon"] == pytest.approx(profile_df.loc[0, "lon"], abs=1e-9)


def test_df_outlier_snap_high_threshold_no_change(profile_df):
    out = df_outlier_snap(profile_df, threshold_m=10_000.0)
    pd.testing.assert_frame_equal(out[["lat", "lon"]], profile_df[["lat", "lon"]])


def test_df_outlier_snap_too_few_points_returns_copy():
    df = pd.DataFrame({"station": ["A"], "lat": [1.0], "lon": [2.0], "elev": [0.0]})
    out = df_outlier_snap(df)
    pd.testing.assert_frame_equal(out, df)


def test_df_elevation_smooth_loess(profile_df):
    out = df_elevation_smooth(profile_df, method="loess", window=2)
    assert out.loc[2, "elev"] < profile_df.loc[2, "elev"]  # spike reduced


def test_df_elevation_smooth_mean(profile_df):
    out = df_elevation_smooth(profile_df, method="mean", window=3)
    assert len(out) == len(profile_df)
    assert "elev" in out


def test_df_elevation_smooth_too_few_points_returns_copy():
    df = pd.DataFrame(
        {
            "station": ["A", "B"],
            "lat": [1.0, 2.0],
            "lon": [3.0, 4.0],
            "elev": [0.0, 1.0],
        }
    )
    out = df_elevation_smooth(df)
    pd.testing.assert_frame_equal(out, df)


def test_df_shift(profile_df):
    out = df_shift(profile_df, delta_lat=0.5, delta_lon=-0.25, delta_elev=10.0)
    assert out["lat"].to_numpy() == pytest.approx(profile_df["lat"].to_numpy() + 0.5)
    assert out["lon"].to_numpy() == pytest.approx(profile_df["lon"].to_numpy() - 0.25)
    assert out["elev"].to_numpy() == pytest.approx(profile_df["elev"].to_numpy() + 10.0)


def test_df_interpolate_missing_nan_only():
    df = pd.DataFrame(
        {
            "station": ["A", "B", "C"],
            "lat": [1.0, np.nan, 3.0],
            "lon": [1.0, np.nan, 3.0],
            "elev": [10.0, np.nan, 30.0],
        }
    )
    out = df_interpolate_missing(df, fill_nan_only=True)
    assert out.loc[1, "lat"] == pytest.approx(2.0)
    assert out.loc[1, "elev"] == pytest.approx(20.0)


def test_df_interpolate_missing_treats_zero_as_missing():
    df = pd.DataFrame(
        {
            "station": ["A", "B", "C"],
            "lat": [1.0, 0.0, 3.0],
            "lon": [1.0, 0.0, 3.0],
            "elev": [10.0, 0.0, 30.0],
        }
    )
    out_default = df_interpolate_missing(df, fill_nan_only=True)
    assert out_default.loc[1, "lat"] == 0.0  # untouched by default

    out_zero_fill = df_interpolate_missing(df, fill_nan_only=False)
    assert out_zero_fill.loc[1, "lat"] == pytest.approx(2.0)


# ------------------------------- dispatcher --------------------------------


def test_apply_coord_correction_to_df_dispatch(profile_df):
    out = apply_coord_correction_to_df(
        "_coord_shift",
        profile_df,
        delta_lat=1.0,
        delta_lon=0.0,
        delta_elev=0.0,
    )
    assert out["lat"].to_numpy() == pytest.approx(profile_df["lat"].to_numpy() + 1.0)


def test_apply_coord_correction_to_df_unknown_name_returns_copy(profile_df):
    out = apply_coord_correction_to_df("_not_a_real_fn", profile_df)
    pd.testing.assert_frame_equal(out, profile_df)
    assert out is not profile_df


@pytest.mark.parametrize(
    "fn_name",
    [
        "_coord_profile_projection",
        "_coord_spacing_regularize",
        "_coord_outlier_snap",
        "_coord_elevation_smooth",
        "_coord_shift",
        "_coord_interpolate_missing",
    ],
)
def test_apply_coord_correction_to_df_all_registered_names(fn_name, profile_df):
    out = apply_coord_correction_to_df(fn_name, profile_df)
    assert isinstance(out, pd.DataFrame)
    assert len(out) == len(profile_df)


# ------------------------------- Sites wrappers ----------------------------


def test_correct_coordinate_shift_sites(profile_sites):
    shifted = correct_coordinate_shift(
        profile_sites, delta_lat=0.001, delta_lon=0.002, delta_elev=5.0
    )
    before = _get_coords_df(profile_sites)
    after = _get_coords_df(shifted)
    assert after["lat"].to_numpy() == pytest.approx(before["lat"].to_numpy() + 0.001)
    assert after["lon"].to_numpy() == pytest.approx(before["lon"].to_numpy() + 0.002)
    assert after["elev"].to_numpy() == pytest.approx(before["elev"].to_numpy() + 5.0)
    # inplace=False (default): original Sites is untouched.
    assert _get_coords_df(profile_sites)["lat"].to_numpy() == pytest.approx(
        before["lat"].to_numpy()
    )


def test_correct_coordinate_shift_inplace(profile_sites):
    before = _get_coords_df(profile_sites)
    result = correct_coordinate_shift(profile_sites, delta_elev=1.0, inplace=True)
    after = _get_coords_df(result)
    assert after["elev"].to_numpy() == pytest.approx(before["elev"].to_numpy() + 1.0)
    # Same object mutated.
    assert _get_coords_df(profile_sites)["elev"].to_numpy() == pytest.approx(
        after["elev"].to_numpy()
    )


def test_correct_profile_projection_sites(profile_sites):
    before = _get_coords_df(profile_sites)
    proj = correct_profile_projection(profile_sites)
    after = _get_coords_df(proj)
    assert after["lon"].std() < before["lon"].std()
    # Original left untouched.
    assert _get_coords_df(profile_sites)["lon"].to_numpy() == pytest.approx(
        before["lon"].to_numpy()
    )


def test_correct_spacing_regularize_sites(profile_sites):
    before = _get_coords_df(profile_sites)
    reg = correct_spacing_regularize(profile_sites, spacing_m=300.0)
    after = _get_coords_df(reg)
    assert len(after) == len(before)


def test_correct_outlier_snap_sites(profile_sites):
    before = _get_coords_df(profile_sites)
    snapped = correct_outlier_snap(profile_sites, threshold_m=200.0)
    after = _get_coords_df(snapped)
    outlier_row = before.index[before["station"] == "S02"][0]
    assert after.loc[outlier_row, "lon"] != pytest.approx(
        before.loc[outlier_row, "lon"]
    )


def test_correct_outlier_snap_high_threshold_returns_unchanged_sites(
    profile_sites,
):
    before = _get_coords_df(profile_sites)
    result = correct_outlier_snap(profile_sites, threshold_m=10_000.0)
    after = _get_coords_df(result)
    assert after["lon"].to_numpy() == pytest.approx(before["lon"].to_numpy())


def test_correct_elevation_smooth_sites(profile_sites):
    before = _get_coords_df(profile_sites)
    smoothed = correct_elevation_smooth(profile_sites, method="loess", window=2)
    after = _get_coords_df(smoothed)
    outlier_row = before.index[before["station"] == "S02"][0]
    assert after.loc[outlier_row, "elev"] < before.loc[outlier_row, "elev"]


def test_correct_interpolate_missing_sites(tmp_path):
    stations = [
        ("S00", 34.00, -118.34, 490.0),
        ("S01", 34.01, -118.34, 500.0),
        ("S02", 0.0, 0.0, 0.0),  # missing GPS fix
        ("S03", 34.03, -118.34, 510.0),
        ("S04", 34.04, -118.34, 520.0),
    ]
    edis = [_make_edi(tmp_path, *row) for row in stations]
    sites = Sites(edis)
    filled = correct_interpolate_missing(sites, fill_nan_only=False)
    after = _get_coords_df(filled)
    row = after.index[after["station"] == "S02"][0]
    assert after.loc[row, "lat"] == pytest.approx(34.02, abs=1e-6)
    assert after.loc[row, "elev"] == pytest.approx(505.0, abs=1e-6)


def test_apply_coords_df_to_sites_mutates_in_place(profile_sites):
    coords_df = _get_coords_df(profile_sites)
    coords_df["lat"] = coords_df["lat"] + 2.0
    apply_coords_df_to_sites(profile_sites, coords_df)
    after = _get_coords_df(profile_sites)
    assert after["lat"].to_numpy() == pytest.approx(coords_df["lat"].to_numpy())


def test_apply_coords_df_to_sites_skips_non_finite(profile_sites):
    before = _get_coords_df(profile_sites)
    coords_df = before.copy()
    coords_df.loc[0, "lat"] = np.nan
    coords_df["lon"] = coords_df["lon"] + 1.0
    apply_coords_df_to_sites(profile_sites, coords_df)
    after = _get_coords_df(profile_sites)
    # Row 0 had a non-finite lat -> left untouched.
    assert after.loc[0, "lon"] == pytest.approx(before.loc[0, "lon"])
    # Other rows updated.
    assert after.loc[1, "lon"] == pytest.approx(before.loc[1, "lon"] + 1.0)


def test_apply_coords_df_to_sites_skips_unknown_station(profile_sites):
    before = _get_coords_df(profile_sites)
    # coords_df only knows about a station absent from the Sites container.
    coords_df = pd.DataFrame(
        {"station": ["NOT_A_SITE"], "lat": [1.0], "lon": [2.0], "elev": [3.0]}
    )
    apply_coords_df_to_sites(profile_sites, coords_df)
    after = _get_coords_df(profile_sites)
    assert after["lat"].to_numpy() == pytest.approx(before["lat"].to_numpy())


# --------------------- single/degenerate-station edge cases ---------------


@pytest.fixture()
def single_site(tmp_path: Path) -> Sites:
    edi = _make_edi(tmp_path, "SOLO", 34.0, -118.0, 100.0)
    return Sites([edi])


def test_correct_profile_projection_single_site_is_noop(single_site):
    before = _get_coords_df(single_site)
    out = correct_profile_projection(single_site)
    after = _get_coords_df(out)
    assert after["lon"].to_numpy() == pytest.approx(before["lon"].to_numpy())


def test_correct_spacing_regularize_single_site_is_noop(single_site):
    before = _get_coords_df(single_site)
    out = correct_spacing_regularize(single_site)
    after = _get_coords_df(out)
    assert after["lat"].to_numpy() == pytest.approx(before["lat"].to_numpy())


def test_correct_outlier_snap_single_site_is_noop(single_site):
    before = _get_coords_df(single_site)
    out = correct_outlier_snap(single_site)
    after = _get_coords_df(out)
    assert after["lat"].to_numpy() == pytest.approx(before["lat"].to_numpy())


def test_correct_elevation_smooth_too_few_sites_is_noop(tmp_path):
    edis = [
        _make_edi(tmp_path, "A", 34.0, -118.0, 100.0),
        _make_edi(tmp_path, "B", 34.01, -118.0, 110.0),
    ]
    sites = Sites(edis)
    before = _get_coords_df(sites)
    out = correct_elevation_smooth(sites)
    after = _get_coords_df(out)
    assert after["elev"].to_numpy() == pytest.approx(before["elev"].to_numpy())


def test_correct_coordinate_shift_all_missing_coords_is_noop(monkeypatch, tmp_path):
    # correct_coordinate_shift builds its coords dict by filtering out
    # rows with nan lat/lon; when every row is nan the dict is empty and
    # the function must return the input sites unchanged. Force that via
    # _get_coords_df since a real EDI header cannot easily produce nan
    # lat/lon (LAT/LONG default to 0.0 rather than missing).
    import pycsamt.gis.coord_correction as cc

    edi = _make_edi(tmp_path, "A", 0.0, 0.0, 0.0)
    sites = Sites([edi])
    all_nan_df = pd.DataFrame(
        {"station": ["A"], "lat": [np.nan], "lon": [np.nan], "elev": [0.0]}
    )
    monkeypatch.setattr(cc, "_get_coords_df", lambda _sites: all_nan_df)

    result = cc.correct_coordinate_shift(sites, delta_lat=1.0)
    assert result is sites


def test_correct_spacing_regularize_strict_spacing_sites(profile_sites):
    before = _get_coords_df(profile_sites)
    out = correct_spacing_regularize(
        profile_sites, spacing_m=50.0, preserve_extent=False
    )
    after = _get_coords_df(out)
    assert len(after) == len(before)


def test_correct_elevation_smooth_mean_method_sites(profile_sites):
    out = correct_elevation_smooth(profile_sites, method="mean", window=3)
    after = _get_coords_df(out)
    assert len(after) == 5


def test_correct_spacing_regularize_pads_when_arange_too_short_sites(
    profile_sites,
):
    out = correct_spacing_regularize(
        profile_sites, spacing_m=10_000.0, preserve_extent=False
    )
    after = _get_coords_df(out)
    assert len(after) == 5


def test_correct_elevation_smooth_mean_scipy_import_error_fallback_sites(
    monkeypatch, profile_sites
):
    import sys

    monkeypatch.setitem(sys.modules, "scipy.ndimage", None)
    out = correct_elevation_smooth(profile_sites, method="mean", window=3)
    after = _get_coords_df(out)
    assert len(after) == 5
    assert after["elev"].notna().all()


def test_df_spacing_regularize_pads_when_arange_too_short(profile_df):
    # A spacing far larger than the profile extent leaves fewer than n
    # points from np.arange(...), exercising the linspace padding fallback.
    out = df_spacing_regularize(profile_df, spacing_m=10_000.0, preserve_extent=False)
    assert len(out) == len(profile_df)


def test_df_elevation_smooth_mean_scipy_import_error_fallback(monkeypatch, profile_df):
    import sys

    monkeypatch.setitem(sys.modules, "scipy.ndimage", None)
    out = df_elevation_smooth(profile_df, method="mean", window=3)
    assert len(out) == len(profile_df)
    assert out["elev"].notna().all()
