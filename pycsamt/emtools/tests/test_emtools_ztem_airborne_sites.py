"""Tests for pycsamt.emtools.ztem accepting AirborneSites (step 2).

Exercises the same code path that lets ztem.py's functions read
genuine tipper-only ZTEM EMTF-XML -- no impedance channel at all --
which pycsamt.site.base.Site/Sites cannot represent. See
pycsamt.airborne.site's module docstring for why, and
pycsamt.emtools._core.ensure_any_sites for the dispatcher under test
indirectly here.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.airborne import AirborneSites, NavigationTrack
from pycsamt.airborne.site import ensure_asites
from pycsamt.airborne.ztem import build_ztem_line
from pycsamt.emtools.ztem import (
    mask_outside_ztem_band,
    phase_rotate_table,
    plot_ztem_band_mask_psection,
    plot_ztem_divergence_profile,
    plot_ztem_divergence_psection,
    plot_ztem_divergence_psection_grid,
    plot_ztem_flight_lines,
    plot_ztem_map,
    plot_ztem_phase_rotation_profile,
    plot_ztem_tipper_profile,
    total_divergence_table,
    ztem_crossover_diagnostics,
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

_ZTEM_FREQS = np.array([30.0, 45.0, 90.0, 180.0, 360.0, 720.0])
_DATA_ROOT = Path(__file__).resolve().parents[3] / "data"
_GOLD_SPRINGS = _DATA_ROOT / "ZTEM" / "gold_springs_nv"


def _crossover_tipper(x: float, freqs: np.ndarray) -> np.ndarray:
    t = np.zeros((freqs.size, 2), dtype=complex)
    shape = (x / 150.0) * np.exp(-((x / 200.0) ** 2))
    for k, f in enumerate(freqs):
        decay = 1.0 / np.sqrt(f / 30.0)
        t[k, 0] = shape * decay + 0.05j * decay
        t[k, 1] = 0.3 * shape * decay - 0.02j * decay
    return t


def _airborne_line(n_stations: int = 9, spacing: float = 100.0):
    x0 = -0.5 * spacing * (n_stations - 1)
    sample_ids = tuple(f"S{i:02d}" for i in range(n_stations))
    east = np.array([x0 + i * spacing for i in range(n_stations)])
    nav = NavigationTrack(
        sample_ids=sample_ids, easting=east, northing=np.zeros(n_stations),
    )
    tipper = np.stack(
        [_crossover_tipper(x, _ZTEM_FREQS) for x in east], axis=0
    )
    return build_ztem_line(
        "L001", nav, tipper, frequency=_ZTEM_FREQS,
    )


def _airborne_sites() -> AirborneSites:
    return AirborneSites.from_line(_airborne_line(), technology="ztem")


def _multiline_airborne_sites(
    n_lines: int = 3, n_stations: int = 7, spacing_deg: float = 0.001,
    line_spacing_deg: float = 0.01,
) -> AirborneSites:
    """A small multi-line synthetic block, with real lat/lon coords
    (unlike ``_airborne_line``'s bare easting/northing), for the
    map-view functions that read ``AirborneSite.coords`` directly."""
    items = []
    lon0 = -114.0
    for li in range(n_lines):
        lat = 38.0 + li * line_spacing_deg
        x0 = -0.5 * spacing_deg * (n_stations - 1)
        lon = np.array([lon0 + x0 + i * spacing_deg for i in range(n_stations)])
        sample_ids = tuple(f"L{li}_{i:02d}" for i in range(n_stations))
        nav = NavigationTrack(
            sample_ids=sample_ids,
            latitude=np.full(n_stations, lat),
            longitude=lon,
        )
        x_m = (lon - lon.mean()) * 111_000.0
        tipper = np.stack(
            [_crossover_tipper(x, _ZTEM_FREQS) for x in x_m], axis=0
        )
        line = build_ztem_line(
            f"L{li}", nav, tipper, frequency=_ZTEM_FREQS,
        )
        items.extend(
            AirborneSites.from_line(line, technology="ztem").as_list()
        )
    return AirborneSites(items)


# ─────────────────────────────────────────────────────────────────────────
# total_divergence_table / phase_rotate_table (read-only tables)
# ─────────────────────────────────────────────────────────────────────────


def test_total_divergence_table_accepts_airborne_sites():
    df = total_divergence_table(_airborne_sites())
    assert not df.empty
    assert df["divergence_real"].abs().max() > 0.0


def test_phase_rotate_table_accepts_airborne_sites():
    df = phase_rotate_table(
        _airborne_sites(), frequency_hz=30.0, n_resample=51,
    )
    assert not df.empty


def test_total_divergence_table_accepts_bare_directory_path(tmp_path):
    _airborne_sites().write_xml(tmp_path)
    df = total_divergence_table(str(tmp_path))
    assert not df.empty


# ─────────────────────────────────────────────────────────────────────────
# mask_outside_ztem_band (the pipeline function)
# ─────────────────────────────────────────────────────────────────────────


def test_mask_outside_ztem_band_returns_airborne_sites():
    out = mask_outside_ztem_band(
        _airborne_sites(), band_hz=(40.0, 200.0),
    )
    assert isinstance(out, AirborneSites)
    df = total_divergence_table(out)
    oob = df[(df["freq_hz"] < 40.0) | (df["freq_hz"] > 200.0)]
    assert oob.empty


def test_mask_outside_ztem_band_drop_rejected_for_airborne_sites():
    with pytest.raises(ValueError):
        mask_outside_ztem_band(
            _airborne_sites(), band_hz=(40.0, 200.0), action="drop",
        )


def test_mask_outside_ztem_band_inplace_true_mutates():
    asites = _airborne_sites()
    before = asites[0].tipper.copy()
    out = mask_outside_ztem_band(
        asites, band_hz=(40.0, 200.0), inplace=True,
    )
    assert out is asites
    assert not np.array_equal(before, asites[0].tipper)


def test_mask_outside_ztem_band_inplace_false_preserves_input():
    asites = _airborne_sites()
    snapshot = asites[0].tipper.copy()
    out = mask_outside_ztem_band(
        asites, band_hz=(40.0, 200.0), inplace=False,
    )
    assert out is not asites
    assert np.array_equal(snapshot, asites[0].tipper, equal_nan=True)


# ─────────────────────────────────────────────────────────────────────────
# Plots
# ─────────────────────────────────────────────────────────────────────────


def test_plot_ztem_divergence_profile_accepts_airborne_sites():
    ax = plot_ztem_divergence_profile(_airborne_sites())
    assert isinstance(ax, plt.Axes)
    plt.close("all")


def test_plot_ztem_divergence_psection_accepts_airborne_sites():
    ax = plot_ztem_divergence_psection(_airborne_sites())
    assert isinstance(ax, plt.Axes)
    plt.close("all")


def test_plot_ztem_phase_rotation_profile_accepts_airborne_sites():
    ax = plot_ztem_phase_rotation_profile(
        _airborne_sites(), frequency_hz=30.0,
    )
    assert isinstance(ax, plt.Axes)
    plt.close("all")


def test_plot_ztem_band_mask_psection_accepts_airborne_sites():
    fig = plot_ztem_band_mask_psection(
        _airborne_sites(), band_hz=(40.0, 200.0),
    )
    assert isinstance(fig, plt.Figure)
    assert len(fig.axes) >= 2  # 2 panels + their colorbars
    plt.close("all")


# ─────────────────────────────────────────────────────────────────────────
# Map view: plot_ztem_flight_lines / plot_ztem_map
# ─────────────────────────────────────────────────────────────────────────


def test_plot_ztem_flight_lines_detects_all_lines():
    ax = plot_ztem_flight_lines(_multiline_airborne_sites(n_lines=3))
    assert isinstance(ax, plt.Axes)
    assert "3 lines" in ax.get_title()
    plt.close("all")


def test_plot_ztem_flight_lines_no_coords():
    ax = plot_ztem_flight_lines(_airborne_sites())  # no lat/lon set
    assert isinstance(ax, plt.Axes)
    assert "no station coordinates" in ax.texts[0].get_text()
    plt.close("all")


def test_plot_ztem_map_tipper_quantity():
    ax = plot_ztem_map(
        _multiline_airborne_sites(n_lines=3), quantity="tipper",
        frequency_hz=30.0,
    )
    assert isinstance(ax, plt.Axes)
    plt.close("all")


def test_plot_ztem_map_divergence_quantity():
    ax = plot_ztem_map(
        _multiline_airborne_sites(n_lines=3), quantity="divergence",
        frequency_hz=30.0,
    )
    assert isinstance(ax, plt.Axes)
    plt.close("all")


def test_plot_ztem_map_invalid_quantity_raises():
    with pytest.raises(ValueError):
        plot_ztem_map(_multiline_airborne_sites(), quantity="bogus")


def test_plot_ztem_map_accepts_bare_directory_path(tmp_path):
    _multiline_airborne_sites(n_lines=3).write_xml(tmp_path)
    ax = plot_ztem_map(str(tmp_path), quantity="tipper", frequency_hz=30.0)
    assert isinstance(ax, plt.Axes)
    plt.close("all")


def test_plot_ztem_divergence_psection_grid_auto_creates_axes():
    fig = plot_ztem_divergence_psection_grid(
        _multiline_airborne_sites(n_lines=5), max_lines=4, n_cols=2,
    )
    assert isinstance(fig, plt.Figure)
    plt.close("all")


def test_plot_ztem_divergence_psection_grid_accepts_existing_axes():
    fig, axes = plt.subplots(2, 2)
    out = plot_ztem_divergence_psection_grid(
        _multiline_airborne_sites(n_lines=5), max_lines=4, n_cols=2,
        axes=axes,
    )
    assert out is fig
    plt.close("all")


def test_plot_ztem_divergence_psection_grid_too_few_axes_raises():
    fig, axes = plt.subplots(1, 2)
    with pytest.raises(ValueError):
        plot_ztem_divergence_psection_grid(
            _multiline_airborne_sites(n_lines=5), max_lines=4, n_cols=2,
            axes=axes,
        )
    plt.close("all")


# ─────────────────────────────────────────────────────────────────────────
# Real synthetic data (skipped if not present locally)
# ─────────────────────────────────────────────────────────────────────────


@pytest.mark.skipif(
    not _GOLD_SPRINGS.exists(), reason="synthetic ZTEM sample data not found"
)
def test_real_synthetic_gold_springs_end_to_end():
    # gold_springs_nv is a real 7-line, 105-station block survey;
    # total_divergence_table/plot_ztem_divergence_profile assume a
    # single flight line (see their own docstring warnings), so
    # select just one line before calling them.
    all_sites = ensure_asites(str(_GOLD_SPRINGS))
    assert len(all_sites) == 105
    line1 = all_sites.select(
        predicate=lambda s: (
            s.emtf.metadata.get("notes", {}).get("ZTEM", {}).get("LineId")
            == "gold_springs_nv_L1"
        )
    )
    assert len(line1) == 15

    df = total_divergence_table(line1)
    assert not df.empty
    masked = mask_outside_ztem_band(
        str(_GOLD_SPRINGS), band_hz=(40.0, 200.0),
    )
    assert isinstance(masked, AirborneSites)
    assert len(masked) == 105
    ax = plot_ztem_divergence_profile(line1)
    assert isinstance(ax, plt.Axes)
    plt.close("all")


@pytest.mark.skipif(
    not _GOLD_SPRINGS.exists(), reason="synthetic ZTEM sample data not found"
)
def test_real_synthetic_gold_springs_flight_lines_and_map():
    ax1 = plot_ztem_flight_lines(str(_GOLD_SPRINGS))
    assert isinstance(ax1, plt.Axes)
    assert "7 lines" in ax1.get_title()
    plt.close("all")

    ax2 = plot_ztem_map(
        str(_GOLD_SPRINGS), quantity="tipper", frequency_hz=90.0,
    )
    assert isinstance(ax2, plt.Axes)
    plt.close("all")

    ax3 = plot_ztem_map(
        str(_GOLD_SPRINGS), quantity="divergence", frequency_hz=90.0,
    )
    assert isinstance(ax3, plt.Axes)
    plt.close("all")


@pytest.mark.skipif(
    not _GOLD_SPRINGS.exists(), reason="synthetic ZTEM sample data not found"
)
def test_real_synthetic_gold_springs_crossover_diagnostics():
    all_sites = ensure_asites(str(_GOLD_SPRINGS))
    line1 = all_sites.select(
        predicate=lambda s: (
            s.emtf.metadata.get("notes", {}).get("ZTEM", {}).get("LineId")
            == "gold_springs_nv_L1"
        )
    )
    diag = ztem_crossover_diagnostics(line1, frequency_hz=90.0)
    assert diag["peak_to_peak_real"] > 0.0
    ax = plot_ztem_tipper_profile(line1, frequency_hz=90.0)
    assert isinstance(ax, plt.Axes)
    plt.close("all")
