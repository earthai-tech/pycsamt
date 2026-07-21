# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.topo.overlay — terrain rendering helpers.

All matplotlib operations run through the Agg non-interactive backend so
no display is required.
"""

from __future__ import annotations

import glob
import os

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.topo.config import TopoConfig, reset_topo
from pycsamt.topo.overlay import (
    add_station_labels,
    draw_topo_section,
    draw_topo_strip,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _reset_topo():
    yield
    reset_topo()
    plt.close("all")


def _simple_inputs(n: int = 8):
    """Synthetic profile with n stations."""
    chain = np.linspace(0.0, 5.0, n)  # km
    elev = 200.0 + 50.0 * np.sin(np.pi * chain / 5.0)  # m, smooth hill
    names = [f"S{i:02d}" for i in range(n)]
    return chain, elev, names


def _make_section_axes():
    fig, ax = plt.subplots()
    ax.set_xlim(0, 5)
    ax.set_ylim(0.0, 0.5)  # km scale — y is absolute elevation
    return fig, ax


def _make_ps_axes():
    """Return (fig, ax) for a pseudosection-style plot."""
    fig, ax = plt.subplots()
    ax.imshow(
        np.random.default_rng(0).random((20, 8)),
        aspect="auto",
        extent=[0, 8, 0, 1],
    )
    return fig, ax


# ---------------------------------------------------------------------------
# draw_topo_section
# ---------------------------------------------------------------------------


class TestDrawTopoSection:
    def test_returns_none(self):
        chain, elev, _ = _simple_inputs()
        _, ax = _make_section_axes()
        result = draw_topo_section(ax, chain, elev)
        assert result is None

    def test_adds_fill_patch(self):
        chain, elev, _ = _simple_inputs()
        _, ax = _make_section_axes()
        n_before = len(ax.collections) + len(ax.patches) + len(ax.lines)
        draw_topo_section(ax, chain, elev)
        n_after = len(ax.collections) + len(ax.patches) + len(ax.lines)
        assert n_after > n_before

    def test_empty_chain_does_not_raise(self):
        _, ax = _make_section_axes()
        draw_topo_section(ax, np.array([]), np.array([]))

    def test_single_station_does_not_raise(self):
        _, ax = _make_section_axes()
        draw_topo_section(ax, np.array([2.0]), np.array([150.0]))

    def test_with_station_names(self):
        chain, elev, names = _simple_inputs()
        _, ax = _make_section_axes()
        draw_topo_section(ax, chain, elev, station_names=names)
        # station labels are Text artists
        texts = [t for t in ax.texts]
        assert len(texts) == len(names)

    def test_without_station_names_no_text(self):
        chain, elev, _ = _simple_inputs()
        _, ax = _make_section_axes()
        draw_topo_section(ax, chain, elev, station_names=None)
        assert len(ax.texts) == 0

    def test_dark_mode_true(self):
        chain, elev, names = _simple_inputs()
        _, ax = _make_section_axes()
        draw_topo_section(ax, chain, elev, names, dark=True)

    def test_light_mode(self):
        chain, elev, names = _simple_inputs()
        _, ax = _make_section_axes()
        draw_topo_section(ax, chain, elev, names, dark=False)

    def test_custom_cfg_fill_color(self):
        chain, elev, _ = _simple_inputs()
        _, ax = _make_section_axes()
        cfg = TopoConfig()
        cfg.configure(fill_color="#ff0000", fill_alpha=0.9)
        draw_topo_section(ax, chain, elev, cfg=cfg)

    def test_no_fill_when_clip_off(self):
        chain, elev, _ = _simple_inputs()
        _, ax = _make_section_axes()
        n_patches_before = len(ax.patches)
        cfg = TopoConfig()
        cfg.configure(clip_below_surface=False)
        draw_topo_section(ax, chain, elev, cfg=cfg)
        # Terrain fill uses ax.fill() → ax.patches; no patch added when disabled.
        # Station markers use ax.scatter() → ax.collections (always added).
        assert len(ax.patches) == n_patches_before

    def test_no_surface_line_when_disabled(self):
        """With show_surface_line=False, no ax.lines are added; pins use scatter."""
        chain, elev, _ = _simple_inputs()
        _, ax = _make_section_axes()
        lines_before = len(ax.lines)
        collections_before = len(ax.collections)
        cfg = TopoConfig()
        cfg.configure(show_surface_line=False, clip_below_surface=False)
        draw_topo_section(ax, chain, elev, cfg=cfg)
        # No new lines (surface line disabled, pins now use scatter not plot)
        assert len(ax.lines) == lines_before
        # Station-pin scatter is the only new collection
        assert len(ax.collections) == collections_before + 1

    def test_exaggeration_applied(self):
        """With exaggeration=2, the surface line y-peak should be ~doubled."""
        chain = np.array([0.0, 1.0, 2.0])
        elev_m = np.array([100.0, 200.0, 100.0])  # max = 200 m = 0.2 km
        _, ax = _make_section_axes()
        cfg = TopoConfig()
        cfg.configure(
            exaggeration=2.0,
            show_surface_line=True,
            clip_below_surface=False,
            station_pins_at_surface=False,
        )
        draw_topo_section(ax, chain, elev_m, cfg=cfg)
        # surface line is lines[0]; station pins are now a scatter (collection)
        surface_line = ax.lines[0]
        assert surface_line.get_ydata().max() == pytest.approx(0.4, rel=0.01)

    def test_custom_station_x(self):
        chain, elev, names = _simple_inputs(n=4)
        _, ax = _make_section_axes()
        sx = chain + 0.1  # shifted x positions for markers
        draw_topo_section(ax, chain, elev, names, station_x_km=sx)

    def test_real_willy_l18(self):
        """Full draw with L18PLT real elevation and chainage."""
        from pycsamt.seg.edi import EDIFile
        from pycsamt.topo.extract import (
            extract_chainage,
            extract_elevation,
            extract_station_names,
        )

        base = os.path.join(
            os.path.dirname(__file__),
            "..",
            "..",
            "..",
            "data",
            "AMT",
            "WILLY_DATA",
            "L18PLT",
        )
        paths = sorted(glob.glob(os.path.join(base, "*.edi")))
        if not paths:
            pytest.skip("WILLY L18PLT data not found")
        edis = [EDIFile(p) for p in paths]
        chain = extract_chainage(edis)
        elev = extract_elevation(edis)
        names = extract_station_names(edis)

        _, ax = _make_section_axes()
        ax.set_xlim(chain[0], chain[-1])
        ax.set_ylim(0.03, 0.25)  # km scale (30–250 m)
        draw_topo_section(ax, chain, elev, names)

        fig = ax.get_figure()
        fig.canvas.draw()  # forces rendering — catches any artist error


# ---------------------------------------------------------------------------
# draw_topo_strip
# ---------------------------------------------------------------------------


class TestDrawTopoStrip:
    def test_returns_axes(self):
        chain, elev, _ = _simple_inputs()
        fig, ax = _make_ps_axes()
        ax_strip = draw_topo_strip(fig, ax, chain, elev)
        assert ax_strip is not None
        assert hasattr(ax_strip, "get_xlim")

    def test_strip_is_new_axes(self):
        chain, elev, _ = _simple_inputs()
        fig, ax = _make_ps_axes()
        n_ax_before = len(fig.axes)
        ax_strip = draw_topo_strip(fig, ax, chain, elev)
        assert len(fig.axes) == n_ax_before + 1
        assert ax_strip is not ax

    def test_main_axes_shrunk(self):
        chain, elev, _ = _simple_inputs()
        fig, ax = _make_ps_axes()
        pos_before = ax.get_position().height
        draw_topo_strip(fig, ax, chain, elev)
        pos_after = ax.get_position().height
        assert pos_after < pos_before

    def test_empty_chain_returns_none(self):
        fig, ax = _make_ps_axes()
        result = draw_topo_strip(fig, ax, np.array([]), np.array([]))
        assert result is None

    def test_with_station_names(self):
        chain, elev, names = _simple_inputs()
        fig, ax = _make_ps_axes()
        ax_strip = draw_topo_strip(fig, ax, chain, elev, station_names=names)
        texts = list(ax_strip.texts)
        assert len(texts) == len(names)

    def test_dark_mode(self):
        chain, elev, _ = _simple_inputs()
        fig, ax = _make_ps_axes()
        ax_strip = draw_topo_strip(fig, ax, chain, elev, dark=True)
        assert ax_strip.get_facecolor() is not None

    def test_light_mode(self):
        chain, elev, _ = _simple_inputs()
        fig, ax = _make_ps_axes()
        draw_topo_strip(fig, ax, chain, elev, dark=False)

    def test_strip_height_ratio_respected(self):
        chain, elev, _ = _simple_inputs()
        fig, ax = _make_ps_axes()
        cfg = TopoConfig()
        cfg.configure(strip_height_ratio=0.30)
        pos_before = ax.get_position().height
        draw_topo_strip(fig, ax, chain, elev, cfg=cfg)
        pos_after = ax.get_position().height
        # main axes should be ~70% of original
        assert pos_after == pytest.approx(pos_before * 0.70, rel=0.02)

    def test_no_strip_when_disabled(self):
        chain, elev, _ = _simple_inputs()
        fig, ax = _make_ps_axes()
        cfg = TopoConfig()
        cfg.configure(show_topo_strip=False)
        # draw_topo_strip always renders — show_topo_strip is the caller-level gate
        # but we can verify the function still returns an axes object
        ax_strip = draw_topo_strip(fig, ax, chain, elev, cfg=cfg)
        assert ax_strip is not None

    def test_real_willy_l22(self):
        """Full strip draw with L22PLT real data."""
        from pycsamt.seg.edi import EDIFile
        from pycsamt.topo.extract import (
            extract_chainage,
            extract_elevation,
        )

        base = os.path.join(
            os.path.dirname(__file__),
            "..",
            "..",
            "..",
            "data",
            "AMT",
            "WILLY_DATA",
            "L22PLT",
        )
        paths = sorted(glob.glob(os.path.join(base, "*.edi")))
        if not paths:
            pytest.skip("WILLY L22PLT data not found")
        edis = [EDIFile(p) for p in paths]
        chain = extract_chainage(edis)
        elev = extract_elevation(edis)

        fig, ax = _make_ps_axes()
        ax_strip = draw_topo_strip(fig, ax, chain, elev)
        fig.canvas.draw()  # final render check
        assert ax_strip is not None


# ---------------------------------------------------------------------------
# add_station_labels
# ---------------------------------------------------------------------------


class TestAddStationLabels:
    def test_adds_correct_number_of_texts(self):
        _, ax = _make_section_axes()
        x = np.linspace(0, 5, 6)
        y = np.ones(6) * 0.2
        names = [f"S{i}" for i in range(6)]
        add_station_labels(ax, x, y, names, color="#ffffff")
        assert len(ax.texts) == 6

    def test_no_labels_for_empty_names(self):
        _, ax = _make_section_axes()
        add_station_labels(
            ax, np.array([]), np.array([]), [], color="#ffffff"
        )
        assert len(ax.texts) == 0

    def test_label_positions(self):
        _, ax = _make_section_axes()
        x = np.array([0.0, 1.0, 2.0])
        y = np.array([0.1, 0.2, 0.15])
        names = ["A", "B", "C"]
        add_station_labels(ax, x, y, names, color="#aaaaaa", offset_km=0.01)
        texts = ax.texts
        xs = [t.get_position()[0] for t in texts]
        np.testing.assert_allclose(xs, x, rtol=1e-9)

    def test_rotation_applied(self):
        _, ax = _make_section_axes()
        x = np.array([1.0])
        y = np.array([0.1])
        add_station_labels(ax, x, y, ["X"], color="#fff", rotation=45)
        assert ax.texts[0].get_rotation() == pytest.approx(45.0)

    def test_fontsize(self):
        _, ax = _make_section_axes()
        add_station_labels(
            ax,
            np.array([1.0]),
            np.array([0.1]),
            ["X"],
            color="#fff",
            fontsize=9,
        )
        assert ax.texts[0].get_fontsize() == pytest.approx(9.0)


# ---------------------------------------------------------------------------
# marker_style override (draw_topo_section / draw_topo_strip)
# ---------------------------------------------------------------------------


class TestMarkerStyleOverride:
    def test_draw_topo_section_default_uses_global_marker(self):
        from pycsamt.api.station import PYCSAMT_STATION_RENDERING

        chain, elev, names = _simple_inputs()
        _, ax = _make_section_axes()
        draw_topo_section(ax, chain, elev, names)
        scatter = ax.collections[0]
        expected = PYCSAMT_STATION_RENDERING.inversion.marker.facecolor
        import matplotlib.colors as mcolors

        assert tuple(scatter.get_facecolor()[0]) == mcolors.to_rgba(expected)

    def test_draw_topo_section_marker_style_override(self):
        from pycsamt.api.station import StationMarkerStyle

        chain, elev, names = _simple_inputs()
        _, ax = _make_section_axes()
        custom = StationMarkerStyle(facecolor="white", edgecolor="black")
        draw_topo_section(ax, chain, elev, names, marker_style=custom)
        scatter = ax.collections[0]
        assert tuple(scatter.get_facecolor()[0]) == (1.0, 1.0, 1.0, 1.0)
        assert tuple(scatter.get_edgecolor()[0]) == (0.0, 0.0, 0.0, 1.0)

    def test_draw_topo_strip_marker_style_override(self):
        from pycsamt.api.station import StationMarkerStyle

        chain, elev, names = _simple_inputs()
        fig, ax = _make_ps_axes()
        custom = StationMarkerStyle(facecolor="white", edgecolor="black")
        ax_strip = draw_topo_strip(fig, ax, chain, elev, names, marker_style=custom)
        # collections[0] is the fill_between ground polygon; the scatter is last.
        scatter = ax_strip.collections[-1]
        assert tuple(scatter.get_facecolor()[0]) == (1.0, 1.0, 1.0, 1.0)

    def test_draw_topo_strip_light_mode_background_is_transparent_by_default(self):
        chain, elev, names = _simple_inputs()
        fig, ax = _make_ps_axes()
        ax_strip = draw_topo_strip(fig, ax, chain, elev, names, dark=False)
        r, g, b, a = ax_strip.get_facecolor()
        assert a == pytest.approx(0.0)

    def test_draw_topo_strip_facecolor_override(self):
        chain, elev, names = _simple_inputs()
        fig, ax = _make_ps_axes()
        ax_strip = draw_topo_strip(
            fig, ax, chain, elev, names, dark=False, facecolor="#eff1f5"
        )
        import matplotlib.colors as mcolors

        assert ax_strip.get_facecolor() == mcolors.to_rgba("#eff1f5")

    def test_draw_topo_strip_dark_mode_background_unchanged(self):
        chain, elev, names = _simple_inputs()
        fig, ax = _make_ps_axes()
        ax_strip = draw_topo_strip(fig, ax, chain, elev, names, dark=True)
        import matplotlib.colors as mcolors

        assert ax_strip.get_facecolor() == mcolors.to_rgba("#181825")
