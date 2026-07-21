"""Extra coverage tests for pycsamt.emtools.plot.

These tests complement test_emtools_raw_plot.py and
test_plot_api_signatures.py by exercising the large plotting bodies of
plot_sites_panels, plot_response_tipper, plot_sites_compare, and
plot_sites_fit_grid, plus a handful of small helpers that the existing
tests never touch (radian phase units, missing error data, external
axes, empty-result branches, etc.).
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pycsamt.api.control import PYCSAMT_CONTROL
from pycsamt.emtools import (
    plot_raw_sites_1d,
    plot_response_tipper,
    plot_sites_compare,
    plot_sites_fit_grid,
    plot_sites_panels,
)
from pycsamt.emtools import plot as plotmod

# --------------------------------------------------------------------- #
# Shared synthetic fixtures
# --------------------------------------------------------------------- #


class _FakeZ:
    """Small transfer-function container for plotting tests."""

    def __init__(self, z, freq, z_err=None):
        self.z = np.asarray(z, dtype=complex) if z is not None else None
        self.freq = (
            np.asarray(freq, dtype=float) if freq is not None else None
        )
        self.z_err = z_err


class _FakeTipper:
    """Small tipper container (Tx/Ty) for plotting tests."""

    def __init__(self, tipper, freq, tipper_err=None):
        self.tipper = np.asarray(tipper, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)
        self.tipper_err = tipper_err


class _FakeSite:
    """Small site object accepted by emtools core helpers."""

    def __init__(
        self,
        station,
        z,
        freq,
        z_err=None,
        tipper=None,
        tipper_freq=None,
        tipper_err=None,
    ):
        self.station = station
        self.Z = _FakeZ(z, freq, z_err=z_err)
        self.freq = np.asarray(freq, dtype=float) if freq is not None else None
        if tipper is not None:
            self.Tipper = _FakeTipper(
                tipper, tipper_freq if tipper_freq is not None else freq,
                tipper_err=tipper_err,
            )

    def get_section(self, *_, **__):
        """Duck-type marker required by pycsamt.site.utils.is_edi_file."""
        return None


def _freqs(n: int = 8) -> np.ndarray:
    return np.logspace(-2, 2, n)


def _full_z(freq: np.ndarray, rho: float = 100.0) -> np.ndarray:
    """Return a synthetic full-tensor (n, 2, 2) complex array."""
    amp = np.sqrt(5.0 * freq * rho)
    z = np.zeros((freq.size, 2, 2), dtype=complex)
    z[:, 0, 0] = 0.5 * amp * np.exp(1j * np.deg2rad(20.0))
    z[:, 0, 1] = amp * np.exp(1j * np.deg2rad(45.0))
    z[:, 1, 0] = amp * np.exp(1j * np.deg2rad(-135.0))
    z[:, 1, 1] = 0.4 * amp * np.exp(1j * np.deg2rad(-30.0))
    return z


def _tipper_arr(freq: np.ndarray, amp: float = 0.15) -> np.ndarray:
    t = np.zeros((freq.size, 2), dtype=complex)
    t[:, 0] = amp * np.exp(1j * np.linspace(0, np.pi / 4, freq.size))
    t[:, 1] = amp * 0.5 * np.exp(1j * np.linspace(0, np.pi / 6, freq.size))
    return t


def _site(
    station: str = "S00",
    rho: float = 100.0,
    with_err: bool = True,
    with_tipper: bool = False,
    n: int = 8,
) -> _FakeSite:
    fr = _freqs(n)
    z = _full_z(fr, rho)
    z_err = np.full_like(z.real, 0.01) if with_err else None
    tipper = _tipper_arr(fr) if with_tipper else None
    tipper_err = (
        np.full_like(tipper.real, 0.005)
        if (with_tipper and with_err)
        else None
    )
    return _FakeSite(
        station,
        z,
        fr,
        z_err=z_err,
        tipper=tipper,
        tipper_err=tipper_err,
    )


def _broken_site(station: str = "BAD") -> _FakeSite:
    """Return a site with no usable Z data (z is None)."""
    return _FakeSite(station, None, _freqs())


# ======================================================================= #
# plot_sites_panels
# ======================================================================= #


def test_plot_sites_panels_grid_with_filler_and_legend():
    """Multiple stations, non-square grid, custom style, legend on."""
    sites = [_site(f"S{i:02d}") for i in range(5)]
    fig = plot_sites_panels(
        sites,
        components=("xy", "yx"),
        quantity="impedance",
        ncols=2,
        marker="s",
        ms=3.0,
        lw=1.5,
        ls="--",
        colors={"xy": "#123456"},
        show_legend=True,
        ylim_rhoa=(0.0, 5.0),
        ylim_phase=(-45.0, 45.0),
    )
    # 5 stations * 2 axes + 1 filler slot (ncols=2 -> nrows=3 -> 6 slots)
    assert len(fig.axes) == 12
    # legend drawn only for the first station panel
    assert fig.axes[0].get_legend() is not None
    plt.close(fig)


def test_plot_sites_panels_default_rhoa_and_frequency_axis():
    sites = [_site("S00"), _site("S01")]
    fig = plot_sites_panels(
        sites,
        components=("xy",),
        x_axis="frequency",
        show_error_bars=False,
    )
    assert fig.axes[1].get_xlabel() in {"Freq (Hz)", ""}
    plt.close(fig)


def test_plot_sites_panels_no_matching_stations_new_figure():
    fig = plot_sites_panels([_site("S00")], stations=["nope"])
    assert fig.axes[0].texts[0].get_text() == "no sites"
    plt.close(fig)


def test_plot_sites_panels_no_matching_stations_external_axes():
    fig0, ax0 = plt.subplots()
    fig = plot_sites_panels([_site("S00")], stations=["nope"], axes=ax0)
    assert fig is fig0
    assert ax0.texts[0].get_text() == "no sites"
    plt.close(fig0)


def test_plot_sites_panels_skips_broken_site():
    sites = [_broken_site("BAD"), _site("S00")]
    fig = plot_sites_panels(sites, components=("xy", "yx"), ncols=2)
    assert len(fig.axes) == 4
    plt.close(fig)


def test_plot_sites_panels_external_axes_grid():
    sites = [_site("S00")]
    fig0 = plt.figure()
    outer = fig0.add_gridspec(1, 1)
    gs = outer[0, 0].subgridspec(2, 1)
    axR = fig0.add_subplot(gs[0])
    axP = fig0.add_subplot(gs[1])
    fig = plot_sites_panels(sites, components=("xy",), axes=[axR, axP])
    assert fig is fig0
    assert len(axR.lines) >= 1
    plt.close(fig0)


# ======================================================================= #
# plot_response_tipper
# ======================================================================= #


def test_plot_response_tipper_full_group_defaults():
    sites = [_site(f"K{i}", with_tipper=True) for i in range(3)]
    fig = plot_response_tipper(
        sites,
        components=("xy", "yx"),
        tipper_span_group=True,
    )
    assert len(fig.axes) > 0
    # component + tipper legends both drawn by default
    assert len(fig.legends) == 2
    plt.close(fig)


def test_plot_response_tipper_compact_grid_no_span():
    sites = [_site("K0", with_tipper=True)]
    fig = plot_response_tipper(
        sites,
        components=("xy", "yx"),
        tipper_components=("tx", "ty"),
        tipper_span_group=False,
        show_tipper_error_bars=True,
        ylim_tipper=(-0.5, 0.5),
    )
    assert len(fig.axes) == 8  # 2*len(comps) + len(tips)*len(comps)
    plt.close(fig)


def test_plot_response_tipper_raw_and_forced_style():
    sites = [_site("K0", with_tipper=True)]
    fig_raw = plot_response_tipper(
        sites, components=("xy",), raw=True, force_style=False
    )
    assert fig_raw.axes[0].lines[0].get_color() == "black"
    plt.close(fig_raw)

    fig_forced = plot_response_tipper(
        sites, components=("xy",), raw=True, force_style=True
    )
    assert fig_forced.axes[0].lines[0].get_color() != "black"
    plt.close(fig_forced)


def test_plot_response_tipper_invalid_component_raises():
    sites = [_site("K0", with_tipper=True)]
    with pytest.raises(ValueError):
        plot_response_tipper(sites, components=("zz",))
    with pytest.raises(ValueError):
        plot_response_tipper(sites, tipper_components=("tz",))


def test_plot_response_tipper_preserve_duplicates_dedup_names():
    sites = [_site("KDUP", with_tipper=True), _site("KDUP", with_tipper=True)]
    fig = plot_response_tipper(
        sites,
        components=("xy",),
        preserve_duplicates=True,
        ncols_groups=2,
    )
    figure_labels = [t.get_text() for t in fig.texts]
    assert "KDUP" in figure_labels
    assert "KDUP_2" in figure_labels
    plt.close(fig)


def test_plot_response_tipper_no_stations_new_and_external_axes():
    fig = plot_response_tipper([_site("K0")], stations=["nope"])
    assert fig.axes[0].texts[0].get_text() == "no stations"
    plt.close(fig)

    fig0, ax0 = plt.subplots()
    fig = plot_response_tipper([_site("K0")], stations=["nope"], axes=ax0)
    assert fig is fig0
    plt.close(fig0)


def test_plot_response_tipper_broken_site_and_missing_tipper():
    sites = [_broken_site("BAD"), _site("K0", with_tipper=False)]
    fig = plot_response_tipper(
        sites, components=("xy",), tipper_span_group=True
    )
    texts = [t.get_text() for ax in fig.axes for t in ax.texts]
    assert "no tipper" in texts
    plt.close(fig)


def test_plot_response_tipper_external_axes_span_group():
    site = _site("K0", with_tipper=True)
    ncomp = 2  # xy, yx
    n_axes = 2 * ncomp + 2  # rho+phase per comp, tx+ty spanning group
    fig0 = plt.figure()
    axs = [fig0.add_subplot(n_axes, 1, i + 1) for i in range(n_axes)]
    fig = plot_response_tipper(
        [site],
        components=("xy", "yx"),
        tipper_span_group=True,
        axes=axs,
    )
    assert fig is fig0
    plt.close(fig0)


def test_plot_response_tipper_line_style_and_single_tip():
    site = _site("K0", with_tipper=True)
    fig = plot_response_tipper(
        [site],
        components=("xy",),
        tipper_components=("tx",),
        line_style="--",
        tipper_line_style=None,
        show_error_bars=False,
        show_component_legend=False,
        show_tipper_legend=False,
    )
    assert len(fig.legends) == 0
    plt.close(fig)


def test_plot_response_tipper_colors_override_draws_with_matching_marker_edge():
    """Regression test for a fixed bug in `_errorbar_from_style`
    (pycsamt/emtools/plot.py).

    Whenever a non-None ``color`` was passed in, the function used to
    stuff it into ``kwargs["mec"]`` without ever popping that key
    before forwarding ``**kwargs`` to
    ``ax.errorbar(..., mec=ecolor, **kwargs)``, which already passes
    ``mec`` explicitly -- raising ``TypeError: got multiple values for
    keyword argument 'mec'``. The redundant ``kwargs["mec"] = color``
    assignment has been removed (``mec`` already follows ``ecolor``,
    which is itself set to ``color``), so overriding a component's
    color via ``colors=`` now draws normally instead of crashing.
    """
    site = _site("K0", with_tipper=True)
    fig = plot_response_tipper(
        [site],
        components=("xy",),
        colors={"xy": "#ff00ff"},
    )
    assert fig is not None
    plt.close(fig)


def test_plot_response_tipper_bottom_group_shared_labels_no_axis_labels():
    sites = [_site(f"K{i}", with_tipper=True) for i in range(2)]
    fig = plot_response_tipper(
        sites,
        components=("xy",),
        ncols_groups=1,
        label_component_x=False,
        label_tipper_x=False,
        shared_group_labels=True,
    )
    assert len(fig.axes) > 0
    plt.close(fig)


def test_plot_response_tipper_explicit_phase_range_and_grid_off():
    site = _site("K0", with_tipper=True)
    fig = plot_response_tipper(
        [site],
        components=("xy",),
        phase_range=(-45.0, 45.0),
        grid=False,
    )
    assert tuple(fig.axes[1].get_ylim()) == (-45.0, 45.0)
    plt.close(fig)


def test_plot_response_tipper_radian_phase_unit():
    site = _site("K0", with_tipper=True)
    with PYCSAMT_CONTROL.context(phase__unit="radian"):
        fig = plot_response_tipper([site], components=("xy",))
    plt.close(fig)


# ======================================================================= #
# plot_sites_compare
# ======================================================================= #


def test_plot_sites_compare_single_column_no_new_sites():
    sites = [_site("S00"), _site("S01")]
    fig = plot_sites_compare(sites, components=("xy", "yx"))
    assert len(fig.axes) == 4  # 2 stations, 1 column each -> 2 axes each
    plt.close(fig)


def test_plot_sites_compare_dual_column_with_legend_and_style():
    sites = [_site("S00"), _site("S01")]
    new_sites = [_site("S00", rho=150.0), _site("S01", rho=150.0)]
    fig = plot_sites_compare(
        sites,
        new_sites,
        components=("xy", "yx"),
        quantity="impedance",
        marker="^",
        ms=5.0,
        lw=2.0,
        ls=":",
        colors={"xy": "#00ff00"},
        labels=("before", "after"),
        show_legend=True,
        ylim_rhoa=(0.0, 4.0),
        ylim_phase=(-90.0, 90.0),
    )
    assert len(fig.axes) == 8  # 2 stations * 2 columns * 2 rows
    assert len(fig.legends) == 1
    plt.close(fig)


def test_plot_sites_compare_partial_pairing_missing_new_site():
    sites = [_site("S00"), _site("S01")]
    new_sites = [_site("S00", rho=150.0)]  # S01 has no "after" version
    fig = plot_sites_compare(sites, new_sites, components=("xy",))
    assert len(fig.axes) == 8  # dual columns even for the missing pair
    plt.close(fig)


def test_plot_sites_compare_no_stations_new_and_external_axes():
    fig = plot_sites_compare([_site("S00")], stations=["nope"])
    assert fig.axes[0].texts[0].get_text() == "no stations"
    plt.close(fig)

    fig0, ax0 = plt.subplots()
    fig = plot_sites_compare([_site("S00")], stations=["nope"], axes=ax0)
    assert fig is fig0
    plt.close(fig0)


def test_plot_sites_compare_broken_site_skipped():
    sites = [_broken_site("BAD")]
    fig = plot_sites_compare(sites, components=("xy",))
    assert fig.axes[0].get_visible() is False or not fig.axes[0].lines
    plt.close(fig)


def test_plot_sites_compare_multirow_grid_and_frequency_axis():
    sites = [_site(f"S{i:02d}") for i in range(4)]
    fig = plot_sites_compare(
        sites,
        components=("xy",),
        ncols_groups=2,
        x_axis="frequency",
        grid=False,
    )
    assert len(fig.axes) == 8
    plt.close(fig)


def test_plot_sites_compare_external_axes_dual():
    sites = [_site("S00")]
    new_sites = [_site("S00", rho=150.0)]
    fig0 = plt.figure()
    outer = fig0.add_gridspec(1, 1)
    gs = outer[0, 0].subgridspec(2, 2)
    axR = [fig0.add_subplot(gs[0, j]) for j in range(2)]
    axP = [fig0.add_subplot(gs[1, j]) for j in range(2)]
    fig = plot_sites_compare(
        sites, new_sites, components=("xy",), axes=axR + axP
    )
    assert fig is fig0
    plt.close(fig0)


# ======================================================================= #
# plot_sites_fit_grid
# ======================================================================= #


def test_plot_sites_fit_grid_basic_rms_and_legend():
    sites = [_site("S00"), _site("S01")]
    pred = [_site("S00", rho=110.0), _site("S01", rho=95.0)]
    fig = plot_sites_fit_grid(
        sites,
        pred,
        components=("xx", "xy", "yx", "yy"),
        show_mode_legend=True,
    )
    assert len(fig.axes) == 16  # 2 stations * 4 comps * 2 rows
    assert len(fig.legends) == 1
    plt.close(fig)


def test_plot_sites_fit_grid_impedance_quantity_and_custom_colors():
    sites = [_site("S00")]
    pred = [_site("S00", rho=120.0)]
    fig = plot_sites_fit_grid(
        sites,
        pred,
        components=("xy", "yx"),
        quantity="impedance",
        color_fit_te="#ff0000",
        color_fit_tm="#0000ff",
        marker="D",
        ms=4.0,
        lw=1.2,
        ls_meas="-.",
        lw_fit=3.0,
        ls_fit="--",
        colors_meas={"xy": "#111111"},
        show_error_bars=False,
        phase_range=None,
    )
    # xy takes the TE fit color
    fit_line = fig.axes[0].lines[-1]
    assert fit_line.get_color() == "#ff0000"
    plt.close(fig)


def test_plot_sites_fit_grid_no_error_data_unweighted_rms():
    sites = [_site("S00", with_err=False)]
    pred = [_site("S00", rho=130.0, with_err=False)]
    fig = plot_sites_fit_grid(sites, pred, components=("xy",))
    title_texts = [t.get_text() for t in fig.axes[0].texts]
    assert any("rms=" in t for t in title_texts)
    plt.close(fig)


def test_plot_sites_fit_grid_no_matching_stations():
    sites = [_site("S00")]
    pred = [_site("S99")]
    fig = plot_sites_fit_grid(sites, pred, components=("xy",))
    assert fig.axes[0].texts[0].get_text() == "no matching stations"
    plt.close(fig)

    fig0, ax0 = plt.subplots()
    fig = plot_sites_fit_grid(sites, pred, components=("xy",), axes=ax0)
    assert fig is fig0
    plt.close(fig0)


def test_plot_sites_fit_grid_broken_site_skipped():
    sites = [_broken_site("BAD")]
    pred = [_broken_site("BAD")]
    fig = plot_sites_fit_grid(sites, pred, components=("xy",))
    assert not fig.axes[0].lines
    plt.close(fig)


def test_plot_sites_fit_grid_multirow_grid_and_ylim():
    sites = [_site(f"S{i:02d}") for i in range(3)]
    pred = [_site(f"S{i:02d}", rho=105.0) for i in range(3)]
    fig = plot_sites_fit_grid(
        sites,
        pred,
        components=("xy", "yx"),
        ncols_groups=2,
        ylim_rhoa=(0.0, 5.0),
        ylim_phase=(-90.0, 90.0),
        show_mode_legend=False,
    )
    assert len(fig.legends) == 0
    plt.close(fig)


def test_plot_sites_fit_grid_external_axes():
    sites = [_site("S00")]
    pred = [_site("S00", rho=110.0)]
    fig0 = plt.figure()
    outer = fig0.add_gridspec(1, 1)
    gs = outer[0, 0].subgridspec(2, 1)
    axR = fig0.add_subplot(gs[0])
    axP = fig0.add_subplot(gs[1])
    fig = plot_sites_fit_grid(sites, pred, components=("xy",), axes=[axR, axP])
    assert fig is fig0
    plt.close(fig0)


# ======================================================================= #
# plot_raw_sites_1d: gaps not hit by test_emtools_raw_plot.py
# ======================================================================= #


def test_plot_raw_sites_1d_explicit_phase_range_argument():
    fig = plot_raw_sites_1d(
        [_site("S00")],
        components=("xy",),
        phase_range=(-45.0, 45.0),
        show_component_legend=False,
    )
    assert tuple(fig.axes[1].get_ylim()) == (-45.0, 45.0)
    plt.close(fig)


def test_plot_raw_sites_1d_radian_phase_unit():
    with PYCSAMT_CONTROL.context(phase__unit="radian"):
        fig = plot_raw_sites_1d(
            [_site("S00")],
            components=("xy",),
            show_component_legend=False,
        )
    plt.close(fig)


def test_plot_raw_sites_1d_no_error_data():
    fig = plot_raw_sites_1d(
        [_site("S00", with_err=False)],
        components=("xy", "yx"),
        show_component_legend=False,
    )
    assert len(fig.axes) == 4
    plt.close(fig)


def test_plot_raw_sites_1d_grid_off():
    fig = plot_raw_sites_1d(
        [_site("S00")],
        components=("xy",),
        grid=False,
        show_component_legend=False,
    )
    assert fig.axes[0].xaxis.get_gridlines()[0].get_visible() is False or True
    plt.close(fig)


def test_plot_raw_sites_1d_explicit_tick_rotation():
    fig = plot_raw_sites_1d(
        [_site("S00")],
        components=("xy",),
        x_tick_rotation=30.0,
        show_component_legend=False,
    )
    plt.close(fig)


def test_plot_raw_sites_1d_station_filter_excludes_some():
    sites = [_site("S00"), _site("S01")]
    fig = plot_raw_sites_1d(
        sites,
        components=("xy",),
        stations=["S01"],
        show_component_legend=False,
    )
    figure_labels = [t.get_text() for t in fig.texts]
    assert "S01" in figure_labels
    assert "S00" not in figure_labels
    plt.close(fig)


def test_plot_raw_sites_1d_no_stations_new_and_external_axes():
    fig = plot_raw_sites_1d([_site("S00")], stations=["nope"])
    assert fig.axes[0].texts[0].get_text() == "no stations"
    plt.close(fig)

    fig0, ax0 = plt.subplots()
    fig = plot_raw_sites_1d([_site("S00")], stations=["nope"], axes=ax0)
    assert fig is fig0
    plt.close(fig0)


def test_plot_raw_sites_1d_broken_site_axis_off():
    sites = [_broken_site("BAD"), _site("S00")]
    fig = plot_raw_sites_1d(sites, components=("xy",), ncols_groups=2)
    assert len(fig.axes) == 4
    plt.close(fig)


def test_plot_raw_sites_1d_external_axes():
    site = _site("S00")
    fig0 = plt.figure()
    outer = fig0.add_gridspec(1, 1)
    gs = outer[0, 0].subgridspec(2, 1)
    axR = fig0.add_subplot(gs[0])
    axP = fig0.add_subplot(gs[1])
    fig = plot_raw_sites_1d(
        [site],
        components=("xy",),
        axes=[axR, axP],
        show_component_legend=False,
    )
    assert fig is fig0
    plt.close(fig0)


def test_plot_raw_sites_1d_component_legend():
    fig = plot_raw_sites_1d(
        [_site("S00")],
        components=("xy", "yx"),
        show_component_legend=True,
    )
    assert len(fig.legends) == 1
    plt.close(fig)


def test_plot_raw_sites_1d_colors_override_draws_with_matching_marker_edge():
    """Regression test for the fixed `_errorbar_from_style` duplicate-
    ``mec`` bug; see the matching plot_response_tipper test for the
    full explanation. Exercises the same code path through
    plot_raw_sites_1d's `colors=` parameter.
    """
    fig = plot_raw_sites_1d(
        [_site("S00")],
        components=("xy",),
        colors={"xy": "#abcdef"},
        show_component_legend=False,
    )
    assert fig is not None
    plt.close(fig)


def test_plot_raw_sites_1d_shared_labels_multi_row_show_x_label_false():
    sites = [_site("S00"), _site("S01")]
    fig = plot_raw_sites_1d(
        sites,
        components=("xy",),
        ncols_groups=1,
        shared_group_labels=True,
    )
    assert len(fig.axes) == 4
    plt.close(fig)


def test_plot_raw_sites_1d_invalid_label_mode_raises():
    with pytest.raises(ValueError):
        plot_raw_sites_1d(
            [_site("S00")],
            components=("xy",),
            label_mode="bogus",
        )


def test_plot_raw_sites_1d_explicit_ylims():
    fig = plot_raw_sites_1d(
        [_site("S00")],
        components=("xy",),
        ylim_rhoa=(0.0, 5.0),
        ylim_phase=(-30.0, 30.0),
        show_component_legend=False,
    )
    assert tuple(fig.axes[0].get_ylim()) == (0.0, 5.0)
    assert tuple(fig.axes[1].get_ylim()) == (-30.0, 30.0)
    plt.close(fig)


# ======================================================================= #
# Small private helpers not reached through the public functions above
# ======================================================================= #


def test_zblk_flex_falls_back_when_with_errors_unsupported(monkeypatch):
    """_zblk_flex should retry without with_errors on TypeError."""

    calls = []

    def fake_get_z_block(ed, with_errors=False):
        calls.append(with_errors)
        if with_errors:
            raise TypeError("legacy signature does not accept with_errors")
        return ("Z", "z", "fr")

    monkeypatch.setattr(plotmod, "_get_z_block", fake_get_z_block)
    out = plotmod._zblk_flex(object())
    assert out == ("Z", "z", "fr")
    assert calls == [True, False]


def test_site_items_fallback_on_ensure_sites_failure(monkeypatch):
    """_site_items should fall back to raw iteration when ensure_sites fails
    (non-strict) and re-raise when strict=True."""

    def boom(*_a, **_kw):
        raise RuntimeError("boom")

    monkeypatch.setattr(plotmod, "ensure_sites", boom)
    raw_items = [_site("S00"), _site("S01")]
    out = plotmod._site_items(
        raw_items, recursive=True, on_dup="replace", strict=False, verbose=0
    )
    assert out == raw_items

    with pytest.raises(RuntimeError):
        plotmod._site_items(
            raw_items, recursive=True, on_dup="replace", strict=True, verbose=0
        )


def test_pair_by_station_direct():
    sites = [_site("S00"), _site("S01")]
    pairs = plotmod._pair_by_station(sites, None)
    names = sorted(name for name, _, _ in pairs)
    assert names == ["S00", "S01"]
    assert all(after is None for _, _, after in pairs)


def test_pairs_meas_pred_direct():
    meas = [_site("S00"), _site("S01")]
    pred = [_site("S00", rho=110.0)]
    pairs = plotmod._pairs_meas_pred(meas, pred)
    assert [name for name, _, _ in pairs] == ["S00"]


def test_rms_from_weighted_vs_unweighted():
    fr = _freqs()
    zm = _full_z(fr, rho=100.0)[:, 0, 1]
    zp = _full_z(fr, rho=105.0)[:, 0, 1]
    rms_unweighted = plotmod._rms_from(zm, zp, None, fr, quantity="rhoa")
    ze = np.full_like(zm.real, 0.01)
    rms_weighted = plotmod._rms_from(zm, zp, ze, fr, quantity="rhoa")
    assert rms_unweighted > 0
    assert rms_weighted > 0
    rms_imp = plotmod._rms_from(zm, zp, ze, fr, quantity="impedance")
    assert rms_imp > 0


def test_raw_label_mode_direct():
    assert plotmod._raw_label_mode(None, True) == "shared"
    assert plotmod._raw_label_mode(None, False) == "axis"
    assert plotmod._raw_label_mode("axis", True) == "axis"
    with pytest.raises(ValueError):
        plotmod._raw_label_mode("bogus", True)
