# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Coverage tests for :mod:`pycsamt.forward.plot`.

Exercises the 1-D/2-D/3-D plotting helpers with small, cheap, *real*
forward-model runs (few frequencies, small grids) rather than mocks,
consistent with ``test_forward_em1d.py``.  Matplotlib uses the ``Agg``
backend so the whole file runs headless.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from pycsamt.api.control import PYCSAMT_CONTROL
from pycsamt.forward.em1d import MT1DForward
from pycsamt.forward.em2d import MT2DForward
from pycsamt.forward.em3d import MT3DForward
from pycsamt.forward.grid2d import Grid2D
from pycsamt.forward.grid3d import Grid3D
from pycsamt.forward.plot import (
    _add_colorbar,
    _phase_label,
    _phase_vals,
    _rho_label,
    _rho_vals,
    _spine_style,
    _x_label,
    _x_log,
    _x_vals,
    plot_model_1d,
    plot_model_2d,
    plot_model_3d,
    plot_pseudosection_2d,
    plot_response_1d,
    plot_response_and_model_1d,
    plot_response_map_3d,
    plot_response_profiles,
    plot_response_section_3d,
    plot_tensor_components_3d,
)
from pycsamt.forward.synthetic import LayeredModel


@pytest.fixture(autouse=True)
def _close_figures():
    """Avoid figure leakage across tests."""
    yield
    plt.close("all")


# ─────────────────────────────────────────────────────────────────────────────
# Shared cheap fixtures
# ─────────────────────────────────────────────────────────────────────────────


@pytest.fixture
def freqs1d():
    return np.logspace(-2, 3, 12)


@pytest.fixture
def model1d():
    return LayeredModel(resistivity=[100.0, 10.0, 500.0], thickness=[200.0, 400.0])


@pytest.fixture
def resp1d(freqs1d, model1d):
    return MT1DForward(freqs1d).run(model1d)


@pytest.fixture
def small_grid2d():
    return Grid2D.halfspace(
        rho=100.0, nx=8, nz=6, x_max=1600.0, z_max=800.0, n_pad=3, n_stations=5
    )


@pytest.fixture
def anomaly_grid2d():
    return Grid2D.with_anomaly(
        bg_rho=100.0,
        anomaly_rho=5.0,
        anomaly_bounds=(400.0, 1200.0, 100.0, 400.0),
        nx=10,
        nz=8,
        x_max=2000.0,
        z_max=1000.0,
        n_pad=2,
        n_stations=6,
    )


@pytest.fixture
def resp2d(anomaly_grid2d):
    freqs = np.logspace(-1, 2, 4)
    return MT2DForward(freqs, anomaly_grid2d, verbose=False).run()


@pytest.fixture
def small_grid3d():
    return Grid3D.halfspace(
        rho=100.0,
        nx=6,
        ny=6,
        nz=5,
        x_max=1200.0,
        y_max=1200.0,
        z_max=600.0,
        n_pad=2,
        nx_stations=3,
        ny_stations=3,
    )


@pytest.fixture
def resp3d(small_grid3d):
    freqs = np.logspace(-1, 2, 3)
    return MT3DForward(freqs, small_grid3d, verbose=False).run()


# ─────────────────────────────────────────────────────────────────────────────
# Private helpers
# ─────────────────────────────────────────────────────────────────────────────


def test_x_vals_label_log_default():
    freqs = np.array([1.0, 10.0, 100.0])
    x = _x_vals(freqs)
    # default view is log10_period -> log10(1/f)
    assert np.allclose(x, np.log10(1.0 / freqs))
    assert _x_label() == r"$\log_{10}T$ (s)"
    assert _x_log() is False


def test_x_vals_label_log_period_view():
    freqs = np.array([1.0, 10.0, 100.0])
    with PYCSAMT_CONTROL.context(x__view="period"):
        x = _x_vals(freqs)
        assert np.allclose(x, 1.0 / freqs)
        assert _x_label() == "Period (s)"
        assert _x_log() is True


def test_rho_vals_label_log_and_linear():
    rho = np.array([1.0, 10.0, 100.0])
    assert np.allclose(_rho_vals(rho), np.log10(rho))
    assert "log" in _rho_label().lower() or r"\log" in _rho_label()
    with PYCSAMT_CONTROL.context(rho__view="linear"):
        assert np.allclose(_rho_vals(rho), rho)
        assert _rho_label() == r"$\rho_a$ ($\Omega\,\mathrm{m}$)"


def test_phase_vals_label_default():
    phase = np.array([10.0, -20.0, 190.0])
    out = _phase_vals(phase)
    # default wraps into (-180, 180]
    assert np.all(out > -180.0) and np.all(out <= 180.0)
    assert _phase_label() == r"Phase ($^\circ$)"


def test_spine_style_sets_grid_and_axisbelow():
    _, ax = plt.subplots()
    _spine_style(ax)
    assert ax.get_axisbelow() is True
    # grid lines were switched on (both x and y major gridlines visible)
    assert all(gl.get_visible() for gl in ax.xaxis.get_gridlines())
    assert all(gl.get_visible() for gl in ax.yaxis.get_gridlines())


def test_add_colorbar_sets_label():
    fig, ax = plt.subplots()
    data = np.array([[1.0, 2.0], [3.0, 4.0]])
    im = ax.imshow(data)
    n_axes_before = len(fig.axes)
    _add_colorbar(fig, ax, im, "my label", fontsize=10)
    assert len(fig.axes) == n_axes_before + 1
    cb_ax = fig.axes[-1]
    assert cb_ax.get_ylabel() == "my label"


# ─────────────────────────────────────────────────────────────────────────────
# plot_response_1d
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_response_1d_basic(resp1d):
    axs = plot_response_1d(resp1d)
    assert isinstance(axs, np.ndarray) and axs.shape == (2,)
    ax_r, ax_p = axs
    assert isinstance(ax_r, Axes) and isinstance(ax_p, Axes)
    # Only TE is drawn since a 1-D ForwardResponse has no rho_a_tm.
    assert len(ax_r.lines) == 1
    assert len(ax_p.lines) == 1
    assert ax_r.get_legend() is not None


def test_plot_response_1d_axes_provided(resp1d):
    fig, axs_in = plt.subplots(2, 1)
    axs_out = plot_response_1d(resp1d, axes=axs_in)
    assert axs_out[0].get_figure() is fig
    assert axs_out[0] is axs_in[0]


def test_plot_response_1d_show_te_false(resp1d):
    axs = plot_response_1d(resp1d, show_te=False)
    ax_r, ax_p = axs
    # plot_tm stays True (default) but response has no TM data -> no lines
    assert len(ax_r.lines) == 0
    assert len(ax_p.lines) == 0


def test_plot_response_1d_modes_te_only(resp1d):
    axs = plot_response_1d(resp1d, modes="te")
    ax_r, _ = axs
    assert len(ax_r.lines) == 1


def test_plot_response_1d_modes_single_element_sequence(resp1d):
    axs = plot_response_1d(resp1d, modes=["te"])
    ax_r, _ = axs
    assert len(ax_r.lines) == 1


def test_plot_response_1d_modes_multi_element_sequence_equivalent_to_both(
    resp1d,
):
    # A 1-D response has no rho_a_tm, so "tm" never draws a curve regardless
    # -- but modes=["te", "tm"] must still select "te" (regression test for
    # a bug where join(["te", "tm"]) == "te,tm" matched neither "te"/"tm"/
    # "both" and silently plotted nothing).
    axs_seq = plot_response_1d(resp1d, modes=["te", "tm"])
    axs_both = plot_response_1d(resp1d, modes="both")
    assert len(axs_seq[0].lines) == len(axs_both[0].lines) == 1


def test_plot_response_1d_explicit_style_overrides(resp1d):
    axs = plot_response_1d(
        resp1d,
        color_te="green",
        color_tm="purple",
        lw=3.0,
        marker_te="d",
        marker_tm="x",
        ms=10.0,
        label_te="Custom TE",
        label_tm="Custom TM",
    )
    ax_r, _ = axs
    line = ax_r.lines[0]
    assert line.get_color() == "green"
    assert line.get_linewidth() == 3.0
    assert line.get_marker() == "d"
    assert line.get_markersize() == 10.0
    assert line.get_label() == "Custom TE"


def test_plot_response_1d_no_lines_no_legend_when_nothing_shown(resp1d):
    # modes="te" forces plot_tm False regardless of show_tm; show_te=False
    # forces plot_te False too -> neither branch draws, no legend created.
    axs = plot_response_1d(resp1d, modes="te", show_te=False)
    ax_r, ax_p = axs
    assert len(ax_r.lines) == 0
    assert ax_r.get_legend() is None


def test_plot_response_1d_title_and_xlog(resp1d):
    with PYCSAMT_CONTROL.context(x__view="frequency"):
        axs = plot_response_1d(resp1d, title="My sounding")
    fig = axs[0].get_figure()
    assert fig._suptitle is not None
    assert fig._suptitle.get_text() == "My sounding"
    assert axs[0].get_xscale() == "log"


def test_plot_response_1d_two_station_response_with_tm(freqs1d):
    """Simulate a multi-station response exposing rho_a_tm/phase_tm.

    plot_response_1d documents support for responses whose rho_a/phase
    are 2-D (n_freqs, n_stations); this is not produced anywhere by the
    real 1-D/2-D solvers, so we hand-build one from two real MT1D runs
    stacked as columns, to exercise that branch with real physics values.
    """
    model_a = LayeredModel(resistivity=[100.0], thickness=[])
    model_b = LayeredModel(resistivity=[10.0, 300.0], thickness=[250.0])
    resp_a = MT1DForward(freqs1d).run(model_a)
    resp_b = MT1DForward(freqs1d).run(model_b)

    @dataclass
    class _MultiStationResponse:
        freqs: np.ndarray
        rho_a: np.ndarray
        phase: np.ndarray
        rho_a_tm: np.ndarray
        phase_tm: np.ndarray

    fake = _MultiStationResponse(
        freqs=freqs1d,
        rho_a=np.column_stack([resp_a.rho_a, resp_b.rho_a]),
        phase=np.column_stack([resp_a.phase, resp_b.phase]),
        rho_a_tm=np.column_stack([resp_b.rho_a, resp_a.rho_a]),
        phase_tm=np.column_stack([resp_b.phase, resp_a.phase]),
    )

    axs = plot_response_1d(fake, modes="both")
    ax_r, ax_p = axs
    assert len(ax_r.lines) == 2
    assert len(ax_p.lines) == 2
    # TE line (col 0) must match resp_a's rho_a in display (log10) space
    te_line = ax_r.lines[0]
    assert np.allclose(te_line.get_ydata(), np.log10(resp_a.rho_a))


# ─────────────────────────────────────────────────────────────────────────────
# plot_model_1d
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_model_1d_single_model(model1d):
    ax = plot_model_1d(model1d)
    assert isinstance(ax, Axes)
    assert len(ax.lines) == 1
    assert ax.get_legend() is None  # n == 1 -> no legend
    assert ax.get_xscale() == "log"  # log_rho default True


def test_plot_model_1d_multiple_models_with_labels():
    m1 = LayeredModel(resistivity=[100.0, 10.0], thickness=[300.0], name="A")
    m2 = LayeredModel(resistivity=[50.0, 500.0], thickness=[150.0], name="B")
    ax = plot_model_1d([m1, m2], labels=["model-1", "model-2"], title="two models")
    assert len(ax.lines) == 2
    leg = ax.get_legend()
    assert leg is not None
    texts = [t.get_text() for t in leg.get_texts()]
    assert texts == ["model-1", "model-2"]
    assert ax.get_title() == "two models"


def test_plot_model_1d_default_labels_use_model_name():
    m1 = LayeredModel(resistivity=[100.0, 10.0], thickness=[300.0], name="Aquifer")
    m2 = LayeredModel(resistivity=[50.0, 500.0], thickness=[150.0], name="")
    ax = plot_model_1d([m1, m2])
    leg = ax.get_legend()
    texts = [t.get_text() for t in leg.get_texts()]
    assert texts == ["Aquifer", "model 2"]


def test_plot_model_1d_linear_scale_and_depth_max(model1d):
    ax = plot_model_1d(model1d, log_rho=False, depth_max=5000.0)
    assert ax.get_xscale() == "linear"
    ys = ax.lines[0].get_ydata()
    assert ys.max() == pytest.approx(5000.0)


def test_plot_model_1d_ax_provided(model1d):
    _, ax_in = plt.subplots()
    ax_out = plot_model_1d(model1d, ax=ax_in)
    assert ax_out is ax_in


def test_plot_model_1d_explicit_lw_alpha(model1d):
    ax = plot_model_1d(model1d, lw=4.0, alpha=0.3)
    line = ax.lines[0]
    assert line.get_linewidth() == 4.0
    assert line.get_alpha() == 0.3


# ─────────────────────────────────────────────────────────────────────────────
# plot_response_and_model_1d
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_response_and_model_1d_with_model(resp1d, model1d):
    fig = plot_response_and_model_1d(resp1d, model1d, title="Combined")
    assert isinstance(fig, Figure)
    assert len(fig.axes) == 3
    ax_m, ax_r, ax_p = fig.axes
    assert ax_m.get_title() == "Earth model"
    assert "resistivity" in ax_r.get_title().lower()
    assert "phase" in ax_p.get_title().lower()
    assert fig._suptitle.get_text() == "Combined"


def test_plot_response_and_model_1d_without_model(resp1d):
    fig = plot_response_and_model_1d(resp1d)
    assert len(fig.axes) == 2


def test_plot_response_and_model_1d_custom_gridspec(resp1d, model1d):
    fig = plot_response_and_model_1d(
        resp1d, model1d, gridspec_kw={"width_ratios": [1, 1, 1]}
    )
    assert len(fig.axes) == 3


def test_plot_response_and_model_1d_2d_response_and_xlog(freqs1d, model1d):
    """Exercise the rho_a.ndim > 1 branch and the log-x-scale branch."""
    resp = MT1DForward(freqs1d).run(model1d)

    @dataclass
    class _MultiStationResponse:
        freqs: np.ndarray
        rho_a: np.ndarray
        phase: np.ndarray

    fake = _MultiStationResponse(
        freqs=freqs1d,
        rho_a=np.column_stack([resp.rho_a, resp.rho_a * 2.0]),
        phase=np.column_stack([resp.phase, resp.phase]),
    )
    with PYCSAMT_CONTROL.context(x__view="period"):
        fig = plot_response_and_model_1d(fake, model1d)
    ax_m, ax_r, ax_p = fig.axes
    assert ax_r.get_xscale() == "log"
    line = ax_r.lines[0]
    assert np.allclose(line.get_ydata(), np.log10(resp.rho_a))


# ─────────────────────────────────────────────────────────────────────────────
# plot_model_2d
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_model_2d_basic(small_grid2d):
    ax = plot_model_2d(small_grid2d)
    assert isinstance(ax, Axes)
    assert len(ax.collections) >= 1  # pcolormesh
    assert ax.get_title() == "2-D resistivity model"
    # colorbar was attached (extra axes on the figure)
    assert len(ax.get_figure().axes) == 2


def test_plot_model_2d_linear_scale_no_clip(anomaly_grid2d):
    ax = plot_model_2d(anomaly_grid2d, log_scale=False, clip_core=False)
    pc = ax.collections[0]
    data = pc.get_array()
    # linear scale -> raw resistivities, including the 5 ohm.m anomaly
    assert np.isclose(data.min(), 5.0, atol=1e-6) or data.min() < 100.0


def test_plot_model_2d_vmin_vmax_and_no_stations(small_grid2d):
    ax = plot_model_2d(small_grid2d, vmin=1.0, vmax=3.0, show_stations=False)
    pc = ax.collections[0]
    assert pc.get_clim() == (1.0, 3.0)


def test_plot_model_2d_custom_title_and_ax(small_grid2d):
    _, ax_in = plt.subplots()
    ax = plot_model_2d(small_grid2d, title="custom title", ax=ax_in)
    assert ax is ax_in
    assert ax.get_title() == "custom title"


# ─────────────────────────────────────────────────────────────────────────────
# plot_pseudosection_2d
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_pseudosection_2d_rho_a_te(resp2d):
    ax = plot_pseudosection_2d(resp2d, mode="te", quantity="rho_a")
    assert isinstance(ax, Axes)
    assert len(ax.collections) >= 1
    assert "resistivity" in ax.get_title().lower()
    assert "TE" in ax.get_title()


def test_plot_pseudosection_2d_phase_tm(resp2d):
    ax = plot_pseudosection_2d(resp2d, mode="tm", quantity="phase")
    assert "Phase" in ax.get_title()
    assert "TM" in ax.get_title()


def test_plot_pseudosection_2d_contours_and_custom_cmap(resp2d):
    ax = plot_pseudosection_2d(resp2d, n_contours=3, cmap="viridis")
    # contour lines add LineCollections beyond the base pcolormesh
    assert len(ax.collections) >= 2


def test_plot_pseudosection_2d_no_stations(resp2d):
    ax = plot_pseudosection_2d(resp2d, show_stations=False)
    assert ax.get_xlabel() == "Distance  (m)"


def test_plot_pseudosection_2d_freqs_ascending_inverts_yaxis(anomaly_grid2d):
    freqs = np.logspace(-1, 2, 4)  # ascending frequency -> descending period
    resp = MT2DForward(freqs, anomaly_grid2d, verbose=False).run()
    ax = plot_pseudosection_2d(resp)
    assert ax.yaxis_inverted()


def test_plot_pseudosection_2d_ax_provided(resp2d):
    _, ax_in = plt.subplots()
    ax = plot_pseudosection_2d(resp2d, ax=ax_in)
    assert ax is ax_in


def test_plot_pseudosection_2d_freqs_descending_no_invert(anomaly_grid2d):
    freqs = np.logspace(2, -1, 4)  # descending frequency -> ascending period
    resp = MT2DForward(freqs, anomaly_grid2d, verbose=False).run()
    ax = plot_pseudosection_2d(resp)
    assert not ax.yaxis_inverted()


# ─────────────────────────────────────────────────────────────────────────────
# plot_response_profiles
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_response_profiles_default(resp2d):
    ax = plot_response_profiles(resp2d, n_freqs_shown=3)
    assert len(ax.lines) == 3
    assert ax.get_legend() is not None


def test_plot_response_profiles_explicit_indices(resp2d):
    ax = plot_response_profiles(resp2d, freq_indices=[0, 1])
    assert len(ax.lines) == 2


def test_plot_response_profiles_phase_tm(resp2d):
    ax = plot_response_profiles(resp2d, mode="tm", quantity="phase", n_freqs_shown=2)
    assert "TM" in ax.get_ylabel()
    assert len(ax.lines) == 2


def test_plot_response_profiles_ax_provided_and_title(resp2d):
    _, ax_in = plt.subplots()
    ax = plot_response_profiles(resp2d, ax=ax_in, title="lateral")
    assert ax is ax_in
    assert ax.get_title() == "lateral"


def test_plot_response_profiles_explicit_lw_alpha(resp2d):
    ax = plot_response_profiles(resp2d, lw=3.5, alpha=0.4, n_freqs_shown=2)
    for line in ax.lines:
        assert line.get_linewidth() == 3.5
        assert line.get_alpha() == 0.4


# ─────────────────────────────────────────────────────────────────────────────
# plot_model_3d
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_model_3d_basic(small_grid3d):
    axs = plot_model_3d(small_grid3d)
    assert isinstance(axs, np.ndarray) and axs.shape == (3,)
    for ax in axs:
        assert len(ax.collections) >= 1
    fig = axs[0].get_figure()
    assert fig._suptitle.get_text() == small_grid3d.name


def test_plot_model_3d_no_clip_no_stations_linear(small_grid3d):
    axs = plot_model_3d(
        small_grid3d, clip_core=False, show_stations=False, log_scale=False,
        title="linear view",
    )
    assert axs.shape == (3,)
    fig = axs[0].get_figure()
    assert fig._suptitle.get_text() == "linear view"


def test_plot_model_3d_stations_scatter(small_grid3d):
    axs = plot_model_3d(small_grid3d, show_stations=True)
    ax_xy = axs[2]
    scatters = [c for c in ax_xy.collections if c.get_label() == "stations"]
    assert len(scatters) == 1
    assert scatters[0].get_offsets().shape[0] == small_grid3d.n_stations


def test_plot_model_3d_vmin_vmax(small_grid3d):
    axs = plot_model_3d(small_grid3d, vmin=1.0, vmax=3.0)
    pc = axs[0].collections[0]
    assert pc.get_clim() == (1.0, 3.0)


# ─────────────────────────────────────────────────────────────────────────────
# plot_response_map_3d
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_response_map_3d_basic(resp3d, small_grid3d):
    ax = plot_response_map_3d(resp3d)
    assert isinstance(ax, Axes)
    scatter = ax.collections[0]
    assert scatter.get_offsets().shape[0] == small_grid3d.n_stations
    assert "Z_XY" in ax.get_title()


def test_plot_response_map_3d_phase_component_and_labels_off(resp3d):
    ax = plot_response_map_3d(
        resp3d, component="yx", quantity="phase", show_labels=False
    )
    assert "Z_YX" in ax.get_title()
    assert len(ax.texts) == 0


def test_plot_response_map_3d_labels_on(resp3d, small_grid3d):
    ax = plot_response_map_3d(resp3d, show_labels=True)
    assert len(ax.texts) == small_grid3d.n_stations


def test_plot_response_map_3d_ax_provided_and_title(resp3d):
    _, ax_in = plt.subplots()
    ax = plot_response_map_3d(resp3d, ax=ax_in, title="custom map")
    assert ax is ax_in
    assert ax.get_title() == "custom map"


def test_plot_response_map_3d_vmin_vmax(resp3d):
    ax = plot_response_map_3d(resp3d, vmin=-1.0, vmax=5.0)
    assert ax.collections[0].get_clim() == (-1.0, 5.0)


def test_plot_response_map_3d_explicit_cmap(resp3d):
    ax = plot_response_map_3d(resp3d, cmap="plasma")
    assert ax.collections[0].get_cmap().name == "plasma"


# ─────────────────────────────────────────────────────────────────────────────
# plot_response_section_3d
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_response_section_3d_default_midrow(resp3d):
    ax = plot_response_section_3d(resp3d)
    assert isinstance(ax, Axes)
    assert len(ax.collections) >= 1
    assert "Z_XY" in ax.get_title()


def test_plot_response_section_3d_explicit_yrow_and_phase(resp3d):
    ax = plot_response_section_3d(resp3d, y_row=0, component="yy", quantity="phase")
    assert "Z_YY" in ax.get_title()


def test_plot_response_section_3d_contours(resp3d):
    ax = plot_response_section_3d(resp3d, n_contours=2)
    assert len(ax.collections) >= 2


def test_plot_response_section_3d_no_stations(resp3d):
    ax = plot_response_section_3d(resp3d, show_stations=False)
    assert ax.get_xlabel() == "x  (m)"


def test_plot_response_section_3d_ax_provided(resp3d):
    _, ax_in = plt.subplots()
    ax = plot_response_section_3d(resp3d, ax=ax_in)
    assert ax is ax_in


def test_plot_response_section_3d_explicit_cmap(resp3d):
    ax = plot_response_section_3d(resp3d, cmap="plasma")
    assert ax.collections[0].get_cmap().name == "plasma"


def test_plot_response_section_3d_single_frequency(small_grid3d):
    # A single frequency collapses the y_log array to length 1, exercising
    # the "else" branches for both the y-edge construction and the
    # invert-y-axis comparison (y_log[0] == y_log[-1] -> no inversion).
    resp = MT3DForward(np.array([50.0]), small_grid3d, verbose=False).run()
    ax = plot_response_section_3d(resp)
    assert not ax.yaxis_inverted()


def test_plot_response_section_3d_single_station_per_row():
    # A single x-station per y-row collapses st_x to length 1, exercising
    # the "else" branch of the x-edge construction.
    grid = Grid3D.halfspace(
        rho=100.0,
        nx=6,
        ny=6,
        nz=5,
        x_max=1200.0,
        y_max=1200.0,
        z_max=600.0,
        n_pad=2,
        nx_stations=1,
        ny_stations=2,
    )
    freqs = np.logspace(-1, 2, 3)
    resp = MT3DForward(freqs, grid, verbose=False).run()
    ax = plot_response_section_3d(resp)
    assert isinstance(ax, Axes)


# ─────────────────────────────────────────────────────────────────────────────
# plot_tensor_components_3d
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_tensor_components_3d_basic(resp3d):
    axs = plot_tensor_components_3d(resp3d)
    assert isinstance(axs, np.ndarray) and axs.shape == (2, 2)
    for ax in axs.ravel():
        assert len(ax.collections) >= 1
    fig = axs[0, 0].get_figure()
    assert "Apparent resistivity" in fig._suptitle.get_text()


def test_plot_tensor_components_3d_phase_and_title(resp3d):
    axs = plot_tensor_components_3d(resp3d, quantity="phase", title="Full tensor")
    fig = axs[0, 0].get_figure()
    assert fig._suptitle.get_text() == "Full tensor"
