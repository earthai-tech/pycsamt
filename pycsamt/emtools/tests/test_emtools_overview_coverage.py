"""Focused branch coverage for the response-overview renderer."""

from __future__ import annotations

from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.colors import Normalize

from pycsamt.emtools import overview


def _impedance(n: int = 4) -> np.ndarray:
    z = np.empty((n, 2, 2), dtype=complex)
    for k in range(n):
        z[k] = np.array(
            [[1.0 + 0.2j, 2.0 + 2.5j], [-2.5 - 2.0j, 0.8 + 0.3j]]
        ) * (k + 1)
    return z


def _patch_overview_inputs(monkeypatch, *, tensor=True, tipper=True):
    freq = np.array([100.0, 10.0, 1.0, 0.1])
    z = _impedance()
    ze = np.full(z.shape, 0.05 + 0.05j)
    ed = SimpleNamespace(edi=object())
    monkeypatch.setattr(overview, "ensure_sites", lambda *a, **k: object())
    monkeypatch.setattr(overview, "_pick_station", lambda *a, **k: ("S01", ed))
    monkeypatch.setattr(overview, "_zblk_flex", lambda obj: (None, z, freq, ze))
    tip = np.array(
        [[0.1 + 0.02j, 0.2 - 0.03j], [0.0j, 0.0j],
         [0.2 - 0.1j, -0.1 + 0.2j], [np.nan + 0j, 0.1j]]
    )
    monkeypatch.setattr(
        overview,
        "_get_t_block",
        lambda *a, **k: (None, tip if tipper else None, freq if tipper else None),
    )
    table = pd.DataFrame(
        {
            "freq": freq,
            "s1": [1.0, 1.2, 0.8, 1.1],
            "s2": [0.5, 0.2, 0.7, 0.4],
            "theta": [0.0, 25.0, 50.0, 80.0],
            "skew": [-8.0, 0.0, 3.0, 12.0],
            "phimin_deg": [20.0, 30.0, 40.0, 50.0],
        }
    )
    monkeypatch.setattr(
        overview,
        "build_phase_tensor_table",
        lambda *a, **k: table.copy() if tensor else table.iloc[:0].copy(),
    )
    return freq, z, ze


def test_overview_drawing_helpers_cover_filters_and_styles():
    assert overview._x_offset(10.0, 1.0, True) == pytest.approx(100.0)
    assert overview._x_offset(10.0, 1.0, False) == pytest.approx(11.0)

    handles, labels = overview._legend_handles(
        ("xy", "yx"),
        lambda comp: overview._component_style(comp, False, False),
        {"xy": "black"},
    )
    assert labels == ["$Z_{xy}$", "$Z_{yx}$"]
    assert handles[0].get_color() == "black"

    fig, (ax1, ax2) = plt.subplots(2, 1)
    tips = overview._draw_induction_arrow_row(
        ax1,
        np.array([1.0, 10.0, 100.0, np.nan]),
        np.array([1.0, 0.0, 1e-12, 1.0]),
        np.array([1.0, 1.0, 0.0, 1.0]),
        color="red",
        tilt_decades=0.2,
        dy_scale=2.0,
        lw=1.0,
        mutation_scale=5.0,
        log_x=True,
    )
    assert len(tips) == 2

    ax2.set_xlim(0.0, 4.0)
    overview._draw_pt_ellipse_row(
        ax2,
        np.array([1.0, 2.0, 3.0]),
        np.array([1.0, 1.0, np.nan]),
        np.array([0.5, 0.1, 0.2]),
        np.array([0.0, 45.0, 10.0]),
        np.array([-2.0, 2.0, 0.0]),
        np.array([0.0, 8.0, 0.0]),
        cmap=plt.get_cmap("viridis"),
        norm=Normalize(-2.0, 2.0),
        scale=1.0,
        min_aspect=0.2,
        edgecolor="k",
        linewidth=0.5,
        alpha=0.8,
        skew_threshold=5.0,
        mark_3d=True,
        cells_per_decade=4.0,
        log_x=False,
    )
    assert len(ax2.patches) == 2
    assert ax2.patches[1].get_linewidth() == pytest.approx(1.5)
    plt.close(fig)


def test_plot_response_overview_full_vertical_layout(monkeypatch):
    _patch_overview_inputs(monkeypatch)
    fig = overview.plot_response_overview(
        object(),
        cbar_orientation="vertical",
        c_by="skew",
        clim=(-10.0, 10.0),
        colors={"xy": "navy"},
        phase_range=(-90.0, 90.0),
        show_phase_error_bars=True,
        ylim_arrows=(-0.5, 0.5),
        ellipse_colorbar_label="beta",
    )
    assert fig._suptitle.get_text() == "S01"
    assert len(fig.axes) == 7
    assert fig.axes[-1].get_ylabel() == "beta"
    plt.close(fig)


def test_plot_response_overview_external_axes_and_placeholders(monkeypatch):
    _patch_overview_inputs(monkeypatch, tensor=False, tipper=False)
    fig, axes = plt.subplots(4, 1)
    result = overview.plot_response_overview(
        object(),
        axes=axes,
        show_diag=False,
        x_view="log10_period",
        log_log_rho=False,
        show_component_legend=False,
        show_ellipse_colorbar=False,
        grid=False,
    )
    assert result is fig
    assert any(t.get_text() == "no tipper" for t in axes[2].texts)
    assert any(t.get_text() == "no phase tensor data" for t in axes[3].texts)
    plt.close(fig)


def test_plot_response_overview_disabled_rows_and_validation(monkeypatch):
    _patch_overview_inputs(monkeypatch)
    fig = overview.plot_response_overview(
        object(), show_diag=False, show_arrows=False, show_ellipses=False
    )
    assert not fig.axes[2].axison
    assert not fig.axes[3].axison
    plt.close(fig)

    with pytest.raises(ValueError, match="cbar_orientation"):
        overview.plot_response_overview(object(), cbar_orientation="diagonal")


def test_plot_response_overview_no_impedance(monkeypatch):
    _patch_overview_inputs(monkeypatch)
    monkeypatch.setattr(overview, "_zblk_flex", lambda obj: (None, None, None))
    fig = overview.plot_response_overview(object(), show_diag=False)
    assert all(ax.texts[0].get_text() == "no impedance data" for ax in fig.axes[:4])
    plt.close(fig)
