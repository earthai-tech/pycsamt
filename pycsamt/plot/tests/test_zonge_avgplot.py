# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
Tests for pycsamt.zonge.AVGPlot plotting methods.

These are **smoke tests**: they exercise each public plotting
API and assert basic return types (Figure/Axes) without
validating visual content. We use the Agg backend so tests can
run headless on CI.
"""
from __future__ import annotations

from typing import List
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import pytest

matplotlib.use("Agg")  # headless backend


# ----------------------------- fixtures ----------------------------- #

@pytest.fixture(scope="session")
def avg_plotter(modern_data_file: Path):
    """Load AMTAVG and construct an AVGPlot instance."""
    try:
        from pycsamt.zonge import AMTAVG # type: ignore
        from pycsamt.plot.zonge import AVGPlot  # type: ignore
    except Exception as exc:  # pragma: no cover
        pytest.skip(f"Import failed: {exc}")

    avg = AMTAVG.from_file(str(modern_data_file))
    plotter = AVGPlot(avg)
    return plotter


@pytest.fixture(scope="session")
def stn_file_any(data_path: Path) -> Path:
    """
    Return a usable STN file (prefer K2.stn; fallback to K1.stn).
    """
    k2 = data_path / "K2.stn"
    k1 = data_path / "K1.stn"
    if k2.exists():
        return k2
    if k1.exists():
        return k1
    pytest.skip("No STN file found (K2.stn/K1.stn).")


@pytest.fixture()
def ensure_topography(avg_plotter, stn_file_any: Path):
    """
    Attach topography if not present; return the plotter.
    """
    if not getattr(avg_plotter.avg, "topo", None):
        avg_plotter.avg.add_topography(stn_file=str(stn_file_any))
    return avg_plotter


@pytest.fixture()
def first_station(avg_plotter) -> float:
    """Get the first station id present in the AVG dataframe."""
    df = avg_plotter.avg.df
    return float(df["station"].iloc[0])


@pytest.fixture()
def subset_frequencies(avg_plotter) -> List[float]:
    """Pick a few representative frequencies for faster tests."""
    freqs = (
        avg_plotter.avg.df["freq"].dropna().unique()
    )
    freqs.sort()
    return list(map(float, freqs[:3]))


@pytest.fixture()
def corrected_df(avg_plotter) -> pd.DataFrame:
    """
    Create a synthetic 'corrected' dataframe compatible with
    plotting APIs (slight perturbations).
    """
    base = avg_plotter.avg.df.copy()
    if "rho" in base:
        base["rho"] = base["rho"] * 0.95
    if "phase" in base:
        base["phase"] = base["phase"] * 0.98
    if "emag" in base:
        base["emag"] = base["emag"] * 1.02
    if "hmag" in base:
        base["hmag"] = base["hmag"] * 0.99
    return base


# ------------------------------ helpers ----------------------------- #

def _assert_fig_ax(fig, ax_like) -> None:
    assert hasattr(fig, "savefig")
    # ax_like may be Axes or ndarray of Axes
    if isinstance(ax_like, np.ndarray):
        assert ax_like.size > 0
        assert all(
            hasattr(a, "plot") for a in ax_like.ravel()
        )
    else:
        assert hasattr(ax_like, "plot")


def _close(fig) -> None:
    try:
        plt.close(fig)
    except Exception:
        pass


# ------------------------------- tests ------------------------------ #

def test_plot_sounding_curves_default(avg_plotter):
    fig, axes = avg_plotter.plot_sounding_curves()
    _assert_fig_ax(fig, axes)
    _close(fig)


def test_plot_sounding_curves_single_station(
    avg_plotter, first_station, corrected_df
):
    fig, axes = avg_plotter.plot_sounding_curves(
        station_id=first_station,
        corrected_df=corrected_df,
        x_axis="frequency",
        phase_in_degrees=True,
    )
    _assert_fig_ax(fig, axes)
    _close(fig)


def test_plot_pseudosection_default(avg_plotter):
    fig, ax = avg_plotter.plot_pseudosection(
        var="rho", comp="ExHy"
    )
    _assert_fig_ax(fig, ax)
    _close(fig)


def test_plot_pseudosections_multi(avg_plotter):
    fig, axes = avg_plotter.plot_pseudosections(
        vars=["rho", "phase"],
        comp="ExHy",
        max_cols=2,
    )
    _assert_fig_ax(fig, axes)
    _close(fig)


def test_plot_tensor_soundings_default(avg_plotter):
    fig, axes = avg_plotter.plot_tensor_soundings(
        tensor="rho", x_axis="frequency", todeg=True
    )
    _assert_fig_ax(fig, axes)
    _close(fig)


def test_plot_tensor_soundings_single_station(
    avg_plotter, first_station
):
    fig, axes = avg_plotter.plot_tensor_soundings(
        station_id=first_station, tensor="z", show_fit=True
    )
    _assert_fig_ax(fig, axes)
    _close(fig)


def test_plot_phase_tensor_ellipses(
    avg_plotter, first_station, subset_frequencies
):
    fig, axes = avg_plotter.plot_phase_tensor_ellipses(
        station_id=first_station,
        frequencies=subset_frequencies,
    )
    _assert_fig_ax(fig, axes)
    _close(fig)


def test_plot_remediation(avg_plotter, first_station, corrected_df):
    fig, ax = avg_plotter.plot_remediation(
        station_id=first_station,
        corrected_df=corrected_df,
    )
    _assert_fig_ax(fig, ax)
    _close(fig)


def test_plot_station_single(avg_plotter, first_station):
    fig = avg_plotter.plot_station(
        station_id=first_station,
        x_axis="frequency",
        todeg=True,
    )
    assert hasattr(fig, "savefig")
    _close(fig)


def test_plot_profile_basic(ensure_topography):
    plt = ensure_topography
    fig, ax = plt.plot_profile(
        station_names="mask",  # fully masked top axis
        right_axis="mask",     # fully hidden right axis
    )
    _assert_fig_ax(fig, ax)
    _close(fig)


def test_plot_profile_with_names(ensure_topography):
    plt = ensure_topography
    names = plt.avg.topo.stations.astype(str)
    fig, ax = plt.plot_profile(
        station_names=list(names),
        tick_step=5,
        right_axis="show",
        top_label="Station",
    )
    _assert_fig_ax(fig, ax)
    _close(fig)


def test_plot_location_map_default(ensure_topography):
    plt = ensure_topography
    fig, ax = plt.plot_location_map(
        kind="scatter", n_levels=8
    )
    _assert_fig_ax(fig, ax)
    _close(fig)


def test_plot_strike_default(avg_plotter):
    fig, ax = avg_plotter.plot_strike(num_bins=24)
    _assert_fig_ax(fig, ax)
    _close(fig)
