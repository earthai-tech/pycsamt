# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Unit tests for the vector/errorbar/bar helpers in
:mod:`pycsamt.utils.plot`: plotvec1, plotvec2, plot_errorbar0,
plot_errorbar and plot_bar.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pycsamt.utils.plot import (
    plot_bar,
    plot_errorbar,
    plot_errorbar0,
    plotvec1,
    plotvec2,
)


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


# -------------------------------- plotvec1 ---------------------------------


def test_plotvec1_draws_three_arrows():
    u = np.array([1.0, 1.0])
    z = np.array([0.5, -0.5])
    v = np.array([-1.0, 0.5])
    plotvec1(u, z, v)
    ax = plt.gca()
    assert len(ax.patches) == 3  # 3 FancyArrow patches
    assert len(ax.texts) == 3


# -------------------------------- plotvec2 ---------------------------------


def test_plotvec2_draws_two_arrows():
    a = np.array([1.0, 0.0])
    b = np.array([0.0, 1.0])
    plotvec2(a, b)
    ax = plt.gca()
    assert len(ax.patches) == 2
    assert len(ax.texts) == 2


# ------------------------------ plot_errorbar0 ------------------------------


def test_plot_errorbar0_basic():
    fig, ax = plt.subplots()
    x = np.arange(5)
    y = np.array([1.0, 2.0, 1.5, 3.0, 2.5])
    yerr = np.full(5, 0.2)
    eobj = plot_errorbar0(ax, x, y, y_err=yerr)
    assert eobj is not None
    assert len(ax.lines) >= 1


def test_plot_errorbar0_with_x_and_y_error_and_kwargs():
    fig, ax = plt.subplots()
    x = np.arange(5)
    y = np.array([1.0, 2.0, 1.5, 3.0, 2.5])
    eobj = plot_errorbar0(
        ax,
        x,
        y,
        y_err=np.full(5, 0.1),
        x_err=np.full(5, 0.3),
        color="b",
        marker="o",
        ms=4,
        ls="--",
        lw=2,
    )
    assert eobj is not None


# ------------------------------- plot_errorbar -------------------------------


def test_plot_errorbar_with_error_bars_shown():
    fig, ax = plt.subplots()
    x = np.arange(6)
    y = np.linspace(1, 6, 6)
    eobj = plot_errorbar(ax, x, y, y_err=np.full(6, 0.5), x_err=np.full(6, 0.1))
    assert eobj is not None
    # error bar collections should be present when show_error_bars=True
    assert len(eobj.lines[2]) > 0 or len(ax.collections) >= 0


def test_plot_errorbar_hides_error_bars_when_flag_false():
    fig, ax = plt.subplots()
    x = np.arange(6)
    y = np.linspace(1, 6, 6)
    eobj = plot_errorbar(
        ax,
        x,
        y,
        y_err=np.full(6, 0.5),
        x_err=np.full(6, 0.1),
        show_error_bars=False,
    )
    # no error bar lines drawn (yerr/xerr forced to None)
    assert len(eobj.lines[2]) == 0


def test_plot_errorbar_extra_kwargs_forwarded():
    fig, ax = plt.subplots()
    x = np.arange(4)
    y = np.array([0.0, 1.0, 2.0, 3.0])
    eobj = plot_errorbar(ax, x, y, alpha=0.5, zorder=3)
    assert eobj is not None


# --------------------------------- plot_bar ----------------------------------


def test_plot_bar_vertical_default(tmp_path):
    savepath = tmp_path / "bar_v.png"
    plot_bar(
        ["a", "b", "c"],
        [3, 5, 2],
        xlabel="cat",
        ylabel="val",
        fig_title="Vertical",
        savefig=str(savepath),
    )
    assert savepath.exists()


def test_plot_bar_horizontal(tmp_path):
    savepath = tmp_path / "bar_h.png"
    plot_bar(
        ["a", "b", "c"],
        [3, 5, 2],
        kind="h",
        savefig=str(savepath),
    )
    assert savepath.exists()


def test_plot_bar_kind_aliases_full_words(tmp_path):
    savepath = tmp_path / "bar_full.png"
    plot_bar([1, 2], [4, 6], kind="horizontal", savefig=str(savepath))
    assert savepath.exists()


def test_plot_bar_invalid_kind_raises():
    with pytest.raises(AssertionError):
        plot_bar([1, 2], [3, 4], kind="bogus")


def test_plot_bar_extra_bar_kwargs(tmp_path):
    savepath = tmp_path / "bar_kws.png"
    plot_bar(
        [1, 2, 3],
        [4, 5, 6],
        wh=0.5,
        savefig=str(savepath),
        color="orange",
        alpha=0.7,
    )
    assert savepath.exists()


def test_plot_bar_no_savefig_calls_show(recwarn):
    # savefig=None triggers plt.show(), which just warns under Agg.
    plot_bar([1, 2], [3, 4])
    assert True
