# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Unit tests for the confidence-interval/ellipse and text-annotation
helpers in :mod:`pycsamt.utils.plot`: plot_confidence,
plot_confidence_ellipse, confidence_ellipse and plot_text.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from pycsamt.utils.plot import (
    confidence_ellipse,
    plot_confidence,
    plot_confidence_ellipse,
    plot_text,
)

try:
    import seaborn  # noqa: F401

    _HAS_SNS = True
except ImportError:
    _HAS_SNS = False


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


# ------------------------------ confidence_ellipse ---------------------------


def test_confidence_ellipse_adds_patch_to_axes():
    rng = np.random.default_rng(0)
    x = rng.normal(size=200)
    y = 0.5 * x + rng.normal(scale=0.5, size=200)
    fig, ax = plt.subplots()
    n_before = len(ax.patches)
    patch = confidence_ellipse(x, y, ax, n_std=2, edgecolor="red")
    assert len(ax.patches) == n_before + 1
    assert patch in ax.patches


def test_confidence_ellipse_mismatched_sizes_raise():
    fig, ax = plt.subplots()
    with pytest.raises(ValueError):
        confidence_ellipse(np.arange(5), np.arange(4), ax)


def test_confidence_ellipse_multiple_std_levels():
    rng = np.random.default_rng(1)
    x = rng.normal(size=100)
    y = rng.normal(size=100)
    fig, ax = plt.subplots()
    for n_std in (1, 2, 3):
        confidence_ellipse(x, y, ax, n_std=n_std)
    assert len(ax.patches) == 3


# ---------------------------- plot_confidence_ellipse -------------------------


def test_plot_confidence_ellipse_runs_and_draws_three_ellipses():
    rng = np.random.default_rng(2)
    x = rng.normal(size=150)
    y = rng.normal(size=150)
    plot_confidence_ellipse(x, y)
    ax = plt.gcf().axes[0]
    assert len(ax.patches) == 3
    assert ax.get_title() == "Different standard deviations"


# -------------------------------- plot_confidence -----------------------------


@pytest.mark.skipif(not _HAS_SNS, reason="seaborn not installed")
def test_plot_confidence_line_kind_with_seaborn():
    x = np.arange(20)
    y = np.arange(20) + np.random.default_rng(3).normal(scale=0.5, size=20)
    ax = plot_confidence(x=x, y=y, kind="line")
    assert ax is not None


@pytest.mark.skipif(not _HAS_SNS, reason="seaborn not installed")
def test_plot_confidence_reg_kind_with_seaborn():
    x = np.arange(20)
    y = 2 * x + np.random.default_rng(4).normal(scale=1.0, size=20)
    ax = plot_confidence(x=x, y=y, kind="reg")
    assert ax is not None


def test_plot_confidence_bootstrap_kind_prints_interval(capsys):
    y = list(np.random.default_rng(5).normal(loc=10, scale=2, size=40))
    ax = plot_confidence(y=y, kind="bootstrap", b_samples=50)
    assert ax is None
    out = capsys.readouterr().out
    assert "confidence interval" in out


def test_plot_confidence_bootstrap_requires_y():
    with pytest.raises(ValueError):
        plot_confidence(x=np.arange(5), kind="bootstrap")


# ---------------------------------- plot_text ----------------------------------


def test_plot_text_basic_with_matching_text():
    x = np.arange(5)
    y = np.arange(5) * 2
    text = [f"S{i}" for i in range(5)]
    ax = plot_text(x, y, text=text)
    assert len(ax.texts) == 5
    # one scatter marker per non-empty text label
    assert len(ax.collections) == 5


def test_plot_text_coerce_generates_basenames():
    x = np.arange(3)
    y = np.arange(3)
    ax = plot_text(x, y, coerce=True, basename="P")
    assert len(ax.texts) == 3


def test_plot_text_show_line_and_legend():
    x = np.arange(4)
    y = np.array([0.0, 1.0, 0.5, 2.0])
    text = ["a", "b", "c", "d"]
    ax = plot_text(
        x,
        y,
        text=text,
        show_line=True,
        show_leg=True,
        linelabel="line1",
        markerlabel="marker1",
        lcolor="blue",
        mcolor="green",
        color="black",
    )
    assert ax.get_legend() is not None
    assert len(ax.lines) == 1


def test_plot_text_step_blanks_some_labels():
    x = np.arange(6)
    y = np.arange(6)
    text = [f"T{i}" for i in range(6)]
    ax = plot_text(x, y, text=text, step=2)
    non_empty = [t for t in ax.texts if t.get_text() != ""]
    assert len(non_empty) < 6


def test_plot_text_missing_text_without_coerce_raises():
    x = np.arange(3)
    y = np.arange(3)
    with pytest.raises(TypeError):
        plot_text(x, y)


def test_plot_text_length_mismatch_without_coerce_raises():
    x = np.arange(4)
    y = np.arange(4)
    with pytest.raises(ValueError):
        plot_text(x, y, text=["only", "two"])


def test_plot_text_none_xy_raises():
    with pytest.raises(TypeError):
        plot_text(None, None, text=["a"])


def test_plot_text_series_inputs_set_labels_from_name():
    df = pd.DataFrame({"east": [0, 1, 2, 3], "north": [0, 1, 2, 3]})
    ax = plot_text(df["east"], df["north"], text=["a", "b", "c", "d"])
    assert ax.get_xlabel() == "east"
    assert ax.get_ylabel() == "north"


def test_plot_text_string_column_names_with_data():
    df = pd.DataFrame({"east": [0, 1, 2, 3], "north": [0, 1, 2, 3]})
    ax = plot_text("east", "north", data=df, text=["a", "b", "c", "d"])
    assert ax.get_xlabel() == "east"
    assert ax.get_ylabel() == "north"


def test_plot_text_external_axes():
    fig, ax0 = plt.subplots()
    x = np.arange(3)
    y = np.arange(3)
    ax = plot_text(x, y, text=["a", "b", "c"], ax=ax0)
    assert ax is ax0
