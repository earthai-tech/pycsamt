"""Tests for TDEM plotting helpers."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pycsamt.api.control import PYCSAMT_CONTROL
from pycsamt.api.style import MultilineStyle
from pycsamt.tdem import (
    TEMAVG,
    PlotDecayCurve,
    PlotElevationProfile,
    PlotGateProfile,
    PlotSurveyMap,
    PlotSurveyOverview,
    PlotTEMAVGSection,
    PlotTEMDashboard,
    PlotTEMZSection,
    PlotTransformedRho,
    StationTickConfig,
    TDEMPlotStyle,
    TEMSounding,
    TEMZPlot,
    plot_decay,
    plot_elevation_profile,
    plot_gate_profile,
    plot_survey_map,
    plot_survey_overview,
    plot_tem_dashboard,
    plot_tem_z_section,
    plot_temavg_section,
    plot_transformed_rho,
    read_temavg_survey,
)

DATA_DIR = Path(__file__).parents[3] / "data" / "TEMAVG" / "JIANGSU"
AVG_FILE = DATA_DIR / "TEM100.AVG"
Z_FILE = DATA_DIR / "TEM100.Z"


def _sounding() -> TEMSounding:
    """Return a small synthetic sounding for plotting."""
    t = np.logspace(-5, -3, 8)
    data = 1e-6 * t**-1.5
    return TEMSounding(
        t,
        data,
        current=8.0,
        tx_area=1e4,
        station_name="S01",
    )


def test_plot_decay_curve_class_and_function():
    """Decay plot helpers should return Matplotlib axes."""
    ax = PlotDecayCurve(_sounding()).plot()
    ax2 = plot_decay(_sounding())

    assert ax.get_xscale() == "log"
    assert ax.get_yscale() == "log"
    assert ax2.get_xscale() == "log"
    plt.close(ax.figure)
    plt.close(ax2.figure)


def test_plot_decay_uses_package_multiline_style():
    """TDEM line colors should accept the shared multiline style API."""
    style = TDEMPlotStyle(
        multiline=MultilineStyle(
            mode="cycle",
            cycle_palette=["#123456"],
        ),
    )
    ax = PlotDecayCurve([_sounding(), _sounding()], style=style).plot()

    assert ax.lines[0].get_color() == "#123456"
    assert ax.lines[1].get_color() == "#123456"
    plt.close(ax.figure)


def test_plot_transformed_rho_returns_figure_by_default():
    """Transformed rho helper should create a two-panel figure."""
    fig = PlotTransformedRho(_sounding()).plot()
    fig2 = plot_transformed_rho(_sounding(), show_phase=False)

    assert len(fig.axes) == 2
    assert (
        fig.axes[0].get_position().height > fig.axes[1].get_position().height
    )
    assert len(fig2.axes) == 1
    plt.close(fig)
    plt.close(fig2)


def test_plot_transformed_rho_can_use_global_control():
    """Transformed rho plots should optionally follow api.control."""
    with PYCSAMT_CONTROL.context(
        rho__view="log10",
        x__view="log10_period",
        phase__range=(-90.0, 90.0),
    ):
        fig = PlotTransformedRho(_sounding(), use_control=True).plot()

    assert r"\log_{10}\rho_a" in fig.axes[0].get_ylabel()
    assert r"\log_{10}T" in fig.axes[-1].get_xlabel()
    plt.close(fig)


@pytest.mark.skipif(
    not AVG_FILE.exists(),
    reason="TEMAVG sample data not available",
)
def test_plot_temavg_section_class_and_function():
    """TEMAVG pseudo-section helpers should draw a color mesh."""
    avg = TEMAVG.read(AVG_FILE)
    ticks = StationTickConfig(every=5, rotation=30, fontsize=7)
    ax = PlotTEMAVGSection(avg, station_ticks=ticks).plot()
    ax2 = plot_temavg_section(avg, value="magnitude", absolute=True)

    assert ax.collections
    assert ax2.collections
    assert ax.get_xticks().size > 0
    assert ax.xaxis.get_label_position() == "top"
    plt.close(ax.figure)
    plt.close(ax2.figure)


@pytest.mark.skipif(
    not Z_FILE.exists(),
    reason="TEMAVG Z sample data not available",
)
def test_plot_tem_z_section_class_and_function():
    """ZPLOT pseudo-section helpers should draw a color mesh."""
    zplot = TEMZPlot.read(Z_FILE)
    ax = PlotTEMZSection(
        zplot,
        station_ticks=StationTickConfig(every="auto", max_ticks=8),
    ).plot()
    ax2 = plot_tem_z_section(zplot)

    assert ax.collections
    assert ax2.collections
    assert ax.xaxis.get_label_position() == "top"
    plt.close(ax.figure)
    plt.close(ax2.figure)


@pytest.mark.skipif(
    not AVG_FILE.exists(),
    reason="TEMAVG sample data not available",
)
def test_plot_survey_map_with_coordinate_rows(tmp_path):
    """Survey maps should plot coordinate tables."""
    coord_path = tmp_path / "coordinates.csv"
    coord_path.write_text(
        "\n".join(
            [
                "Profile,point,Gauss Coordinate,,Relative coordinate,,H(m),",
                ",,X(m),Y(m),X(m),Y(m),,",
                "100,100,4291789.0,19510112.0,100.0,100.0,1102.0,",
                "100,120,4291809.0,19510120.0,120.0,100.5,1103.0,",
                "100,140,4291829.0,19510140.0,140.0,101.0,1101.0,",
                "200,100,4291889.0,19510212.0,100.0,200.0,1095.0,",
            ]
        ),
        encoding="utf-8",
    )
    survey = read_temavg_survey(
        DATA_DIR,
        pattern="TEM100.AVG",
        coordinate_file=coord_path,
    )

    ax = PlotSurveyMap(survey).plot()
    ax2 = plot_survey_map(
        survey.coordinates,
        contour=True,
        title=None,
    )

    assert ax.collections
    assert ax2.collections
    assert ax2.get_title() == ""
    ax3 = PlotElevationProfile(
        survey.coordinates,
        profiles=100,
        x="distance",
    ).plot()
    ax4 = plot_elevation_profile(
        survey.coordinates,
        profiles=[100, 200],
        x="point",
        title=None,
    )
    assert ax3.lines
    assert ax4.lines
    assert ax3.xaxis.get_label_position() == "top"
    fig = PlotSurveyOverview(
        survey.coordinates,
        profile=100,
        title=None,
    ).plot()
    fig2 = plot_survey_overview(
        survey.coordinates,
        profile=100,
        title=None,
    )
    assert len(fig.axes) >= 2
    assert len(fig2.axes) >= 2
    plt.close(ax.figure)
    plt.close(ax2.figure)
    plt.close(ax3.figure)
    plt.close(ax4.figure)
    plt.close(fig)
    plt.close(fig2)


@pytest.mark.skipif(
    not AVG_FILE.exists(),
    reason="TEMAVG sample data not available",
)
def test_plot_gate_profile_class_and_function():
    """Gate-profile helpers should draw selected windows."""
    avg = TEMAVG.read(AVG_FILE)
    ticks = StationTickConfig(every=10, rotation=0)
    ax = PlotGateProfile(
        avg,
        windows=[1, 5, 10],
        station_ticks=ticks,
    ).plot()
    ax2 = plot_gate_profile(avg, windows=[1, 8], log_y=False)

    assert len(ax.lines) == 3
    assert len(ax2.lines) == 2
    assert ax.get_xticks().size > 0
    assert ax.xaxis.get_label_position() == "top"
    plt.close(ax.figure)
    plt.close(ax2.figure)


@pytest.mark.skipif(
    not (AVG_FILE.exists() and Z_FILE.exists()),
    reason="TEMAVG sample data not available",
)
def test_plot_tem_dashboard_class_and_function():
    """Dashboard helper should combine decay, transform, and sections."""
    avg = TEMAVG.read(AVG_FILE)
    zplot = TEMZPlot.read(Z_FILE)
    survey = read_temavg_survey(DATA_DIR, pattern="TEM100.AVG")
    soundings = survey.to_soundings(stems=["TEM100"])

    fig = PlotTEMDashboard(avg, zplot, soundings).plot()
    fig2 = plot_tem_dashboard(avg, zplot, soundings[:4])

    assert len(fig.axes) >= 4
    assert len(fig2.axes) >= 4
    plt.close(fig)
    plt.close(fig2)
