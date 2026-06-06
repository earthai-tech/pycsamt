"""Tests for global section-view controls."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from pycsamt.api import (
    PYCSAMT_SECTION,
    SectionAxisStyle,
    SectionColorbarStyle,
    SectionFigureStyle,
    SectionStyle,
    configure_section,
    reset_section,
)


def test_section_presets_are_available():
    """Section presets should expose consistent view defaults."""
    reset_section()

    pseudo = PYCSAMT_SECTION.style_for("pseudosection")
    inversion = PYCSAMT_SECTION.style_for("inversion")
    dynamic = PYCSAMT_SECTION.style_for("dynamic")

    assert pseudo.axis.y_direction == "down"
    assert pseudo.station_preset == "pseudosection"
    assert inversion.station_preset == "inversion"
    assert dynamic.figure.figsize == "dynamic"


def test_section_configure_and_context_restore():
    """Nested section configuration should restore after context."""
    reset_section()
    original = PYCSAMT_SECTION.pseudosection.figure.figsize

    configure_section(pseudosection__figure__figsize=(11.0, 4.0))

    assert PYCSAMT_SECTION.pseudosection.figure.figsize == (11.0, 4.0)
    with PYCSAMT_SECTION.context(
        pseudosection__figure__figsize=(6.0, 2.5),
    ):
        assert PYCSAMT_SECTION.pseudosection.figure.figsize == (6.0, 2.5)

    assert PYCSAMT_SECTION.pseudosection.figure.figsize == (11.0, 4.0)
    reset_section()
    assert PYCSAMT_SECTION.pseudosection.figure.figsize == original


def test_section_style_applies_axis_and_stations():
    """Section styles should apply axes and station rendering."""
    style = SectionStyle(
        figure=SectionFigureStyle(figsize=(4.0, 2.0)),
        axis=SectionAxisStyle(
            xlabel="Station",
            ylabel="Depth (m)",
            y_direction="down",
            station_side="top",
        ),
        colorbar=SectionColorbarStyle(max_ticks=4),
    )
    fig, ax = plt.subplots(figsize=style.figure.figsize)

    style.apply_axis(ax)
    style.apply_stations(ax, [0, 1, 2], labels=["S0", "S1", "S2"])

    assert ax.yaxis_inverted()
    assert ax.get_ylabel() == "Depth (m)"
    assert ax.xaxis.get_label_position() == "top"
    plt.close(fig)


def test_dynamic_section_resolves_data_aware_size():
    """Dynamic sections should size from stations and y samples."""
    reset_section()
    style = PYCSAMT_SECTION.style_for("dynamic")

    small = style.figsize_for(n_stations=4, n_y=5, labels=["S1"])
    large = style.figsize_for(
        n_stations=80,
        n_y=100,
        labels=[f"Station_{i:03d}" for i in range(80)],
    )

    assert small[0] >= style.figure.min_width
    assert small[1] >= style.figure.min_height
    assert large[0] > small[0]
    assert large[1] > small[1]
    assert large[0] <= style.figure.max_width
    assert large[1] <= style.figure.max_height


def test_section_colorbar_uses_smart_ticks():
    """Section colorbars should honor max tick configuration."""
    style = SectionStyle(
        colorbar=SectionColorbarStyle(size="4%", max_ticks=3),
    )
    fig, ax = plt.subplots()
    image = ax.imshow(np.arange(9).reshape(3, 3))

    cbar = style.add_colorbar(image, ax, label="Value")

    assert cbar is not None
    assert cbar.ax.get_ylabel() == "Value"
    assert len(cbar.get_ticks()) <= 5
    plt.close(fig)
