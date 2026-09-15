"""Tests for global section-view controls."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

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


def test_figure_style_resolve_non_dynamic_and_no_colorbar():
    style = SectionFigureStyle(figsize=(5.0, 2.5))
    assert style.resolve() == (5.0, 2.5)

    dynamic = SectionFigureStyle(figsize="dynamic")
    with_cbar = dynamic.resolve(n_stations=40, n_y=5, colorbar=True)
    without_cbar = dynamic.resolve(n_stations=40, n_y=5, colorbar=False)
    assert with_cbar[0] > without_cbar[0]


def test_axis_style_rejects_bad_y_direction_and_aspect():
    ax_style = SectionAxisStyle(y_direction="sideways")
    fig, ax = plt.subplots()
    with pytest.raises(ValueError, match="y_direction must be one of"):
        ax_style.apply(ax)
    plt.close(fig)

    bad_aspect = SectionAxisStyle(aspect="weird")
    fig, ax = plt.subplots()
    with pytest.raises(ValueError, match="aspect must be"):
        bad_aspect.apply(ax)
    plt.close(fig)


def test_axis_style_applies_label_pads_and_skips_title():
    style = SectionAxisStyle(
        xlabel="X", ylabel="Y", xlabel_pad=5.0, ylabel_pad=6.0, title=False,
    )
    fig, ax = plt.subplots()
    style.apply(ax, title="Should be skipped")
    assert ax.get_xlabel() == "X"
    assert ax.get_ylabel() == "Y"
    assert ax.get_title() == ""
    plt.close(fig)


def test_axis_style_sets_title_when_enabled_and_skips_empty_xlabel():
    style = SectionAxisStyle(xlabel="", ylabel=None, title=True)
    fig, ax = plt.subplots()
    style.apply(ax, title="My Title")
    assert ax.get_title() == "My Title"
    assert ax.get_xlabel() == ""
    plt.close(fig)


def test_axis_style_y_direction_none_and_aspect_none_skip_both():
    style = SectionAxisStyle(y_direction="none", aspect=None)
    fig, ax = plt.subplots()
    was_inverted = ax.yaxis_inverted()
    style.apply(ax)
    assert ax.yaxis_inverted() == was_inverted  # untouched
    plt.close(fig)


def test_axis_style_up_direction_un_inverts_when_needed():
    style = SectionAxisStyle(y_direction="up")
    fig, ax = plt.subplots()
    ax.invert_yaxis()  # start inverted
    style.apply(ax)
    assert not ax.yaxis_inverted()
    plt.close(fig)


def test_colorbar_style_add_returns_none_when_disabled():
    style = SectionColorbarStyle(show=False)
    fig, ax = plt.subplots()
    image = ax.imshow(np.arange(4).reshape(2, 2))
    assert style.add(image, ax) is None
    plt.close(fig)


def test_section_style_topo_active_delegates_to_topo_config(monkeypatch):
    style = SectionStyle()
    import pycsamt.topo.config as topo_config

    monkeypatch.setattr(
        topo_config.PYCSAMT_TOPO, "is_active_for", lambda y_type: True
    )
    assert style.topo_active() is True


def test_section_style_copy_overrides_without_mutating_original():
    style = SectionStyle()
    copied = style.copy(station_preset="inversion")
    assert copied.station_preset == "inversion"
    assert style.station_preset == "pseudosection"


def test_figsize_for_without_labels_uses_default_label_len():
    reset_section()
    style = PYCSAMT_SECTION.style_for("dynamic")
    size = style.figsize_for(n_stations=4, n_y=5, labels=None)
    assert len(size) == 2


def test_style_for_unknown_preset_raises():
    with pytest.raises(ValueError, match="section preset must be one of"):
        PYCSAMT_SECTION.style_for("bogus")


def test_context_with_preset_reverts_after_block():
    reset_section()
    original = PYCSAMT_SECTION.pseudosection.station_preset
    with PYCSAMT_SECTION.context("inversion"):
        assert PYCSAMT_SECTION.pseudosection.station_preset == "inversion"
    assert PYCSAMT_SECTION.pseudosection.station_preset == original


def test_use_preset_copies_into_pseudosection_slot():
    reset_section()
    PYCSAMT_SECTION.use_preset("inversion")
    assert PYCSAMT_SECTION.pseudosection.station_preset == "inversion"
    reset_section()


def test_summary_and_repr_list_every_preset():
    text = PYCSAMT_SECTION.summary()
    assert "PyCSAMTSection" in text
    assert "pseudosection:" in text
    assert repr(PYCSAMT_SECTION) == text
