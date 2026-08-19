"""Tests for the reusable Occam1D plotting-style API."""

import pytest

from pycsamt.api.occam1d import (
    PYCSAMT_OCCAM1D,
    Occam1DPlotStyle,
    resolve_occam1d_style,
)


def test_default_hides_model_legend_and_exposes_all_artists():
    style = resolve_occam1d_style()

    assert style.model_legend is False
    assert style.response_legend is True
    assert style.observed.marker == "o"
    assert style.predicted.visible is True
    assert style.iteration.visible is True
    assert style.target.visible is True


def test_copy_applies_nested_overrides_without_mutating_source():
    original = Occam1DPlotStyle()
    changed = original.copy(
        observed__marker="x",
        predicted__color="purple",
        model__linewidth=3.5,
        target__visible=False,
    )

    assert changed.observed.marker == "x"
    assert changed.predicted.color == "purple"
    assert changed.model.linewidth == 3.5
    assert changed.target.visible is False
    assert original.observed.marker == "o"
    assert original.target.visible is True


def test_context_restores_live_default():
    before = PYCSAMT_OCCAM1D.default.observed.marker
    with PYCSAMT_OCCAM1D.context(observed__marker="^"):
        assert PYCSAMT_OCCAM1D.default.observed.marker == "^"
    assert PYCSAMT_OCCAM1D.default.observed.marker == before


def test_unknown_preset_and_override_fail_clearly():
    with pytest.raises(ValueError, match="Unknown Occam1D style preset"):
        resolve_occam1d_style("missing")
    with pytest.raises(ValueError, match="style path"):
        Occam1DPlotStyle().copy(model__unknown=1)
