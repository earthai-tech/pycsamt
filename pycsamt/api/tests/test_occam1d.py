"""Tests for the reusable Occam1D plotting-style API."""

import pytest

from pycsamt.api.occam1d import (
    PYCSAMT_OCCAM1D,
    Occam1DArtistStyle,
    Occam1DPlotStyle,
    configure_occam1d_style,
    reset_occam1d_style,
    resolve_occam1d_style,
    use_occam1d_style,
)


@pytest.fixture(autouse=True)
def _restore_global_occam1d_style():
    yield
    PYCSAMT_OCCAM1D.reset()


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
    with pytest.raises(ValueError, match="style path"):
        Occam1DPlotStyle().copy(unknown_top_level=1)


def test_artist_style_kwargs_omits_none_values():
    artist = Occam1DArtistStyle(marker=None, markerfacecolor=None)
    kw = artist.kwargs()
    assert "marker" not in kw
    assert "color" in kw


def test_artist_style_kwargs_excludes_marker_group_when_requested():
    artist = Occam1DArtistStyle(marker="o", markersize=6.0)
    kw = artist.kwargs(include_marker=False)
    assert "marker" not in kw
    assert "markersize" not in kw
    assert "color" in kw


@pytest.mark.parametrize("preset", ["pycsamt", "publication", "minimal", "default"])
def test_style_for_every_named_preset(preset):
    style = PYCSAMT_OCCAM1D.style_for(preset)
    assert isinstance(style, Occam1DPlotStyle)


def test_publication_and_minimal_presets_apply_their_overrides():
    pub = PYCSAMT_OCCAM1D.style_for("publication")
    assert pub.observed.markerfacecolor == "white"
    assert pub.predicted.color == "0.15"

    minimal = PYCSAMT_OCCAM1D.style_for("minimal")
    assert minimal.response_legend is False
    assert minimal.roughness.visible is False


def test_use_replaces_live_default_and_returns_a_copy():
    returned = PYCSAMT_OCCAM1D.use("publication")
    assert PYCSAMT_OCCAM1D.default.observed.markerfacecolor == "white"
    assert returned is not PYCSAMT_OCCAM1D.default


def test_reset_restores_pycsamt_house_style():
    PYCSAMT_OCCAM1D.use("publication")
    PYCSAMT_OCCAM1D.reset()
    assert PYCSAMT_OCCAM1D.default.observed.markerfacecolor == "black"


def test_context_with_preset_reverts_after_block():
    PYCSAMT_OCCAM1D.reset()
    with PYCSAMT_OCCAM1D.context("minimal") as live:
        assert live.response_legend is False
        assert PYCSAMT_OCCAM1D.default.response_legend is False
    assert PYCSAMT_OCCAM1D.default.response_legend is True


def test_resolve_occam1d_style_accepts_preset_name_and_instance():
    from_preset = resolve_occam1d_style("minimal")
    assert from_preset.response_legend is False

    custom = Occam1DPlotStyle().copy(model__color="teal")
    from_instance = resolve_occam1d_style(custom)
    assert from_instance.model.color == "teal"
    assert from_instance is not custom


def test_resolve_occam1d_style_rejects_bad_type():
    with pytest.raises(TypeError, match="style must be None"):
        resolve_occam1d_style(42)


def test_module_level_configure_use_reset():
    use_occam1d_style("minimal")
    assert PYCSAMT_OCCAM1D.default.response_legend is False
    configure_occam1d_style(convergence_legend=False)
    assert PYCSAMT_OCCAM1D.default.convergence_legend is False
    reset_occam1d_style()
    assert PYCSAMT_OCCAM1D.default.response_legend is True


def test_validate_rejects_bool_and_non_numeric_and_negative_values():
    with pytest.raises(TypeError, match="must be a real number"):
        Occam1DPlotStyle().copy(model__linewidth=True)
    with pytest.raises(TypeError, match="must be a real number"):
        Occam1DPlotStyle().copy(model__linewidth="not-a-number")
    with pytest.raises(ValueError, match="cannot be negative"):
        Occam1DPlotStyle().copy(model__linewidth=-1.0)
    with pytest.raises(TypeError, match="must be a bool"):
        Occam1DPlotStyle().copy(model_legend="not-a-bool")
