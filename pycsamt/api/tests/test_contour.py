"""Tests for package-wide contour configuration."""

from __future__ import annotations

import pytest

from pycsamt.api.contour import (
    PYCSAMT_CONTOUR,
    ContourStyle,
    configure_contour,
    reset_contour,
    use_contour,
)


def teardown_function():
    reset_contour()


def test_review_style_is_the_default():
    style = PYCSAMT_CONTOUR.default
    assert style.enabled is True
    assert style.levels == 7
    assert style.linewidths == 0.8
    assert style.alpha == 0.8


def test_configure_and_resolve_overrides():
    PYCSAMT_CONTOUR.configure(linewidths=1.2, colors="white")
    enabled, kwargs, labels = PYCSAMT_CONTOUR.resolve(
        True, {"levels": 4}
    )
    assert enabled is True
    assert kwargs["levels"] == 4
    assert kwargs["linewidths"] == 1.2
    assert kwargs["colors"] == "white"
    assert labels == {}


def test_context_restores_style():
    before = PYCSAMT_CONTOUR.default.linewidths
    with PYCSAMT_CONTOUR.context("publication", linewidths=1.5):
        assert PYCSAMT_CONTOUR.default.labels is True
        assert PYCSAMT_CONTOUR.default.linewidths == 1.5
    assert PYCSAMT_CONTOUR.default.labels is False
    assert PYCSAMT_CONTOUR.default.linewidths == before


def test_off_preset_disables_default():
    PYCSAMT_CONTOUR.use_preset("off")
    enabled, _, _ = PYCSAMT_CONTOUR.resolve()
    assert enabled is False


def test_contour_style_label_kwargs():
    style = ContourStyle(label_fmt="%.1f", label_fontsize=9.0, label_inline=False)
    kw = style.label_kwargs()
    assert kw == {"fmt": "%.1f", "fontsize": 9.0, "inline": False}


def test_style_for_unknown_preset_raises():
    with pytest.raises(ValueError, match="contour preset must be one of"):
        PYCSAMT_CONTOUR.style_for("bogus")


def test_use_preset_default_is_a_noop():
    PYCSAMT_CONTOUR.configure(linewidths=42.0)
    PYCSAMT_CONTOUR.use_preset("default")  # must not overwrite the tweak
    assert PYCSAMT_CONTOUR.default.linewidths == 42.0


def test_configure_preset_prefixed_path_and_unknown_key():
    PYCSAMT_CONTOUR.configure(review__linewidths=2.5)
    assert PYCSAMT_CONTOUR.review.linewidths == 2.5
    with pytest.raises(AttributeError, match="unknown contour setting"):
        PYCSAMT_CONTOUR.configure(bogus_setting=1)


def test_context_without_preset_or_kwargs_is_a_noop_block():
    before = PYCSAMT_CONTOUR.default.linewidths
    with PYCSAMT_CONTOUR.context():
        assert PYCSAMT_CONTOUR.default.linewidths == before
    assert PYCSAMT_CONTOUR.default.linewidths == before


def test_summary_and_repr_describe_the_live_style():
    text = PYCSAMT_CONTOUR.summary()
    assert "PyCSAMTContour" in text
    assert "enabled" in text
    assert repr(PYCSAMT_CONTOUR) == text


def test_module_level_configure_and_use_contour():
    use_contour("publication")
    assert PYCSAMT_CONTOUR.default.labels is True
    configure_contour(linewidths=3.3)
    assert PYCSAMT_CONTOUR.default.linewidths == 3.3
