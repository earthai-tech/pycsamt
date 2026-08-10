"""Tests for package-wide contour configuration."""

from __future__ import annotations

from pycsamt.api.contour import PYCSAMT_CONTOUR, reset_contour


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
