"""Tests for the package-wide site ordering policy."""

from __future__ import annotations

from pathlib import Path

import pytest

from pycsamt.api import (
    PYCSAMT_ORDERING,
    SiteOrderingConfig,
    configure_ordering,
    reset_ordering,
)
from pycsamt.emtools._core import ensure_sites


@pytest.fixture(autouse=True)
def _reset_global_ordering():
    reset_ordering()
    yield
    reset_ordering()


@pytest.fixture
def l22() -> Path:
    path = Path(__file__).parents[3] / "data" / "AMT" / "WILLY_DATA" / "L22PLT"
    if not path.exists():
        pytest.skip("L22PLT integration data is not available")
    return path


def test_default_policy_is_auto() -> None:
    assert PYCSAMT_ORDERING.mode == "auto"
    assert PYCSAMT_ORDERING.min_linearity == pytest.approx(0.95)


def test_configure_ordering_changes_singleton() -> None:
    returned = configure_ordering(mode="input", min_linearity=0.99)
    assert returned is PYCSAMT_ORDERING
    assert PYCSAMT_ORDERING.mode == "input"
    assert PYCSAMT_ORDERING.min_linearity == pytest.approx(0.99)


def test_aliases_are_normalized() -> None:
    configure_ordering(mode="spatial")
    assert PYCSAMT_ORDERING.mode == "chainage"
    configure_ordering(mode="name")
    assert PYCSAMT_ORDERING.mode == "station"


def test_invalid_configuration_is_atomic() -> None:
    before = PYCSAMT_ORDERING.clone()
    with pytest.raises(ValueError, match="mode must be one of"):
        configure_ordering(mode="mystery", min_linearity=0.2)
    assert PYCSAMT_ORDERING == before


@pytest.mark.parametrize(
    "key,value",
    [
        ("min_linearity", -0.1),
        ("max_cross_track_ratio", 1.1),
        ("min_coordinate_fraction", 2.0),
    ],
)
def test_thresholds_must_be_probabilities(key: str, value: float) -> None:
    with pytest.raises(ValueError, match="between 0 and 1"):
        configure_ordering(**{key: value})


def test_unknown_configuration_key_is_rejected() -> None:
    with pytest.raises(AttributeError, match="Unknown ordering config"):
        configure_ordering(unknown=True)


def test_context_temporarily_overrides_and_restores() -> None:
    configure_ordering(mode="auto", min_linearity=0.95)
    with PYCSAMT_ORDERING.context(mode="input", min_linearity=0.8):
        assert PYCSAMT_ORDERING.mode == "input"
        assert PYCSAMT_ORDERING.min_linearity == pytest.approx(0.8)
    assert PYCSAMT_ORDERING.mode == "auto"
    assert PYCSAMT_ORDERING.min_linearity == pytest.approx(0.95)


def test_context_restores_after_exception() -> None:
    with pytest.raises(RuntimeError):
        with PYCSAMT_ORDERING.context(mode="station"):
            raise RuntimeError("stop")
    assert PYCSAMT_ORDERING.mode == "auto"


def test_reset_restores_all_defaults() -> None:
    configure_ordering(
        mode="input",
        min_linearity=0.7,
        max_cross_track_ratio=0.4,
        min_coordinate_fraction=0.2,
    )
    reset_ordering()
    assert PYCSAMT_ORDERING == SiteOrderingConfig()


def test_ensure_sites_uses_global_policy_when_unspecified(l22: Path) -> None:
    configure_ordering(mode="input")
    sites = ensure_sites(l22)
    assert sites.ordering["requested"] == "input"
    assert [site.name for site in sites][:3] == [
        "22-013VF",
        "22-025AF",
        "22-10U",
    ]


def test_explicit_per_call_policy_overrides_global(l22: Path) -> None:
    configure_ordering(mode="input")
    sites = ensure_sites(l22, order_by="auto")
    assert sites.ordering["applied"] == "chainage"
    assert [site.name for site in sites][:3] == ["22-1BF", "22-2VF", "22-3U"]


def test_sites_ordered_without_argument_uses_global_policy(l22: Path) -> None:
    raw = ensure_sites(l22, order_by="input")
    configure_ordering(mode="auto")
    ordered = raw.ordered()
    assert ordered.ordering["applied"] == "chainage"
    assert [site.name for site in ordered][-1] == "22-025AF"


def test_global_thresholds_control_auto_acceptance(l22: Path) -> None:
    configure_ordering(mode="auto", min_linearity=1.0)
    sites = ensure_sites(l22)
    assert sites.ordering["applied"] == "input"
    assert "credible straight line" in sites.ordering["reason"]
