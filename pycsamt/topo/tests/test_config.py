# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.topo.config — TopoConfig singleton and helpers."""

import contextlib

import pytest

import pycsamt.api.topo as api_topo
from pycsamt.topo.config import (
    PYCSAMT_TOPO,
    TopoConfig,
    configure_topo,
    reset_topo,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _reset_after_each():
    """Always reset the global singleton between tests."""
    yield
    reset_topo()


# ---------------------------------------------------------------------------
# Default values
# ---------------------------------------------------------------------------


class TestTopoConfigDefaults:
    def test_enabled_false(self):
        cfg = TopoConfig()
        assert cfg.enabled is False

    def test_source_sites(self):
        assert TopoConfig().source == "sites"

    def test_exaggeration_one(self):
        assert TopoConfig().exaggeration == 1.0

    def test_interp_method_linear(self):
        assert TopoConfig().interp_method == "linear"

    def test_fill_color(self):
        assert TopoConfig().fill_color == "#a89070"

    def test_fill_alpha(self):
        assert TopoConfig().fill_alpha == pytest.approx(0.40)

    def test_line_color(self):
        assert TopoConfig().line_color == "#6b4e2a"

    def test_line_width(self):
        assert TopoConfig().line_width == pytest.approx(1.2)

    def test_show_surface_line_true(self):
        assert TopoConfig().show_surface_line is True

    def test_clip_below_surface_true(self):
        assert TopoConfig().clip_below_surface is True

    def test_station_pins_at_surface_true(self):
        assert TopoConfig().station_pins_at_surface is True

    def test_show_topo_strip_true(self):
        assert TopoConfig().show_topo_strip is True

    def test_strip_height_ratio(self):
        assert TopoConfig().strip_height_ratio == pytest.approx(0.18)

    def test_elev_array_none(self):
        assert TopoConfig().elev_array is None

    def test_elev_file_none(self):
        assert TopoConfig().elev_file is None


# ---------------------------------------------------------------------------
# configure() method
# ---------------------------------------------------------------------------


class TestConfigure:
    def test_enable(self):
        cfg = TopoConfig()
        cfg.configure(enabled=True)
        assert cfg.enabled is True

    def test_partial_update_leaves_others(self):
        cfg = TopoConfig()
        cfg.configure(exaggeration=3.5)
        assert cfg.exaggeration == pytest.approx(3.5)
        assert cfg.enabled is False  # unchanged

    def test_update_multiple_fields(self):
        cfg = TopoConfig()
        cfg.configure(enabled=True, exaggeration=2.0, fill_alpha=0.6)
        assert cfg.enabled is True
        assert cfg.exaggeration == pytest.approx(2.0)
        assert cfg.fill_alpha == pytest.approx(0.6)

    def test_update_string_field(self):
        cfg = TopoConfig()
        cfg.configure(interp_method="cubic")
        assert cfg.interp_method == "cubic"

    def test_unknown_kwarg_raises(self):
        cfg = TopoConfig()
        with pytest.raises((TypeError, AttributeError, ValueError)):
            cfg.configure(nonexistent_field=99)

    def test_configure_returns_self(self):
        cfg = TopoConfig()
        ret = cfg.configure(enabled=True)
        assert ret is cfg


# ---------------------------------------------------------------------------
# reset() method
# ---------------------------------------------------------------------------


class TestReset:
    def test_reset_enabled(self):
        cfg = TopoConfig()
        cfg.configure(enabled=True)
        cfg.reset()
        assert cfg.enabled is False

    def test_reset_exaggeration(self):
        cfg = TopoConfig()
        cfg.configure(exaggeration=5.0)
        cfg.reset()
        assert cfg.exaggeration == pytest.approx(1.0)

    def test_reset_multiple_fields(self):
        cfg = TopoConfig()
        cfg.configure(enabled=True, fill_alpha=0.9, exaggeration=4.0)
        cfg.reset()
        assert cfg.enabled is False
        assert cfg.fill_alpha == pytest.approx(0.40)
        assert cfg.exaggeration == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# context() context manager
# ---------------------------------------------------------------------------


class TestContextManager:
    def test_context_restores_enabled(self):
        cfg = TopoConfig()
        cfg.configure(enabled=True)
        with cfg.context(enabled=False):
            assert cfg.enabled is False
        assert cfg.enabled is True

    def test_context_restores_exaggeration(self):
        cfg = TopoConfig()
        with cfg.context(exaggeration=7.0):
            assert cfg.exaggeration == pytest.approx(7.0)
        assert cfg.exaggeration == pytest.approx(1.0)

    def test_context_restores_on_exception(self):
        cfg = TopoConfig()
        with contextlib.suppress(RuntimeError):
            with cfg.context(enabled=True):
                raise RuntimeError("boom")
        assert cfg.enabled is False

    def test_context_multiple_fields(self):
        cfg = TopoConfig()
        with cfg.context(enabled=True, exaggeration=2.5, fill_alpha=0.8):
            assert cfg.enabled is True
            assert cfg.exaggeration == pytest.approx(2.5)
        assert cfg.enabled is False
        assert cfg.exaggeration == pytest.approx(1.0)

    def test_nested_contexts(self):
        cfg = TopoConfig()
        with cfg.context(exaggeration=2.0):
            with cfg.context(exaggeration=4.0):
                assert cfg.exaggeration == pytest.approx(4.0)
            assert cfg.exaggeration == pytest.approx(2.0)
        assert cfg.exaggeration == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# clone() method
# ---------------------------------------------------------------------------


class TestClone:
    def test_clone_independent(self):
        cfg = TopoConfig()
        cfg.configure(enabled=True, exaggeration=3.0)
        clone = cfg.clone()
        clone.configure(exaggeration=9.0)
        assert cfg.exaggeration == pytest.approx(3.0)

    def test_clone_copies_values(self):
        cfg = TopoConfig()
        cfg.configure(enabled=True, exaggeration=2.5, fill_alpha=0.7)
        clone = cfg.clone()
        assert clone.enabled is True
        assert clone.exaggeration == pytest.approx(2.5)
        assert clone.fill_alpha == pytest.approx(0.7)


# ---------------------------------------------------------------------------
# Global singleton functions
# ---------------------------------------------------------------------------


class TestGlobalSingleton:
    def test_configure_topo_updates_singleton(self):
        configure_topo(enabled=True, exaggeration=2.0)
        assert PYCSAMT_TOPO.enabled is True
        assert PYCSAMT_TOPO.exaggeration == pytest.approx(2.0)

    def test_reset_topo_restores_defaults(self):
        configure_topo(enabled=True, exaggeration=5.0)
        reset_topo()
        assert PYCSAMT_TOPO.enabled is False
        assert PYCSAMT_TOPO.exaggeration == pytest.approx(1.0)

    def test_singleton_identity(self):
        from pycsamt.topo.config import PYCSAMT_TOPO as T2

        assert PYCSAMT_TOPO is T2


# ---------------------------------------------------------------------------
# API re-exports
# ---------------------------------------------------------------------------


class TestApiReExports:
    def test_api_topo_config_class(self):
        assert api_topo.TopoConfig is TopoConfig

    def test_api_topo_singleton(self):
        assert api_topo.PYCSAMT_TOPO is PYCSAMT_TOPO

    def test_api_topo_configure_fn(self):
        api_topo.configure_topo(exaggeration=3.3)
        assert PYCSAMT_TOPO.exaggeration == pytest.approx(3.3)

    def test_api_topo_reset_fn(self):
        api_topo.configure_topo(enabled=True)
        api_topo.reset_topo()
        assert PYCSAMT_TOPO.enabled is False
