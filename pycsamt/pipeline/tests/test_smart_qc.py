"""Tests for pycsamt.pipeline._smart_qc — data-driven QC gating.

Uses the same minimal EDI-like test-double pattern already established in
pycsamt/emtools/tests/test_emtools_fieldzone.py (a small ``_FakeZ``/
``_FakeSite`` pair that ``_get_z_block``/``_get_t_block`` understand),
self-contained in this file per this test suite's convention.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

from pycsamt.pipeline._smart_qc import (
    _any_multicomponent,
    _any_tipper,
    induction_multiperiod_map_smart,
    induction_section_smart,
    phase_tensor_smart,
    response_tipper_smart,
    tipper_components_smart,
    tipper_hodograms_smart,
)

# ─────────────────────────────────────────────────────────────────────────────
# Test doubles
# ─────────────────────────────────────────────────────────────────────────────


class _FakeZ:
    def __init__(self, z, freq):
        self.z = np.asarray(z, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)


class _FakeTipper:
    def __init__(self, t, freq):
        self.tipper = np.asarray(t, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)


class _FakeSite:
    """Minimal EDI-like object understood by _get_z_block/_get_t_block."""

    def __init__(self, station, z, freq, tipper=None):
        self.station = station
        self.Z = _FakeZ(z, freq)
        self.freq = np.asarray(freq, dtype=float)
        if tipper is not None:
            self.Tipper = _FakeTipper(tipper, freq)

    def get_section(self, *_, **__):  # required by is_edi_file() duck-type
        return None


def _make_z(n_freq=6, rho=100.0, *, only=None):
    """Synthetic (n,2,2) impedance tensor; *only* forces single-component."""
    freqs = np.logspace(-2, 2, n_freq)
    z_abs = np.sqrt(5.0 * freqs * rho)
    z = np.zeros((n_freq, 2, 2), dtype=complex)
    val = z_abs * (1.0 + 1j) / np.sqrt(2)
    if only in (None, "xy"):
        z[:, 0, 1] = val
    if only in (None, "yx"):
        z[:, 1, 0] = -val
    return z, freqs


def _multicomponent_site(station="S01") -> _FakeSite:
    z, fr = _make_z()
    return _FakeSite(station, z, fr)


def _single_component_site(station="S01") -> _FakeSite:
    z, fr = _make_z(only="xy")
    return _FakeSite(station, z, fr)


def _with_tipper_site(station="S01") -> _FakeSite:
    z, fr = _make_z()
    t = np.full((len(fr), 2), 0.1 + 0.05j, dtype=complex)
    return _FakeSite(station, z, fr, tipper=t)


def _no_tipper_site(station="S01") -> _FakeSite:
    z, fr = _make_z()
    return _FakeSite(station, z, fr)


# ─────────────────────────────────────────────────────────────────────────────
# Detection helpers
# ─────────────────────────────────────────────────────────────────────────────


class TestAnyMulticomponent:
    def test_true_when_both_components_present(self):
        assert _any_multicomponent([_multicomponent_site()]) is True

    def test_false_when_single_component_only(self):
        assert _any_multicomponent([_single_component_site()]) is False

    def test_true_if_any_one_of_several_is_multicomponent(self):
        sites = [_single_component_site("A"), _multicomponent_site("B")]
        assert _any_multicomponent(sites) is True

    def test_false_for_empty_collection(self):
        assert _any_multicomponent([]) is False


class TestAnyTipper:
    def test_true_when_tipper_present(self):
        assert _any_tipper([_with_tipper_site()]) is True

    def test_false_when_tipper_absent(self):
        assert _any_tipper([_no_tipper_site()]) is False

    def test_true_if_any_one_of_several_has_tipper(self):
        sites = [_no_tipper_site("A"), _with_tipper_site("B")]
        assert _any_tipper(sites) is True

    def test_false_for_empty_collection(self):
        assert _any_tipper([]) is False


# ─────────────────────────────────────────────────────────────────────────────
# Gated wrappers
# ─────────────────────────────────────────────────────────────────────────────


class TestPhaseTensorSmart:
    def test_returns_none_for_single_component(self):
        assert phase_tensor_smart([_single_component_site()]) is None

    def test_returns_figure_for_multicomponent(self):
        result = phase_tensor_smart([_multicomponent_site()])
        assert result is not None


class TestTipperGatedWrappers:
    @pytest.mark.parametrize(
        "fn",
        [
            induction_multiperiod_map_smart,
            induction_section_smart,
            response_tipper_smart,
            tipper_components_smart,
            tipper_hodograms_smart,
        ],
    )
    def test_returns_none_without_tipper(self, fn):
        assert fn([_no_tipper_site()]) is None

    @pytest.mark.parametrize(
        "fn",
        [
            induction_multiperiod_map_smart,
            induction_section_smart,
            response_tipper_smart,
            tipper_components_smart,
            tipper_hodograms_smart,
        ],
    )
    def test_returns_something_with_tipper(self, fn):
        result = fn([_with_tipper_site()])
        assert result is not None
