"""Tests for pycsamt.pipeline._preview — deterministic raw-vs-processed preview."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

from pycsamt.pipeline._preview import (
    pick_preview_stations,
    plot_processed_preview,
    plot_raw_preview,
)

# ─────────────────────────────────────────────────────────────────────────────
# Test doubles
# ─────────────────────────────────────────────────────────────────────────────


class _FakeZ:
    def __init__(self, z, freq):
        self.z = np.asarray(z, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)


class _FakeSite:
    def __init__(self, station):
        freqs = np.logspace(-2, 2, 6)
        z_abs = np.sqrt(5.0 * freqs * 100.0)
        z = np.zeros((6, 2, 2), dtype=complex)
        z[:, 0, 1] = z_abs * (1.0 + 1j) / np.sqrt(2)
        z[:, 1, 0] = -z_abs * (1.0 + 1j) / np.sqrt(2)
        self.station = station
        self.Z = _FakeZ(z, freqs)
        self.freq = freqs

    def get_section(self, *_, **__):
        return None


def _sites(n: int) -> list[_FakeSite]:
    return [_FakeSite(f"S{i:02d}") for i in range(n)]


# ─────────────────────────────────────────────────────────────────────────────
# pick_preview_stations
# ─────────────────────────────────────────────────────────────────────────────


class TestPickPreviewStations:
    def test_returns_n_stations(self):
        picked = pick_preview_stations(_sites(10), n=3, random_state=0)
        assert len(picked) == 3

    def test_caps_at_available_stations(self):
        picked = pick_preview_stations(_sites(2), n=3, random_state=0)
        assert len(picked) == 2

    def test_empty_collection_returns_empty(self):
        assert pick_preview_stations([], n=3, random_state=0) == []

    def test_deterministic_for_same_seed(self):
        sites = _sites(10)
        a = pick_preview_stations(sites, n=3, random_state=42)
        b = pick_preview_stations(sites, n=3, random_state=42)
        assert a == b

    def test_different_seeds_can_differ(self):
        sites = _sites(20)
        a = pick_preview_stations(sites, n=3, random_state=0)
        b = pick_preview_stations(sites, n=3, random_state=1)
        assert a != b

    def test_picked_names_are_real_station_names(self):
        sites = _sites(5)
        picked = pick_preview_stations(sites, n=3, random_state=0)
        real_names = {s.station for s in sites}
        assert set(picked).issubset(real_names)

    def test_default_n_is_three(self):
        picked = pick_preview_stations(_sites(10), random_state=0)
        assert len(picked) == 3


# ─────────────────────────────────────────────────────────────────────────────
# plot_raw_preview / plot_processed_preview
# ─────────────────────────────────────────────────────────────────────────────


class TestPreviewPlots:
    def test_raw_preview_returns_something(self):
        assert plot_raw_preview(_sites(5)) is not None

    def test_processed_preview_returns_something(self):
        assert plot_processed_preview(_sites(5)) is not None

    def test_empty_sites_does_not_raise(self):
        # plot_raw_sites_1d already handles "no stations" gracefully.
        result = plot_raw_preview([])
        assert result is not None

    def test_before_and_after_use_same_default_stations(self):
        # Same seed, same station set => same station subset selected
        # independently by each call, matching the documented preview
        # contract (no cross-step state channel).
        sites = _sites(10)
        before = pick_preview_stations(sites, random_state=0)
        after = pick_preview_stations(sites, random_state=0)
        assert before == after
