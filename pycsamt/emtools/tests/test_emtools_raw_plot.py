"""Tests for raw 1-D station plotting."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from pycsamt.api.control import PYCSAMT_CONTROL
from pycsamt.emtools import plot_raw_sites_1d


class _FakeZ:
    """Small transfer-function container for plotting tests."""

    def __init__(self, z, freq, z_err=None):
        self.z = np.asarray(z, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)
        self.z_err = z_err


class _FakeSite:
    """Small site object accepted by emtools core helpers."""

    def __init__(self, station, z, freq, z_err=None):
        self.station = station
        self.Z = _FakeZ(z, freq, z_err=z_err)
        self.freq = np.asarray(freq, dtype=float)


def _freqs(n: int = 6) -> np.ndarray:
    """Return a compact frequency vector."""
    return np.logspace(-2, 2, n)


def _site(station: str = "S00", rho: float = 100.0) -> _FakeSite:
    """Return a synthetic full-tensor site."""
    freq = _freqs()
    amp = np.sqrt(5.0 * freq * rho)
    z = np.zeros((freq.size, 2, 2), dtype=complex)
    z[:, 0, 0] = 0.5 * amp * np.exp(1j * np.deg2rad(20.0))
    z[:, 0, 1] = amp * np.exp(1j * np.deg2rad(45.0))
    z[:, 1, 0] = amp * np.exp(1j * np.deg2rad(-135.0))
    z[:, 1, 1] = 0.4 * amp * np.exp(1j * np.deg2rad(-30.0))
    z_err = np.full_like(z.real, 0.01)
    return _FakeSite(station, z, freq, z_err=z_err)


def test_plot_raw_sites_1d_uses_black_raw_style():
    """Raw 1-D panels should use the global raw-data style."""
    fig = plot_raw_sites_1d(
        [_site("S00")],
        components=("xy", "yx"),
        show_component_legend=False,
    )

    assert len(fig.axes) == 4
    assert fig.axes[0].lines[0].get_color() == "black"
    assert fig.axes[1].lines[0].get_color() == "black"
    plt.close(fig)


def test_plot_raw_sites_1d_can_force_component_style():
    """Processed/component style can be forced explicitly."""
    fig = plot_raw_sites_1d(
        [_site("S00")],
        components=("xy",),
        raw=True,
        force_style=True,
        show_component_legend=False,
    )

    assert fig.axes[0].lines[0].get_color() != "black"
    plt.close(fig)


def test_plot_raw_sites_1d_uses_control_context():
    """Global control should set axis labels and phase range."""
    with PYCSAMT_CONTROL.context(
        phase__range=(-90.0, 90.0),
        x__view="frequency",
        rho__view="linear",
    ):
        fig = plot_raw_sites_1d(
            [_site("S00")],
            components=("xy",),
            show_component_legend=False,
        )

    text = [item.get_text() for item in fig.texts]
    assert r"$\rho_a$ ($\Omega\,\mathrm{m}$)" in text
    assert "Freq (Hz)" in text
    assert tuple(fig.axes[1].get_ylim()) == (-90.0, 90.0)
    assert fig.axes[1].get_xscale() == "log"
    plt.close(fig)


def test_plot_raw_sites_1d_uses_one_shared_x_label_per_group():
    """Repeated component x labels should be replaced by one group label."""
    fig = plot_raw_sites_1d(
        [_site("S00")],
        components=("xx", "xy", "yx", "yy"),
        show_component_legend=False,
    )

    axis_xlabels = [ax.get_xlabel() for ax in fig.axes]
    figure_labels = [item.get_text() for item in fig.texts]

    assert axis_xlabels.count(r"$\log_{10}T$ (s)") == 0
    assert figure_labels.count(r"$\log_{10}T$ (s)") == 1
    assert figure_labels.count("S00") == 1
    visible_tick_labels = [
        label.get_text()
        for label in fig.axes[4].get_xticklabels()
        if label.get_visible() and label.get_text()
    ]
    assert visible_tick_labels
    plt.close(fig)


def test_plot_raw_sites_1d_axis_label_mode_keeps_old_behavior():
    """Axis label mode should keep repeated labels with rotated ticks."""
    fig = plot_raw_sites_1d(
        [_site("S00")],
        components=("xx", "xy"),
        label_mode="axis",
        show_component_legend=False,
    )

    axis_xlabels = [ax.get_xlabel() for ax in fig.axes]
    rotations = [
        label.get_rotation()
        for label in fig.axes[2].get_xticklabels()
        if label.get_visible()
    ]

    assert axis_xlabels.count(r"$\log_{10}T$ (s)") == 2
    assert fig.axes[0].get_title() == "S00"
    assert rotations and all(rotation == 90.0 for rotation in rotations)
    plt.close(fig)
