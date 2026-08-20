"""Tests for pycsamt.emtools.ztem"""

from __future__ import annotations

import numpy as np
import pytest

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.emtools.ztem import (
    mask_outside_ztem_band,
    phase_rotate_table,
    plot_ztem_band_mask_psection,
    plot_ztem_divergence_profile,
    plot_ztem_divergence_psection,
    plot_ztem_phase_rotation_profile,
    plot_ztem_tipper_profile,
    total_divergence_table,
    ztem_crossover_diagnostics,
)

# ─────────────────────────────────────────────────────────────────────────
# Shared fixture (mirrors pycsamt/emtools/tests/test_emtools_afmag.py)
# ─────────────────────────────────────────────────────────────────────────


class _FakeZ:
    def __init__(self, z, freq):
        self.z = np.asarray(z, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)


class _FakeTipper:
    def __init__(self, tipper, freq):
        self.tipper = np.asarray(tipper, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)


class _FakeSite:
    def __init__(
        self, station, z, freq, *, tipper=None, east=None, north=None
    ):
        self.station = station
        # ensure_sites/to_sites only recognizes an item as EDI-like once
        # it has a ``Z`` attribute, even for a tipper-only test.
        self.Z = _FakeZ(z, freq)
        self.freq = np.asarray(freq, dtype=float)
        if tipper is not None:
            self.Tipper = _FakeTipper(tipper, freq)
        if east is not None:
            self.east = float(east)
            self.north = float(north)

    def get_section(self, *_, **__):
        return None


_ZTEM_FREQS = np.array([30.0, 45.0, 90.0, 180.0, 360.0, 720.0])


def _tipper_crossover(x: float, freqs: np.ndarray) -> np.ndarray:
    """Synthetic Tzx/Tzy with a real-part crossover over a localized
    conductor at x=0 -- decays back to ~0 at the profile edges (unlike
    a saturating tanh) so the FFT-based Hilbert transform used by
    ``phase_rotate_table`` is not dominated by wrap-around edge
    artifacts.
    """
    t = np.zeros((freqs.size, 2), dtype=complex)
    shape = (x / 150.0) * np.exp(-((x / 200.0) ** 2))
    for k, f in enumerate(freqs):
        decay = 1.0 / np.sqrt(f / 30.0)
        t[k, 0] = shape * decay + 0.05j * decay
        t[k, 1] = 0.3 * shape * decay - 0.02j * decay
    return t


def _site(name: str, x: float, *, with_tipper: bool = True) -> _FakeSite:
    fr = _ZTEM_FREQS
    z = np.zeros((fr.size, 2, 2), dtype=complex)
    tip = _tipper_crossover(x, fr) if with_tipper else None
    return _FakeSite(name, z, fr, tipper=tip, east=x, north=0.0)


def _profile(n_sites: int = 9, spacing: float = 100.0) -> list:
    x0 = -0.5 * spacing * (n_sites - 1)
    return [
        _site(f"S{i:02d}", x0 + i * spacing) for i in range(n_sites)
    ]


# ─────────────────────────────────────────────────────────────────────────
# total_divergence_table
# ─────────────────────────────────────────────────────────────────────────


class TestTotalDivergenceTable:
    def test_table_shape_and_columns(self):
        sites = _profile(6)
        df = total_divergence_table(sites)
        assert not df.empty
        expected_cols = {
            "station_a",
            "station_b",
            "x_m",
            "dx_m",
            "freq_hz",
            "period_s",
            "divergence_real",
            "divergence_imag",
            "divergence_abs",
        }
        assert expected_cols.issubset(df.columns)
        # 6 stations -> 5 adjacent pairs, times 6 frequencies
        assert len(df) == 5 * _ZTEM_FREQS.size

    def test_peaks_near_the_contact(self):
        # A tanh-shaped real crossover at x=0 has its steepest slope
        # (largest |derivative|) at the pair straddling x=0.
        sites = _profile(9, spacing=100.0)
        df = total_divergence_table(sites, component="tzx")
        f0 = float(df["freq_hz"].iloc[0])
        sub = df[df["freq_hz"] == f0].reset_index(drop=True)
        peak_idx = sub["divergence_real"].abs().idxmax()
        assert sub.loc[peak_idx, "station_a"] in {"S03", "S04"}

    def test_invalid_component_raises(self):
        with pytest.raises(ValueError):
            total_divergence_table(_profile(3), component="bogus")

    def test_single_station_returns_empty(self):
        df = total_divergence_table(_profile(1))
        assert df.empty

    def test_no_tipper_returns_empty(self):
        sites = [_site("S00", 0.0, with_tipper=False)]
        df = total_divergence_table(sites)
        assert df.empty

    def test_tzy_component(self):
        sites = _profile(5)
        df = total_divergence_table(sites, component="tzy")
        assert not df.empty


# ─────────────────────────────────────────────────────────────────────────
# phase_rotate_table
# ─────────────────────────────────────────────────────────────────────────


class TestPhaseRotateTable:
    def test_table_shape_and_columns(self):
        sites = _profile(9)
        df = phase_rotate_table(sites, frequency_hz=30.0)
        assert not df.empty
        expected_cols = {
            "x_m",
            "nearest_station",
            "freq_hz",
            "period_s",
            "raw",
            "rotated",
            "envelope",
        }
        assert expected_cols.issubset(df.columns)
        assert np.allclose(df["freq_hz"], 30.0)

    def test_crossover_becomes_a_peak(self):
        # The raw in-phase crossover is ~0 at the contact and the
        # Hilbert-rotated response should peak there instead. x_m is
        # chainage from the first station, so the contact (at physical
        # east=0) sits at chainage = -x0.
        n_sites, spacing = 11, 100.0
        x0 = -0.5 * spacing * (n_sites - 1)
        sites = _profile(n_sites, spacing=spacing)
        df = phase_rotate_table(
            sites, frequency_hz=30.0, n_resample=201,
        )
        contact_chainage = -x0
        i0 = int(np.argmin(np.abs(df["x_m"] - contact_chainage)))
        assert abs(df["raw"].iloc[i0]) < 0.15
        peak_at_center = np.abs(df["rotated"].iloc[i0]) > 0.3 * np.max(
            np.abs(df["rotated"])
        )
        assert peak_at_center

    def test_both_frequency_and_period_raises(self):
        with pytest.raises(ValueError):
            phase_rotate_table(
                _profile(6), frequency_hz=30.0, period_s=1.0,
            )

    def test_invalid_component_raises(self):
        with pytest.raises(ValueError):
            phase_rotate_table(_profile(6), component="bogus")

    def test_invalid_part_raises(self):
        with pytest.raises(ValueError):
            phase_rotate_table(_profile(6), part="bogus")

    def test_too_few_stations_returns_empty(self):
        df = phase_rotate_table(_profile(2))
        assert df.empty


# ─────────────────────────────────────────────────────────────────────────
# mask_outside_ztem_band
# ─────────────────────────────────────────────────────────────────────────


class TestMaskOutsideZtemBand:
    def test_masks_out_of_band_frequencies(self):
        sites = _profile(3)
        out = mask_outside_ztem_band(sites, band_hz=(40.0, 200.0))
        df = total_divergence_table(out)
        out_of_band = df[
            (df["freq_hz"] < 40.0) | (df["freq_hz"] > 200.0)
        ]
        assert out_of_band.empty

    def test_default_band_uses_ztem_system_spec(self):
        sites = _profile(3)
        out = mask_outside_ztem_band(sites)
        df = total_divergence_table(out)
        assert (df["freq_hz"] >= 22.0).all()
        assert (df["freq_hz"] <= 720.0).all()

    def test_band_hz_and_system_spec_mutually_exclusive(self):
        from pycsamt.airborne.ztem import ZTEMSystemSpec

        with pytest.raises(ValueError):
            mask_outside_ztem_band(
                _profile(2),
                band_hz=(30.0, 100.0),
                system_spec=ZTEMSystemSpec(),
            )

    def test_invalid_action_raises(self):
        with pytest.raises(ValueError):
            mask_outside_ztem_band(_profile(2), action="bogus")

    def test_invalid_system_spec_type_raises(self):
        with pytest.raises(TypeError):
            mask_outside_ztem_band(_profile(2), system_spec="not-a-spec")

    def test_drop_action_removes_rows(self):
        sites = _profile(3)
        out = mask_outside_ztem_band(
            sites, band_hz=(40.0, 200.0), action="drop",
        )
        df = total_divergence_table(out)
        assert not (
            (df["freq_hz"] < 40.0) | (df["freq_hz"] > 200.0)
        ).any()

    def test_inplace_false_leaves_input_untouched(self):
        sites = _profile(2)
        before = total_divergence_table(sites)
        mask_outside_ztem_band(
            sites, band_hz=(40.0, 200.0), inplace=False,
        )
        after = total_divergence_table(sites)
        assert np.allclose(
            before["divergence_real"].to_numpy(),
            after["divergence_real"].to_numpy(),
        )


# ─────────────────────────────────────────────────────────────────────────
# Plots
# ─────────────────────────────────────────────────────────────────────────


class TestZtemCrossoverDiagnostics:
    def test_finds_the_synthetic_crossover_near_the_target(self):
        sites = _profile(9, spacing=100.0)
        diag = ztem_crossover_diagnostics(sites, frequency_hz=30.0)
        profile = diag["profile"]
        center_m = float(profile["position_m"].median())
        assert diag["freq_hz"] == pytest.approx(30.0)
        # the real part crosses zero exactly at the profile's centre
        # station (_tipper_crossover's odd shape function); the
        # imaginary part is a constant positive offset at this single
        # frequency (no sign change at all), so it has no crossover.
        assert abs(diag["crossover_real_m"] - center_m) < 1e-6
        assert np.isnan(diag["crossover_imag_m"])
        assert diag["peak_to_peak_real"] > 0.0
        assert not profile.empty

    def test_invalid_component_raises(self):
        with pytest.raises(ValueError):
            ztem_crossover_diagnostics(_profile(6), component="bogus")

    def test_single_station_raises(self):
        with pytest.raises(ValueError):
            ztem_crossover_diagnostics(_profile(1))

    def test_no_tipper_raises(self):
        sites = [_site("S00", 0.0, with_tipper=False)]
        with pytest.raises(ValueError):
            ztem_crossover_diagnostics(sites)


class TestPlots:
    def test_plot_ztem_divergence_profile(self):
        ax = plot_ztem_divergence_profile(_profile(6))
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_plot_ztem_divergence_profile_external_axes(self):
        _, ax_in = plt.subplots()
        ax = plot_ztem_divergence_profile(_profile(5), ax=ax_in)
        assert ax is ax_in
        plt.close("all")

    def test_plot_ztem_divergence_profile_no_tipper(self):
        sites = [_site("S00", 0.0, with_tipper=False)]
        ax = plot_ztem_divergence_profile(sites)
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_plot_ztem_divergence_psection(self):
        ax = plot_ztem_divergence_psection(_profile(7))
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_plot_ztem_phase_rotation_profile(self):
        ax = plot_ztem_phase_rotation_profile(
            _profile(9), frequency_hz=30.0,
        )
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_plot_ztem_band_mask_psection(self):
        fig = plot_ztem_band_mask_psection(
            _profile(6), band_hz=(40.0, 200.0),
        )
        assert isinstance(fig, plt.Figure)
        plt.close("all")

    def test_plot_ztem_tipper_profile(self):
        ax = plot_ztem_tipper_profile(_profile(9), frequency_hz=30.0)
        assert isinstance(ax, plt.Axes)
        assert "in-phase/quadrature" in ax.get_title()
        plt.close("all")

    def test_plot_ztem_tipper_profile_no_tipper(self):
        sites = [_site("S00", 0.0, with_tipper=False)]
        ax = plot_ztem_tipper_profile(sites)
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_plot_ztem_divergence_psection_grid_and_contour(self):
        ax = plot_ztem_divergence_psection(
            _profile(9), show_grid=True, show_contour=True,
        )
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_plot_ztem_divergence_psection_overlays_can_be_disabled(self):
        ax = plot_ztem_divergence_psection(
            _profile(9), show_grid=False, show_contour=False,
        )
        assert isinstance(ax, plt.Axes)
        plt.close("all")
