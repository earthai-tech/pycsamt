"""
Tests for:
  - correct_static_shift (remove_noise.py) — Torres-Verdín & Bostick (1992) AMA
  - plot_confidence_profile (qc.py) — Kouadio et al. (2024) Fig. 3 style
"""

from __future__ import annotations

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection

from pycsamt.api import (
    APIFrame,
    APIResult,
    configure_api_view,
    reset_api_view,
)
from pycsamt.compat.numpy import trapz as _trapz
from pycsamt.emtools.dimensionality import (
    pre2d_inversion_assessment,
)
from pycsamt.emtools.frequency import (
    FrequencyEditResult,
    drop_low_confidence_frequencies,
    edit_frequencies_by_confidence,
    frequency_edit_decision_table,
    frequency_edit_report,
    mask_low_confidence_frequencies,
    plot_frequency_edit_decisions,
    plot_frequency_edit_summary,
    recover_low_confidence_frequencies,
)
from pycsamt.emtools.qc import (
    confidence_ratio,
    frequency_confidence_table,
    plot_confidence_band_summary,
    plot_confidence_profile,
    plot_frequency_confidence_psection,
    plot_station_confidence_dashboard,
    plot_station_confidence_spectrum,
    station_confidence_table,
)
from pycsamt.emtools.remove_noise import (
    EMAPFilterResult,
    apply_emap_filter,
    confidence_gated_emap_filter,
    correct_static_shift,
    drop_freqs_manual,
    emap_filter_report,
    emi_mitigation_report,
    fixed_length_moving_average,
    plot_emap_filter_profile,
    plot_emap_filter_psection,
    trimmed_moving_average,
)

# ─────────────────────────────────────────────────────────────────────────────
# Shared fake-site infrastructure
# ─────────────────────────────────────────────────────────────────────────────


class _FakeZ:
    def __init__(self, z, freq):
        self.z = np.array(z, dtype=complex)
        self.freq = np.array(freq, dtype=float)

    def compute_resistivity_phase(self):
        """Mimics the real Z class's marker attribute so tests exercise
        the same "strict container" code path as production (see the
        _set_array_field NaN-write bug this fixture previously masked)."""
        return None


class _FakeSite:
    def __init__(self, station, z, freq, east=None, north=None):
        self.station = station
        self.Z = _FakeZ(z, freq)
        self.freq = np.array(freq, dtype=float)
        if east is not None:
            self.east = float(east)
            self.north = float(north)

    def get_section(self, *_, **__):
        return None


def _freqs(n: int = 10, f_lo: float = 1.0, f_hi: float = 1e4) -> np.ndarray:
    return np.logspace(np.log10(f_lo), np.log10(f_hi), n)


def _make_z(freqs: np.ndarray, rho: float = 100.0) -> np.ndarray:
    """1-D isotropic Z: Z_xy = -Z_yx, |Z| = sqrt(5 f ρ), diag = 0."""
    z_abs = np.sqrt(5.0 * freqs * rho)
    z = np.zeros((freqs.size, 2, 2), dtype=complex)
    z[:, 0, 1] = z_abs * (1 + 1j) / np.sqrt(2)
    z[:, 1, 0] = -z_abs * (1 + 1j) / np.sqrt(2)
    return z


def _rho_from_z(freqs: np.ndarray, rho: float = 100.0) -> np.ndarray:
    """True ρ_a (det) from the synthetic Z above."""
    z_abs_sq = 5.0 * freqs * rho
    return 0.2 * z_abs_sq / freqs


def _site(
    station: str, rho: float = 100.0, n: int = 10, east=None, north=None
) -> _FakeSite:
    fr = _freqs(n)
    return _FakeSite(station, _make_z(fr, rho), fr, east=east, north=north)


def _shifted_site(
    station: str,
    rho: float = 100.0,
    ss: float = 4.0,
    n: int = 10,
    east=None,
    north=None,
) -> _FakeSite:
    """Z with static shift C in ρ_a: Z_shifted = Z × sqrt(C)."""
    fr = _freqs(n)
    z = _make_z(fr, rho)
    z_s = z * np.sqrt(ss)  # ρ_a_shifted = ss × ρ_a_true
    return _FakeSite(station, z_s, fr, east=east, north=north)


def _line_sites(n: int = 5, rho: float = 100.0, spacing: float = 200.0) -> list:
    """n stations along a line with identical ρ_a (no static shift)."""
    return [_site(f"S{i:02d}", rho, east=i * spacing, north=0.0) for i in range(n)]


def test_confidence_ratio_weighted_formula_and_error():
    scores = {
        "coverage": 1.0,
        "uncertainty": 0.8,
        "offdiag": np.nan,
        "diagonal": 0.6,
    }
    weights = {
        "coverage": 0.5,
        "uncertainty": 0.25,
        "diagonal": 0.25,
    }

    cr, err = confidence_ratio(
        scores,
        weights=weights,
        n_freq=12,
        return_error=True,
    )

    assert cr == pytest.approx(0.85)
    assert err == pytest.approx(np.std([1.0, 0.8, 0.6]))


def test_confidence_ratio_single_score_binomial_error():
    cr, err = confidence_ratio(
        {"coverage": 0.75},
        weights={"coverage": 1.0},
        n_freq=12,
        return_error=True,
    )

    assert cr == pytest.approx(0.75)
    assert err == pytest.approx(np.sqrt(0.75 * 0.25 / 12.0))


# ─────────────────────────────────────────────────────────────────────────────
# correct_static_shift — basic contract
# ─────────────────────────────────────────────────────────────────────────────


def test_emi_mitigation_report_records_remote_reference_and_harmonics():
    fr = np.array([10.0, 50.0, 100.0, 250.0])
    site = _FakeSite("S0", _make_z(fr), fr)
    report = emi_mitigation_report(
        [site],
        remote_reference_attempted=False,
        remote_reference_reason="no independent remote site was acquired",
        mains_hz=50.0,
        n_harm=3,
        tol_hz=0.01,
        notch_mode="interp",
    )

    assert list(report["station"]) == ["S0"]
    assert bool(report.loc[0, "remote_reference_attempted"]) is False
    assert bool(report.loc[0, "remote_reference_available"]) is False
    assert report.loc[0, "harmonic_z_samples"] == 2
    assert "notch_powerline" in report.loc[0, "applied_measures"]


def test_pre2d_inversion_assessment_records_strike_and_gb_status(monkeypatch):
    import pandas as pd

    from pycsamt.emtools import dimensionality as dim_mod

    site = _site("S0")
    dim = pd.DataFrame(
        {
            "station": ["S0", "S0", "S0", "S0"],
            "period": [0.1, 1.0, 10.0, 100.0],
            "dim": [0, 1, 2, 2],
            "beta_abs": [1.0, 2.0, 8.0, 10.0],
            "ellipt_abs": [0.05, 0.3, 0.4, 0.5],
        }
    )
    strike = pd.DataFrame(
        {
            "station": ["S0"],
            "ang": [35.0],
            "iqr": [12.0],
            "lo": [0.1],
            "hi": [100.0],
            "n": [4],
        }
    )
    curve = pd.DataFrame(
        {
            "station": ["S0", "S0", "S0", "S0"],
            "period": [0.1, 1.0, 10.0, 100.0],
            "ang": [30.0, 33.0, 40.0, 42.0],
        }
    )

    monkeypatch.setattr(dim_mod, "ensure_sites", lambda sites, **_: sites)
    monkeypatch.setattr(dim_mod, "classify_dimensionality", lambda *_, **__: dim)
    monkeypatch.setattr(dim_mod, "estimate_strike_sweep", lambda *_, **__: strike)
    monkeypatch.setattr(
        dim_mod, "estimate_strike_phase_tensor", lambda *_, **__: strike
    )
    monkeypatch.setattr(dim_mod, "estimate_strike_consensus", lambda *_, **__: strike)
    monkeypatch.setattr(dim_mod, "strike_curve_sweep", lambda *_, **__: curve)

    report = pre2d_inversion_assessment(
        [site],
        band=(0.1, 100.0),
        rotation_applied=True,
        groom_bailey_attempted=False,
    )

    assert list(report["station"]) == ["S0"]
    assert report.loc[0, "frac_3d"] == pytest.approx(0.5)
    assert report.loc[0, "strike_consensus_deg"] == pytest.approx(35.0)
    assert report.loc[0, "rotation_angle_deg"] == pytest.approx(35.0)
    assert bool(report.loc[0, "rotated_to_strike"]) is True
    assert bool(report.loc[0, "groom_bailey_applied"]) is False
    assert "Groom-Bailey" in report.loc[0, "groom_bailey_reason"]


class TestCorrectStaticShiftContract:
    def test_returns_sites_object(self):
        from pycsamt.site.base import Sites

        sites = _line_sites(5)
        result = correct_static_shift(sites)
        assert isinstance(result, Sites)

    def test_returns_same_number_of_sites(self):
        sites = _line_sites(4)
        result = correct_static_shift(sites)
        count = sum(1 for _ in result)
        assert count == 4

    def test_inplace_false_returns_new_object(self):
        sites = _line_sites(3)
        result = correct_static_shift(sites, inplace=False)
        assert result is not sites

    def test_z_shape_preserved(self):
        from pycsamt.emtools._core import _iter_items

        sites = _line_sites(3)
        result = correct_static_shift(sites)
        for ed in _iter_items(result):
            z, fr = _z_freq(ed)
            assert z is not None
            assert z.ndim == 3
            assert z.shape[1:] == (2, 2)

    def test_freq_count_preserved(self):
        from pycsamt.emtools._core import _iter_items

        sites = _line_sites(3, rho=100.0)
        result = correct_static_shift(sites)
        for ed in _iter_items(result):
            z, fr = _z_freq(ed)
            assert z is not None
            assert z.shape[0] == fr.size

    def test_single_site_no_change(self):
        """With only one station, Hanning window covers it alone → C=1 → no change."""
        from pycsamt.emtools._core import _iter_items

        site = _site("S0", rho=100.0)
        result = correct_static_shift([site], window_m=1500.0)
        for ed in _iter_items(result):
            z, _ = _z_freq(ed)
            assert z is not None
            np.testing.assert_allclose(np.abs(z), np.abs(_make_z(_freqs())), rtol=1e-10)

    def test_uniform_profile_no_change(self):
        """Uniform ρ_a across all stations → spatial mean = individual → C=1."""
        from pycsamt.emtools._core import _iter_items

        sites = _line_sites(5, rho=200.0, spacing=300.0)
        z_before = [
            _z_freq(ed)[0].copy() for ed in _iter_items(ensure_sites_local(sites))
        ]
        result = correct_static_shift(sites, window_m=2000.0)
        z_after = [_z_freq(ed)[0] for ed in _iter_items(result)]
        for zb, za in zip(z_before, z_after):
            np.testing.assert_allclose(np.abs(za), np.abs(zb), rtol=1e-6)

    def test_comp_det_default(self):
        sites = _line_sites(3)
        result = correct_static_shift(sites, comp="det")
        assert result is not None

    def test_comp_xy(self):
        sites = _line_sites(3)
        result = correct_static_shift(sites, comp="xy")
        assert result is not None

    def test_comp_yx(self):
        sites = _line_sites(3)
        result = correct_static_shift(sites, comp="yx")
        assert result is not None

    def test_empty_list_returns_sites(self):
        """Empty input doesn't crash (returns empty Sites)."""
        from pycsamt.site.base import Sites

        result = correct_static_shift([])
        assert isinstance(result, Sites)


def ensure_sites_local(sites):
    from pycsamt.emtools._core import ensure_sites

    return ensure_sites(sites, recursive=True, strict=False)


def _z_freq(ed):
    """Return (z, freq) from either a _FakeSite or a Sites-wrapped Site."""
    z = getattr(ed, "z", None)
    fr = getattr(ed, "freq", None)
    if isinstance(z, np.ndarray) and z.ndim == 3 and isinstance(fr, np.ndarray):
        return z, fr
    Z_obj = getattr(ed, "Z", None)
    if Z_obj is not None:
        z2 = getattr(Z_obj, "z", None)
        fr2 = getattr(Z_obj, "freq", None)
        if isinstance(z2, np.ndarray) and isinstance(fr2, np.ndarray):
            return z2, fr2
    return None, None


# ─────────────────────────────────────────────────────────────────────────────
# correct_static_shift — correction direction
# ─────────────────────────────────────────────────────────────────────────────


class TestCorrectStaticShiftDirection:
    def _build_shifted_line(self, n=5, shift_idx=2, ss=4.0):
        """n stations; one outlier with static shift ss, rest at rho=100."""
        sites = []
        for i in range(n):
            rho_local = 100.0
            if i == shift_idx:
                site = _shifted_site(
                    f"S{i:02d}",
                    rho=rho_local,
                    ss=ss,
                    east=i * 300.0,
                    north=0.0,
                )
            else:
                site = _site(f"S{i:02d}", rho=rho_local, east=i * 300.0, north=0.0)
            sites.append(site)
        return sites

    def test_shifted_station_z_magnitude_reduced(self):
        """Shifted station has larger |Z|; after correction |Z| should decrease."""
        from pycsamt.emtools._core import _iter_items

        sites = self._build_shifted_line(n=5, shift_idx=2, ss=4.0)
        # Capture before-snapshots (copy!) before correction modifies in-place
        S0 = ensure_sites_local(sites)
        mag_before = [
            float(np.mean(np.abs(_z_freq(ed)[0][:, 0, 1]))) for ed in _iter_items(S0)
        ]
        result = correct_static_shift(sites, window_m=3000.0)
        mag_after = [
            float(np.mean(np.abs(_z_freq(ed)[0][:, 0, 1])))
            for ed in _iter_items(result)
        ]
        assert mag_after[2] < mag_before[2], "shifted station: |Z| should decrease"

    def test_phase_unchanged(self):
        """Static shift is a real scalar — it should NOT change the impedance phase."""
        from pycsamt.emtools._core import _iter_items

        sites = self._build_shifted_line(n=5, shift_idx=2, ss=4.0)
        S0 = ensure_sites_local(sites)
        # Capture phase before (static shift doesn't touch phase)
        ph_before = [np.angle(_z_freq(ed)[0][:, 0, 1]).copy() for ed in _iter_items(S0)]
        result = correct_static_shift(sites, window_m=3000.0)
        ph_after = [np.angle(_z_freq(ed)[0][:, 0, 1]) for ed in _iter_items(result)]
        for pb, pa in zip(ph_before, ph_after):
            np.testing.assert_allclose(pa, pb, atol=1e-10)

    def test_correction_factor_positive(self):
        """C = sqrt(ρ_smooth / ρ_obs) must always be > 0 for off-diagonal Z."""
        from pycsamt.emtools._core import _iter_items

        sites = self._build_shifted_line(n=5, shift_idx=1, ss=9.0)
        result = correct_static_shift(sites, window_m=2400.0)
        for ed in _iter_items(result):
            z, _ = _z_freq(ed)
            assert np.all(np.isfinite(z))
            assert np.all(np.abs(z[:, 0, 1]) > 0)  # off-diagonal > 0
            assert np.all(np.abs(z[:, 1, 0]) > 0)


# ─────────────────────────────────────────────────────────────────────────────
# EMAP-style FLMA/TMA spatial filters
# ─────────────────────────────────────────────────────────────────────────────


class TestEMAPSpatialFilters:
    def _outlier_line(self, n=5, outlier_idx=2, outlier_rho=1600.0):
        sites = []
        for i in range(n):
            rho = outlier_rho if i == outlier_idx else 100.0
            sites.append(
                _site(
                    f"S{i:02d}",
                    rho=rho,
                    n=6,
                    east=i * 200.0,
                    north=0.0,
                )
            )
        return sites

    def _outlier_line_with_partial_row(self):
        sites = self._outlier_line()
        sites[2].Z.z[2, 0, 0] = np.nan
        return sites

    def _xy_mean(self, sites, station_index):
        from pycsamt.emtools._core import _iter_items

        items = list(_iter_items(sites))
        z, _ = _z_freq(items[station_index])
        return float(np.mean(np.abs(z[:, 0, 1])))

    def test_fixed_length_moving_average_reduces_center_outlier(self):
        sites = self._outlier_line()
        before = self._xy_mean(sites, 2)

        out = fixed_length_moving_average(
            sites,
            window=3,
            component="offdiag",
        )
        after = self._xy_mean(out, 2)

        assert after < before

    def test_trimmed_moving_average_suppresses_isolated_outlier(self):
        sites = self._outlier_line()
        before = self._xy_mean(sites, 2)

        out = trimmed_moving_average(
            sites,
            window=5,
            component="xy",
        )
        after = self._xy_mean(out, 2)

        assert after < before

    def test_apply_emap_filter_dispatches_flma(self):
        sites = self._outlier_line()
        out = apply_emap_filter(
            sites,
            method="flma",
            window=3,
            component="xy",
        )

        assert out is not None

    def test_apply_emap_filter_rejects_bad_method(self):
        sites = self._outlier_line()

        with pytest.raises(ValueError, match="method"):
            apply_emap_filter(sites, method="unknown")

    def test_confidence_gated_emap_filter_returns_result(self):
        before = self._outlier_line_with_partial_row()
        source = self._outlier_line_with_partial_row()

        result = confidence_gated_emap_filter(
            source,
            before_sites=before,
            method="flma",
            confidence_method="presence",
            component="xy",
            ci_hi=0.90,
            ci_lo=0.50,
            window=3,
        )

        assert isinstance(result, EMAPFilterResult)
        assert not result.report.empty
        assert not result.decisions.empty
        assert result.n_blended >= 1

    def test_confidence_gated_emap_filter_preserves_high_confidence_rows(
        self,
    ):
        before = self._outlier_line()
        source = self._outlier_line()

        result = confidence_gated_emap_filter(
            source,
            before_sites=before,
            method="flma",
            confidence_method="presence",
            component="xy",
            ci_hi=0.90,
            ci_lo=0.50,
            window=3,
        )

        assert result.n_preserved > 0
        assert result.n_blended == 0
        assert result.n_filtered == 0

    def test_emap_filter_report_columns(self):
        before = self._outlier_line()
        source = self._outlier_line()
        after = fixed_length_moving_average(
            source,
            window=3,
            component="xy",
        )

        report = emap_filter_report(
            before,
            after,
            component="xy",
            frequency_hz=10.0,
        )

        expected = {
            "station",
            "component",
            "n_matched_freq",
            "median_delta_log10_abs_z",
            "rms_delta_log10_abs_z",
            "before_log10_abs_z",
            "after_log10_abs_z",
        }
        assert expected.issubset(report.columns)
        assert len(report) == 5

    def test_plot_emap_filter_profile_uses_station_axis(self):
        before = self._outlier_line()
        source = self._outlier_line()
        after = trimmed_moving_average(
            source,
            window=5,
            component="xy",
        )

        ax = plot_emap_filter_profile(
            before,
            after,
            method="tma",
            component="xy",
            frequency_hz=10.0,
        )

        assert ax.xaxis.get_label_position() == "top"
        assert ax.get_legend() is not None
        plt.close("all")

    def test_plot_emap_filter_psection_returns_three_panels(self):
        before = self._outlier_line()
        source = self._outlier_line()
        after = fixed_length_moving_average(
            source,
            window=3,
            component="xy",
        )

        fig = plot_emap_filter_psection(
            before,
            after,
            method="flma",
            component="xy",
        )

        assert isinstance(fig, plt.Figure)
        assert len(fig.axes) >= 3
        assert fig.axes[0].xaxis.get_label_position() == "top"
        plt.close("all")


# ─────────────────────────────────────────────────────────────────────────────
# correct_static_shift — Hanning weight helper
# ─────────────────────────────────────────────────────────────────────────────


class TestHanningWeights:
    def _w(self, dx, w_H):
        from pycsamt.emtools.remove_noise import (
            _hanning_weights,
        )

        return _hanning_weights(np.asarray(dx, dtype=float), w_H)

    def test_centre_is_max(self):
        w = self._w(np.linspace(-500, 500, 101), 1000.0)
        assert w[50] == pytest.approx(w.max(), rel=1e-6)

    def test_zero_at_half_window(self):
        w = self._w(np.array([-500.0, 500.0]), 1000.0)
        np.testing.assert_allclose(w, 0.0, atol=1e-12)

    def test_zero_outside_window(self):
        w = self._w(np.array([-600.0, 600.0]), 1000.0)
        np.testing.assert_allclose(w, 0.0, atol=1e-12)

    def test_nonnegative(self):
        dx = np.linspace(-800, 800, 200)
        w = self._w(dx, 1000.0)
        assert np.all(w >= 0.0)

    def test_integral_approx_one(self):
        dx = np.linspace(-499.9, 499.9, 10000)
        w = self._w(dx, 1000.0)
        val = float(_trapz(w, dx))
        assert val == pytest.approx(1.0, abs=0.01)

    def test_symmetry(self):
        dx = np.linspace(-400, 400, 201)
        w = self._w(dx, 1000.0)
        np.testing.assert_allclose(w, w[::-1], atol=1e-12)


# ─────────────────────────────────────────────────────────────────────────────
# correct_static_shift — _station_positions helper
# ─────────────────────────────────────────────────────────────────────────────


class TestStationPositions:
    def test_fallback_spacing(self):
        from pycsamt.emtools._core import _station_positions

        eds = [_site(f"S{i}") for i in range(4)]
        pos = _station_positions(eds, spacing_m=300.0)
        np.testing.assert_allclose(pos, [0, 300, 600, 900])

    def test_with_east_north(self):
        from pycsamt.emtools._core import _station_positions

        eds = [
            _site("S0", east=0.0, north=0.0),
            _site("S1", east=500.0, north=0.0),
            _site("S2", east=1000.0, north=0.0),
        ]
        pos = _station_positions(eds, spacing_m=200.0)
        assert pos[0] == pytest.approx(0.0, abs=1.0)
        assert pos[1] == pytest.approx(500.0, abs=1.0)
        assert pos[2] == pytest.approx(1000.0, abs=1.0)

    def test_monotone_on_line(self):
        from pycsamt.emtools._core import _station_positions

        eds = [_site(f"S{i}", east=float(i * 400), north=0.0) for i in range(5)]
        pos = _station_positions(eds, spacing_m=200.0)
        assert np.all(np.diff(pos) > 0)

    def test_length_matches_input(self):
        from pycsamt.emtools._core import _station_positions

        eds = [_site(f"S{i}") for i in range(7)]
        pos = _station_positions(eds)
        assert len(pos) == 7

    def test_single_site_fallback(self):
        from pycsamt.emtools._core import _station_positions

        eds = [_site("S0")]
        pos = _station_positions(eds, spacing_m=250.0)
        assert pos[0] == pytest.approx(0.0)

    def test_latlon_fallback_matches_real_span(self):
        """Regression for GH report: stations without east/north attrs
        (the normal case for real EDI-backed Site objects, which only
        expose lat/lon) must not fall back to index*spacing_m — that
        silently stretched a ~2050 m AMT line out to ~8000 m."""
        from types import SimpleNamespace

        from pycsamt.emtools._core import _station_positions

        n = 41
        span_m = 2050.0
        lat0 = 10.0
        lats = [lat0 + i * (span_m / (n - 1)) / 111_000.0 for i in range(n)]

        def _latlon_site(i, lat):
            head = SimpleNamespace(lat=lat, long=120.0, elev=0.0)
            edi = SimpleNamespace(
                station=f"S{i}", get_section=lambda *_a, **_k: head
            )
            return SimpleNamespace(edi=edi)

        eds = [_latlon_site(i, lat) for i, lat in enumerate(lats)]
        pos = _station_positions(eds, spacing_m=200.0)

        # rel=1e-2: the synthetic fixture spaces latitudes using a flat
        # 111 km/deg constant, while the real UTM projection accounts for
        # Earth's ellipticity, so a small (~0.3%) difference is expected.
        assert pos.max() == pytest.approx(span_m, rel=1e-2)
        assert pos.max() != pytest.approx((n - 1) * 200.0)

    def test_latlon_offsets_are_cached_on_the_edi_object(self):
        from pycsamt.emtools._core import (
            _GEO_OFFSET_CACHE_ATTR,
            _station_latlon_offsets,
        )

        eds = [
            _site("S0", east=None, north=None),
            _site("S1", east=None, north=None),
        ]
        # attach lat/lon via get_section, like a real EDI head would.
        from types import SimpleNamespace

        for ed, lat in zip(eds, (10.0, 10.01)):
            ed.get_section = lambda *_a, _h=SimpleNamespace(
                lat=lat, long=120.0, elev=0.0
            ), **_k: _h

        first = _station_latlon_offsets(eds)
        assert all(o is not None for o in first)
        cached = getattr(eds[0], _GEO_OFFSET_CACHE_ATTR, None)
        assert cached is not None

        # second call reuses the cached (lat, lon)-keyed offset unchanged.
        second = _station_latlon_offsets(eds)
        assert second == first

    def test_force_spacing_bypasses_coordinates(self):
        from pycsamt.emtools._core import _station_positions

        eds = [
            _site("S0", east=0.0, north=0.0),
            _site("S1", east=500.0, north=0.0),
            _site("S2", east=1000.0, north=0.0),
        ]
        pos = _station_positions(eds, spacing_m=100.0, force_spacing=True)
        np.testing.assert_allclose(pos, [0.0, 100.0, 200.0])


# ─────────────────────────────────────────────────────────────────────────────
# plot_confidence_profile — return type and axes labels
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotConfidenceProfile:
    def _sites_all_good(self, n=5):
        """All stations with 100 % valid Z → CI = 1.0."""
        return [_site(f"S{i:02d}", east=i * 300.0, north=0.0) for i in range(n)]

    def _sites_mixed(self):
        """Three stations: good, some NaN, mostly NaN."""
        fr = _freqs(10)

        z_good = _make_z(fr, 100.0)

        z_partial = _make_z(fr, 100.0).copy()
        z_partial[:4] = np.nan  # 4 of 10 bad → CI = 0.6

        z_bad = _make_z(fr, 100.0).copy()
        z_bad[:8] = np.nan  # 8 of 10 bad → CI = 0.2

        return [
            _FakeSite("S_good", z_good, fr, east=0.0, north=0.0),
            _FakeSite("S_partial", z_partial, fr, east=300.0, north=0.0),
            _FakeSite("S_bad", z_bad, fr, east=600.0, north=0.0),
        ]

    def test_returns_axes(self):
        sites = self._sites_all_good()
        ax = plot_confidence_profile(sites)
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_accepts_existing_ax(self):
        fig, ax_in = plt.subplots()
        sites = self._sites_all_good()
        ax_out = plot_confidence_profile(sites, ax=ax_in)
        assert ax_out is ax_in
        plt.close("all")

    def test_ylim_includes_thresholds(self):
        sites = self._sites_all_good()
        ax = plot_confidence_profile(sites)
        ylo, yhi = ax.get_ylim()
        assert ylo <= 0.85
        assert yhi >= 1.0
        plt.close("all")

    def test_has_two_threshold_lines(self):
        sites = self._sites_all_good()
        ax = plot_confidence_profile(sites)
        hlines = [l for l in ax.get_lines() if l.get_linestyle() == "--"]
        assert len(hlines) >= 2
        plt.close("all")

    def test_xlabel_contains_distance(self):
        sites = self._sites_all_good()
        ax = plot_confidence_profile(sites)
        assert (
            "distance" in ax.get_xlabel().lower() or "dist" in ax.get_xlabel().lower()
        )
        plt.close("all")

    def test_scatter_points_count(self):
        sites = self._sites_all_good(n=4)
        ax = plot_confidence_profile(sites)
        colls = [coll for coll in ax.collections if isinstance(coll, PathCollection)]
        n_pts = sum(len(c.get_offsets()) for c in colls)
        assert n_pts == 4
        plt.close("all")

    def test_mixed_ci_colors(self):
        """green / pink / red for good / partial / bad stations."""
        sites = self._sites_mixed()
        ax = plot_confidence_profile(sites, ci_hi=0.95, ci_lo=0.50)
        colls = [c for c in ax.collections if len(c.get_offsets()) > 0]
        assert len(colls) > 0
        plt.close("all")

    def test_custom_thresholds(self):
        sites = self._sites_all_good()
        ax = plot_confidence_profile(sites, ci_hi=0.80, ci_lo=0.40)
        legend_text = " ".join(t.get_text() for t in ax.get_legend().get_texts())
        assert "0.80" in legend_text or "0.8" in legend_text
        assert "0.40" in legend_text or "0.4" in legend_text
        plt.close("all")

    def test_empty_sites_returns_axes(self):
        ax = plot_confidence_profile([])
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_figsize_kwarg(self):
        sites = self._sites_all_good(3)
        ax = plot_confidence_profile(sites, figsize=(6.0, 3.0))
        w, h = ax.get_figure().get_size_inches()
        assert w == pytest.approx(6.0)
        assert h == pytest.approx(3.0)
        plt.close("all")

    def test_no_coordinates_fallback(self):
        """Sites without east/north use integer spacing — should not crash."""
        sites = [_site(f"S{i}") for i in range(4)]
        ax = plot_confidence_profile(sites, spacing_m=150.0)
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_ci_values_between_0_and_1(self):
        sites = self._sites_mixed()
        ax = plot_confidence_profile(sites)
        for coll in ax.collections:
            for _x, y in coll.get_offsets():
                assert 0.0 <= y <= 1.0
        plt.close("all")

    def test_composite_confidence_table_columns(self):
        sites = self._sites_mixed()
        tb = station_confidence_table(sites, method="composite")

        expected = {
            "station",
            "distance_m",
            "confidence",
            "coverage",
            "confidence_err",
            "offdiag",
            "diagonal",
            "phase",
            "spatial",
        }
        assert expected.issubset(tb.columns)
        assert tb["confidence"].between(0.0, 1.0).all()
        assert tb["confidence_err"].between(0.0, 1.0).all()

    def test_station_confidence_table_api_flag(self):
        sites = self._sites_mixed()
        plain = station_confidence_table(sites, method="composite", api=False)
        view = station_confidence_table(sites, method="composite", api=True)

        assert isinstance(view, APIFrame)
        assert view.kind == "emtools.qc.station_confidence"
        assert view.df.equals(plain)

    def test_station_confidence_table_obeys_disabled_api_view(self):
        """With global disabled, api=None returns raw pandas; api=True still forces APIFrame."""
        sites = self._sites_mixed()
        try:
            configure_api_view(backend=False)

            out_default = station_confidence_table(sites, method="composite")
            assert not isinstance(out_default, APIFrame)
            assert "confidence" in out_default.columns

            out_forced = station_confidence_table(sites, method="composite", api=True)
            assert isinstance(out_forced, APIFrame)
        finally:
            reset_api_view()

    def test_composite_confidence_penalizes_bad_station(self):
        sites = self._sites_mixed()
        tb = station_confidence_table(sites, method="composite")
        good = float(tb.loc[tb.station == "S_good", "confidence"].iloc[0])
        bad = float(tb.loc[tb.station == "S_bad", "confidence"].iloc[0])

        assert bad < good

    def test_plot_confidence_profile_composite_method(self):
        sites = self._sites_mixed()
        ax = plot_confidence_profile(sites, method="composite")

        assert "composite" in ax.get_title().lower()
        assert ax.get_legend() is not None
        plt.close("all")

    def test_plot_confidence_profile_full_shade_mode(self):
        sites = self._sites_mixed()
        ax = plot_confidence_profile(
            sites,
            ci_hi=0.95,
            ci_lo=0.50,
            shade_mode="full",
        )

        assert len(ax.patches) > 0
        plt.close("all")

    def _sites_all_bad(self, n=25):
        """n stations all below the default ci_lo, to exercise low-point
        label decluttering."""
        fr = _freqs(10)
        z_bad = _make_z(fr, 100.0).copy()
        z_bad[:8] = np.nan  # 8 of 10 bad -> CI = 0.2 < ci_lo
        return [
            _FakeSite(f"S{i:02d}", z_bad, fr, east=i * 300.0, north=0.0)
            for i in range(n)
        ]

    def test_annotate_low_auto_thins_many_low_points(self):
        """Regression: with most/all stations flagged low, the default
        annotate_low=True must not label every single one -- that clutters
        and duplicates the top station-name axis (GH #76 follow-up)."""
        sites = self._sites_all_bad(n=25)
        ax = plot_confidence_profile(sites)
        assert 0 < len(ax.texts) < 25
        plt.close("all")

    def test_annotate_low_step_1_labels_every_low_point(self):
        sites = self._sites_all_bad(n=25)
        ax = plot_confidence_profile(sites, annotate_low_step=1)
        assert len(ax.texts) == 25
        plt.close("all")

    def test_annotate_low_false_disables_low_point_labels(self):
        sites = self._sites_all_bad(n=25)
        ax = plot_confidence_profile(sites, annotate_low=False)
        assert len(ax.texts) == 0
        plt.close("all")

    def test_frequency_confidence_table_columns(self):
        sites = self._sites_mixed()
        tb = frequency_confidence_table(sites, method="composite")

        expected = {
            "station",
            "frequency_hz",
            "log10_period",
            "confidence",
            "confidence_err",
            "coverage",
            "offdiag",
            "diagonal",
            "phase",
            "spatial",
            "flags",
        }
        assert expected.issubset(tb.columns)
        assert tb["confidence"].between(0.0, 1.0).all()

    def test_frequency_confidence_table_api_flag(self):
        sites = self._sites_mixed()
        plain = frequency_confidence_table(sites, method="presence", api=False)
        view = frequency_confidence_table(sites, method="presence", api=True)

        assert isinstance(view, APIFrame)
        assert view.kind == "emtools.qc.frequency_confidence"
        assert view.df.equals(plain)

    def test_frequency_confidence_table_has_one_row_per_sample(self):
        sites = self._sites_mixed()
        tb = frequency_confidence_table(sites, method="presence")

        assert len(tb) == 30

    def test_plot_frequency_confidence_psection(self):
        sites = self._sites_mixed()
        ax = plot_frequency_confidence_psection(sites, method="composite")

        assert ax.images
        assert ax.xaxis.get_label_position() == "top"
        assert "frequency confidence" in ax.get_title().lower()
        plt.close("all")

    def test_plot_station_confidence_spectrum(self):
        sites = self._sites_mixed()
        ax = plot_station_confidence_spectrum(
            sites,
            station="S_partial",
            method="composite",
        )

        assert "S_partial" in ax.get_title()
        assert ax.get_legend() is not None
        plt.close("all")

    def test_plot_station_confidence_dashboard(self):
        sites = self._sites_mixed()
        fig = plot_station_confidence_dashboard(
            sites,
            station="S_partial",
            method="composite",
        )

        assert isinstance(fig, plt.Figure)
        assert len(fig.axes) == 6
        plt.close("all")

    def test_plot_confidence_band_summary(self):
        sites = self._sites_mixed()
        ax = plot_confidence_band_summary(sites, method="composite")

        assert "period-band" in ax.get_title().lower()
        assert ax.get_legend() is not None
        plt.close("all")


class TestConfidenceFrequencyEditing:
    def _site_with_recoverable_rows(self):
        fr = _freqs(6)
        z = _make_z(fr, 100.0)
        z[2, 0, 0] = np.nan
        z[4] = np.nan
        return [_FakeSite("S_edit", z, fr, east=0.0, north=0.0)]

    def _first_z_block(self, sites):
        from pycsamt.emtools._core import (
            _get_z_block,
            _iter_items,
        )

        ed = next(_iter_items(sites))
        return _get_z_block(ed)

    def test_mask_low_confidence_keeps_grid(self):
        sites = self._site_with_recoverable_rows()
        out = mask_low_confidence_frequencies(
            sites,
            method="presence",
            threshold=0.80,
        )
        _, z, fr = self._first_z_block(out)

        assert fr.size == 6
        assert np.isnan(z[2]).any()
        assert np.isnan(z[4]).all()

    def test_drop_low_confidence_removes_bad_rows(self):
        sites = self._site_with_recoverable_rows()
        out = drop_low_confidence_frequencies(
            sites,
            method="presence",
            threshold=0.80,
        )
        _, z, fr = self._first_z_block(out)

        assert fr.size == 4
        assert z.shape[0] == 4

    def test_recover_low_confidence_interpolates_recoverable_rows(self):
        sites = self._site_with_recoverable_rows()
        out = recover_low_confidence_frequencies(
            sites,
            method="presence",
            ci_hi=0.90,
            ci_lo=0.50,
            reject="mask",
        )
        _, z, fr = self._first_z_block(out)

        assert fr.size == 6
        assert np.isfinite(z[2]).all()
        assert np.isnan(z[4]).all()

    def test_recover_low_confidence_can_drop_rejected_rows(self):
        sites = self._site_with_recoverable_rows()
        out = recover_low_confidence_frequencies(
            sites,
            method="presence",
            ci_hi=0.90,
            ci_lo=0.50,
            reject="drop",
        )
        _, z, fr = self._first_z_block(out)

        assert fr.size == 5
        assert z.shape[0] == 5

    def test_frequency_edit_report_counts_changes(self):
        before = self._site_with_recoverable_rows()
        sites = self._site_with_recoverable_rows()
        out = recover_low_confidence_frequencies(
            sites,
            method="presence",
            ci_hi=0.90,
            ci_lo=0.50,
            reject="drop",
        )
        report = frequency_edit_report(before, out, method="presence")

        row = report.iloc[0]
        assert int(row["n_freq_before"]) == 6
        assert int(row["n_freq_after"]) == 5
        assert int(row["n_dropped"]) == 1

    def test_frequency_edit_report_api_flag(self):
        before = self._site_with_recoverable_rows()
        sites = self._site_with_recoverable_rows()
        out = recover_low_confidence_frequencies(
            sites,
            method="presence",
            ci_hi=0.90,
            ci_lo=0.50,
            reject="drop",
        )
        plain = frequency_edit_report(before, out, method="presence", api=False)
        view = frequency_edit_report(before, out, method="presence", api=True)

        assert isinstance(view, APIFrame)
        assert view.kind == "emtools.frequency.edit_report"
        assert view.df.equals(plain)

    def test_frequency_edit_decision_table_marks_recovery_and_drop(self):
        before = self._site_with_recoverable_rows()
        sites = self._site_with_recoverable_rows()
        out = recover_low_confidence_frequencies(
            sites,
            method="presence",
            ci_hi=0.90,
            ci_lo=0.50,
            reject="drop",
        )
        decisions = frequency_edit_decision_table(
            before,
            out,
            method="presence",
        )

        actions = set(decisions["action"])
        assert "recovered" in actions
        assert "dropped" in actions

    def test_frequency_edit_decision_table_api_flag(self):
        before = self._site_with_recoverable_rows()
        sites = self._site_with_recoverable_rows()
        out = recover_low_confidence_frequencies(
            sites,
            method="presence",
            ci_hi=0.90,
            ci_lo=0.50,
            reject="drop",
        )
        plain = frequency_edit_decision_table(before, out, method="presence", api=False)
        view = frequency_edit_decision_table(
            before,
            out,
            method="presence",
            api=True,
        )

        assert isinstance(view, APIFrame)
        assert view.kind == "emtools.frequency.edit_decisions"
        assert view.df.equals(plain)

    def test_edit_frequencies_by_confidence_returns_result(self):
        before = self._site_with_recoverable_rows()
        sites = self._site_with_recoverable_rows()
        result = edit_frequencies_by_confidence(
            sites,
            before_sites=before,
            mode="recover",
            method="presence",
            ci_hi=0.90,
            ci_lo=0.50,
            reject="drop",
            api=False,
        )

        assert isinstance(result, FrequencyEditResult)
        assert not result.report.empty
        assert not result.decisions.empty
        assert result.n_dropped == 1
        assert result.n_recovered >= 1

    def test_edit_frequencies_by_confidence_api_result(self):
        before = self._site_with_recoverable_rows()
        sites = self._site_with_recoverable_rows()
        result = edit_frequencies_by_confidence(
            sites,
            before_sites=before,
            mode="recover",
            method="presence",
            ci_hi=0.90,
            ci_lo=0.50,
            reject="drop",
            api=True,
        )

        assert isinstance(result, APIResult)
        assert isinstance(result.report, APIFrame)
        assert isinstance(result.decisions, APIFrame)
        assert result.kind == "emtools.frequency.edit"
        assert result.n_dropped == 1
        assert result.n_recovered >= 1
        assert "dropped" in set(result.decisions.df["action"])

    def test_edit_frequencies_by_confidence_api_result_can_be_disabled(self):
        """With global disabled: api=None returns FrequencyEditResult; api=True still forces APIResult."""
        before = self._site_with_recoverable_rows()
        sites = self._site_with_recoverable_rows()
        try:
            configure_api_view(backend=False)

            # Default (api=None) respects global disabled → raw result.
            result_default = edit_frequencies_by_confidence(
                sites,
                before_sites=before,
                mode="recover",
                method="presence",
                ci_hi=0.90,
                ci_lo=0.50,
                reject="drop",
            )
            assert isinstance(result_default, FrequencyEditResult)
            assert not isinstance(result_default, APIResult)
            assert result_default.n_dropped == 1

            # Explicit api=True overrides global disabled → APIResult.
            result_forced = edit_frequencies_by_confidence(
                sites,
                before_sites=before,
                mode="recover",
                method="presence",
                ci_hi=0.90,
                ci_lo=0.50,
                reject="drop",
                api=True,
            )
            assert isinstance(result_forced, APIResult)
        finally:
            reset_api_view()

    def test_edit_frequencies_by_confidence_rejects_bad_mode(self):
        sites = self._site_with_recoverable_rows()

        with pytest.raises(ValueError, match="mode"):
            edit_frequencies_by_confidence(
                sites,
                mode="unknown",
                method="presence",
            )

    def test_plot_frequency_edit_summary(self):
        before = self._site_with_recoverable_rows()
        sites = self._site_with_recoverable_rows()
        out = mask_low_confidence_frequencies(
            sites,
            method="presence",
            threshold=0.80,
        )
        ax = plot_frequency_edit_summary(before, out, method="presence")

        assert "frequency edit" in ax.get_title().lower()
        assert ax.xaxis.get_label_position() == "top"
        assert ax.get_legend() is not None
        plt.close("all")

    def test_plot_frequency_edit_decisions(self):
        before = self._site_with_recoverable_rows()
        sites = self._site_with_recoverable_rows()
        out = recover_low_confidence_frequencies(
            sites,
            method="presence",
            ci_hi=0.90,
            ci_lo=0.50,
            reject="drop",
        )
        ax = plot_frequency_edit_decisions(before, out, method="presence")

        assert ax.images
        assert "decisions" in ax.get_title().lower()
        assert ax.xaxis.get_label_position() == "top"
        plt.close("all")


class TestManualFrequencyDrop:
    def test_drop_freqs_manual_removes_rows_without_injecting_nan(self):
        site = _site("S00", n=6)
        drop = float(site.freq[2])

        out = drop_freqs_manual([site], drop_freqs=[drop])
        edited = list(out)[0]
        raw = getattr(edited, "edi", edited)
        Z = raw.Z

        assert Z.freq.size == 5
        assert Z.z.shape == (5, 2, 2)
        assert not np.any(np.isclose(Z.freq, drop))
        assert np.isfinite(Z.z).all()


class TestApplyEachInplaceSemantics:
    """Regression: ``inplace=False`` must not mutate the caller's sites.

    The underlying ``Z.z`` was previously scaled in place regardless of
    the ``inplace`` flag (``_apply_each`` rebuilt a new container around
    the *same* mutated items). That corrupted the "before" data, so
    before/after comparison plots showed identical images.
    """

    def _shifted_line(self, n=5, shift_idx=2, ss=4.0):
        sites = []
        for i in range(n):
            if i == shift_idx:
                sites.append(
                    _shifted_site(
                        f"S{i:02d}",
                        rho=100.0,
                        ss=ss,
                        east=i * 300.0,
                        north=0.0,
                    )
                )
            else:
                sites.append(_site(f"S{i:02d}", rho=100.0, east=i * 300.0, north=0.0))
        return sites

    def test_inplace_false_preserves_original_and_changes_copy(self):
        from pycsamt.emtools.ss import correct_ss_ama

        sites = ensure_sites_local(self._shifted_line())
        items_before = list(sites)
        z_orig = [_z_freq(ed)[0].copy() for ed in items_before]

        corrected = correct_ss_ama(sites, half_window=3, inplace=False)

        # original is untouched
        z_after_orig = [_z_freq(ed)[0] for ed in list(sites)]
        for a, b in zip(z_orig, z_after_orig):
            assert np.allclose(
                a, b, equal_nan=True
            ), "inplace=False mutated the caller's sites"

        # the returned copy actually differs from the original
        z_corr = [_z_freq(ed)[0] for ed in list(corrected)]
        changed = any(
            not np.allclose(o, c, equal_nan=True) for o, c in zip(z_orig, z_corr)
        )
        assert changed, "correction produced no change in the copy"

    def test_inplace_true_mutates_original(self):
        from pycsamt.emtools.ss import correct_ss_ama

        sites = ensure_sites_local(self._shifted_line())
        z_orig = [_z_freq(ed)[0].copy() for ed in list(sites)]

        correct_ss_ama(sites, half_window=3, inplace=True)

        z_now = [_z_freq(ed)[0] for ed in list(sites)]
        changed = any(
            not np.allclose(o, c, equal_nan=True) for o, c in zip(z_orig, z_now)
        )
        assert changed, "inplace=True did not mutate the original"
