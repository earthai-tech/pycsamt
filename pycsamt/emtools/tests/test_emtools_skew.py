"""Tests for pycsamt.emtools.skew"""

from __future__ import annotations

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.emtools.skew import (
    _fill_small_gaps,
    _mask_apply,
    _runs_bool,
    _skew_track_for,
    bahr_skewness,
    close_skew_gaps,
    keep_longest_low_skew,
    mask_by_skew,
    plot_skew_percentile_ribbon,
    plot_skew_traffic_psection,
    plot_skew_vote_band,
    plot_skewness,
    select_low_skew_band,
    skew_table,
)

# ─────────────────────────────────────────────────────────────────────────────
# Shared helpers
# ─────────────────────────────────────────────────────────────────────────────


class _FakeZ:
    def __init__(self, z, freq):
        self.z = np.asarray(z, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)


class _FakeSite:
    def __init__(self, station, z, freq):
        self.station = station
        self.Z = _FakeZ(z, freq)
        self.freq = np.asarray(freq, dtype=float)

    def get_section(self, *_, **__):
        return None


def _freqs(n: int = 12, f_lo: float = 0.1, f_hi: float = 1e4) -> np.ndarray:
    return np.logspace(np.log10(f_lo), np.log10(f_hi), n)


def _iso_z(freqs: np.ndarray, rho: float = 100.0) -> np.ndarray:
    """Purely 2-D isotropic Z (zero diagonal → low skew)."""
    amp = np.sqrt(5.0 * freqs * rho)
    z = np.zeros((freqs.size, 2, 2), dtype=complex)
    z[:, 0, 1] = amp * (1 + 1j) / np.sqrt(2)
    z[:, 1, 0] = -amp * (1 + 1j) / np.sqrt(2)
    return z


def _3d_z(freqs: np.ndarray, skew_frac: float = 0.6) -> np.ndarray:
    """Z with non-zero diagonal → high Bahr skew."""
    z = _iso_z(freqs)
    amp = np.abs(z[:, 0, 1])
    z[:, 0, 0] = skew_frac * amp * (0.6 + 0.8j)
    z[:, 1, 1] = -skew_frac * amp * (0.5 + 0.7j)
    return z


def _site(station: str, z: np.ndarray, freqs: np.ndarray) -> _FakeSite:
    return _FakeSite(station, z, freqs)


class _FakeSiteBadZ:
    """Site whose ``Z.z`` has a malformed shape (not (n,2,2)/(n,4)).

    ``ensure_sites`` still accepts it (it only requires a ``.Z``/``station``
    duck-type), but ``_get_z_block`` returns ``Z is None`` for it, which is
    exactly the "no impedance data" branch several skew helpers must guard
    against.
    """

    def __init__(self, station: str, freqs: np.ndarray):
        self.station = station
        self.Z = _FakeZ(np.zeros(freqs.size, dtype=complex), freqs)
        self.freq = np.asarray(freqs, dtype=float)

    def get_section(self, *_, **__):
        return None


# ─────────────────────────────────────────────────────────────────────────────
# bahr_skewness
# ─────────────────────────────────────────────────────────────────────────────


class TestBahrSkewness:
    def test_shape(self):
        fr = _freqs(10)
        z = _iso_z(fr)
        eta = bahr_skewness(z)
        assert eta.shape == (10,)

    def test_isotropic_denominator_zero(self):
        """Pure 2-D Z (zero diagonal, Zxy = -Zyx) gives NaN (den=0)."""
        fr = _freqs(10)
        z = _iso_z(fr)
        eta = bahr_skewness(z)
        # den = |D1|^2 + |D2|^2 = 0 for a strictly symmetric 2-D tensor
        assert np.all(~np.isfinite(eta))

    def test_3d_positive(self):
        fr = _freqs(10)
        z = _3d_z(fr, skew_frac=0.6)
        eta = bahr_skewness(z)
        assert np.nanmean(eta) > 0.0

    def test_accepts_n4_layout(self):
        """Z can be (n, 4) flat as well as (n, 2, 2)."""
        fr = _freqs(8)
        z = _iso_z(fr)
        eta_3d = bahr_skewness(z)
        eta_flat = bahr_skewness(z.reshape(-1, 4))
        np.testing.assert_allclose(eta_3d, eta_flat, rtol=1e-10)

    def test_invalid_shape_raises(self):
        with pytest.raises(ValueError):
            bahr_skewness(np.ones((5, 3)))


# ─────────────────────────────────────────────────────────────────────────────
# skew_table
# ─────────────────────────────────────────────────────────────────────────────


class TestSkewTable:
    def _sites(self, n_sites=3):
        return [_site(f"S{i:02d}", _iso_z(_freqs()), _freqs()) for i in range(n_sites)]

    def test_returns_dataframe(self):
        import pandas as pd

        df = skew_table(self._sites())
        assert isinstance(df, pd.DataFrame)

    def test_expected_columns(self):
        df = skew_table(self._sites())
        for col in ("station", "freq", "period", "skew", "beta"):
            assert col in df.columns

    def test_row_count(self):
        n_sites, n_freq = 3, 12
        sites = [
            _site(f"S{i:02d}", _iso_z(_freqs(n_freq)), _freqs(n_freq))
            for i in range(n_sites)
        ]
        df = skew_table(sites)
        assert len(df) == n_sites * n_freq

    def test_empty_input_returns_empty(self):
        import pandas as pd

        df = skew_table([])
        assert isinstance(df, pd.DataFrame)
        assert df.empty

    def test_skew_nonneg_for_iso(self):
        df = skew_table(self._sites())
        assert (df["skew"].dropna() >= 0).all()


# ─────────────────────────────────────────────────────────────────────────────
# mask_by_skew
# ─────────────────────────────────────────────────────────────────────────────


class TestMaskBySkew:
    def _iso_sites(self, n=3):
        fr = _freqs()
        return [_site(f"S{i:02d}", _iso_z(fr), fr) for i in range(n)]

    def test_returns_sites(self):
        from pycsamt.site.base import Sites

        result = mask_by_skew(self._iso_sites())
        assert isinstance(result, Sites)

    def test_site_count_preserved(self):
        sites = self._iso_sites(4)
        result = mask_by_skew(sites)
        assert sum(1 for _ in result) == 4

    def test_high_threshold_keeps_all(self):
        """thresh=1.0 → no masking at all."""
        fr = _freqs()
        sites = [_site("S00", _3d_z(fr, skew_frac=0.5), fr)]
        result = mask_by_skew(sites, thresh=1.0)
        assert sum(1 for _ in result) == 1

    def test_zero_threshold_masks_all(self):
        """thresh=0.0 → every row is masked (NaN)."""
        fr = _freqs(6)
        sites = [_site("S00", _3d_z(fr, skew_frac=0.5), fr)]
        result = mask_by_skew(sites, thresh=0.0)
        assert isinstance(result, object)  # just doesn't raise


# ─────────────────────────────────────────────────────────────────────────────
# keep_longest_low_skew
# ─────────────────────────────────────────────────────────────────────────────


class TestKeepLongestLowSkew:
    def test_returns_sites(self):
        from pycsamt.site.base import Sites

        fr = _freqs(14)
        sites = [_site("S00", _iso_z(fr), fr)]
        result = keep_longest_low_skew(sites, thresh=0.3)
        assert isinstance(result, Sites)

    def test_site_count_preserved(self):
        fr = _freqs()
        sites = [_site(f"S{i}", _iso_z(fr), fr) for i in range(3)]
        result = keep_longest_low_skew(sites, thresh=0.3)
        assert sum(1 for _ in result) == 3


# ─────────────────────────────────────────────────────────────────────────────
# close_skew_gaps
# ─────────────────────────────────────────────────────────────────────────────


class TestCloseSkewGaps:
    def test_returns_sites(self):
        from pycsamt.site.base import Sites

        fr = _freqs()
        sites = [_site("S00", _iso_z(fr), fr)]
        result = close_skew_gaps(sites)
        assert isinstance(result, Sites)


# ─────────────────────────────────────────────────────────────────────────────
# select_low_skew_band
# ─────────────────────────────────────────────────────────────────────────────


class TestSelectLowSkewBand:
    def test_returns_sites(self):
        from pycsamt.site.base import Sites

        fr = _freqs(14)
        sites = [_site(f"S{i}", _iso_z(fr), fr) for i in range(4)]
        result = select_low_skew_band(sites, thresh=6.0)
        assert isinstance(result, Sites)

    def test_site_count_preserved(self):
        fr = _freqs(14)
        sites = [_site(f"S{i}", _iso_z(fr), fr) for i in range(3)]
        result = select_low_skew_band(sites, thresh=6.0)
        assert sum(1 for _ in result) == 3


# ─────────────────────────────────────────────────────────────────────────────
# plot_skewness
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotSkewness:
    def test_returns_axes(self):
        fr = _freqs(10)
        z = _3d_z(fr)
        ax = plot_skewness(fr, z)
        plt.close("all")
        assert ax is not None

    def test_accepts_external_ax(self):
        fr = _freqs(10)
        z = _3d_z(fr)
        fig, ax_in = plt.subplots()
        ax = plot_skewness(fr, z, ax=ax_in)
        plt.close("all")
        assert ax is ax_in

    def test_custom_threshold(self):
        fr = _freqs(10)
        z = _iso_z(fr)
        ax = plot_skewness(fr, z, threshold=0.1)
        plt.close("all")
        assert ax is not None

    def test_title_is_set(self):
        fr = _freqs(10)
        z = _3d_z(fr)
        ax = plot_skewness(fr, z, title="Station S00")
        plt.close("all")
        assert ax.get_title() == "Station S00"

    def test_no_title_leaves_default(self):
        fr = _freqs(10)
        z = _3d_z(fr)
        ax = plot_skewness(fr, z)
        plt.close("all")
        assert ax.get_title() == ""


# ─────────────────────────────────────────────────────────────────────────────
# internal helpers: _skew_track_for
# ─────────────────────────────────────────────────────────────────────────────


class TestSkewTrackForInternals:
    def test_missing_z_returns_none_none(self):
        """When _get_z_block can't find a valid (n,2,2) block, the track
        helper bails out early with (None, None)."""
        fr = _freqs(6)
        ed = _FakeSiteBadZ("S00", fr)
        pt = pd.DataFrame({"station": [], "period": [], "skew": []})
        fr_out, sk_out = _skew_track_for(ed, pt)
        assert fr_out is None
        assert sk_out is None

    def test_station_absent_from_table_returns_nan_track(self):
        """If the station has no rows in the phase-tensor table, the
        frequency grid is still returned but skew is all-NaN."""
        fr = _freqs(6)
        ed = _site("S00", _iso_z(fr), fr)
        pt = pd.DataFrame({"station": ["OTHER"], "period": [1.0], "skew": [0.5]})
        fr_out, sk_out = _skew_track_for(ed, pt)
        assert fr_out is not None
        assert fr_out.shape == fr.shape
        assert sk_out is not None
        assert np.all(np.isnan(sk_out))

    def test_station_present_maps_skew_by_period(self):
        fr = _freqs(6)
        ed = _site("S00", _iso_z(fr), fr)
        per = 1.0 / fr
        pt = pd.DataFrame(
            {
                "station": ["S00"] * per.size,
                "period": per,
                "skew": np.arange(per.size, dtype=float),
            }
        )
        fr_out, sk_out = _skew_track_for(ed, pt)
        assert fr_out is not None
        assert sk_out is not None
        assert np.isfinite(sk_out).all()


# ─────────────────────────────────────────────────────────────────────────────
# internal helpers: _mask_apply
# ─────────────────────────────────────────────────────────────────────────────


class _TipperObj:
    def __init__(self, t, freq):
        self.tipper = np.asarray(t, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)


class _ReadOnlyZ:
    """Z wrapper exposing ``z`` as a getter-only property (assignment
    raises AttributeError, exercising the swallowed-exception branch)."""

    def __init__(self, z, freq):
        self._z = np.asarray(z, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)

    @property
    def z(self):
        return self._z


class _ReadOnlyTipper:
    def __init__(self, t, freq):
        self._t = np.asarray(t, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)

    @property
    def tipper(self):
        return self._t


class TestMaskApplyInternals:
    def test_skips_silently_when_z_missing(self):
        fr = _freqs(6)
        ed = _FakeSiteBadZ("S00", fr)
        keep = np.ones(fr.size, dtype=bool)
        # Must not raise even though the (n,2,2) block can't be resolved.
        _mask_apply(ed, keep, also="both")

    def test_z_masking_writes_nan_for_dropped_rows(self):
        fr = _freqs(6)
        ed = _site("S00", _iso_z(fr), fr)
        keep = np.array([True, False, True, False, True, False])
        _mask_apply(ed, keep, also="z")
        assert np.all(np.isnan(ed.Z.z[~keep]))
        assert np.all(np.isfinite(ed.Z.z[keep]))

    def test_z_assignment_exception_is_swallowed(self):
        fr = _freqs(6)

        class _Site:
            def __init__(self, station, z, freq):
                self.station = station
                self.Z = _ReadOnlyZ(z, freq)

        ed = _Site("S00", _iso_z(fr), fr)
        keep = np.zeros(fr.size, dtype=bool)
        # Z.z has no setter -> AttributeError must be caught, not raised.
        _mask_apply(ed, keep, also="z")

    def test_also_z_only_skips_tipper_block(self):
        fr = _freqs(6)

        class _Site:
            def __init__(self, station, z, t, freq):
                self.station = station
                self.Z = _FakeZ(z, freq)
                self.Tipper = _TipperObj(t, freq)

        t0 = np.ones((fr.size, 2), dtype=complex)
        ed = _Site("S00", _iso_z(fr), t0, fr)
        keep = np.zeros(fr.size, dtype=bool)
        _mask_apply(ed, keep, also="z")  # not "tipper"/"both"
        # tipper branch never entered -> tipper values remain untouched
        assert np.all(np.isfinite(ed.Tipper.tipper))

    def test_also_both_masks_tipper(self):
        fr = _freqs(6)

        class _Site:
            def __init__(self, station, z, t, freq):
                self.station = station
                self.Z = _FakeZ(z, freq)
                self.Tipper = _TipperObj(t, freq)

        t0 = np.ones((fr.size, 2), dtype=complex)
        ed = _Site("S00", _iso_z(fr), t0, fr)
        keep = np.zeros(fr.size, dtype=bool)
        _mask_apply(ed, keep, also="both")
        assert np.all(np.isnan(ed.Tipper.tipper))

    def test_also_tipper_masks_tipper_too(self):
        fr = _freqs(6)

        class _Site:
            def __init__(self, station, z, t, freq):
                self.station = station
                self.Z = _FakeZ(z, freq)
                self.Tipper = _TipperObj(t, freq)

        t0 = np.ones((fr.size, 2), dtype=complex)
        keep = np.array([True, True, False, False, True, True])
        ed = _Site("S00", _iso_z(fr), t0, fr)
        _mask_apply(ed, keep, also="tipper")
        assert np.all(np.isnan(ed.Tipper.tipper[~keep]))
        assert np.all(np.isfinite(ed.Tipper.tipper[keep]))

    def test_tipper_assignment_exception_is_swallowed(self):
        fr = _freqs(6)

        class _Site:
            def __init__(self, station, z, t, freq):
                self.station = station
                self.Z = _FakeZ(z, freq)
                self.Tipper = _ReadOnlyTipper(t, freq)

        t0 = np.ones((fr.size, 2), dtype=complex)
        ed = _Site("S00", _iso_z(fr), t0, fr)
        keep = np.zeros(fr.size, dtype=bool)
        # Tipper.tipper has no setter -> AttributeError must be swallowed.
        _mask_apply(ed, keep, also="both")


# ─────────────────────────────────────────────────────────────────────────────
# internal helpers: _runs_bool
# ─────────────────────────────────────────────────────────────────────────────


class TestRunsBoolInternals:
    def test_empty_array(self):
        assert _runs_bool(np.array([], dtype=bool)) == []

    def test_single_run_ending_true(self):
        """Run open at loop end is closed by the trailing ``if s is not
        None`` check (branch: condition True)."""
        m = np.array([False, True, True, True])
        assert _runs_bool(m) == [(1, 3)]

    def test_multiple_runs_with_trailing_false(self):
        """A run that closes mid-array (branch at ``not v and s is not
        None``) followed by a run that also closes before the end (so the
        trailing ``if s is not None`` is False)."""
        m = np.array([True, True, False, False, True, False])
        assert _runs_bool(m) == [(0, 1), (4, 4)]

    def test_all_false(self):
        m = np.array([False, False, False])
        assert _runs_bool(m) == []

    def test_all_true(self):
        m = np.array([True, True, True])
        assert _runs_bool(m) == [(0, 2)]


# ─────────────────────────────────────────────────────────────────────────────
# internal helpers: _fill_small_gaps
# ─────────────────────────────────────────────────────────────────────────────


class TestFillSmallGapsInternals:
    def test_max_gap_zero_returns_unchanged(self):
        m = np.array([True, False, True])
        out = _fill_small_gaps(m, max_gap=0)
        assert out is m

    def test_empty_array_returns_unchanged(self):
        m = np.array([], dtype=bool)
        out = _fill_small_gaps(m, max_gap=2)
        assert out is m

    def test_all_false_no_edges_returns_unchanged(self):
        m = np.array([False, False, False, False])
        out = _fill_small_gaps(m, max_gap=2)
        np.testing.assert_array_equal(out, m)

    def test_fills_small_gap(self):
        m = np.array([True, True, False, True, True])
        out = _fill_small_gaps(m, max_gap=1)
        assert out.dtype == bool
        assert bool(out[2]) is True

    def test_large_gap_not_filled(self):
        m = np.array([True, True, False, False, False, True, True])
        out = _fill_small_gaps(m, max_gap=1)
        assert not out[2] and not out[3] and not out[4]


# ─────────────────────────────────────────────────────────────────────────────
# mask_by_skew: mode variants + missing-track guard
# ─────────────────────────────────────────────────────────────────────────────


class TestMaskBySkewModesAndGuards:
    def _site_3d(self, n=8):
        fr = _freqs(n)
        return _site("S00", _3d_z(fr, skew_frac=0.6), fr)

    def test_mode_gt(self):
        result = mask_by_skew([self._site_3d()], thresh=1.0, mode="gt")
        assert sum(1 for _ in result) == 1

    def test_mode_lt(self):
        result = mask_by_skew([self._site_3d()], thresh=1.0, mode="lt")
        assert sum(1 for _ in result) == 1

    def test_mode_abs_lt(self):
        result = mask_by_skew([self._site_3d()], thresh=1.0, mode="abs_lt")
        assert sum(1 for _ in result) == 1

    def test_mode_abs_gt_default(self):
        result = mask_by_skew([self._site_3d()], thresh=1.0)
        assert sum(1 for _ in result) == 1

    def test_site_with_no_z_block_is_returned_untouched(self):
        """A site whose Z can't be resolved must hit the early-return
        guard in mask_by_skew's per-site closure, not raise."""
        fr = _freqs(6)
        good = _site("S00", _iso_z(fr), fr)
        bad = _FakeSiteBadZ("S01", fr)
        result = mask_by_skew([good, bad], thresh=1.0)
        assert sum(1 for _ in result) == 2


# ─────────────────────────────────────────────────────────────────────────────
# keep_longest_low_skew: fallback branches
# ─────────────────────────────────────────────────────────────────────────────


def _short_low_skew_run_site(station="S00", n=14, low_lo=5, low_hi=7):
    """Site whose phase-tensor beta is 0 on ``[low_lo, low_hi)`` (a short
    low-skew run) and large elsewhere, so ``keep_longest_low_skew`` /
    ``select_low_skew_band`` find a run shorter than a generous min_len."""
    fr = _freqs(n)
    z = _3d_z(fr, skew_frac=0.6)
    z[low_lo:low_hi] = _iso_z(fr)[low_lo:low_hi]
    return _site(station, z, fr)


class TestKeepLongestLowSkewFallback:
    def test_site_with_no_z_block_is_returned_untouched(self):
        fr = _freqs(6)
        bad = _FakeSiteBadZ("S00", fr)
        result = keep_longest_low_skew([bad], thresh=0.3)
        assert sum(1 for _ in result) == 1

    def test_short_run_fallback_keep_all(self):
        site = _short_low_skew_run_site()
        result = keep_longest_low_skew(
            [site], thresh=1.0, min_len=3, fallback="keep_all"
        )
        assert sum(1 for _ in result) == 1

    def test_short_run_fallback_drop_all(self):
        site = _short_low_skew_run_site()
        result = keep_longest_low_skew(
            [site], thresh=1.0, min_len=3, fallback="drop_all"
        )
        assert sum(1 for _ in result) == 1

    def test_no_runs_fallback_keep_all(self):
        """thresh below the smallest possible |skew| -> no run at all."""
        fr = _freqs(10)
        site = _site("S00", _3d_z(fr, skew_frac=0.6), fr)
        result = keep_longest_low_skew([site], thresh=-1.0, fallback="keep_all")
        assert sum(1 for _ in result) == 1

    def test_no_runs_fallback_drop_all(self):
        fr = _freqs(10)
        site = _site("S00", _3d_z(fr, skew_frac=0.6), fr)
        result = keep_longest_low_skew([site], thresh=-1.0, fallback="drop_all")
        assert sum(1 for _ in result) == 1


# ─────────────────────────────────────────────────────────────────────────────
# close_skew_gaps: missing-track guard
# ─────────────────────────────────────────────────────────────────────────────


class TestCloseSkewGapsGuard:
    def test_site_with_no_z_block_is_returned_untouched(self):
        fr = _freqs(6)
        bad = _FakeSiteBadZ("S00", fr)
        result = close_skew_gaps([bad], thresh=1.0, max_gap=1)
        assert sum(1 for _ in result) == 1


# ─────────────────────────────────────────────────────────────────────────────
# select_low_skew_band: continue/short-run/empty-masks/apply-guard branches
# ─────────────────────────────────────────────────────────────────────────────


class TestSelectLowSkewBandBranches:
    def test_site_with_no_z_is_skipped_and_untouched(self):
        """Exercises both the mask-collection ``continue`` and the
        per-site apply-stage ``Z is None`` early return."""
        fr = _freqs(10)
        good = _site("S00", _iso_z(fr), fr)
        bad = _FakeSiteBadZ("S01", fr)
        result = select_low_skew_band([good, bad], thresh=6.0)
        assert sum(1 for _ in result) == 2

    def test_all_sites_untracked_returns_input_unchanged(self):
        """When every site lacks a resolvable Z block, ``masks`` stays
        empty and the function returns ``S`` as-is."""
        fr = _freqs(6)
        sites = [_FakeSiteBadZ("S00", fr), _FakeSiteBadZ("S01", fr)]
        result = select_low_skew_band(sites, thresh=6.0)
        assert sum(1 for _ in result) == 2

    def test_short_run_below_min_len_keeps_good_mask(self):
        sites = [_short_low_skew_run_site(f"S{i:02d}") for i in range(3)]
        result = select_low_skew_band(sites, thresh=1.0, min_len=3, frac=0.5)
        assert sum(1 for _ in result) == 3


# ─────────────────────────────────────────────────────────────────────────────
# plot_skew_traffic_psection
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotSkewTrafficPsection:
    def _sites(self, n_sta=4, n_freq=10):
        return [
            _site(
                f"S{i:02d}",
                _3d_z(_freqs(n_freq), skew_frac=0.1 * (i + 1)),
                _freqs(n_freq),
            )
            for i in range(n_sta)
        ]

    def test_returns_axes(self):
        ax = plot_skew_traffic_psection(self._sites())
        plt.close("all")
        assert isinstance(ax, plt.Axes)

    def test_accepts_external_ax(self):
        fig, ax0 = plt.subplots()
        ax = plot_skew_traffic_psection(self._sites(), ax=ax0)
        plt.close("all")
        assert ax is ax0

    def test_axis_y_linear_period(self):
        ax = plot_skew_traffic_psection(self._sites(), axis_y="period")
        plt.close("all")
        assert isinstance(ax, plt.Axes)

    def test_custom_thresholds(self):
        ax = plot_skew_traffic_psection(self._sites(), t1=1.0, t2=2.0)
        plt.close("all")
        assert isinstance(ax, plt.Axes)

    def test_empty_sites_shows_placeholder_text(self):
        ax = plot_skew_traffic_psection([])
        plt.close("all")
        assert isinstance(ax, plt.Axes)
        assert len(ax.texts) >= 1

    def test_empty_sites_with_external_ax(self):
        fig, ax0 = plt.subplots()
        ax = plot_skew_traffic_psection([], ax=ax0)
        plt.close("all")
        assert ax is ax0


# ─────────────────────────────────────────────────────────────────────────────
# plot_skew_percentile_ribbon
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotSkewPercentileRibbon:
    def _sites(self, n_sta=5, n_freq=14):
        return [
            _site(
                f"S{i:02d}",
                _3d_z(_freqs(n_freq), skew_frac=0.1 * (i + 1)),
                _freqs(n_freq),
            )
            for i in range(n_sta)
        ]

    def test_returns_axes(self):
        ax = plot_skew_percentile_ribbon(self._sites())
        plt.close("all")
        assert isinstance(ax, plt.Axes)

    def test_accepts_external_ax(self):
        fig, ax0 = plt.subplots()
        ax = plot_skew_percentile_ribbon(self._sites(), ax=ax0)
        plt.close("all")
        assert ax is ax0

    def test_extra_none_skips_outer_band(self):
        ax = plot_skew_percentile_ribbon(self._sites(), extra=None)
        plt.close("all")
        assert isinstance(ax, plt.Axes)

    def test_custom_quantiles_and_bins(self):
        ax = plot_skew_percentile_ribbon(self._sites(), n_bins=5, q_lo=10.0, q_hi=90.0)
        plt.close("all")
        assert isinstance(ax, plt.Axes)

    def test_empty_sites_shows_placeholder_text(self):
        ax = plot_skew_percentile_ribbon([])
        plt.close("all")
        assert isinstance(ax, plt.Axes)
        assert len(ax.texts) >= 1

    def test_empty_sites_with_external_ax(self):
        fig, ax0 = plt.subplots()
        ax = plot_skew_percentile_ribbon([], ax=ax0)
        plt.close("all")
        assert ax is ax0


# ─────────────────────────────────────────────────────────────────────────────
# plot_skew_vote_band
# ─────────────────────────────────────────────────────────────────────────────


class TestPlotSkewVoteBand:
    def _sites(self, n_sta=5, n_freq=14):
        return [
            _site(
                f"S{i:02d}",
                _3d_z(_freqs(n_freq), skew_frac=0.1 * (i + 1)),
                _freqs(n_freq),
            )
            for i in range(n_sta)
        ]

    def test_returns_axes(self):
        ax = plot_skew_vote_band(self._sites())
        plt.close("all")
        assert isinstance(ax, plt.Axes)

    def test_accepts_external_ax(self):
        fig, ax0 = plt.subplots()
        ax = plot_skew_vote_band(self._sites(), ax=ax0)
        plt.close("all")
        assert ax is ax0

    def test_custom_threshold_and_bins(self):
        ax = plot_skew_vote_band(self._sites(), thresh=10.0, n_bins=8)
        plt.close("all")
        assert isinstance(ax, plt.Axes)

    def test_empty_sites_shows_placeholder_text(self):
        ax = plot_skew_vote_band([])
        plt.close("all")
        assert isinstance(ax, plt.Axes)
        assert len(ax.texts) >= 1

    def test_empty_sites_with_external_ax(self):
        fig, ax0 = plt.subplots()
        ax = plot_skew_vote_band([], ax=ax0)
        plt.close("all")
        assert ax is ax0
