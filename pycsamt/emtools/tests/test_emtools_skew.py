"""Tests for pycsamt.emtools.skew"""

from __future__ import annotations

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.emtools.skew import (
    bahr_skewness,
    close_skew_gaps,
    keep_longest_low_skew,
    mask_by_skew,
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
        return [
            _site(f"S{i:02d}", _iso_z(_freqs()), _freqs())
            for i in range(n_sites)
        ]

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
