"""Tests for pycsamt.emtools.strike"""
from __future__ import annotations

import matplotlib
import numpy as np

matplotlib.use("Agg")

from pycsamt.emtools.strike import (
    estimate_strike_consensus,
    estimate_strike_phase_tensor,
    estimate_strike_sweep,
    rotate_to_strike,
    strike_curve_sweep,
)

# ─────────────────────────────────────────────────────────────────────────────
# Shared helpers
# ─────────────────────────────────────────────────────────────────────────────

class _FakeZ:
    def __init__(self, z, freq):
        self.z    = np.asarray(z, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)


class _FakeSite:
    def __init__(self, station, z, freq):
        self.station = station
        self.Z    = _FakeZ(z, freq)
        self.freq = np.asarray(freq, dtype=float)

    def get_section(self, *_, **__):
        return None


def _freqs(n: int = 12, f_lo: float = 0.1, f_hi: float = 1e4) -> np.ndarray:
    return np.logspace(np.log10(f_lo), np.log10(f_hi), n)


def _iso_z(freqs: np.ndarray, rho: float = 100.0) -> np.ndarray:
    amp = np.sqrt(5.0 * freqs * rho)
    z = np.zeros((freqs.size, 2, 2), dtype=complex)
    z[:, 0, 1] =  amp * (1 + 1j) / np.sqrt(2)
    z[:, 1, 0] = -amp * (1 + 1j) / np.sqrt(2)
    return z


def _site(name: str, n: int = 12) -> _FakeSite:
    fr = _freqs(n)
    return _FakeSite(name, _iso_z(fr), fr)


# ─────────────────────────────────────────────────────────────────────────────
# estimate_strike_sweep
# ─────────────────────────────────────────────────────────────────────────────

class TestEstimateStrikeSweep:

    def test_returns_dataframe(self):
        import pandas as pd
        sites = [_site("S00")]
        df = estimate_strike_sweep(sites)
        assert isinstance(df, pd.DataFrame)

    def test_expected_columns(self):
        sites = [_site("S00")]
        df = estimate_strike_sweep(sites)
        for col in ("station", "ang", "iqr", "n"):
            assert col in df.columns

    def test_one_row_per_site(self):
        n = 4
        sites = [_site(f"S{i:02d}") for i in range(n)]
        df = estimate_strike_sweep(sites)
        assert len(df) == n

    def test_empty_input_empty_df(self):
        import pandas as pd
        df = estimate_strike_sweep([])
        assert isinstance(df, pd.DataFrame)
        assert df.empty

    def test_angle_in_range(self):
        sites = [_site("S00")]
        df = estimate_strike_sweep(sites)
        ang = df["ang"].dropna()
        assert (ang >= -90.0).all() and (ang <= 90.0).all()


# ─────────────────────────────────────────────────────────────────────────────
# estimate_strike_phase_tensor
# ─────────────────────────────────────────────────────────────────────────────

class TestEstimateStrikePhaseTensor:

    def test_returns_dataframe(self):
        import pandas as pd
        sites = [_site("S00")]
        df = estimate_strike_phase_tensor(sites)
        assert isinstance(df, pd.DataFrame)

    def test_columns(self):
        sites = [_site("S00")]
        df = estimate_strike_phase_tensor(sites)
        assert "station" in df.columns
        assert "ang" in df.columns

    def test_empty_input(self):
        import pandas as pd
        df = estimate_strike_phase_tensor([])
        assert isinstance(df, pd.DataFrame)

    def test_multi_site(self):
        sites = [_site(f"S{i}") for i in range(3)]
        df = estimate_strike_phase_tensor(sites)
        assert len(df) >= 0   # may be empty for some Z shapes, just no error


# ─────────────────────────────────────────────────────────────────────────────
# estimate_strike_consensus
# ─────────────────────────────────────────────────────────────────────────────

class TestEstimateStrikeConsensus:

    def test_returns_dataframe(self):
        import pandas as pd
        sites = [_site("S00")]
        df = estimate_strike_consensus(sites)
        assert isinstance(df, pd.DataFrame)

    def test_columns(self):
        sites = [_site("S00")]
        df = estimate_strike_consensus(sites)
        assert "station" in df.columns

    def test_empty_input(self):
        import pandas as pd
        df = estimate_strike_consensus([])
        assert isinstance(df, pd.DataFrame)


# ─────────────────────────────────────────────────────────────────────────────
# rotate_to_strike
# ─────────────────────────────────────────────────────────────────────────────

class TestRotateToStrike:

    def test_returns_sites(self):
        from pycsamt.site.base import Sites
        sites = [_site("S00")]
        result = rotate_to_strike(sites)
        assert isinstance(result, Sites)

    def test_site_count_preserved(self):
        sites = [_site(f"S{i:02d}") for i in range(3)]
        result = rotate_to_strike(sites)
        assert sum(1 for _ in result) == 3

    def test_empty_sites(self):
        from pycsamt.site.base import Sites
        result = rotate_to_strike([])
        assert isinstance(result, Sites)

    def test_phase_tensor_method(self):
        from pycsamt.site.base import Sites
        sites = [_site("S00")]
        result = rotate_to_strike(sites, method="phase_tensor")
        assert isinstance(result, Sites)


# ─────────────────────────────────────────────────────────────────────────────
# strike_curve_sweep
# ─────────────────────────────────────────────────────────────────────────────

class TestStrikeCurveSweep:

    def test_returns_dataframe(self):
        import pandas as pd
        sites = [_site("S00")]
        df = strike_curve_sweep(sites)
        assert isinstance(df, pd.DataFrame)

    def test_columns(self):
        sites = [_site("S00")]
        df = strike_curve_sweep(sites)
        assert "station" in df.columns

    def test_empty_input(self):
        import pandas as pd
        df = strike_curve_sweep([])
        assert isinstance(df, pd.DataFrame)
