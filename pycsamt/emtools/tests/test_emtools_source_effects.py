"""
Tests for pycsamt.emtools.source_effects
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.emtools.source_effects import (
    BETA_THRESH_PCT,
    detect_source_overprint,
    overprint_beta,
    plot_overprint_section,
    source_overprint_table,
)

# ─────────────────────────────────────────────────────────────────────────────
# Fixtures
# ─────────────────────────────────────────────────────────────────────────────


class _FakeZ:
    def __init__(self, z, freq):
        self.z = np.asarray(z, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)


class _FakeSite:
    """Minimal EDI-like object (duck-type compatible with is_edi_file)."""

    def __init__(self, station, z, freq, offset=None):
        self.station = station
        self.Z = _FakeZ(z, freq)
        self.freq = np.asarray(freq, dtype=float)
        if offset is not None:
            self.source_offset = float(offset)

    def get_section(self, *_, **__):
        return None


def _make_z(freqs: np.ndarray, rho: float) -> np.ndarray:
    """Z with constant ρ_a = rho; pycsamt convention |Z| = √(5 f ρ)."""
    z_abs = np.sqrt(5.0 * freqs * rho)
    z = np.zeros((freqs.size, 2, 2), dtype=complex)
    z[:, 0, 1] = z_abs * (1 + 1j) / np.sqrt(2)
    z[:, 1, 0] = -z_abs * (1 + 1j) / np.sqrt(2)
    return z


def _site(
    station: str,
    rho: float = 100.0,
    f_lo: float = 0.1,
    f_hi: float = 1e4,
    n: int = 10,
    offset: float = 5000.0,
) -> _FakeSite:
    fr = np.logspace(np.log10(f_lo), np.log10(f_hi), n)
    return _FakeSite(station, _make_z(fr, rho), fr, offset=offset)


def _site_no_offset(station: str, rho: float = 100.0, n: int = 8) -> _FakeSite:
    fr = np.logspace(0, 4, n)
    return _FakeSite(station, _make_z(fr, rho), fr)


# ─────────────────────────────────────────────────────────────────────────────
# overprint_beta — pure math
# ─────────────────────────────────────────────────────────────────────────────


def test_beta_scalar_positive():
    b = overprint_beta(100.0, 10.0, 5000.0)
    assert float(b) >= 0.0


def test_beta_decreases_with_frequency():
    """Higher frequency → farther from source in wavelengths → lower β."""
    b_lo = float(overprint_beta(100.0, 1.0, 5000.0))
    b_hi = float(overprint_beta(100.0, 1000.0, 5000.0))
    assert b_lo > b_hi, f"β at 1 Hz={b_lo:.2f} should exceed β at 1 kHz={b_hi:.2f}"


def test_beta_decreases_with_offset():
    """Larger offset → smaller ground-wave contamination."""
    b_near = float(overprint_beta(100.0, 10.0, 2000.0))
    b_far = float(overprint_beta(100.0, 10.0, 10_000.0))
    assert b_near > b_far, f"near β={b_near:.2f} should exceed far β={b_far:.2f}"


def test_beta_array_shape():
    freqs = np.logspace(0, 4, 6)
    beta = overprint_beta(100.0, freqs, 5000.0)
    assert beta.shape == (6,)


def test_beta_array_monotone():
    freqs = np.logspace(0, 4, 8)
    beta = overprint_beta(100.0, freqs, 5000.0)
    # β should decrease monotonically with frequency
    assert (np.diff(beta) < 0).all()


def test_beta_zero_offset_nan():
    b = overprint_beta(100.0, 10.0, 0.0)
    assert np.isnan(float(b))


def test_beta_negative_offset_nan():
    b = overprint_beta(100.0, 10.0, -500.0)
    assert np.isnan(float(b))


def test_beta_broadcast():
    rhos = np.array([10.0, 100.0, 1000.0])
    result = overprint_beta(rhos, 10.0, 5000.0)
    assert result.shape == (3,)
    assert (result > 0).all()


def test_beta_threshold_constant():
    assert BETA_THRESH_PCT == pytest.approx(3.0)


def test_beta_transition_zone_high():
    """In the transition zone (kr ~ 1), β should exceed 3 %."""
    # kr ≈ 1: δ_B = 356√(ρ/f) = offset → f = ρ(356/offset)²
    # ρ=100, offset=5000 → f ≈ 100 × (356/5000)² ≈ 0.507 Hz
    b = float(overprint_beta(100.0, 0.5, 5000.0))
    assert b > 3.0, f"transition-zone β={b:.2f} should exceed 3 %"


def test_beta_far_field_low():
    """In the far field (kr >> 3), β should be well below 3 %."""
    # kr ≈ 20: large offset + high frequency
    b = float(overprint_beta(100.0, 1000.0, 5000.0))
    assert b < 3.0, f"far-field β={b:.2f} should be below 3 %"


# ─────────────────────────────────────────────────────────────────────────────
# detect_source_overprint — sites-based
# ─────────────────────────────────────────────────────────────────────────────


def test_detect_columns():
    sites = [_site(f"S{i:02d}") for i in range(3)]
    df = detect_source_overprint(sites)
    expected = {
        "station",
        "freq_hz",
        "period_s",
        "offset_m",
        "rho_a_ohmm",
        "kr",
        "beta_pct",
        "overprint_flag",
    }
    assert expected.issubset(df.columns)


def test_detect_row_count():
    n_freq = 8
    sites = [_site("S00", n=n_freq)]
    df = detect_source_overprint(sites)
    assert len(df) == n_freq


def test_detect_empty_input():
    df = detect_source_overprint([])
    assert df.empty


def test_detect_beta_positive():
    sites = [_site("S00")]
    df = detect_source_overprint(sites)
    valid = df["beta_pct"].dropna()
    assert (valid >= 0).all()


def test_detect_kr_positive():
    sites = [_site("S00")]
    df = detect_source_overprint(sites)
    valid = df["kr"].dropna()
    assert (valid > 0).all()


def test_detect_no_offset_nan_beta():
    """Sites without offset should produce NaN beta_pct."""
    sites = [_site_no_offset("S00")]
    df = detect_source_overprint(sites)
    assert df["beta_pct"].isna().all()


def test_detect_offset_dict():
    sites = [_site_no_offset("S00", n=6), _site_no_offset("S01", n=6)]
    df = detect_source_overprint(sites, source_offset={"S00": 3000.0, "S01": 8000.0})
    grp = df.groupby("station")["beta_pct"]
    beta_near = grp.mean()["S00"]
    beta_far = grp.mean()["S01"]
    assert beta_near > beta_far


def test_detect_overprint_flag_dtype():
    sites = [_site("S00")]
    df = detect_source_overprint(sites)
    assert df["overprint_flag"].dtype == bool


def test_detect_multiple_sites():
    sites = [_site(f"S{i:02d}") for i in range(5)]
    df = detect_source_overprint(sites)
    assert df["station"].nunique() == 5


# ─────────────────────────────────────────────────────────────────────────────
# source_overprint_table — per-station summary
# ─────────────────────────────────────────────────────────────────────────────


def test_table_columns():
    sites = [_site(f"S{i:02d}") for i in range(3)]
    df = source_overprint_table(sites)
    expected = {
        "station",
        "n_freq",
        "offset_m",
        "beta_max_pct",
        "beta_mean_pct",
        "n_overprint",
        "overprint_frac",
        "lf_slope",
        "hf_slope",
        "slope_delta",
        "overprint_flag",
    }
    assert expected.issubset(df.columns)


def test_table_one_row_per_station():
    n = 4
    sites = [_site(f"S{i:02d}") for i in range(n)]
    df = source_overprint_table(sites)
    assert len(df) == n


def test_table_empty_input():
    df = source_overprint_table([])
    assert df.empty


def test_table_beta_max_ge_mean():
    sites = [_site("S00", n=10)]
    df = source_overprint_table(sites)
    max_b = df["beta_max_pct"].iloc[0]
    mean_b = df["beta_mean_pct"].iloc[0]
    if np.isfinite(max_b) and np.isfinite(mean_b):
        assert max_b >= mean_b


def test_table_overprint_frac_in_01():
    sites = [_site(f"S{i:02d}") for i in range(3)]
    df = source_overprint_table(sites)
    frac = df["overprint_frac"].dropna()
    assert ((frac >= 0) & (frac <= 1)).all()


def test_table_slope_columns_present():
    sites = [_site("S00", n=12)]
    df = source_overprint_table(sites, f_split=1.0)
    assert "lf_slope" in df.columns
    assert "hf_slope" in df.columns
    assert "slope_delta" in df.columns


# ─────────────────────────────────────────────────────────────────────────────
# plot_overprint_section
# ─────────────────────────────────────────────────────────────────────────────


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_returns_axes():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sites = [_site(f"S{i:02d}") for i in range(3)]
    ax = plot_overprint_section(sites)
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_empty_input():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ax = plot_overprint_section([])
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_existing_axes():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _, ax_in = plt.subplots()
    sites = [_site(f"S{i:02d}") for i in range(3)]
    ax_out = plot_overprint_section(sites, ax=ax_in)
    assert ax_out is ax_in
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_linear_color():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sites = [_site(f"S{i:02d}") for i in range(3)]
    ax = plot_overprint_section(sites, log_color=False)
    assert ax is not None
    plt.close("all")
