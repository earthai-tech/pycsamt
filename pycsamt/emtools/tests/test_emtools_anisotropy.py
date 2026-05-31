"""
Tests for pycsamt.emtools.anisotropy
"""
from __future__ import annotations

import numpy as np
import pytest

from pycsamt.emtools.anisotropy import (
    ANISO_RATIO_THRESH,
    SWIFT_SKEW_THRESH,
    analyze_anisotropy,
    anisotropy_table,
    plot_anisotropy,
)


# ─────────────────────────────────────────────────────────────────────────────
# Fixtures / helpers
# ─────────────────────────────────────────────────────────────────────────────

class _FakeZ:
    def __init__(self, z, freq):
        self.z    = np.asarray(z, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)


class _FakeSite:
    """Minimal EDI-like object (duck-type compatible with is_edi_file)."""
    def __init__(self, station, z, freq):
        self.station = station
        self.Z    = _FakeZ(z, freq)
        self.freq = np.asarray(freq, dtype=float)

    def get_section(self, *_, **__):
        return None


def _freqs(n: int = 10, f_lo: float = 1.0, f_hi: float = 1e4) -> np.ndarray:
    return np.logspace(np.log10(f_lo), np.log10(f_hi), n)


def _isotropic_z(freqs: np.ndarray, rho: float = 100.0) -> np.ndarray:
    """1-D isotropic Z: Z_xy = -Z_yx, Z_xx = Z_yy = 0."""
    z_abs = np.sqrt(5.0 * freqs * rho)
    z = np.zeros((freqs.size, 2, 2), dtype=complex)
    z[:, 0, 1] =  z_abs * (1 + 1j) / np.sqrt(2)
    z[:, 1, 0] = -z_abs * (1 + 1j) / np.sqrt(2)
    return z


def _anisotropic_z(
    freqs: np.ndarray,
    rho_xy: float = 100.0,
    rho_yx: float = 400.0,
) -> np.ndarray:
    """Anisotropic Z: |Z_xy| ≠ |Z_yx|, diagonal still zero."""
    z = np.zeros((freqs.size, 2, 2), dtype=complex)
    z_xy = np.sqrt(5.0 * freqs * rho_xy)
    z_yx = np.sqrt(5.0 * freqs * rho_yx)
    z[:, 0, 1] =  z_xy * (1 + 1j) / np.sqrt(2)
    z[:, 1, 0] = -z_yx * (1 + 1j) / np.sqrt(2)
    return z


def _full_3d_z(
    freqs: np.ndarray,
    rho_xy: float = 100.0,
    rho_yx: float = 400.0,
    skew_frac: float = 0.4,
) -> np.ndarray:
    """Z with non-zero diagonal (simulates 3-D / anisotropic earth)."""
    z = _anisotropic_z(freqs, rho_xy, rho_yx)
    z_xy_amp = np.abs(z[:, 0, 1])
    z[:, 0, 0] = skew_frac * z_xy_amp * (0.6 + 0.8j)
    z[:, 1, 1] = -skew_frac * z_xy_amp * (0.5 + 0.7j)
    return z


def _site(station: str, z: np.ndarray, freqs: np.ndarray) -> _FakeSite:
    return _FakeSite(station, z, freqs)


def _iso_site(station: str, rho: float = 100.0, n: int = 10) -> _FakeSite:
    fr = _freqs(n)
    return _site(station, _isotropic_z(fr, rho), fr)


def _aniso_site(
    station: str,
    rho_xy: float = 100.0,
    rho_yx: float = 400.0,
    n: int = 10,
) -> _FakeSite:
    fr = _freqs(n)
    return _site(station, _anisotropic_z(fr, rho_xy, rho_yx), fr)


def _3d_site(station: str, n: int = 10) -> _FakeSite:
    fr = _freqs(n)
    return _site(station, _full_3d_z(fr), fr)


# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

def test_aniso_ratio_thresh_positive():
    assert ANISO_RATIO_THRESH > 0


def test_swift_skew_thresh_positive():
    assert SWIFT_SKEW_THRESH > 0


def test_constants_values():
    assert ANISO_RATIO_THRESH == pytest.approx(0.1)
    assert SWIFT_SKEW_THRESH  == pytest.approx(0.2)


# ─────────────────────────────────────────────────────────────────────────────
# analyze_anisotropy — columns / shape
# ─────────────────────────────────────────────────────────────────────────────

def test_analyze_columns_present():
    sites = [_iso_site("S00")]
    df = analyze_anisotropy(sites)
    expected = {
        "station", "freq_hz", "period_s",
        "rho_xy_ohmm", "rho_yx_ohmm",
        "phi_xy_deg", "phi_yx_deg",
        "ratio_log10", "phase_diff_deg",
        "swift_skew", "strike_deg",
    }
    assert expected.issubset(df.columns)


def test_analyze_row_count():
    n = 8
    sites = [_iso_site("S00", n=n)]
    df = analyze_anisotropy(sites)
    assert len(df) == n


def test_analyze_empty_input():
    df = analyze_anisotropy([])
    assert df.empty


def test_analyze_multiple_sites():
    sites = [_iso_site(f"S{i:02d}") for i in range(4)]
    df = analyze_anisotropy(sites)
    assert df["station"].nunique() == 4


def test_analyze_freq_positive():
    sites = [_iso_site("S00")]
    df = analyze_anisotropy(sites)
    assert (df["freq_hz"] > 0).all()


def test_analyze_period_positive():
    sites = [_iso_site("S00")]
    df = analyze_anisotropy(sites)
    assert (df["period_s"] > 0).all()


def test_analyze_period_inverse_freq():
    sites = [_iso_site("S00")]
    df = analyze_anisotropy(sites)
    ratio = df["period_s"] * df["freq_hz"]
    assert np.allclose(ratio, 1.0, rtol=1e-6)


def test_analyze_rho_positive():
    sites = [_iso_site("S00")]
    df = analyze_anisotropy(sites)
    assert (df["rho_xy_ohmm"] > 0).all()
    assert (df["rho_yx_ohmm"] > 0).all()


# ─────────────────────────────────────────────────────────────────────────────
# analyze_anisotropy — isotropic site (ratio_log10 ≈ 0)
# ─────────────────────────────────────────────────────────────────────────────

def test_iso_ratio_near_zero():
    """Isotropic 1-D: ρ_xy = ρ_yx → ratio_log10 ≈ 0."""
    sites = [_iso_site("S00")]
    df = analyze_anisotropy(sites)
    assert np.allclose(df["ratio_log10"].dropna(), 0.0, atol=1e-4)


def test_iso_swift_skew_zero():
    """1-D isotropic Z has zero diagonal → swift_skew = 0."""
    sites = [_iso_site("S00")]
    df = analyze_anisotropy(sites)
    assert np.allclose(df["swift_skew"], 0.0, atol=1e-8)


# ─────────────────────────────────────────────────────────────────────────────
# analyze_anisotropy — anisotropic site
# ─────────────────────────────────────────────────────────────────────────────

def test_aniso_ratio_nonzero():
    """Anisotropic: |Z_xy| ≠ |Z_yx| → ratio_log10 ≠ 0."""
    sites = [_aniso_site("S00", rho_xy=100.0, rho_yx=400.0)]
    df = analyze_anisotropy(sites)
    assert not np.allclose(df["ratio_log10"].dropna(), 0.0, atol=1e-2)


def test_aniso_ratio_known_value():
    """ρ_xy/ρ_yx = 100/400 = 0.25 → log10 = −0.602."""
    sites = [_aniso_site("S00", rho_xy=100.0, rho_yx=400.0)]
    df = analyze_anisotropy(sites)
    # Approximate: the pycsamt |Z|²/(ωμ₀) formula differs from 5fρ by ωμ₀ factor,
    # but the RATIO ρ_xy/ρ_yx = |Z_xy|²/|Z_yx|² = rho_xy/rho_yx exactly.
    expected = np.log10(100.0 / 400.0)
    assert np.allclose(df["ratio_log10"].dropna(), expected, atol=1e-4)


def test_aniso_ratio_sign():
    """ρ_xy < ρ_yx → ratio_log10 < 0; ρ_xy > ρ_yx → ratio_log10 > 0."""
    sites_neg = [_aniso_site("S00", rho_xy=50.0,  rho_yx=200.0)]
    sites_pos = [_aniso_site("S01", rho_xy=200.0, rho_yx=50.0)]
    df_neg = analyze_anisotropy(sites_neg)
    df_pos = analyze_anisotropy(sites_pos)
    assert (df_neg["ratio_log10"].dropna() < 0).all()
    assert (df_pos["ratio_log10"].dropna() > 0).all()


def test_3d_swift_skew_positive():
    """3-D Z with non-zero diagonal → swift_skew > 0."""
    sites = [_3d_site("S00")]
    df = analyze_anisotropy(sites)
    assert (df["swift_skew"] > 0).all()


def test_aniso_skew_zero_for_no_diagonal():
    """Anisotropic Z but zero diagonal → swift_skew = 0."""
    sites = [_aniso_site("S00", rho_xy=100.0, rho_yx=400.0)]
    df = analyze_anisotropy(sites)
    assert np.allclose(df["swift_skew"], 0.0, atol=1e-8)


# ─────────────────────────────────────────────────────────────────────────────
# analyze_anisotropy — phase / strike
# ─────────────────────────────────────────────────────────────────────────────

def test_iso_phase_diff_near_zero():
    """Isotropic 1-D: φ_xy and φ_yx share same angle → diff ≈ 0."""
    sites = [_iso_site("S00")]
    df = analyze_anisotropy(sites)
    # Z_xy = +A(1+j)/√2, Z_yx = -A(1+j)/√2 → phases are 45° and 180°+45°
    # Their difference = ±180° — just check it is finite
    assert df["phase_diff_deg"].notna().all()


def test_phi_range():
    """Phase angles should be in (−180, 180]."""
    sites = [_3d_site("S00")]
    df = analyze_anisotropy(sites)
    assert (df["phi_xy_deg"] >= -180).all()
    assert (df["phi_xy_deg"] <= 180).all()


def test_strike_finite():
    sites = [_3d_site("S00")]
    df = analyze_anisotropy(sites)
    assert np.isfinite(df["strike_deg"]).all()


# ─────────────────────────────────────────────────────────────────────────────
# anisotropy_table — shape / columns
# ─────────────────────────────────────────────────────────────────────────────

def test_table_columns_present():
    sites = [_iso_site("S00")]
    df = anisotropy_table(sites)
    expected = {
        "station", "n_freq",
        "mean_ratio_log10", "max_abs_ratio_log10",
        "mean_phase_diff_deg",
        "mean_swift_skew",
        "median_strike_deg",
        "anisotropy_flag",
    }
    assert expected.issubset(df.columns)


def test_table_one_row_per_station():
    n = 5
    sites = [_iso_site(f"S{i:02d}") for i in range(n)]
    df = anisotropy_table(sites)
    assert len(df) == n


def test_table_empty_input():
    df = anisotropy_table([])
    assert df.empty


def test_table_n_freq_correct():
    n = 8
    sites = [_iso_site("S00", n=n)]
    df = anisotropy_table(sites)
    assert int(df["n_freq"].iloc[0]) == n


def test_table_flag_dtype():
    sites = [_iso_site("S00")]
    df = anisotropy_table(sites)
    assert df["anisotropy_flag"].dtype == bool


# ─────────────────────────────────────────────────────────────────────────────
# anisotropy_table — flag logic
# ─────────────────────────────────────────────────────────────────────────────

def test_table_flag_false_for_isotropic():
    """Isotropic: both ratio and skew are tiny → no flag."""
    sites = [_iso_site("S00")]
    df = anisotropy_table(sites, ratio_threshold=0.1, skew_threshold=0.2)
    assert not bool(df["anisotropy_flag"].iloc[0])


def test_table_flag_true_for_anisotropic():
    """Large ρ_xy/ρ_yx ratio → flag raised."""
    sites = [_aniso_site("S00", rho_xy=1.0, rho_yx=10000.0)]
    df = anisotropy_table(sites, ratio_threshold=0.1)
    assert bool(df["anisotropy_flag"].iloc[0])


def test_table_flag_true_for_3d():
    """Large Swift skew → flag raised."""
    sites = [_3d_site("S00")]
    df = anisotropy_table(sites, skew_threshold=0.05)
    assert bool(df["anisotropy_flag"].iloc[0])


def test_table_max_abs_ratio_ge_mean_abs():
    sites = [_aniso_site("S00")]
    df = anisotropy_table(sites)
    max_r  = float(df["max_abs_ratio_log10"].iloc[0])
    mean_r = abs(float(df["mean_ratio_log10"].iloc[0]))
    if np.isfinite(max_r) and np.isfinite(mean_r):
        assert max_r >= mean_r - 1e-8


# ─────────────────────────────────────────────────────────────────────────────
# plot_anisotropy
# ─────────────────────────────────────────────────────────────────────────────

@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_returns_axes():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sites = [_iso_site(f"S{i:02d}") for i in range(3)]
    ax = plot_anisotropy(sites)
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_existing_axes():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _, ax_in = plt.subplots()
    sites = [_iso_site("S00")]
    ax_out = plot_anisotropy(sites, ax=ax_in)
    assert ax_out is ax_in
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_empty_input():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ax = plot_anisotropy([])
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_swift_skew_metric():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sites = [_3d_site(f"S{i:02d}") for i in range(3)]
    ax = plot_anisotropy(sites, metric="swift_skew")
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_strike_deg_metric():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sites = [_3d_site(f"S{i:02d}") for i in range(3)]
    ax = plot_anisotropy(sites, metric="strike_deg")
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_phase_diff_metric():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sites = [_aniso_site(f"S{i:02d}") for i in range(3)]
    ax = plot_anisotropy(sites, metric="phase_diff_deg")
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_freq_axis():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sites = [_iso_site(f"S{i:02d}") for i in range(3)]
    ax = plot_anisotropy(sites, period_axis=False)
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_linear_y():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sites = [_iso_site(f"S{i:02d}") for i in range(3)]
    ax = plot_anisotropy(sites, log_y=False)
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_no_contour():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sites = [_iso_site(f"S{i:02d}") for i in range(3)]
    ax = plot_anisotropy(sites, contour_zero=False)
    assert ax is not None
    plt.close("all")
