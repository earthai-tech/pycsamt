"""
Tests for pycsamt.emtools.diag
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.emtools.diag import (
    COVERAGE_THRESH,
    coverage_score,
    coverage_table,
    plot_polar_coverage,
    plot_polar_errors,
    plot_width_drift,
    rho_coverage,
    rho_error_stats,
)

# ─────────────────────────────────────────────────────────────────────────────
# Fixtures / helpers
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


def _freqs(n: int = 10, f_lo: float = 1.0, f_hi: float = 1e4) -> np.ndarray:
    return np.logspace(np.log10(f_lo), np.log10(f_hi), n)


def _make_z(freqs: np.ndarray, rho: float = 100.0) -> np.ndarray:
    """1-D isotropic Z: |Z| = sqrt(5 f ρ), Z_xx = Z_yy = 0."""
    z_abs = np.sqrt(5.0 * freqs * rho)
    z = np.zeros((freqs.size, 2, 2), dtype=complex)
    z[:, 0, 1] = z_abs * (1 + 1j) / np.sqrt(2)
    z[:, 1, 0] = -z_abs * (1 + 1j) / np.sqrt(2)
    return z


def _rho_obs_from_z(freqs: np.ndarray, rho: float = 100.0) -> np.ndarray:
    """Compute expected rho_a (xy) matching the Cagniard formula used in diag."""
    z_abs_sq = 5.0 * freqs * rho  # |Z_xy|²
    return 0.2 * z_abs_sq / freqs


def _site(station: str, rho: float = 100.0, n: int = 10) -> _FakeSite:
    fr = _freqs(n)
    return _FakeSite(station, _make_z(fr, rho), fr)


def _wide_bounds(station: str, n: int = 10) -> tuple:
    """Bounds that always cover (lo=1, hi=1e15)."""
    return {station: np.full(n, 1.0)}, {station: np.full(n, 1e15)}


def _tight_above(station: str, freqs: np.ndarray, rho: float = 100.0) -> tuple:
    """Bounds placed above rho_obs → 0 % coverage."""
    r = _rho_obs_from_z(freqs, rho)
    return {station: r * 1.5}, {station: r * 3.0}


def _covering_bounds(station: str, freqs: np.ndarray, rho: float = 100.0) -> tuple:
    """Bounds that bracket rho_obs → 100 % coverage."""
    r = _rho_obs_from_z(freqs, rho)
    return {station: r * 0.5}, {station: r * 2.0}


# ─────────────────────────────────────────────────────────────────────────────
# COVERAGE_THRESH constant
# ─────────────────────────────────────────────────────────────────────────────


def test_coverage_thresh_value():
    assert COVERAGE_THRESH == pytest.approx(0.9)


def test_coverage_thresh_in_01():
    assert 0.0 < COVERAGE_THRESH <= 1.0


# ─────────────────────────────────────────────────────────────────────────────
# coverage_score — pure math
# ─────────────────────────────────────────────────────────────────────────────


def test_cs_all_covered():
    y = np.array([1.0, 2.0, 3.0])
    assert coverage_score(y, 0.0, 4.0) == pytest.approx(1.0)


def test_cs_none_covered():
    y = np.array([1.0, 2.0, 3.0])
    assert coverage_score(y, 5.0, 10.0) == pytest.approx(0.0)


def test_cs_partial():
    y = np.array([1.0, 5.0, 9.0])
    lo = np.array([0.0, 0.0, 0.0])
    hi = np.array([3.0, 3.0, 3.0])
    assert coverage_score(y, lo, hi) == pytest.approx(1.0 / 3.0)


def test_cs_exact_boundary():
    """Boundary values are covered (inclusive)."""
    assert coverage_score(1.0, 1.0, 1.0) == pytest.approx(1.0)


def test_cs_returns_float():
    assert isinstance(coverage_score([1.0], [0.5], [2.0]), float)


def test_cs_empty():
    score = coverage_score([], [], [])
    assert np.isnan(score) or score == pytest.approx(0.0) or score == pytest.approx(1.0)


# ─────────────────────────────────────────────────────────────────────────────
# rho_coverage — columns / shape
# ─────────────────────────────────────────────────────────────────────────────


def test_rc_columns():
    s = _site("S00")
    lo, hi = _wide_bounds("S00")
    df = rho_coverage([s], lo, hi)
    expected = {
        "station",
        "freq_hz",
        "period_s",
        "rho_obs",
        "q_lo",
        "q_hi",
        "covered",
        "width_pct",
    }
    assert expected.issubset(df.columns)


def test_rc_row_count():
    n = 8
    s = _site("S00", n=n)
    lo, hi = _wide_bounds("S00", n)
    df = rho_coverage([s], lo, hi)
    assert len(df) == n


def test_rc_empty_input():
    df = rho_coverage([], {}, {})
    assert df.empty


def test_rc_covered_dtype():
    s = _site("S00")
    lo, hi = _wide_bounds("S00")
    df = rho_coverage([s], lo, hi)
    assert df["covered"].dtype == bool


def test_rc_width_pct_positive():
    s = _site("S00")
    lo, hi = _wide_bounds("S00")
    df = rho_coverage([s], lo, hi)
    assert (df["width_pct"].dropna() > 0).all()


def test_rc_rho_obs_positive():
    s = _site("S00")
    lo, hi = _wide_bounds("S00")
    df = rho_coverage([s], lo, hi)
    assert (df["rho_obs"] > 0).all()


def test_rc_freq_positive():
    s = _site("S00")
    lo, hi = _wide_bounds("S00")
    df = rho_coverage([s], lo, hi)
    assert (df["freq_hz"] > 0).all()


def test_rc_period_inverse_freq():
    s = _site("S00")
    lo, hi = _wide_bounds("S00")
    df = rho_coverage([s], lo, hi)
    assert np.allclose(df["period_s"] * df["freq_hz"], 1.0, rtol=1e-6)


def test_rc_multiple_sites():
    sites = [_site(f"S{i:02d}") for i in range(3)]
    lo = {f"S{i:02d}": np.full(10, 1.0) for i in range(3)}
    hi = {f"S{i:02d}": np.full(10, 1e15) for i in range(3)}
    df = rho_coverage(sites, lo, hi)
    assert df["station"].nunique() == 3


def test_rc_missing_station_skipped():
    s = _site("S00")
    # lo/hi have no key for "S00"
    df = rho_coverage([s], {"OTHER": np.ones(10)}, {"OTHER": np.ones(10) * 1e15})
    assert df.empty


# ─────────────────────────────────────────────────────────────────────────────
# rho_coverage — coverage correctness
# ─────────────────────────────────────────────────────────────────────────────


def test_rc_all_covered_wide_bounds():
    """Wide bounds [1, 1e15] cover every ρ_obs."""
    s = _site("S00")
    lo, hi = _wide_bounds("S00")
    df = rho_coverage([s], lo, hi)
    assert df["covered"].all()


def test_rc_none_covered_above_obs():
    """Bounds placed above rho_obs → 0 % coverage."""
    fr = _freqs(8)
    s = _FakeSite("S00", _make_z(fr), fr)
    lo, hi = _tight_above("S00", fr)
    df = rho_coverage([s], lo, hi)
    assert not df["covered"].any()


def test_rc_all_covered_bracketing_bounds():
    """Bounds that bracket rho_obs → 100 % coverage."""
    fr = _freqs(8)
    s = _FakeSite("S00", _make_z(fr), fr)
    lo, hi = _covering_bounds("S00", fr)
    df = rho_coverage([s], lo, hi)
    assert df["covered"].all()


def test_rc_scalar_bounds_broadcast():
    """Scalar q_lo=1, q_hi=1e15 covers every ρ_obs (always positive)."""
    s = _site("S00")
    df = rho_coverage([s], 1.0, 1e15)
    assert df["covered"].all()


# ─────────────────────────────────────────────────────────────────────────────
# rho_error_stats — columns / shape
# ─────────────────────────────────────────────────────────────────────────────


def test_re_columns():
    fr = _freqs(8)
    s = _FakeSite("S00", _make_z(fr), fr)
    r = _rho_obs_from_z(fr)
    df = rho_error_stats([s], {"S00": r * 1.1})
    expected = {
        "station",
        "freq_hz",
        "period_s",
        "rho_obs",
        "rho_pred",
        "rel_err_pct",
        "abs_err_pct",
    }
    assert expected.issubset(df.columns)


def test_re_row_count():
    n = 8
    fr = _freqs(n)
    s = _FakeSite("S00", _make_z(fr), fr)
    r = _rho_obs_from_z(fr)
    df = rho_error_stats([s], {"S00": r})
    assert len(df) == n


def test_re_empty_input():
    df = rho_error_stats([], {})
    assert df.empty


def test_re_missing_station_skipped():
    s = _site("S00")
    df = rho_error_stats([s], {"OTHER": np.ones(10)})
    assert df.empty


def test_re_zero_error_perfect_model():
    """Model = obs → rel_err_pct = 0 everywhere."""
    fr = _freqs(8)
    s = _FakeSite("S00", _make_z(fr), fr)
    r = _rho_obs_from_z(fr)
    df = rho_error_stats([s], {"S00": r})
    assert np.allclose(df["rel_err_pct"].dropna(), 0.0, atol=1e-4)


def test_re_positive_error_over_prediction():
    """Model > obs → rel_err_pct > 0."""
    fr = _freqs(8)
    s = _FakeSite("S00", _make_z(fr), fr)
    r = _rho_obs_from_z(fr)
    df = rho_error_stats([s], {"S00": r * 2.0})
    assert (df["rel_err_pct"].dropna() > 0).all()


def test_re_negative_error_under_prediction():
    """Model < obs → rel_err_pct < 0."""
    fr = _freqs(8)
    s = _FakeSite("S00", _make_z(fr), fr)
    r = _rho_obs_from_z(fr)
    df = rho_error_stats([s], {"S00": r * 0.5})
    assert (df["rel_err_pct"].dropna() < 0).all()


def test_re_abs_err_non_negative():
    fr = _freqs(8)
    s = _FakeSite("S00", _make_z(fr), fr)
    r = _rho_obs_from_z(fr)
    df = rho_error_stats([s], {"S00": r * 0.5})
    assert (df["abs_err_pct"].dropna() >= 0).all()


def test_re_rel_err_known_value():
    """50 % over-prediction → rel_err_pct = 50."""
    fr = _freqs(8)
    s = _FakeSite("S00", _make_z(fr), fr)
    r = _rho_obs_from_z(fr)
    df = rho_error_stats([s], {"S00": r * 1.5})
    assert np.allclose(df["rel_err_pct"].dropna(), 50.0, atol=1e-3)


def test_re_multiple_sites():
    sites = [_site(f"S{i:02d}") for i in range(3)]
    fr = _freqs(10)
    r = _rho_obs_from_z(fr)
    pred = {f"S{i:02d}": r for i in range(3)}
    df = rho_error_stats(sites, pred)
    assert df["station"].nunique() == 3


# ─────────────────────────────────────────────────────────────────────────────
# coverage_table — shape / columns
# ─────────────────────────────────────────────────────────────────────────────


def test_ct_columns():
    s = _site("S00")
    lo, hi = _wide_bounds("S00")
    df = coverage_table([s], lo, hi)
    expected = {
        "station",
        "n_freq",
        "empirical_cov",
        "mean_width_pct",
        "calibrated_flag",
    }
    assert expected.issubset(df.columns)


def test_ct_one_row_per_station():
    n = 4
    sites = [_site(f"S{i:02d}") for i in range(n)]
    lo = {f"S{i:02d}": np.full(10, 1.0) for i in range(n)}
    hi = {f"S{i:02d}": np.full(10, 1e15) for i in range(n)}
    df = coverage_table(sites, lo, hi)
    assert len(df) == n


def test_ct_empty_input():
    df = coverage_table([], {}, {})
    assert df.empty


def test_ct_flag_dtype():
    s = _site("S00")
    lo, hi = _wide_bounds("S00")
    df = coverage_table([s], lo, hi)
    assert df["calibrated_flag"].dtype == bool


def test_ct_n_freq_correct():
    n = 8
    s = _site("S00", n=n)
    lo, hi = _wide_bounds("S00", n)
    df = coverage_table([s], lo, hi)
    assert int(df["n_freq"].iloc[0]) == n


def test_ct_empirical_cov_in_01():
    s = _site("S00")
    lo, hi = _wide_bounds("S00")
    df = coverage_table([s], lo, hi)
    assert (df["empirical_cov"] >= 0.0).all()
    assert (df["empirical_cov"] <= 1.0).all()


def test_ct_calibrated_true_full_coverage():
    """100 % coverage ≥ nominal=0.9 → flag True."""
    s = _site("S00")
    lo, hi = _wide_bounds("S00")
    df = coverage_table([s], lo, hi, nominal=0.9)
    assert bool(df["calibrated_flag"].iloc[0])


def test_ct_calibrated_false_zero_coverage():
    """0 % coverage < nominal=0.9 → flag False."""
    fr = _freqs(8)
    s = _FakeSite("S00", _make_z(fr), fr)
    lo, hi = _tight_above("S00", fr)
    df = coverage_table([s], lo, hi, nominal=0.9)
    assert not bool(df["calibrated_flag"].iloc[0])


def test_ct_mean_width_positive():
    s = _site("S00")
    lo, hi = _wide_bounds("S00")
    df = coverage_table([s], lo, hi)
    w = df["mean_width_pct"].iloc[0]
    if np.isfinite(w):
        assert w > 0


# ─────────────────────────────────────────────────────────────────────────────
# plot_polar_coverage
# ─────────────────────────────────────────────────────────────────────────────


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_ppc_returns_axes():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sites = [_site(f"S{i:02d}") for i in range(3)]
    ax = plot_polar_coverage(sites, 1.0, 1e15)
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_ppc_empty_input():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ax = plot_polar_coverage([], {}, {})
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_ppc_existing_axes():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _, ax_in = plt.subplots(subplot_kw={"projection": "polar"})
    sites = [_site("S00")]
    ax_out = plot_polar_coverage(sites, 1.0, 1e15, ax=ax_in)
    assert ax_out is ax_in
    plt.close("all")


# ─────────────────────────────────────────────────────────────────────────────
# plot_polar_errors
# ─────────────────────────────────────────────────────────────────────────────


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_ppe_returns_axes():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fr = _freqs(10)
    r = _rho_obs_from_z(fr)
    sites = [_FakeSite("S00", _make_z(fr), fr)]
    ax = plot_polar_errors(sites, {"S00": r * 1.2})
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_ppe_empty_input():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ax = plot_polar_errors([], {})
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_ppe_existing_axes():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _, ax_in = plt.subplots(subplot_kw={"projection": "polar"})
    fr = _freqs(10)
    r = _rho_obs_from_z(fr)
    sites = [_FakeSite("S00", _make_z(fr), fr)]
    ax_out = plot_polar_errors(sites, {"S00": r}, ax=ax_in)
    assert ax_out is ax_in
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_ppe_multiple_sites():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fr = _freqs(10)
    r = _rho_obs_from_z(fr)
    sites = [_FakeSite(f"S{i:02d}", _make_z(fr), fr) for i in range(4)]
    pred = {f"S{i:02d}": r * (1.0 + 0.1 * i) for i in range(4)}
    ax = plot_polar_errors(sites, pred)
    assert ax is not None
    plt.close("all")


# ─────────────────────────────────────────────────────────────────────────────
# plot_width_drift
# ─────────────────────────────────────────────────────────────────────────────


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_pwd_cartesian():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sites = [_site(f"S{i:02d}") for i in range(3)]
    ax = plot_width_drift(sites, 1.0, 1e15, polar=False)
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_pwd_polar():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sites = [_site(f"S{i:02d}") for i in range(3)]
    ax = plot_width_drift(sites, 1.0, 1e15, polar=True)
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_pwd_empty_input():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ax = plot_width_drift([], {}, {})
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_pwd_existing_axes():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _, ax_in = plt.subplots()
    sites = [_site("S00")]
    ax_out = plot_width_drift(sites, 1.0, 1e15, ax=ax_in)
    assert ax_out is ax_in
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_pwd_yx_component():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sites = [_site("S00")]
    ax = plot_width_drift(sites, 1.0, 1e15, rho_comp="yx")
    assert ax is not None
    plt.close("all")
