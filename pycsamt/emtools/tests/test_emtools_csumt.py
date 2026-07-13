"""
Tests for pycsamt.emtools.csumt
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.emtools.csumt import (
    BOSTICK_CONST,
    F_MAX_CSUMT,
    F_MIN_CSUMT,
    bostick_depth,
    bostick_depth_from_rho,
    depth_coverage_table,
    frequency_for_depth,
    frequency_schedule,
    plot_depth_section,
    vertical_resolution,
    vertical_resolution_pair,
)

# ----------------------------- fixtures ----------------------------------- #


class _FakeZ:
    def __init__(self, z, freq):
        self.z = np.asarray(z, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)


class _FakeSite:
    """Minimal EDI-like object (passes is_edi_file duck-type check)."""

    def __init__(self, station, z, freq):
        self.station = station
        self.Z = _FakeZ(z, freq)
        self.freq = np.asarray(freq, dtype=float)

    def get_section(self, *_, **__):
        return None


def _make_z(freqs: np.ndarray, rho: float) -> np.ndarray:
    """Build Z with constant ρ_a = rho (pycsamt unit: |Z| = sqrt(5 f ρ))."""
    z_abs = np.sqrt(5.0 * freqs * rho)
    z = np.zeros((freqs.size, 2, 2), dtype=complex)
    z[:, 0, 1] = z_abs * (1 + 1j) / np.sqrt(2)
    z[:, 1, 0] = -z_abs * (1 + 1j) / np.sqrt(2)
    return z


def _site(
    station: str,
    rho: float = 100.0,
    f_lo: float = 1e1,
    f_hi: float = 1e4,
    n: int = 10,
) -> _FakeSite:
    fr = np.logspace(np.log10(f_lo), np.log10(f_hi), n)
    return _FakeSite(station, _make_z(fr, rho), fr)


# -------------------------------------------------------------------------- #
# bostick_depth_from_rho  (pure math)                                         #
# -------------------------------------------------------------------------- #


def test_bostick_depth_formula():
    rho, f = 100.0, 1.0
    expected = BOSTICK_CONST * np.sqrt(rho / f)
    assert np.isclose(bostick_depth_from_rho(rho, f), expected)


def test_bostick_depth_scalar():
    d = bostick_depth_from_rho(100.0, 100.0)
    assert d > 0


def test_bostick_depth_array():
    freqs = np.logspace(1, 4, 8)
    d = bostick_depth_from_rho(100.0, freqs)
    assert d.shape == freqs.shape
    assert (d > 0).all()
    assert (np.diff(d) < 0).all(), "depth should decrease as frequency rises"


def test_bostick_depth_higher_rho_deeper():
    d_lo = bostick_depth_from_rho(10.0, 100.0)
    d_hi = bostick_depth_from_rho(1000.0, 100.0)
    assert d_hi > d_lo


def test_bostick_depth_zero_freq_no_crash():
    d = bostick_depth_from_rho(100.0, 0.0)
    assert np.isfinite(d) or np.isinf(d)  # just no exception


# -------------------------------------------------------------------------- #
# vertical_resolution_pair  (pure math)                                       #
# -------------------------------------------------------------------------- #


def test_vr_pair_positive():
    dr = vertical_resolution_pair(100.0, 100.0, 1000.0)
    assert dr > 0


def test_vr_pair_formula():
    rho, f_lo, f_hi = 100.0, 100.0, 1000.0
    expected = (
        BOSTICK_CONST
        * np.sqrt(rho)
        * (1.0 / np.sqrt(f_lo) - 1.0 / np.sqrt(f_hi))
    )
    assert np.isclose(vertical_resolution_pair(rho, f_lo, f_hi), expected)


def test_vr_pair_equal_frequencies():
    assert vertical_resolution_pair(100.0, 500.0, 500.0) == pytest.approx(0.0)


def test_vr_pair_larger_spacing_higher_delta():
    dr_wide = vertical_resolution_pair(100.0, 10.0, 10_000.0)
    dr_narrow = vertical_resolution_pair(100.0, 100.0, 1_000.0)
    assert dr_wide > dr_narrow


def test_vr_pair_higher_rho_higher_delta():
    dr_lo = vertical_resolution_pair(10.0, 100.0, 1000.0)
    dr_hi = vertical_resolution_pair(1000.0, 100.0, 1000.0)
    assert dr_hi > dr_lo


# -------------------------------------------------------------------------- #
# frequency_for_depth  (pure math)                                            #
# -------------------------------------------------------------------------- #


def test_freq_for_depth_formula():
    d, rho = 100.0, 100.0
    f = frequency_for_depth(d, rho)
    # invert: bostick_depth_from_rho(rho, f) should equal d
    assert np.isclose(bostick_depth_from_rho(rho, f), d, rtol=1e-6)


def test_freq_for_depth_deeper_lower_freq():
    f_shallow = frequency_for_depth(10.0, 100.0)
    f_deep = frequency_for_depth(100.0, 100.0)
    assert f_shallow > f_deep


def test_freq_for_depth_array():
    depths = np.array([5.0, 10.0, 50.0, 100.0])
    freqs = frequency_for_depth(depths, 100.0)
    assert freqs.shape == depths.shape
    assert (np.diff(freqs) < 0).all(), "deeper → lower frequency"


def test_bostick_roundtrip():
    """frequency_for_depth ∘ bostick_depth_from_rho ≈ identity."""
    rho = 50.0
    f_orig = np.logspace(1, 5, 10)
    d = bostick_depth_from_rho(rho, f_orig)
    f_back = frequency_for_depth(d, rho)
    assert np.allclose(f_orig, f_back, rtol=1e-6)


# -------------------------------------------------------------------------- #
# frequency_schedule  (pure math)                                             #
# -------------------------------------------------------------------------- #


def test_schedule_returns_array():
    depths = np.array([5.0, 10.0, 20.0, 50.0])
    freqs = frequency_schedule(depths, 100.0)
    assert isinstance(freqs, np.ndarray)


def test_schedule_within_bounds():
    depths = np.linspace(5, 30, 8)
    freqs = frequency_schedule(
        depths, 100.0, f_min=F_MIN_CSUMT, f_max=F_MAX_CSUMT
    )
    assert (freqs >= F_MIN_CSUMT).all()
    assert (freqs <= F_MAX_CSUMT).all()


def test_schedule_sorted():
    freqs = frequency_schedule([5.0, 10.0, 20.0], 100.0)
    assert (np.diff(freqs) >= 0).all()


def test_schedule_empty_out_of_range():
    """Target depths that map to frequencies outside [f_min, f_max] → empty."""
    very_deep = np.array([10_000.0])  # f ≪ 9.6 kHz for any reasonable rho
    freqs = frequency_schedule(
        very_deep, 100.0, f_min=F_MIN_CSUMT, f_max=F_MAX_CSUMT
    )
    assert freqs.size == 0


def test_schedule_as_khz():
    depths = [5.0, 10.0, 20.0]
    freqs_hz = frequency_schedule(depths, 100.0, as_khz=False)
    freqs_khz = frequency_schedule(depths, 100.0, as_khz=True)
    assert np.allclose(freqs_hz / 1e3, freqs_khz, rtol=1e-6)


def test_schedule_min_resolution_inserts_freqs():
    """Requesting fine min_resolution should insert extra frequencies."""
    depths = [5.0, 30.0]
    f_coarse = frequency_schedule(depths, 100.0)
    f_fine = frequency_schedule(depths, 100.0, min_resolution_m=1.0)
    assert f_fine.size >= f_coarse.size


def test_schedule_fill_decades():
    depths = [5.0, 30.0]
    f_plain = frequency_schedule(depths, 100.0, fill_decades=False)
    f_fill = frequency_schedule(
        depths, 100.0, fill_decades=True, per_decade=4
    )
    assert f_fill.size >= f_plain.size


def test_schedule_single_depth():
    freqs = frequency_schedule(10.0, 100.0)
    assert freqs.size <= 1  # may be empty if outside range


# -------------------------------------------------------------------------- #
# bostick_depth  (sites-based)                                                #
# -------------------------------------------------------------------------- #


def test_bostick_depth_sites_columns():
    sites = [_site(f"S{i:02d}") for i in range(3)]
    df = bostick_depth(sites)
    expected = {"station", "freq_hz", "period_s", "rho_a_ohmm", "depth_m"}
    assert expected.issubset(df.columns)


def test_bostick_depth_sites_row_count():
    n_freq = 8
    sites = [_site("S00", n=n_freq)]
    df = bostick_depth(sites)
    assert len(df) == n_freq


def test_bostick_depth_sites_formula_consistency():
    """D from sites should match bostick_depth_from_rho(ρ_a, f)."""
    rho = 100.0
    sites = [_site("S00", rho=rho, n=6)]
    df = bostick_depth(sites)
    expected = BOSTICK_CONST * np.sqrt(df["rho_a_ohmm"] / df["freq_hz"])
    assert np.allclose(df["depth_m"].values, expected.values, rtol=1e-6)


def test_bostick_depth_sites_positive():
    sites = [_site(f"S{i:02d}") for i in range(3)]
    df = bostick_depth(sites)
    assert (df["depth_m"] > 0).all()


def test_bostick_depth_sites_empty():
    df = bostick_depth([])
    assert df.empty


def test_bostick_depth_multiple_sites():
    sites = [_site(f"S{i:02d}", rho=10.0 * (i + 1)) for i in range(4)]
    df = bostick_depth(sites)
    assert len(df["station"].unique()) == 4


# -------------------------------------------------------------------------- #
# vertical_resolution  (sites-based)                                          #
# -------------------------------------------------------------------------- #


def test_vr_sites_columns():
    sites = [_site("S00", n=8)]
    df = vertical_resolution(sites)
    expected = {
        "station",
        "freq_lo_hz",
        "freq_hi_hz",
        "depth_lo_m",
        "depth_hi_m",
        "delta_depth_m",
        "rho_a_ohmm",
    }
    assert expected.issubset(df.columns)


def test_vr_sites_row_count():
    n = 8
    sites = [_site("S00", n=n)]
    df = vertical_resolution(sites)
    assert len(df) == n - 1


def test_vr_sites_positive_delta():
    sites = [_site("S00", n=8)]
    df = vertical_resolution(sites)
    # lower freq → deeper; ΔD should be positive
    assert (df["delta_depth_m"] > 0).all(), df[
        ["freq_lo_hz", "freq_hi_hz", "delta_depth_m"]
    ]


def test_vr_sites_rho_override():
    rho = 200.0
    sites = [_site("S00", n=6)]
    df = vertical_resolution(sites, rho_override=rho)
    assert (df["rho_a_ohmm"] == rho).all()


def test_vr_sites_empty():
    df = vertical_resolution([])
    assert df.empty


# -------------------------------------------------------------------------- #
# depth_coverage_table  (sites-based)                                         #
# -------------------------------------------------------------------------- #


def test_coverage_columns():
    sites = [_site(f"S{i:02d}") for i in range(3)]
    df = depth_coverage_table(sites)
    expected = {
        "station",
        "n_freq",
        "freq_min_hz",
        "freq_max_hz",
        "depth_min_m",
        "depth_max_m",
        "mean_resolution_m",
        "median_resolution_m",
    }
    assert expected.issubset(df.columns)


def test_coverage_one_row_per_station():
    n_sites = 4
    sites = [_site(f"S{i:02d}") for i in range(n_sites)]
    df = depth_coverage_table(sites)
    assert len(df) == n_sites


def test_coverage_depth_ordering():
    """depth_min < depth_max for every station."""
    sites = [_site(f"S{i:02d}", n=8) for i in range(3)]
    df = depth_coverage_table(sites)
    assert (df["depth_min_m"] < df["depth_max_m"]).all()


def test_coverage_empty():
    df = depth_coverage_table([])
    assert df.empty


# -------------------------------------------------------------------------- #
# plot_depth_section  (sites-based)                                           #
# -------------------------------------------------------------------------- #


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_returns_axes():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sites = [_site(f"S{i:02d}") for i in range(3)]
    ax = plot_depth_section(sites)
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_empty_input():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ax = plot_depth_section([])
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_existing_axes():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _, ax_in = plt.subplots()
    sites = [_site(f"S{i:02d}") for i in range(3)]
    ax_out = plot_depth_section(sites, ax=ax_in)
    assert ax_out is ax_in
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_linear_color():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sites = [_site(f"S{i:02d}") for i in range(3)]
    ax = plot_depth_section(sites, log_color=False)
    assert ax is not None
    plt.close("all")
