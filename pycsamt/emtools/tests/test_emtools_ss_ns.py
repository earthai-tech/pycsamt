"""
Tests for the near-surface effect detection extension in
pycsamt.emtools.ss (detect_near_surface, plot_ns_detection).
"""
from __future__ import annotations

import numpy as np
import pytest

from pycsamt.api import APIFrame, reset_api_view
from pycsamt.emtools.ss import (
    detect_near_surface,
    estimate_ss_ama,
    estimate_ss_bilateral,
    estimate_ss_loess,
    estimate_ss_refmedian,
    plot_ns_detection,
)

# ----------------------------- fixtures ----------------------------------- #

class _FakeZ:
    def __init__(self, z, freq):
        self.z    = np.asarray(z, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)


class _FakeSite:
    """Minimal EDI-like object (passes is_edi_file duck-type check)."""
    def __init__(self, station, z, freq):
        self.station = station
        self.Z    = _FakeZ(z, freq)
        self.freq = np.asarray(freq, dtype=float)

    def get_section(self, *_, **__):
        return None


# pycsamt unit convention: ρ_a = 0.2 * |Z|² / f
# so |Z| = sqrt(5 * f * rho_a)
def _make_z(freqs: np.ndarray, rho: float) -> np.ndarray:
    z_abs = np.sqrt(5.0 * freqs * rho)
    z = np.zeros((freqs.size, 2, 2), dtype=complex)
    z[:, 0, 1] =  z_abs * (1 + 1j) / np.sqrt(2)
    z[:, 1, 0] = -z_abs * (1 + 1j) / np.sqrt(2)
    return z


def _clean_site(station: str, rho: float = 100.0, n: int = 12) -> _FakeSite:
    """Flat ρ_a curve — no distortion."""
    fr = np.logspace(-2, 3, n)
    return _FakeSite(station, _make_z(fr, rho), fr)


def _static_site(station: str, rho: float = 100.0,
                 shift: float = 0.4, n: int = 12) -> _FakeSite:
    """Flat ρ_a curve shifted uniformly → static effect."""
    fr = np.logspace(-2, 3, n)
    rho_shifted = rho * (10.0 ** shift)
    return _FakeSite(station, _make_z(fr, rho_shifted), fr)


def _ns_site(station: str, rho: float = 100.0,
             hf_noise: float = 0.6, f_split: float = 1.0,
             n: int = 12) -> _FakeSite:
    """
    ρ_a curve with large high-frequency scatter → near-surface effect.
    HF data are randomly perturbed by hf_noise (in log10 units).
    """
    fr = np.logspace(-2, 3, n)
    rng = np.random.default_rng(42)
    log_rho = np.full(n, np.log10(rho))
    hf = fr >= f_split
    log_rho[hf] += rng.normal(0, hf_noise, size=hf.sum())
    rho_arr = 10.0 ** log_rho
    z = np.zeros((n, 2, 2), dtype=complex)
    for i in range(n):
        z_abs = np.sqrt(5.0 * fr[i] * rho_arr[i])
        z[i, 0, 1] =  z_abs * (1 + 1j) / np.sqrt(2)
        z[i, 1, 0] = -z_abs * (1 + 1j) / np.sqrt(2)
    return _FakeSite(station, z, fr)


def _profile(sites, offsets=None):
    """Build a multi-site profile list."""
    return list(sites)


# -------------------------------------------------------------------------- #
# detect_near_surface — columns and shape                                     #
# -------------------------------------------------------------------------- #

def test_detect_ns_columns():
    reset_api_view()
    sites = [_clean_site(f"S{i:02d}") for i in range(5)]
    df = detect_near_surface(sites)
    expected = {
        "station", "n_hf", "n_lf",
        "sigma_hf", "sigma_lf", "ns_index",
        "slope_hf", "slope_lf", "gradient_delta",
        "ss_delta_log10", "ns_flag", "ss_flag", "distortion_type",
    }
    assert expected.issubset(df.columns)


def test_detect_ns_api_flag_overrides_global():
    reset_api_view()
    sites = [_clean_site(f"S{i:02d}") for i in range(5)]

    plain = detect_near_surface(sites, api=False)
    viewed = detect_near_surface(sites, api=True)

    assert not isinstance(plain, APIFrame)
    assert isinstance(viewed, APIFrame)
    assert viewed.kind == "emtools.ss.near_surface"
    assert viewed.df.equals(plain)


def test_static_shift_estimators_can_return_api_frames_for_empty_input():
    reset_api_view()

    for func, kind in (
        (estimate_ss_ama, "emtools.ss.ama"),
        (estimate_ss_loess, "emtools.ss.loess"),
        (estimate_ss_bilateral, "emtools.ss.bilateral"),
        (estimate_ss_refmedian, "emtools.ss.refmedian"),
    ):
        viewed = func([], api=True)

        assert isinstance(viewed, APIFrame)
        assert viewed.kind == kind
        assert viewed.df.empty


def test_estimate_ss_ama_single_station_returns_empty_table():
    """A single-station survey yields no AMA neighbours.

    Regression: ``rows`` ends up empty, so ``pd.DataFrame.from_records``
    produced a column-less frame and ``sort_values("station")`` raised
    ``KeyError: 'station'``.  It must instead return an empty,
    correctly-typed table so correction degrades to a no-op.
    """
    tbl = estimate_ss_ama([_clean_site("S1")])
    assert list(tbl.columns) == [
        "station", "delta_log10_rho", "fac_rho", "fac_z", "n_used",
    ]
    assert len(tbl) == 0


def test_detect_ns_row_count():
    n = 6
    sites = [_clean_site(f"S{i:02d}") for i in range(n)]
    df = detect_near_surface(sites)
    assert len(df) == n


def test_detect_ns_empty_input():
    df = detect_near_surface([])
    assert df.empty


def test_detect_ns_single_site():
    """Single site: no AMA neighbors → residuals are NaN → still produces row."""
    df = detect_near_surface([_clean_site("S00")])
    assert len(df) == 1


# -------------------------------------------------------------------------- #
# Zone labels                                                                 #
# -------------------------------------------------------------------------- #

def test_detect_ns_clean_profile():
    """Flat ρ_a for all sites → all 'clean' or 'static' (near-zero shifts)."""
    sites = [_clean_site(f"S{i:02d}") for i in range(7)]
    df = detect_near_surface(sites, ns_threshold=2.0, ss_threshold=0.1)
    types = set(df["distortion_type"].unique())
    assert types.issubset({"clean", "static"}), types


def test_detect_ns_flags_bool_dtype():
    sites = [_clean_site(f"S{i:02d}") for i in range(4)]
    df = detect_near_surface(sites)
    assert df["ns_flag"].dtype == bool
    assert df["ss_flag"].dtype == bool


def test_detect_ns_type_valid_values():
    sites = [_clean_site(f"S{i:02d}") for i in range(5)]
    df = detect_near_surface(sites)
    valid = {"clean", "static", "near_surface", "mixed"}
    assert set(df["distortion_type"].unique()).issubset(valid)


def test_detect_ns_near_surface_detected():
    """
    Profile where one site has large HF scatter → should be flagged as
    'near_surface' or 'mixed'.
    """
    sites = (
        [_clean_site(f"S{i:02d}") for i in range(4)]
        + [_ns_site("NS", hf_noise=1.2, f_split=1.0)]
        + [_clean_site(f"S{i:02d}") for i in range(5, 8)]
    )
    df = detect_near_surface(
        sites,
        f_split=1.0,
        ns_threshold=1.5,
        ss_threshold=0.5,
    )
    ns_row = df[df["station"] == "NS"]
    assert len(ns_row) == 1
    assert ns_row["ns_flag"].iloc[0] or ns_row["ns_index"].iloc[0] > 1.0, \
        ns_row[["ns_index", "distortion_type"]]


def test_detect_ns_static_detected():
    """
    One site uniformly shifted → 'static' or 'mixed'.
    """
    sites = (
        [_clean_site(f"S{i:02d}") for i in range(3)]
        + [_static_site("SS", shift=0.5)]
        + [_clean_site(f"S{i:02d}") for i in range(4, 7)]
    )
    df = detect_near_surface(
        sites,
        ns_threshold=2.0,
        ss_threshold=0.15,
    )
    ss_row = df[df["station"] == "SS"]
    assert len(ss_row) == 1
    # The shifted site should have a non-zero ss_delta
    assert abs(ss_row["ss_delta_log10"].iloc[0]) > 0.05, \
        ss_row[["ss_delta_log10", "distortion_type"]]


# -------------------------------------------------------------------------- #
# NS index properties                                                          #
# -------------------------------------------------------------------------- #

def test_ns_index_positive():
    sites = [_clean_site(f"S{i:02d}") for i in range(5)]
    df = detect_near_surface(sites)
    valid = df["ns_index"].dropna()
    assert (valid >= 0).all()


def test_ns_index_larger_for_ns_site():
    """NS site should have higher η than a clean site."""
    clean = [_clean_site(f"C{i:02d}") for i in range(4)]
    ns    = [_ns_site("NS", hf_noise=1.5)]
    rest  = [_clean_site(f"R{i:02d}") for i in range(5, 8)]
    df = detect_near_surface(clean + ns + rest, f_split=1.0)
    row_clean = df[df["station"].str.startswith("C")]["ns_index"].median()
    row_ns    = df[df["station"] == "NS"]["ns_index"].iloc[0]
    if np.isfinite(row_ns) and np.isfinite(row_clean):
        assert row_ns >= row_clean


# -------------------------------------------------------------------------- #
# Gradient delta                                                               #
# -------------------------------------------------------------------------- #

def test_gradient_delta_non_negative():
    sites = [_clean_site(f"S{i:02d}") for i in range(5)]
    df = detect_near_surface(sites)
    valid = df["gradient_delta"].dropna()
    assert (valid >= 0).all()


def test_gradient_delta_present():
    sites = [_clean_site(f"S{i:02d}") for i in range(5)]
    df = detect_near_surface(sites)
    assert "gradient_delta" in df.columns
    assert df["gradient_delta"].notna().any()


# -------------------------------------------------------------------------- #
# n_hf / n_lf counts                                                          #
# -------------------------------------------------------------------------- #

def test_hf_lf_counts_sum_to_total():
    sites = [_clean_site(f"S{i:02d}", n=10) for i in range(4)]
    df = detect_near_surface(sites, f_split=1.0)
    assert ((df["n_hf"] + df["n_lf"]) <= 10).all()


def test_f_split_affects_counts():
    sites = [_clean_site(f"S{i:02d}") for i in range(4)]
    df_lo = detect_near_surface(sites, f_split=0.01)   # nearly all HF
    df_hi = detect_near_surface(sites, f_split=1e4)    # nearly all LF
    assert df_lo["n_hf"].sum() >= df_hi["n_hf"].sum()


# -------------------------------------------------------------------------- #
# plot_ns_detection                                                            #
# -------------------------------------------------------------------------- #

@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ns_returns_axes():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sites = [_clean_site(f"S{i:02d}") for i in range(5)]
    ax = plot_ns_detection(sites)
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ns_empty():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ax = plot_ns_detection([])
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ns_accepts_existing_axes():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax_in = plt.subplots()
    sites = [_clean_site(f"S{i:02d}") for i in range(4)]
    ax_out = plot_ns_detection(sites, ax=ax_in, show_ss=False)
    assert ax_out is ax_in
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_ns_no_ss_overlay():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sites = [_clean_site(f"S{i:02d}") for i in range(4)]
    ax = plot_ns_detection(sites, show_ss=False)
    assert ax is not None
    plt.close("all")
