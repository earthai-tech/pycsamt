"""
Tests for pycsamt.emtools.gradient_imaging
  - rho_spatial_gradient
  - rho_frequency_gradient
  - rho_joint_gradient
  - plot_gradient_section
"""
from __future__ import annotations

import numpy as np
import pytest

from pycsamt.emtools.gradient_imaging import (
    rho_spatial_gradient,
    rho_frequency_gradient,
    rho_joint_gradient,
    plot_gradient_section,
)

# ─────────────────────────────────────────────────────────────────────────────
# Fixtures
# ─────────────────────────────────────────────────────────────────────────────


class _FakeZ:
    def __init__(self, z, freq):
        self.z    = np.asarray(z, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)


class _FakeSite:
    """Minimal EDI-like object accepted by ensure_sites."""
    def __init__(self, station, z, freq, east=None, north=None):
        self.station = station
        self.Z    = _FakeZ(z, freq)
        self.freq = np.asarray(freq, dtype=float)
        if east is not None:
            self.east = float(east)
        if north is not None:
            self.north = float(north)

    def get_section(self, *_, **__):
        return None


def _make_z(freqs: np.ndarray, rho: float) -> np.ndarray:
    """Synthetic Z: constant ρ_a = rho, phase = 45°."""
    z_abs = np.sqrt(5.0 * freqs * rho)
    z = np.zeros((freqs.size, 2, 2), dtype=complex)
    z[:, 0, 1] =  z_abs * (1 + 1j) / np.sqrt(2)
    z[:, 1, 0] = -z_abs * (1 + 1j) / np.sqrt(2)
    return z


def _uniform_sites(
    n: int = 5,
    rho: float = 100.0,
    n_freq: int = 6,
    f_lo: float = 1.0,
    f_hi: float = 1e3,
    spacing_m: float = 200.0,
) -> list:
    """N collinear sites with identical ρ_a = rho on a uniform grid."""
    freqs = np.logspace(np.log10(f_lo), np.log10(f_hi), n_freq)
    sites = []
    for i in range(n):
        s = _FakeSite(
            f"S{i:02d}",
            _make_z(freqs, rho),
            freqs,
            east=float(i) * spacing_m,
            north=0.0,
        )
        sites.append(s)
    return sites


def _heterogeneous_sites(n: int = 5, n_freq: int = 6) -> list:
    """Sites with linearly varying ρ_a (100 to 500 Ω·m)."""
    freqs = np.logspace(0, 3, n_freq)
    rho_vals = np.linspace(100.0, 500.0, n)
    sites = []
    for i in range(n):
        s = _FakeSite(
            f"S{i:02d}",
            _make_z(freqs, rho_vals[i]),
            freqs,
            east=float(i) * 200.0,
            north=0.0,
        )
        sites.append(s)
    return sites


# ─────────────────────────────────────────────────────────────────────────────
# rho_spatial_gradient — DataFrame shape and columns
# ─────────────────────────────────────────────────────────────────────────────

def test_spatial_columns():
    sites = _uniform_sites()
    df = rho_spatial_gradient(sites)
    expected = {
        "station_a", "station_b", "x_m", "dx_m",
        "freq_hz", "period_s", "depth_m",
        "rho_a_ohmm", "delta_rho_x",
    }
    assert expected.issubset(df.columns)


def test_spatial_row_count():
    n, f = 5, 6
    sites = _uniform_sites(n=n, n_freq=f)
    df = rho_spatial_gradient(sites)
    assert len(df) == (n - 1) * f


def test_spatial_two_sites():
    sites = _uniform_sites(n=2, n_freq=4)
    df = rho_spatial_gradient(sites)
    assert len(df) == 4                 # 1 pair × 4 freqs


def test_spatial_one_site_empty():
    """Only one station → no pairs → empty DataFrame."""
    sites = _uniform_sites(n=1)
    df = rho_spatial_gradient(sites)
    assert df.empty


def test_spatial_empty_input():
    df = rho_spatial_gradient([])
    assert df.empty


def test_spatial_invalid_comp():
    sites = _uniform_sites()
    with pytest.raises(ValueError):
        rho_spatial_gradient(sites, comp="zz")


# ─────────────────────────────────────────────────────────────────────────────
# rho_spatial_gradient — values
# ─────────────────────────────────────────────────────────────────────────────

def test_spatial_uniform_rho_zero_gradient():
    """Uniform ρ_a → Δρ_a^x = 0 everywhere."""
    sites = _uniform_sites(rho=100.0)
    df = rho_spatial_gradient(sites)
    assert np.allclose(df["delta_rho_x"].values, 0.0, atol=1e-8)


def test_spatial_gradient_sign():
    """Increasing ρ_a along the line → positive Δρ_a^x."""
    sites = _heterogeneous_sites()
    df = rho_spatial_gradient(sites)
    # Each pair (j, j+1) has rho(j+1) > rho(j), so gradient > 0
    assert (df["delta_rho_x"] > 0).all()


def test_spatial_gradient_magnitude():
    """Δρ_a^x = ρ_a(j+1) - ρ_a(j) at any frequency."""
    n, n_freq = 3, 4
    rho_vals = [50.0, 200.0, 80.0]
    freqs = np.logspace(0, 2, n_freq)
    sites = [
        _FakeSite(f"S{i:02d}", _make_z(freqs, rho_vals[i]), freqs,
                  east=float(i) * 200.0, north=0.0)
        for i in range(n)
    ]
    df = rho_spatial_gradient(sites)
    # Pair (S00, S01) at any frequency
    pair01 = df[df.station_a == "S00"]
    expected = rho_vals[1] - rho_vals[0]  # ~ 150 Ω·m
    assert np.allclose(pair01["delta_rho_x"].values, expected, rtol=1e-5)


def test_spatial_rho_a_is_mean_of_pair():
    """rho_a_ohmm column is the mean of the two adjacent stations."""
    n, n_freq = 2, 3
    freqs = np.logspace(0, 2, n_freq)
    sites = [
        _FakeSite("A", _make_z(freqs, 100.0), freqs, east=0.0, north=0.0),
        _FakeSite("B", _make_z(freqs, 300.0), freqs, east=200.0, north=0.0),
    ]
    df = rho_spatial_gradient(sites)
    expected_mean = 200.0  # (100 + 300) / 2
    assert np.allclose(df["rho_a_ohmm"].values, expected_mean, rtol=1e-5)


def test_spatial_depth_positive():
    sites = _uniform_sites()
    df = rho_spatial_gradient(sites)
    assert (df["depth_m"] > 0).all()


def test_spatial_period_is_inverse_freq():
    sites = _uniform_sites(n=2, n_freq=4)
    df = rho_spatial_gradient(sites)
    assert np.allclose(df["period_s"].values, 1.0 / df["freq_hz"].values)


def test_spatial_comp_xy():
    sites = _uniform_sites()
    df = rho_spatial_gradient(sites, comp="xy")
    assert len(df) > 0


def test_spatial_comp_yx():
    sites = _uniform_sites()
    df = rho_spatial_gradient(sites, comp="yx")
    assert len(df) > 0


# ─────────────────────────────────────────────────────────────────────────────
# rho_frequency_gradient — DataFrame shape and columns
# ─────────────────────────────────────────────────────────────────────────────

def test_frequency_columns():
    sites = _uniform_sites()
    df = rho_frequency_gradient(sites)
    expected = {
        "station", "x_m",
        "freq_hz", "period_s", "depth_m",
        "rho_a_ohmm", "delta_rho_z",
    }
    assert expected.issubset(df.columns)


def test_frequency_row_count():
    n, f = 4, 7
    sites = _uniform_sites(n=n, n_freq=f)
    df = rho_frequency_gradient(sites)
    assert len(df) == n * (f - 1)


def test_frequency_one_freq_empty():
    """Only one frequency → no pairs → empty DataFrame."""
    freqs = np.array([1.0])
    s = _FakeSite("S00", _make_z(freqs, 100.0), freqs)
    df = rho_frequency_gradient([s])
    assert df.empty


def test_frequency_empty_input():
    df = rho_frequency_gradient([])
    assert df.empty


def test_frequency_invalid_comp():
    sites = _uniform_sites()
    with pytest.raises(ValueError):
        rho_frequency_gradient(sites, comp="xy2")


# ─────────────────────────────────────────────────────────────────────────────
# rho_frequency_gradient — values
# ─────────────────────────────────────────────────────────────────────────────

def test_frequency_uniform_rho_zero_gradient():
    """Constant ρ_a → Δρ_a^z = 0 everywhere."""
    sites = _uniform_sites(rho=100.0)
    df = rho_frequency_gradient(sites)
    assert np.allclose(df["delta_rho_z"].values, 0.0, atol=1e-8)


def test_frequency_gradient_sign_with_slope():
    """ρ_a increasing with frequency → positive Δρ_a^z."""
    n_freq = 6
    freqs  = np.logspace(0, 3, n_freq)
    # Make rho_a increase with frequency: rho_a = 0.2*|Z|²/f = rho → rho = f*constant
    # Use artificial Z where |Z|² ∝ f²  i.e. rho_a = 0.2*f
    rho_a_vals = 0.2 * freqs  # increases with f
    z_abs = np.sqrt(5.0 * freqs * rho_a_vals)
    z = np.zeros((n_freq, 2, 2), dtype=complex)
    z[:, 0, 1] =  z_abs * (1 + 1j) / np.sqrt(2)
    z[:, 1, 0] = -z_abs * (1 + 1j) / np.sqrt(2)
    s = _FakeSite("S00", z, freqs)
    df = rho_frequency_gradient([s])
    assert (df["delta_rho_z"] > 0).all()


def test_frequency_gradient_magnitude():
    """Δρ_a^z = ρ_a(f_k) - ρ_a(f_{k-1}) for each adjacent pair."""
    n_freq = 4
    freqs  = np.logspace(0, 2, n_freq)
    rho_base = 100.0
    s = _FakeSite("S00", _make_z(freqs, rho_base), freqs)
    df = rho_frequency_gradient([s])
    # For constant rho_a, all differences are 0
    assert np.allclose(df["delta_rho_z"].values, 0.0, atol=1e-8)


def test_frequency_depth_positive():
    sites = _uniform_sites()
    df = rho_frequency_gradient(sites)
    assert (df["depth_m"] > 0).all()


def test_frequency_station_count():
    n = 5
    sites = _uniform_sites(n=n, n_freq=4)
    df = rho_frequency_gradient(sites)
    assert df["station"].nunique() == n


# ─────────────────────────────────────────────────────────────────────────────
# rho_joint_gradient — DataFrame shape and columns
# ─────────────────────────────────────────────────────────────────────────────

def test_joint_columns():
    sites = _uniform_sites()
    df = rho_joint_gradient(sites)
    expected = {
        "station_a", "station_b", "x_m", "dx_m",
        "freq_hz", "period_s", "depth_m",
        "delta_rho_zx",
    }
    assert expected.issubset(df.columns)


def test_joint_row_count():
    n, f = 5, 6
    sites = _uniform_sites(n=n, n_freq=f)
    df = rho_joint_gradient(sites)
    assert len(df) == (n - 1) * (f - 1)


def test_joint_one_site_empty():
    sites = _uniform_sites(n=1, n_freq=4)
    df = rho_joint_gradient(sites)
    assert df.empty


def test_joint_one_freq_empty():
    freqs = np.array([1.0])
    sites = [_FakeSite(f"S{i:02d}", _make_z(freqs, 100.0), freqs,
                       east=float(i) * 200.0, north=0.0) for i in range(3)]
    df = rho_joint_gradient(sites)
    assert df.empty


def test_joint_empty_input():
    df = rho_joint_gradient([])
    assert df.empty


def test_joint_invalid_comp():
    sites = _uniform_sites()
    with pytest.raises(ValueError):
        rho_joint_gradient(sites, comp="bad")


# ─────────────────────────────────────────────────────────────────────────────
# rho_joint_gradient — values
# ─────────────────────────────────────────────────────────────────────────────

def test_joint_uniform_rho_zero():
    """Uniform ρ_a over all stations and frequencies → Δρ_a^{zx} = 0."""
    sites = _uniform_sites(rho=100.0, n=5, n_freq=6)
    df = rho_joint_gradient(sites)
    assert np.allclose(df["delta_rho_zx"].values, 0.0, atol=1e-8)


def test_joint_linear_rho_only_spatial_zero():
    """
    ρ_a varies laterally but NOT with frequency → Δρ_a^x = f(x) only.
    Then Δρ_a^{zx} = Δρ_a^x(f_k) - Δρ_a^x(f_{k-1}) = 0.
    """
    n, n_freq = 4, 5
    freqs = np.logspace(0, 2, n_freq)
    rho_vals = [50.0, 100.0, 200.0, 400.0]  # constant in freq, varies in space
    sites = [
        _FakeSite(f"S{i:02d}", _make_z(freqs, rho_vals[i]), freqs,
                  east=float(i) * 200.0, north=0.0)
        for i in range(n)
    ]
    df = rho_joint_gradient(sites)
    assert np.allclose(df["delta_rho_zx"].values, 0.0, atol=1e-6)


def test_joint_nonzero_when_both_vary():
    """
    ρ_a varies with BOTH frequency and position → Δρ_a^{zx} ≠ 0.
    Construct: ρ_a(j, f) = rho_j * freq_k (product of station and frequency factor).
    Then:
      Δρ_a^x(j, f_k) = (rho_{j+1} - rho_j) * f_k
      Δρ_a^{zx}(j, f_k) = (rho_{j+1} - rho_j) * (f_k - f_{k-1}) ≠ 0
    """
    n_freq = 4
    freqs  = np.linspace(1.0, 4.0, n_freq)
    rho_j  = [10.0, 30.0]   # rho values for 2 stations

    def _z_for_rho_f(rho_base, freqs):
        rho_a_arr = rho_base * freqs   # rho varies with frequency
        z_abs = np.sqrt(5.0 * freqs * rho_a_arr)
        z = np.zeros((len(freqs), 2, 2), dtype=complex)
        z[:, 0, 1] =  z_abs * (1 + 1j) / np.sqrt(2)
        z[:, 1, 0] = -z_abs * (1 + 1j) / np.sqrt(2)
        return z

    sites = [
        _FakeSite("A", _z_for_rho_f(rho_j[0], freqs), freqs, east=0.0,   north=0.0),
        _FakeSite("B", _z_for_rho_f(rho_j[1], freqs), freqs, east=200.0, north=0.0),
    ]
    df = rho_joint_gradient(sites)
    assert not np.allclose(df["delta_rho_zx"].values, 0.0, atol=1e-3)


def test_joint_antisymmetric_gradient():
    """
    Δρ_a^{zx}(j,f_k) = [ρ(j+1,f_k)-ρ(j,f_k)] - [ρ(j+1,f_{k-1})-ρ(j,f_{k-1})].
    Verify against a manual calculation.
    """
    freqs = np.array([1.0, 10.0])   # 2 frequencies
    rho = {
        "A": np.array([100.0, 200.0]),  # ρ_a at f=1, f=10
        "B": np.array([150.0, 400.0]),
    }

    def _z_arr(rho_arr, freqs):
        z = np.zeros((len(freqs), 2, 2), dtype=complex)
        for k, (f, r) in enumerate(zip(freqs, rho_arr)):
            z_abs = np.sqrt(5.0 * f * r)
            z[k, 0, 1] =  z_abs * (1 + 1j) / np.sqrt(2)
            z[k, 1, 0] = -z_abs * (1 + 1j) / np.sqrt(2)
        return z

    sites = [
        _FakeSite("A", _z_arr(rho["A"], freqs), freqs, east=0.0,   north=0.0),
        _FakeSite("B", _z_arr(rho["B"], freqs), freqs, east=200.0, north=0.0),
    ]
    df = rho_joint_gradient(sites)
    assert len(df) == 1   # (N-1=1 pairs) × (F-1=1 freq pairs)

    # Manual: Δρ_a^x(f=10) - Δρ_a^x(f=1)
    # Δρ_a^x(f=1)  = ρ(B,1)  - ρ(A,1)  = 150 - 100 = 50
    # Δρ_a^x(f=10) = ρ(B,10) - ρ(A,10) = 400 - 200 = 200
    # Δρ_a^{zx} = 200 - 50 = 150
    expected = 150.0
    assert np.isclose(df["delta_rho_zx"].iloc[0], expected, rtol=1e-5)


def test_joint_depth_positive():
    sites = _uniform_sites(n=3, n_freq=4)
    df = rho_joint_gradient(sites)
    if not df.empty:
        assert (df["depth_m"] > 0).all()


def test_joint_period_inverse_freq():
    sites = _uniform_sites(n=3, n_freq=4)
    df = rho_joint_gradient(sites)
    if not df.empty:
        assert np.allclose(df["period_s"].values, 1.0 / df["freq_hz"].values)


# ─────────────────────────────────────────────────────────────────────────────
# plot_gradient_section
# ─────────────────────────────────────────────────────────────────────────────

@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_joint_returns_axes():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sites = _uniform_sites(n=4, n_freq=5)
    ax = plot_gradient_section(sites, quantity="joint")
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_spatial_returns_axes():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sites = _uniform_sites(n=4, n_freq=5)
    ax = plot_gradient_section(sites, quantity="spatial")
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_frequency_returns_axes():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sites = _uniform_sites(n=4, n_freq=5)
    ax = plot_gradient_section(sites, quantity="frequency")
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_empty_input():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ax = plot_gradient_section([])
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_existing_axes():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _, ax_in = plt.subplots()
    sites = _heterogeneous_sites()
    ax_out = plot_gradient_section(sites, ax=ax_in)
    assert ax_out is ax_in
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_freq_axis():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sites = _uniform_sites(n=3, n_freq=4)
    ax = plot_gradient_section(sites, period_axis=False)
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_custom_vlim():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sites = _heterogeneous_sites()
    ax = plot_gradient_section(sites, vlim=(-50.0, 50.0))
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_invalid_quantity():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sites = _uniform_sites()
    with pytest.raises(ValueError):
        plot_gradient_section(sites, quantity="diagonal")
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_alias_x():
    """quantity='x' is an alias for 'spatial'."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sites = _uniform_sites(n=3, n_freq=4)
    ax = plot_gradient_section(sites, quantity="x")
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_alias_z():
    """quantity='z' is an alias for 'frequency'."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sites = _uniform_sites(n=3, n_freq=4)
    ax = plot_gradient_section(sites, quantity="z")
    assert ax is not None
    plt.close("all")
