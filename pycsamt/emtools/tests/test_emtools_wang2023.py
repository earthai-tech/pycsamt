"""
Tests for Wang & Lin (2023) methods in pycsamt.emtools.source_effects:
  - normalize_response
  - correct_near_field
  - plot_normalized_response
"""
from __future__ import annotations

import numpy as np
import pytest

from pycsamt.emtools.source_effects import (
    normalize_response,
    correct_near_field,
    plot_normalized_response,
)

# ─────────────────────────────────────────────────────────────────────────────
# Shared fixtures
# ─────────────────────────────────────────────────────────────────────────────


class _FakeZ:
    def __init__(self, z, freq):
        self.z    = np.asarray(z, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)


class _FakeSite:
    """Minimal EDI-like object accepted by ensure_sites."""
    def __init__(self, station, z, freq, offset=None):
        self.station = station
        self.Z    = _FakeZ(z, freq)
        self.freq = np.asarray(freq, dtype=float)
        if offset is not None:
            self.source_offset = float(offset)

    def get_section(self, *_, **__):
        return None


def _make_z(freqs: np.ndarray, rho: float) -> np.ndarray:
    """Synthetic Z with constant ρ_a = rho and phase = 45°.

    Convention: ρ_a = 0.2 |Z|² / f  →  |Z| = √(5 f ρ).
    """
    z_abs = np.sqrt(5.0 * freqs * rho)
    z = np.zeros((freqs.size, 2, 2), dtype=complex)
    z[:, 0, 1] =  z_abs * (1 + 1j) / np.sqrt(2)
    z[:, 1, 0] = -z_abs * (1 + 1j) / np.sqrt(2)
    return z


def _site(
    station: str = "S01",
    rho: float = 100.0,
    n: int = 8,
    f_lo: float = 0.1,
    f_hi: float = 1e3,
    offset: float = 5000.0,
) -> _FakeSite:
    fr = np.logspace(np.log10(f_lo), np.log10(f_hi), n)
    return _FakeSite(station, _make_z(fr, rho), fr, offset=offset)


def _site_no_off(
    station: str = "S00",
    rho: float = 100.0,
    n: int = 6,
) -> _FakeSite:
    fr = np.logspace(0, 3, n)
    return _FakeSite(station, _make_z(fr, rho), fr)


# ─────────────────────────────────────────────────────────────────────────────
# normalize_response — columns and shape
# ─────────────────────────────────────────────────────────────────────────────

def test_normalize_response_columns():
    sites = [_site()]
    df = normalize_response(sites)
    expected = {
        "station", "freq_hz", "period_s", "offset_m",
        "rho_a_ohmm", "rho_n", "phi_obs_deg", "phi_ref_deg",
        "phi_diff_deg", "zone", "kr",
    }
    assert expected.issubset(df.columns)


def test_normalize_response_row_count():
    n = 10
    sites = [_site(n=n)]
    df = normalize_response(sites)
    assert len(df) == n


def test_normalize_response_multi_site_rows():
    n = 6
    sites = [_site(f"S{i:02d}", n=n) for i in range(4)]
    df = normalize_response(sites)
    assert len(df) == 4 * n


def test_normalize_response_empty_input():
    df = normalize_response([])
    assert df.empty
    assert set(_normalize_response_required_cols()).issubset(df.columns)


def _normalize_response_required_cols():
    return ["station", "freq_hz", "rho_n", "phi_diff_deg", "zone", "kr"]


# ─────────────────────────────────────────────────────────────────────────────
# normalize_response — values
# ─────────────────────────────────────────────────────────────────────────────

def test_normalize_response_rho_n_positive():
    sites = [_site()]
    df = normalize_response(sites)
    assert (df["rho_n"] > 0).all()


def test_normalize_response_rho_n_at_ref_is_unity():
    """When rho_obs == rho_ref the normalised value should be 1."""
    rho_ref = 150.0
    sites = [_site(rho=rho_ref)]
    df = normalize_response(sites, rho_ref=rho_ref)
    assert np.allclose(df["rho_n"].values, 1.0, rtol=1e-6)


def test_normalize_response_rho_n_doubles_when_obs_is_2x_ref():
    rho_ref = 100.0
    sites = [_site(rho=2 * rho_ref)]
    df = normalize_response(sites, rho_ref=rho_ref)
    assert np.allclose(df["rho_n"].values, 2.0, rtol=1e-6)


def test_normalize_response_phi_diff_zero_for_ideal_z():
    """Zxy has phase 45°; using comp='xy' and phi_ref=45° → phi_diff ≈ 0."""
    sites = [_site()]
    df = normalize_response(sites, phi_ref_deg=45.0, comp="xy")
    assert np.allclose(df["phi_diff_deg"].values, 0.0, atol=1e-8)


def test_normalize_response_phi_ref_stored_in_df():
    phi_ref = 30.0
    sites = [_site()]
    df = normalize_response(sites, phi_ref_deg=phi_ref)
    assert np.allclose(df["phi_ref_deg"].values, phi_ref, atol=1e-8)


def test_normalize_response_phi_diff_correct():
    """phi_diff = phi_obs - phi_ref regardless of phi_ref value."""
    phi_ref = 30.0
    sites = [_site()]
    df = normalize_response(sites, phi_ref_deg=phi_ref)
    expected = df["phi_obs_deg"] - phi_ref
    assert np.allclose(df["phi_diff_deg"].values, expected.values, atol=1e-8)


# ─────────────────────────────────────────────────────────────────────────────
# normalize_response — zone and kr
# ─────────────────────────────────────────────────────────────────────────────

def test_normalize_response_zone_labels_valid():
    sites = [_site(offset=3000.0)]
    df = normalize_response(sites)
    valid_zones = {"far", "transition", "near"}
    assert set(df["zone"].dropna()).issubset(valid_zones)


def test_normalize_response_no_offset_zone_is_none():
    sites = [_site_no_off()]
    df = normalize_response(sites)
    assert df["zone"].isna().all()
    assert df["kr"].isna().all()


def test_normalize_response_far_field_large_offset():
    """Very large offset → all far-field zones for low ρ."""
    sites = [_site(rho=10.0, offset=500_000.0)]
    df = normalize_response(sites)
    assert (df["zone"] == "far").all()


def test_normalize_response_near_field_small_offset():
    """Very small offset → all near-field zones for high ρ."""
    sites = [_site(rho=1000.0, offset=1.0)]
    df = normalize_response(sites)
    assert (df["zone"] == "near").all()


def test_normalize_response_kr_positive_when_offset_given():
    sites = [_site(offset=3000.0)]
    df = normalize_response(sites)
    assert (df["kr"] > 0).all()


# ─────────────────────────────────────────────────────────────────────────────
# normalize_response — comp parameter
# ─────────────────────────────────────────────────────────────────────────────

def test_normalize_response_comp_xy():
    sites = [_site()]
    df = normalize_response(sites, comp="xy")
    assert len(df) > 0
    assert (df["rho_n"] > 0).all()


def test_normalize_response_comp_yx():
    sites = [_site()]
    df = normalize_response(sites, comp="yx")
    assert len(df) > 0


def test_normalize_response_invalid_comp():
    sites = [_site()]
    with pytest.raises(ValueError):
        normalize_response(sites, comp="zz")


def test_normalize_response_dict_offset():
    s1 = _site_no_off("A")
    s2 = _site_no_off("B")
    df = normalize_response([s1, s2], source_offset={"A": 2000.0, "B": 8000.0})
    assert not df[df.station == "A"]["kr"].isna().any()
    assert not df[df.station == "B"]["kr"].isna().any()


# ─────────────────────────────────────────────────────────────────────────────
# correct_near_field — return type and shape
# ─────────────────────────────────────────────────────────────────────────────

def test_correct_near_field_returns_sites():
    from pycsamt.site.base import Sites
    sites = [_site("S01"), _site("S02")]
    result = correct_near_field(sites, source_offset=5000.0)
    assert isinstance(result, Sites)


def test_correct_near_field_station_count():
    n = 4
    sites = [_site(f"S{i:02d}") for i in range(n)]
    result = correct_near_field(sites, source_offset=5000.0)
    assert len(result) == n


def test_correct_near_field_empty_input():
    from pycsamt.site.base import Sites
    result = correct_near_field([], source_offset=1000.0)
    assert isinstance(result, Sites)
    assert len(result) == 0


# ─────────────────────────────────────────────────────────────────────────────
# correct_near_field — correction values
# ─────────────────────────────────────────────────────────────────────────────

def _z_from_sites(result):
    """Extract all Z tensors from a Sites result as a list of arrays."""
    tensors = []
    for ed in result:
        try:
            edi = getattr(ed, "edi", None) or ed
            Z = getattr(edi, "Z", None) or getattr(ed, "Z", None)
            if Z is not None:
                tensors.append(np.asarray(Z.z, dtype=complex).copy())
        except Exception:
            pass
    return tensors


def test_correct_near_field_far_field_near_identity():
    """In the far field (large offset) F ≈ 1 → correction ≈ identity."""
    offset = 500_000.0
    s = _site("S01", rho=100.0, n=8, f_lo=0.1, f_hi=1e3, offset=offset)
    z_before = s.Z.z.copy()

    result = correct_near_field([s], source_offset=offset)
    tensors = _z_from_sites(result)
    assert len(tensors) == 1
    z_after = tensors[0]
    # Far-field correction leaves Z unchanged to within 2 %
    assert np.allclose(z_after, z_before, rtol=0.02), (
        f"max rel diff = {np.max(np.abs(z_after - z_before) / (np.abs(z_before) + 1e-30)):.4f}"
    )


def test_correct_near_field_near_field_modifies_z():
    """In the near field (tiny offset) F ≠ 1 → Z is changed."""
    offset = 1.0
    s = _site("S01", rho=1000.0, n=6, f_lo=0.1, f_hi=1.0, offset=offset)
    z_before = s.Z.z.copy()

    result = correct_near_field([s], source_offset=offset)
    tensors = _z_from_sites(result)
    assert len(tensors) == 1
    z_after = tensors[0]
    assert not np.allclose(z_after, z_before, rtol=0.01), (
        "Near-field correction should change Z"
    )


def test_correct_near_field_inplace_modifies_original():
    """inplace=True should modify Z in-place on the original EDI objects."""
    offset = 1.0
    s = _site("S01", rho=1000.0, n=4, offset=offset)
    # Grab a reference to the underlying Z wrapper
    z_id_before = id(s.Z.z)
    correct_near_field([s], source_offset=offset, inplace=True)
    z_id_after = id(s.Z.z)
    # The object may be replaced (Z.z = new_array), but values should differ
    # from what they would be without correction (which we can check via
    # a fresh uncorrected site with the same rho)
    s2 = _site("S01", rho=1000.0, n=4, offset=offset)
    assert not np.allclose(s.Z.z, s2.Z.z, rtol=0.01), (
        "inplace=True should have modified Z.z"
    )


def test_correct_near_field_dict_offset():
    """Source offset supplied as {station: r} dict."""
    s1 = _site_no_off("A")
    s2 = _site_no_off("B")
    result = correct_near_field(
        [s1, s2],
        source_offset={"A": 1.0, "B": 500_000.0},
    )
    assert len(result) == 2


def test_correct_near_field_z_finite():
    """Corrected Z must not contain NaN or Inf."""
    sites = [_site(f"S{i:02d}", offset=3000.0) for i in range(3)]
    result = correct_near_field(sites, source_offset=3000.0)
    tensors = _z_from_sites(result)
    for z in tensors:
        assert np.isfinite(z).all(), "Z_corrected contains non-finite values"


# ─────────────────────────────────────────────────────────────────────────────
# plot_normalized_response
# ─────────────────────────────────────────────────────────────────────────────

@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_normalized_response_returns_tuple():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sites = [_site(f"S{i:02d}") for i in range(3)]
    result = plot_normalized_response(sites, rho_ref=100.0)
    assert isinstance(result, tuple)
    assert len(result) == 2
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_normalized_response_axes_not_none():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sites = [_site(f"S{i:02d}") for i in range(3)]
    ax1, ax2 = plot_normalized_response(sites, rho_ref=100.0)
    assert ax1 is not None
    assert ax2 is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_normalized_response_empty_input():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ax1, ax2 = plot_normalized_response([])
    assert ax1 is not None
    assert ax2 is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_normalized_response_existing_axes():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _, (a1, a2) = plt.subplots(1, 2)
    sites = [_site(f"S{i:02d}") for i in range(3)]
    r1, r2 = plot_normalized_response(sites, axes=(a1, a2))
    assert r1 is a1
    assert r2 is a2
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_normalized_response_freq_axis():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sites = [_site(f"S{i:02d}") for i in range(2)]
    ax1, ax2 = plot_normalized_response(sites, period_axis=False)
    assert ax1 is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_normalized_response_custom_cmaps():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sites = [_site(f"S{i:02d}") for i in range(2)]
    ax1, ax2 = plot_normalized_response(
        sites, cmap_rho="viridis", cmap_phi="plasma"
    )
    assert ax1 is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_normalized_response_custom_lims():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sites = [_site(f"S{i:02d}") for i in range(2)]
    ax1, ax2 = plot_normalized_response(
        sites, rho_n_lim=(0.5, 2.0), phi_diff_lim=(-30.0, 30.0)
    )
    assert ax1 is not None
    plt.close("all")
