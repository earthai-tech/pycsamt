"""
Tests for pycsamt.emtools.fieldzone
"""
from __future__ import annotations

import numpy as np
import pytest

from pycsamt.emtools.fieldzone import (
    classify_field_zones,
    near_field_factor,
    plot_field_zones,
)

# ----------------------------- fixtures ----------------------------------- #

class _FakeZ:
    """Minimal Z-like object understood by _get_z_block."""
    def __init__(self, z, freq):
        self.z    = np.asarray(z, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)


class _FakeSite:
    """Minimal EDI-like object (passes is_edi_file duck-type check)."""
    def __init__(self, station, z, freq, offset=None):
        self.station = station
        self.Z = _FakeZ(z, freq)
        self.freq = np.asarray(freq, dtype=float)
        if offset is not None:
            self.source_offset = float(offset)

    def get_section(self, *_, **__):   # required by is_edi_file()
        return None


def _make_z(n_freq=6, rho=100.0):
    """
    Build a synthetic impedance tensor (n, 2, 2) with constant ρ_a = rho.

    Uses pycsamt's unit convention: ρ_a = 0.2 * |Z|² / f,
    so |Z| = sqrt(5 * f * ρ).
    """
    freqs = np.logspace(-2, 2, n_freq)          # 0.01 – 100 Hz
    z_abs = np.sqrt(5.0 * freqs * rho)           # gives ρ_a = rho via 0.2*|Z|²/f
    z_full = np.zeros((n_freq, 2, 2), dtype=complex)
    z_full[:, 0, 1] =  z_abs * (1.0 + 1j) / np.sqrt(2)   # Zxy
    z_full[:, 1, 0] = -z_abs * (1.0 + 1j) / np.sqrt(2)   # Zyx
    return z_full, freqs


def _one_site(station="S01", rho=100.0, n_freq=6, offset=None):
    z, fr = _make_z(n_freq=n_freq, rho=rho)
    return _FakeSite(station, z, fr, offset=offset)


# -------------------------------------------------------------------------- #
# classify_field_zones                                                        #
# -------------------------------------------------------------------------- #

def test_classify_returns_dataframe_columns():
    site = _one_site(offset=5000.0)
    df = classify_field_zones([site], source_offset=5000.0)
    expected = {
        "station", "freq_hz", "period_s", "offset_m",
        "rho_a_ohmm", "delta_bostick_m", "kr", "zone",
    }
    assert expected.issubset(df.columns)


def test_classify_row_count():
    n = 8
    site = _one_site(n_freq=n, offset=3000.0)
    df = classify_field_zones([site], source_offset=3000.0)
    assert len(df) == n


def test_classify_two_sites():
    s1 = _one_site("A", offset=1000.0)
    s2 = _one_site("B", offset=8000.0)
    df = classify_field_zones([s1, s2], source_offset={"A": 1000.0, "B": 8000.0})
    assert set(df["station"].unique()) == {"A", "B"}


def test_classify_zone_labels_valid():
    site = _one_site(offset=3000.0)
    df = classify_field_zones([site], source_offset=3000.0)
    assert set(df["zone"].unique()).issubset({"far", "transition", "near"})


def test_classify_far_field_large_offset():
    """Very large offset → all far-field zones."""
    site = _one_site(rho=10.0, n_freq=6, offset=500_000.0)
    df = classify_field_zones([site], source_offset=500_000.0)
    assert (df["zone"] == "far").all(), df[["freq_hz", "kr", "zone"]]


def test_classify_near_field_small_offset():
    """Very small offset → all near-field zones."""
    site = _one_site(rho=1000.0, n_freq=6, offset=1.0)
    df = classify_field_zones([site], source_offset=1.0)
    assert (df["zone"] == "near").all(), df[["freq_hz", "kr", "zone"]]


def test_classify_kr_monotone_with_offset():
    """Higher offset → higher |kr| at same frequency."""
    s1 = _one_site("near", offset=100.0)
    s2 = _one_site("far",  offset=50_000.0)
    df = classify_field_zones([s1, s2],
                              source_offset={"near": 100.0, "far": 50_000.0})
    kr_near = df[df.station == "near"]["kr"].values
    kr_far  = df[df.station == "far"]["kr"].values
    assert (kr_far > kr_near).all()


def test_classify_offset_from_site_attribute():
    """source_offset=None → reads ed.source_offset."""
    site = _one_site(offset=4000.0)   # stores as ed.source_offset
    df = classify_field_zones([site], source_offset=None)
    assert len(df) > 0
    assert (df["offset_m"] == 4000.0).all()


def test_classify_missing_offset_skips_site(recwarn):
    """Site with no offset and source_offset=None → skipped with warning."""
    site = _FakeSite("X", *_make_z())   # no offset set
    df = classify_field_zones([site], source_offset=None, verbose=1)
    assert df.empty


def test_classify_dict_offset_case_insensitive():
    site = _one_site("S01")
    df = classify_field_zones([site], source_offset={"s01": 2000.0})
    assert len(df) > 0


def test_classify_empty_input_returns_empty_df():
    df = classify_field_zones([], source_offset=1000.0)
    assert df.empty


def test_classify_custom_thresholds():
    site = _one_site(rho=100.0, offset=3000.0)
    df_default = classify_field_zones([site], source_offset=3000.0)
    df_strict  = classify_field_zones([site], source_offset=3000.0,
                                      far_threshold=100.0,
                                      near_threshold=50.0)
    # with very high far_threshold nothing should be "far"
    assert "far" not in df_strict["zone"].values


def test_classify_bostick_formula():
    """δ_B = 356 √(ρ/f); verify against hand-calc for one row."""
    rho = 100.0
    freq = 1.0
    offset = 2000.0
    site = _one_site(rho=rho, n_freq=1, offset=offset)
    # override freq to exactly 1 Hz
    site.Z.freq = np.array([freq])
    site.freq   = np.array([freq])
    # pycsamt unit: |Z| = sqrt(5 * f * rho)
    z_abs = np.sqrt(5.0 * freq * rho)
    site.Z.z = np.zeros((1, 2, 2), dtype=complex)
    site.Z.z[0, 0, 1] =  z_abs * (1 + 1j) / np.sqrt(2)
    site.Z.z[0, 1, 0] = -z_abs * (1 + 1j) / np.sqrt(2)
    df = classify_field_zones([site], source_offset=offset)
    expected_db = 356.0 * np.sqrt(rho / freq)
    expected_kr = offset / expected_db
    assert np.isclose(df["delta_bostick_m"].iloc[0], expected_db, rtol=1e-6)
    assert np.isclose(df["kr"].iloc[0], expected_kr, rtol=1e-6)


# -------------------------------------------------------------------------- #
# near_field_factor                                                           #
# -------------------------------------------------------------------------- #

def test_nf_factor_columns():
    site = _one_site(offset=5000.0)
    df = near_field_factor([site], source_offset=5000.0)
    expected = {"station", "freq_hz", "period_s", "offset_m",
                "rho_a_ohmm", "kr", "nf_factor"}
    assert expected.issubset(df.columns)


def test_nf_factor_far_field_near_unity():
    """For large |kr| the correction factor should approach 1."""
    site = _one_site(rho=10.0, n_freq=6, offset=500_000.0)
    df = near_field_factor([site], source_offset=500_000.0)
    # In the far field |F| → 1; allow generous tolerance because ρ_a is
    # itself estimated from the synthetic Z (not exactly rho)
    assert np.allclose(df["nf_factor"].values, 1.0, atol=0.15), \
        df[["kr", "nf_factor"]]


def test_nf_factor_near_field_deviates():
    """For small |kr| the correction factor should be != 1."""
    site = _one_site(rho=1000.0, n_freq=6, offset=1.0)
    df = near_field_factor([site], source_offset=1.0)
    assert not np.allclose(df["nf_factor"].values, 1.0, atol=0.05)


def test_nf_factor_empty_input():
    df = near_field_factor([], source_offset=1000.0)
    assert df.empty


def test_nf_factor_positive():
    site = _one_site(offset=3000.0)
    df = near_field_factor([site], source_offset=3000.0)
    assert (df["nf_factor"] > 0).all()


# -------------------------------------------------------------------------- #
# plot_field_zones                                                            #
# -------------------------------------------------------------------------- #

@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_returns_axes():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    site = _one_site(offset=3000.0)
    ax = plot_field_zones([site], source_offset=3000.0)
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_empty_returns_axes():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ax = plot_field_zones([], source_offset=1000.0)
    assert ax is not None
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_accepts_existing_axes():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax_in = plt.subplots()
    site = _one_site(offset=3000.0)
    ax_out = plot_field_zones([site], source_offset=3000.0, ax=ax_in)
    assert ax_out is ax_in
    plt.close("all")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_no_contour():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    site = _one_site(offset=3000.0)
    ax = plot_field_zones([site], source_offset=3000.0, contour_kr=False)
    assert ax is not None
    plt.close("all")
