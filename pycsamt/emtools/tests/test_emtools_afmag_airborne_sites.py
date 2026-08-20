"""Tests for pycsamt.emtools.afmag accepting AirborneSites.

AFMAG keeps two genuinely different response families
(pycsamt.airborne.afmag), neither tipper-shaped: a bare scalar tilt
(original comparator) and a complex (nf,3,2) interstation tensor
(AirMt). This module exercises the two dedicated functions built for
those real shapes -- airmt_tilt_angles/original_afmag_tilt_table --
plus the existing tipper-shaped functions' graceful degradation on
AirborneSites input that has no tipper transfer function at all.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.airborne import AirborneSites, NavigationTrack
from pycsamt.airborne.afmag import build_airmt_line, build_original_afmag_line
from pycsamt.emtools.afmag import (
    afmag_tilt_angles,
    airmt_tilt_angles,
    flag_motion_susceptible_band,
    original_afmag_conductor_diagnostics,
    original_afmag_tilt_table,
    plot_airmt_tilt_profile,
    plot_airmt_tilt_psection,
    plot_original_afmag_dual_frequency_profile,
    plot_original_afmag_tilt_profile,
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

_DATA_ROOT = Path(__file__).resolve().parents[3] / "data"
_ABITIBI = _DATA_ROOT / "AFMAG" / "abitibi_on"
_YULONG = _DATA_ROOT / "AFMAG" / "yulong_belt_cn"

_AIRMT_FREQS = np.array([25.0, 100.0, 500.0])
_ORIGINAL_FREQS = np.array([150.0, 510.0])


def _interstation_tensor(n_stations: int, x: np.ndarray) -> np.ndarray:
    nf = _AIRMT_FREQS.size
    data = np.zeros((n_stations, nf, 3, 2), dtype=complex)
    for i in range(n_stations):
        for k in range(nf):
            data[i, k, 0, 0] = 1.0
            data[i, k, 1, 1] = 1.0
            shape = (x[i] / 150.0) * np.exp(-((x[i] / 200.0) ** 2))
            data[i, k, 2, 0] = 0.15 * shape + 0.02j * shape
            data[i, k, 2, 1] = 0.05 * shape - 0.01j * shape
    return data


def _airmt_sites(n_stations: int = 9, spacing: float = 100.0) -> AirborneSites:
    x0 = -0.5 * spacing * (n_stations - 1)
    x = np.array([x0 + i * spacing for i in range(n_stations)])
    nav = NavigationTrack(
        sample_ids=tuple(f"AM{i:02d}" for i in range(n_stations)),
        easting=x,
        northing=np.zeros(n_stations),
    )
    tensor = _interstation_tensor(n_stations, x)
    line = build_airmt_line(
        "AIRMT_L001", nav, tensor, frequency=_AIRMT_FREQS,
    )
    return AirborneSites.from_line(line, technology="afmag_airmt")


def _original_afmag_sites(
    n_stations: int = 9, spacing: float = 60.0
) -> AirborneSites:
    x0 = -0.5 * spacing * (n_stations - 1)
    x = np.array([x0 + i * spacing for i in range(n_stations)])
    nav = NavigationTrack(
        sample_ids=tuple(f"OA{i:02d}" for i in range(n_stations)),
        easting=x,
        northing=np.zeros(n_stations),
    )
    tilt = np.zeros((n_stations, _ORIGINAL_FREQS.size))
    for i in range(n_stations):
        shape = (x[i] / 60.0) * np.exp(-((x[i] / 100.0) ** 2))
        tilt[i, :] = 8.0 * shape
    line = build_original_afmag_line(
        "ORIG_L001", nav, tilt, frequency=_ORIGINAL_FREQS,
    )
    return AirborneSites.from_line(line, technology="afmag_original")


# ─────────────────────────────────────────────────────────────────────────
# airmt_tilt_angles
# ─────────────────────────────────────────────────────────────────────────


def test_airmt_tilt_angles_reads_interstation_tensor_hz_row():
    df = airmt_tilt_angles(_airmt_sites())
    assert not df.empty
    expected_cols = {
        "station", "freq", "period",
        "tilt_real_deg", "tilt_real_azimuth_deg",
        "tilt_imag_deg", "tilt_imag_azimuth_deg",
        "tilt_resultant_deg",
    }
    assert expected_cols.issubset(df.columns)
    assert df["tilt_real_deg"].abs().max() > 0.0


def test_airmt_tilt_angles_empty_for_original_afmag_sites():
    # original-comparator sites have no interstation tensor at all.
    df = airmt_tilt_angles(_original_afmag_sites())
    assert df.empty


def test_airmt_tilt_angles_accepts_bare_directory_path(tmp_path):
    _airmt_sites().write_xml(tmp_path)
    df = airmt_tilt_angles(str(tmp_path))
    assert not df.empty


# ─────────────────────────────────────────────────────────────────────────
# original_afmag_tilt_table
# ─────────────────────────────────────────────────────────────────────────


def test_original_afmag_tilt_table_reads_scalar_tilt():
    df = original_afmag_tilt_table(_original_afmag_sites())
    assert not df.empty
    assert {"station", "freq", "period", "tilt_deg"}.issubset(df.columns)
    assert df["tilt_deg"].abs().max() > 0.0


def test_original_afmag_tilt_table_empty_for_airmt_sites():
    # AirMt sites have no afmag_tilt transfer function at all.
    df = original_afmag_tilt_table(_airmt_sites())
    assert df.empty


def test_original_afmag_tilt_table_accepts_bare_directory_path(tmp_path):
    _original_afmag_sites().write_xml(tmp_path)
    df = original_afmag_tilt_table(str(tmp_path))
    assert not df.empty


# ─────────────────────────────────────────────────────────────────────────
# Graceful degradation of the existing tipper-shaped functions
# ─────────────────────────────────────────────────────────────────────────


def test_afmag_tilt_angles_empty_for_airmt_sites():
    # afmag_tilt_angles reads Site.tipper; AirMt has no tipper TF.
    df = afmag_tilt_angles(_airmt_sites())
    assert df.empty


def test_flag_motion_susceptible_band_no_ops_gracefully_on_original_afmag():
    sites = _original_afmag_sites()
    out = flag_motion_susceptible_band(
        sites,
        inclination=60.0,
        declination=5.0,
        roll_amplitude_deg=10.0,
        pitch_amplitude_deg=5.0,
        threshold=1e-6,
    )
    assert isinstance(out, AirborneSites)
    assert len(out) == len(sites)


def test_flag_motion_susceptible_band_drop_rejected_for_airborne_sites():
    with pytest.raises(ValueError):
        flag_motion_susceptible_band(
            _original_afmag_sites(),
            inclination=60.0,
            declination=5.0,
            roll_amplitude_deg=10.0,
            pitch_amplitude_deg=5.0,
            threshold=1e-6,
            action="drop",
        )


# ─────────────────────────────────────────────────────────────────────────
# Plots
# ─────────────────────────────────────────────────────────────────────────


def test_plot_airmt_tilt_profile():
    ax = plot_airmt_tilt_profile(_airmt_sites())
    assert isinstance(ax, plt.Axes)
    plt.close("all")


def test_plot_original_afmag_tilt_profile():
    ax = plot_original_afmag_tilt_profile(_original_afmag_sites())
    assert isinstance(ax, plt.Axes)
    plt.close("all")


def test_plot_airmt_tilt_psection():
    ax = plot_airmt_tilt_psection(_airmt_sites())
    assert isinstance(ax, plt.Axes)
    plt.close("all")


def test_plot_original_afmag_dual_frequency_profile():
    ax = plot_original_afmag_dual_frequency_profile(_original_afmag_sites())
    assert isinstance(ax, plt.Axes)
    assert "peak-to-peak ratio" in ax.get_title()
    plt.close("all")


# ─────────────────────────────────────────────────────────────────────────
# Ward (1959)-style conductor diagnostics
# ─────────────────────────────────────────────────────────────────────────


def test_original_afmag_conductor_diagnostics_finds_the_synthetic_target():
    # _original_afmag_sites() builds a sign-reversing crossover shape
    # centred on the profile's middle station (see that fixture's
    # x-dependent shape function), so the crossover position should
    # land near the profile's own centre, and both frequencies are
    # identical here (same shape function used for both), giving a
    # peak-to-peak ratio of exactly 1.
    diag = original_afmag_conductor_diagnostics(_original_afmag_sites())
    profile = diag["profile"]
    center_m = float(profile["position_m"].median())
    assert diag["freq_low_hz"] == pytest.approx(150.0)
    assert diag["freq_high_hz"] == pytest.approx(510.0)
    assert abs(diag["crossover_low_m"] - center_m) < 60.0
    assert diag["peak_to_peak_ratio"] == pytest.approx(1.0, abs=1e-9)
    assert not profile.empty


def test_original_afmag_conductor_diagnostics_needs_two_frequencies():
    from pycsamt.airborne import NavigationTrack
    from pycsamt.airborne.afmag import build_original_afmag_line

    nav = NavigationTrack(
        sample_ids=("S1", "S2", "S3"),
        easting=np.array([0.0, 60.0, 120.0]),
        northing=np.zeros(3),
    )
    tilt = np.array([[1.0], [2.0], [-1.0]])
    line = build_original_afmag_line(
        "ONEFREQ", nav, tilt, frequency=np.array([150.0]),
    )
    single_freq_sites = AirborneSites.from_line(line)
    with pytest.raises(ValueError):
        original_afmag_conductor_diagnostics(single_freq_sites)


# ─────────────────────────────────────────────────────────────────────────
# Real synthetic data (skipped if not present locally)
# ─────────────────────────────────────────────────────────────────────────


@pytest.mark.skipif(
    not _ABITIBI.exists(), reason="synthetic original-AFMAG data not found"
)
def test_real_synthetic_original_afmag_end_to_end():
    df = original_afmag_tilt_table(str(_ABITIBI))
    assert not df.empty
    ax = plot_original_afmag_tilt_profile(str(_ABITIBI))
    assert isinstance(ax, plt.Axes)
    plt.close("all")
    diag = original_afmag_conductor_diagnostics(str(_ABITIBI))
    assert diag["freq_low_hz"] == pytest.approx(150.0)
    assert diag["freq_high_hz"] == pytest.approx(510.0)
    assert 0.0 < diag["peak_to_peak_ratio"] < 2.0
    ax2 = plot_original_afmag_dual_frequency_profile(str(_ABITIBI))
    assert isinstance(ax2, plt.Axes)
    plt.close("all")
    out = flag_motion_susceptible_band(
        str(_ABITIBI),
        inclination=60.0,
        declination=5.0,
        roll_amplitude_deg=10.0,
        pitch_amplitude_deg=5.0,
        threshold=1e-6,
    )
    assert isinstance(out, AirborneSites)
    assert len(out) == 13


@pytest.mark.skipif(
    not _YULONG.exists(), reason="synthetic AirMt data not found"
)
def test_real_synthetic_airmt_end_to_end():
    df = airmt_tilt_angles(str(_YULONG))
    assert not df.empty
    ax = plot_airmt_tilt_profile(str(_YULONG))
    assert isinstance(ax, plt.Axes)
    plt.close("all")
    ax2 = plot_airmt_tilt_psection(str(_YULONG))
    assert isinstance(ax2, plt.Axes)
    plt.close("all")
