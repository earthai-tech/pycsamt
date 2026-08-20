"""Tests for pycsamt.emtools.mobilemt accepting AirborneSites/paths.

``ensure_mobilemt_dataset`` deliberately keeps returning an
``AirborneEMDataset`` -- this module's functions are organized around
flight-line chainage, not a flat station list (see the module
docstring for why that stays one-directional, unlike ZTEM/AFMAG's
``ensure_any_sites``). What this file covers is the *input* side: an
``AirborneSites``/``AirborneSite`` produced elsewhere in ``emtools``,
or a bare path/directory of MobileMT EMTF-XML, should flow straight
into every public function here without an intermediate
``AirborneEMDataset`` construction step.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.airborne import AirborneEMDataset, NavigationTrack
from pycsamt.airborne.mobilemt import build_mobilemt_line
from pycsamt.airborne.site import AirborneSite, AirborneSites, ensure_asites
from pycsamt.emtools.mobilemt import (
    admittance_determinant_table,
    admittance_skew_table,
    admittance_table,
    ensure_mobilemt_dataset,
    mask_outside_mobilemt_band,
    plot_mobilemt_admittance_profile,
    plot_mobilemt_conductivity_psection,
    plot_mobilemt_skew_profile,
)

_DATA_ROOT = Path(__file__).resolve().parents[3] / "data"
_FLAMMEFJELD = _DATA_ROOT / "mobileMT" / "flammefjeld_greenland"
_TIMISKAMING = _DATA_ROOT / "mobileMT" / "timiskaming_kimberlite_on"

_MMT_FREQS = np.array([25.0, 100.0, 1000.0, 8000.0])


def _admittance_1d(freqs: np.ndarray, rho: float = 100.0) -> np.ndarray:
    amp = np.sqrt(5.0 * freqs * rho)
    z = np.zeros((freqs.size, 2, 2), dtype=complex)
    z[:, 0, 1] = amp * (1 + 1j) / np.sqrt(2)
    z[:, 1, 0] = -amp * (1 + 1j) / np.sqrt(2)
    y2 = np.linalg.inv(z)
    y = np.zeros((freqs.size, 3, 2), dtype=complex)
    y[:, :2, :] = y2
    y[:, 2, 0] = 0.01 * y2[:, 0, 0]
    y[:, 2, 1] = 0.01 * y2[:, 1, 1]
    return y


def _line(line_id: str, n_samples: int = 6, spacing: float = 100.0):
    sample_ids = tuple(f"{line_id}.{i:03d}" for i in range(n_samples))
    easting = np.arange(n_samples, dtype=float) * spacing
    northing = np.zeros(n_samples)
    nav = NavigationTrack(
        sample_ids=sample_ids, easting=easting, northing=northing,
    )
    admittance = np.stack(
        [_admittance_1d(_MMT_FREQS) for _ in range(n_samples)], axis=0,
    )
    return build_mobilemt_line(
        line_id, nav, admittance, frequency=_MMT_FREQS,
    )


def _asites(n_lines: int = 1) -> AirborneSites:
    items = []
    for i in range(n_lines):
        line = _line(f"L{i:03d}")
        items.extend(
            AirborneSites.from_line(line, technology="mobilemt").as_list()
        )
    return AirborneSites(items)


# ─────────────────────────────────────────────────────────────────────────
# ensure_mobilemt_dataset
# ─────────────────────────────────────────────────────────────────────────


def test_ensure_mobilemt_dataset_from_airborne_sites_one_line():
    ds = ensure_mobilemt_dataset(_asites())
    assert isinstance(ds, AirborneEMDataset)
    assert len(list(ds.iter_lines())) == 1
    line = next(ds.iter_lines())
    assert line.n_records == 6


def test_ensure_mobilemt_dataset_from_airborne_sites_regroups_by_line_id():
    ds = ensure_mobilemt_dataset(_asites(n_lines=2))
    line_ids = {line.line_id for line in ds.iter_lines()}
    assert len(line_ids) == 2


def test_ensure_mobilemt_dataset_from_single_airborne_site():
    one = _asites()[0]
    assert isinstance(one, AirborneSite)
    ds = ensure_mobilemt_dataset(one)
    assert isinstance(ds, AirborneEMDataset)
    assert sum(line.n_records for line in ds.iter_lines()) == 1


def test_ensure_mobilemt_dataset_rejects_unsupported_type():
    with pytest.raises(TypeError):
        ensure_mobilemt_dataset(42)


# ─────────────────────────────────────────────────────────────────────────
# Tables
# ─────────────────────────────────────────────────────────────────────────


def test_admittance_table_accepts_airborne_sites():
    df = admittance_table(_asites())
    assert not df.empty
    assert len(df) == 6 * _MMT_FREQS.size


def test_admittance_determinant_table_accepts_airborne_sites():
    df = admittance_determinant_table(_asites())
    assert not df.empty


def test_admittance_skew_table_accepts_airborne_sites():
    df = admittance_skew_table(_asites())
    assert not df.empty


def test_tables_agree_between_dataset_and_airborne_sites_input():
    line = _line("L001")
    ds = AirborneEMDataset(name="s", lines={line.line_id: line})
    asites = AirborneSites.from_line(line, technology="mobilemt")

    df_ds = admittance_table(ds)
    df_as = admittance_table(asites)
    assert len(df_ds) == len(df_as)
    np.testing.assert_allclose(
        df_ds["Yxy_real"].to_numpy(), df_as["Yxy_real"].to_numpy(),
    )


# ─────────────────────────────────────────────────────────────────────────
# mask_outside_mobilemt_band
# ─────────────────────────────────────────────────────────────────────────


def test_mask_outside_mobilemt_band_accepts_airborne_sites():
    out = mask_outside_mobilemt_band(_asites(), band_hz=(50.0, 5000.0))
    assert isinstance(out, AirborneEMDataset)
    line = next(out.iter_lines())
    record = next(line.iter_records())
    tf = record.emtf.get_transfer_function("mobilemt_admittance")
    Y = np.asarray(tf.data)
    assert not np.isfinite(Y[0]).any()  # 25 Hz masked out
    assert np.isfinite(Y[1]).any()  # 100 Hz kept


def test_mask_outside_mobilemt_band_does_not_mutate_source_sites():
    asites = _asites()
    original = np.asarray(
        asites[0].admittance
    ).copy()
    mask_outside_mobilemt_band(asites, band_hz=(50.0, 5000.0))
    np.testing.assert_array_equal(np.asarray(asites[0].admittance), original)


# ─────────────────────────────────────────────────────────────────────────
# Plots
# ─────────────────────────────────────────────────────────────────────────


def test_plot_mobilemt_admittance_profile_accepts_airborne_sites():
    ax = plot_mobilemt_admittance_profile(_asites())
    assert isinstance(ax, plt.Axes)
    plt.close("all")


def test_plot_mobilemt_conductivity_psection_accepts_airborne_sites():
    ax = plot_mobilemt_conductivity_psection(_asites())
    assert isinstance(ax, plt.Axes)
    plt.close("all")


def test_plot_mobilemt_skew_profile_accepts_airborne_sites():
    ax = plot_mobilemt_skew_profile(_asites())
    assert isinstance(ax, plt.Axes)
    plt.close("all")


# ─────────────────────────────────────────────────────────────────────────
# Bare path / directory input
# ─────────────────────────────────────────────────────────────────────────


def test_ensure_mobilemt_dataset_accepts_directory_path(tmp_path):
    _asites().write_xml(tmp_path)
    ds = ensure_mobilemt_dataset(str(tmp_path))
    assert isinstance(ds, AirborneEMDataset)
    assert sum(line.n_records for line in ds.iter_lines()) == 6


def test_admittance_table_accepts_bare_directory_path(tmp_path):
    _asites().write_xml(tmp_path)
    df = admittance_table(str(tmp_path))
    assert not df.empty


# ─────────────────────────────────────────────────────────────────────────
# Real synthetic data (skipped if not present locally)
# ─────────────────────────────────────────────────────────────────────────


@pytest.mark.skipif(
    not _FLAMMEFJELD.exists(), reason="synthetic MobileMT data not found"
)
def test_real_synthetic_flammefjeld_end_to_end_via_path():
    df = admittance_table(str(_FLAMMEFJELD))
    assert not df.empty
    asites = ensure_asites(str(_FLAMMEFJELD))
    df2 = admittance_determinant_table(asites)
    assert not df2.empty
    ax = plot_mobilemt_conductivity_psection(asites)
    assert isinstance(ax, plt.Axes)
    plt.close("all")


@pytest.mark.skipif(
    not _TIMISKAMING.exists(), reason="synthetic MobileMT data not found"
)
def test_real_synthetic_timiskaming_end_to_end_via_asites():
    asites = ensure_asites(str(_TIMISKAMING))
    out = mask_outside_mobilemt_band(asites, band_hz=(50.0, 20000.0))
    assert isinstance(out, AirborneEMDataset)


# ─────────────────────────────────────────────────────────────────────────
# Native apparent-conductivity XML round trip (see
# _recover_native_apparent_conductivity in pycsamt.emtools.mobilemt)
# ─────────────────────────────────────────────────────────────────────────


@pytest.mark.skipif(
    not _FLAMMEFJELD.exists(), reason="synthetic MobileMT data not found"
)
def test_real_synthetic_flammefjeld_native_apparent_conductivity_recovered():
    # AirborneSite.from_xml never populated AirborneEMRecord.fields, so
    # this silently returned all-nan before the read-side recovery fix.
    df = admittance_table(str(_FLAMMEFJELD))
    sigma = df["apparent_conductivity_native_Sm"]
    assert sigma.notna().all()
    assert (sigma > 0.0).all()


def test_native_apparent_conductivity_survives_xml_write_and_read(tmp_path):
    # One vector per sample, one value per frequency in _MMT_FREQS
    # (ascending), matching how the real generator note is shaped.
    sigma = np.array([0.01, 0.02, 0.03, 0.04])
    line = _line("L000", n_samples=1)
    # Mirror scripts/airborne/generate_synthetic_survey_xml.py's own
    # write-side workaround: AirborneEMRecord.fields has no EMTF-XML
    # representation, so the native vector must be stashed as text in
    # metadata["notes"] to survive a write_xml/ensure_asites round trip.
    doc = line.get_record(line.navigation.sample_ids[0]).emtf
    doc.metadata.setdefault("notes", {}).setdefault("MobileMT", {})[
        "ApparentConductivitySm"
    ] = ",".join(f"{v:.10g}" for v in sigma)

    sites = AirborneSites.from_line(line, technology="mobilemt")
    sites.write_xml(tmp_path)

    reloaded = ensure_asites(str(tmp_path))
    df = admittance_table(reloaded).sort_values("freq_hz")
    recovered = df["apparent_conductivity_native_Sm"].to_numpy()
    np.testing.assert_allclose(recovered, sigma)
