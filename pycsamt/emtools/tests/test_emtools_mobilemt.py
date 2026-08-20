"""Tests for pycsamt.emtools.mobilemt"""

from __future__ import annotations

import numpy as np
import pytest

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.airborne import AirborneEMDataset, AirborneEMLine, NavigationTrack
from pycsamt.airborne.mobilemt import (
    MOBILEMT_APPARENT_CONDUCTIVITY_FIELD,
    MobileMTSystemSpec,
    build_mobilemt_line,
)
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

_MMT_FREQS = np.array([25.0, 100.0, 1000.0, 8000.0])


def _admittance_1d(freqs: np.ndarray, rho: float = 100.0) -> np.ndarray:
    """A synthetic uniform-half-space-like admittance, one sample.

    Built by inverting the classical impedance quarter-space form
    used elsewhere in this test suite (see
    ``test_emtools_afmag.py::_iso_z``): a 2x2 impedance
    ``Zxy = -Zyx = sqrt(5*f*rho)*(1+1j)/sqrt(2)`` is matrix-inverted
    to give the admittance ``Y = Z^-1``, then a small Hz row is
    appended to reach the MobileMT (nf, 3, 2) shape.
    """
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


def _line(
    line_id: str,
    n_samples: int = 6,
    *,
    spacing: float = 100.0,
    rho: float = 100.0,
    with_sigma: bool = True,
) -> AirborneEMLine:
    sample_ids = tuple(f"{line_id}.{i:03d}" for i in range(n_samples))
    easting = np.arange(n_samples, dtype=float) * spacing
    northing = np.zeros(n_samples)
    nav = NavigationTrack(
        sample_ids=sample_ids, easting=easting, northing=northing,
    )
    admittance = np.stack(
        [_admittance_1d(_MMT_FREQS, rho=rho) for _ in range(n_samples)],
        axis=0,
    )
    sigma = None
    if with_sigma:
        sigma = np.full((n_samples, _MMT_FREQS.size), 1.0 / rho)
    return build_mobilemt_line(
        line_id,
        nav,
        admittance,
        frequency=_MMT_FREQS,
        apparent_conductivity=sigma,
    )


def _dataset(**kwargs) -> AirborneEMDataset:
    line = _line("L001", **kwargs)
    return AirborneEMDataset(name="survey", lines={line.line_id: line})


# ─────────────────────────────────────────────────────────────────────────
# ensure_mobilemt_dataset
# ─────────────────────────────────────────────────────────────────────────


class TestEnsureMobilemtDataset:
    def test_passes_through_dataset(self):
        ds = _dataset()
        assert ensure_mobilemt_dataset(ds) is ds

    def test_wraps_single_line(self):
        line = _line("L002")
        ds = ensure_mobilemt_dataset(line)
        assert isinstance(ds, AirborneEMDataset)
        assert ds.get_line("L002") is line

    def test_invalid_type_raises(self):
        # An int is not a dataset/line/AirborneSites/AirborneSite/path.
        with pytest.raises(TypeError):
            ensure_mobilemt_dataset(42)

    def test_nonexistent_path_resolves_to_empty_dataset(self):
        # Strings are now accepted as EMTF-XML paths (see
        # test_emtools_mobilemt_airborne_sites.py); a path that
        # resolves to nothing mirrors ensure_asites's own
        # silent-empty-result convention rather than raising.
        ds = ensure_mobilemt_dataset("not-a-dataset")
        assert isinstance(ds, AirborneEMDataset)
        assert sum(line.n_records for line in ds.iter_lines()) == 0


# ─────────────────────────────────────────────────────────────────────────
# admittance_table
# ─────────────────────────────────────────────────────────────────────────


class TestAdmittanceTable:
    def test_table_shape_and_columns(self):
        ds = _dataset(n_samples=5)
        df = admittance_table(ds)
        assert not df.empty
        expected_cols = {
            "line_id", "sample_id", "x_m", "freq_hz", "period_s",
            "Yxx_real", "Yxx_imag", "Yxy_real", "Yxy_imag",
            "Yyx_real", "Yyx_imag", "Yyy_real", "Yyy_imag",
            "Yhzx_real", "Yhzx_imag", "Yhzy_real", "Yhzy_imag",
            "apparent_conductivity_native_Sm",
        }
        assert expected_cols.issubset(df.columns)
        assert len(df) == 5 * _MMT_FREQS.size

    def test_native_conductivity_carried_through(self):
        ds = _dataset(n_samples=3, rho=50.0, with_sigma=True)
        df = admittance_table(ds)
        assert np.allclose(
            df["apparent_conductivity_native_Sm"].to_numpy(), 1.0 / 50.0
        )

    def test_missing_native_conductivity_is_nan(self):
        ds = _dataset(n_samples=3, with_sigma=False)
        df = admittance_table(ds)
        assert df["apparent_conductivity_native_Sm"].isna().all()

    def test_empty_dataset_returns_empty(self):
        ds = AirborneEMDataset(name="empty")
        df = admittance_table(ds)
        assert df.empty


# ─────────────────────────────────────────────────────────────────────────
# admittance_determinant_table
# ─────────────────────────────────────────────────────────────────────────


class TestAdmittanceDeterminantTable:
    def test_recovers_the_known_half_space_resistivity(self):
        # Y was built as the matrix inverse of a known-rho impedance,
        # so the theoretical sigma_a must recover 1/rho.
        rho = 80.0
        ds = _dataset(n_samples=4, rho=rho)
        df = admittance_determinant_table(ds)
        assert not df.empty
        assert np.allclose(
            df["theoretical_sigma_a_Sm"].to_numpy(), 1.0 / rho, rtol=1e-6,
        )
        assert np.allclose(
            df["theoretical_rho_a_ohm_m"].to_numpy(), rho, rtol=1e-6,
        )

    def test_phase_matches_uniform_half_space(self):
        # A uniform half-space has a 45-degree MT determinant phase;
        # by construction Y = Z^-1 here, so the derived phase should
        # match the classical value.
        ds = _dataset(n_samples=2, rho=100.0)
        df = admittance_determinant_table(ds)
        assert np.allclose(
            df["theoretical_phase_deg"].to_numpy(), 45.0, atol=1.0,
        )

    def test_native_conductivity_carried_through(self):
        ds = _dataset(n_samples=2, rho=100.0, with_sigma=True)
        df = admittance_determinant_table(ds)
        assert np.allclose(
            df["apparent_conductivity_native_Sm"].to_numpy(), 0.01,
        )


# ─────────────────────────────────────────────────────────────────────────
# admittance_skew_table
# ─────────────────────────────────────────────────────────────────────────


class TestAdmittanceSkewTable:
    def test_zero_skew_for_ideal_1d_admittance(self):
        # For the synthetic 1D admittance, Yxx=Yyy=0, so skew is 0.
        ds = _dataset(n_samples=3)
        df = admittance_skew_table(ds)
        assert not df.empty
        assert np.allclose(df["skew"].to_numpy(), 0.0, atol=1e-9)

    def test_empty_dataset_returns_empty(self):
        df = admittance_skew_table(AirborneEMDataset(name="empty"))
        assert df.empty


# ─────────────────────────────────────────────────────────────────────────
# mask_outside_mobilemt_band
# ─────────────────────────────────────────────────────────────────────────


class TestMaskOutsideMobilemtBand:
    def test_masks_out_of_band_frequencies(self):
        ds = _dataset(n_samples=3)
        out = mask_outside_mobilemt_band(ds, band_hz=(50.0, 2000.0))
        df = admittance_table(out)
        out_of_band = df[
            (df["freq_hz"] < 50.0) | (df["freq_hz"] > 2000.0)
        ]
        assert out_of_band["Yxx_real"].isna().all()
        in_band = df[(df["freq_hz"] >= 50.0) & (df["freq_hz"] <= 2000.0)]
        assert not in_band["Yxx_real"].isna().any()

    def test_masks_native_conductivity_too(self):
        ds = _dataset(n_samples=2, with_sigma=True)
        out = mask_outside_mobilemt_band(ds, band_hz=(50.0, 2000.0))
        df = admittance_table(out)
        out_of_band = df[
            (df["freq_hz"] < 50.0) | (df["freq_hz"] > 2000.0)
        ]
        assert out_of_band["apparent_conductivity_native_Sm"].isna().all()

    def test_default_band_uses_mobilemt_system_spec(self):
        ds = _dataset(n_samples=2)
        out = mask_outside_mobilemt_band(ds)
        df = admittance_table(out)
        assert (df.dropna(subset=["Yxx_real"])["freq_hz"] >= 19.0).all()
        assert (
            df.dropna(subset=["Yxx_real"])["freq_hz"] <= 26_000.0
        ).all()

    def test_band_hz_and_system_spec_mutually_exclusive(self):
        with pytest.raises(ValueError):
            mask_outside_mobilemt_band(
                _dataset(n_samples=1),
                band_hz=(30.0, 100.0),
                system_spec=MobileMTSystemSpec(),
            )

    def test_invalid_system_spec_type_raises(self):
        with pytest.raises(TypeError):
            mask_outside_mobilemt_band(
                _dataset(n_samples=1), system_spec="not-a-spec",
            )

    def test_inplace_false_leaves_input_untouched(self):
        ds = _dataset(n_samples=2)
        before = admittance_table(ds)
        mask_outside_mobilemt_band(
            ds, band_hz=(50.0, 2000.0), inplace=False,
        )
        after = admittance_table(ds)
        assert np.allclose(
            before["Yxx_real"].to_numpy(), after["Yxx_real"].to_numpy(),
        )


# ─────────────────────────────────────────────────────────────────────────
# Plots
# ─────────────────────────────────────────────────────────────────────────


class TestPlots:
    def test_plot_mobilemt_admittance_profile(self):
        ax = plot_mobilemt_admittance_profile(_dataset(n_samples=5))
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_plot_mobilemt_admittance_profile_external_axes(self):
        _, ax_in = plt.subplots()
        ax = plot_mobilemt_admittance_profile(
            _dataset(n_samples=4), ax=ax_in,
        )
        assert ax is ax_in
        plt.close("all")

    def test_plot_mobilemt_admittance_profile_missing_line_id(self):
        ax = plot_mobilemt_admittance_profile(
            _dataset(n_samples=3), line_id="does-not-exist",
        )
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_plot_mobilemt_conductivity_psection_theoretical(self):
        ax = plot_mobilemt_conductivity_psection(_dataset(n_samples=6))
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_plot_mobilemt_conductivity_psection_native(self):
        ax = plot_mobilemt_conductivity_psection(
            _dataset(n_samples=6, with_sigma=True), source="native",
        )
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_plot_mobilemt_conductivity_psection_invalid_source(self):
        with pytest.raises(ValueError):
            plot_mobilemt_conductivity_psection(
                _dataset(n_samples=3), source="bogus",
            )

    def test_plot_mobilemt_skew_profile(self):
        ax = plot_mobilemt_skew_profile(_dataset(n_samples=5))
        assert isinstance(ax, plt.Axes)
        plt.close("all")
