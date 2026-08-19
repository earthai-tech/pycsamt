"""Tests for pycsamt.emtools.afmag"""

from __future__ import annotations

import numpy as np
import pytest

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.emtools.afmag import (
    afmag_tilt_angles,
    coil_normal_direction,
    correct_motion_induced_noise,
    euler_rotation_matrix,
    flag_motion_susceptible_band,
    geomagnetic_field_direction,
    motion_coupling_angle,
    motion_coupling_cosine,
    motion_susceptibility_table,
    plot_afmag_correction_comparison,
    plot_afmag_tilt_polar,
    plot_afmag_tilt_profile,
    plot_afmag_tilt_psection,
    plot_motion_coupling,
    plot_motion_susceptibility_map,
    simulate_motion_induced_voltage,
)

# ─────────────────────────────────────────────────────────────────────────
# Shared fixture (mirrors pycsamt/emtools/tests/test_emtools_tf.py)
# ─────────────────────────────────────────────────────────────────────────


class _FakeZ:
    def __init__(self, z, freq):
        self.z = np.asarray(z, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)


class _FakeTipper:
    def __init__(self, tipper, freq):
        self.tipper = np.asarray(tipper, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)


class _FakeSite:
    def __init__(
        self, station, z, freq, *, tipper=None, east=None, north=None
    ):
        self.station = station
        # ensure_sites/to_sites only recognizes an item as EDI-like once
        # it has a ``Z`` attribute, even for a tipper-only (AFMAG) test;
        # a dummy Z is therefore required, matching test_emtools_tf.py.
        self.Z = _FakeZ(z, freq)
        self.freq = np.asarray(freq, dtype=float)
        if tipper is not None:
            self.Tipper = _FakeTipper(tipper, freq)
        if east is not None:
            self.east = float(east)
            self.north = float(north)

    def get_section(self, *_, **__):
        return None


def _freqs(n: int = 12, f_lo: float = 1.0, f_hi: float = 1e3) -> np.ndarray:
    return np.logspace(np.log10(f_lo), np.log10(f_hi), n)


def _iso_z(freqs: np.ndarray, rho: float = 100.0) -> np.ndarray:
    amp = np.sqrt(5.0 * freqs * rho)
    z = np.zeros((freqs.size, 2, 2), dtype=complex)
    z[:, 0, 1] = amp * (1 + 1j) / np.sqrt(2)
    z[:, 1, 0] = -amp * (1 + 1j) / np.sqrt(2)
    return z


def _tipper(freqs: np.ndarray, amp: float = 0.15) -> np.ndarray:
    t = np.zeros((freqs.size, 2), dtype=complex)
    t[:, 0] = amp * np.exp(1j * np.linspace(0, np.pi / 4, freqs.size))
    t[:, 1] = amp * 0.5 * np.exp(1j * np.linspace(0, np.pi / 6, freqs.size))
    return t


def _site(
    name: str, *, with_tipper: bool = True, east=None, north=None
) -> _FakeSite:
    fr = _freqs()
    tip = _tipper(fr) if with_tipper else None
    return _FakeSite(name, _iso_z(fr), fr, tipper=tip, east=east, north=north)


def _profile(n_sites: int = 4) -> list:
    return [
        _site(f"S{i:02d}", east=i * 200.0, north=0.0)
        for i in range(n_sites)
    ]


# ─────────────────────────────────────────────────────────────────────────
# Section 1 — motion-coupling physics
# ─────────────────────────────────────────────────────────────────────────


class TestEulerRotationMatrix:
    def test_identity(self):
        r = euler_rotation_matrix(0.0, 0.0, 0.0)
        assert np.allclose(r, np.eye(3))

    def test_orthogonal_and_proper(self):
        r = euler_rotation_matrix(35.0, -20.0, 60.0)
        assert np.allclose(r.T, np.linalg.inv(r), atol=1e-10)
        assert np.isclose(np.linalg.det(r), 1.0, atol=1e-10)

    def test_batched_shape(self):
        yaw = np.linspace(0, 90, 10)
        r = euler_rotation_matrix(yaw, 0.0, 0.0)
        assert r.shape == (10, 3, 3)

    def test_roll_90_matches_hand_derivation(self):
        # yaw=pitch=0, roll=90deg: R = [[1,0,0],[0,0,-1],[0,1,0]]
        r = euler_rotation_matrix(0.0, 0.0, 90.0)
        expected = np.array(
            [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]]
        )
        assert np.allclose(r, expected, atol=1e-10)


class TestGeomagneticFieldDirection:
    def test_unit_vector(self):
        b = geomagnetic_field_direction(45.0, 30.0)
        assert np.isclose(np.linalg.norm(b), 1.0)

    def test_vertical_field_at_pole(self):
        b = geomagnetic_field_direction(90.0, 0.0)
        assert np.allclose(b, [0.0, 0.0, -1.0], atol=1e-10)

    def test_horizontal_field_at_equator(self):
        b = geomagnetic_field_direction(0.0, 0.0)
        assert np.allclose(b, [0.0, 1.0, 0.0], atol=1e-10)


class TestCoilNormalDirection:
    def test_identity_points_up(self):
        n = coil_normal_direction(0.0, 0.0, 0.0)
        assert np.allclose(n, [0.0, 0.0, 1.0])

    def test_unit_vector_for_arbitrary_attitude(self):
        n = coil_normal_direction(12.0, -8.0, 40.0)
        assert np.isclose(np.linalg.norm(n), 1.0)


class TestMotionCouplingCosine:
    def test_bounded(self):
        yaw = np.linspace(-180, 180, 50)
        pitch = np.linspace(-89, 89, 50)
        roll = np.linspace(-89, 89, 50)
        c = motion_coupling_cosine(yaw, pitch, roll, 45.0, 10.0)
        assert np.all(c >= -1.0 - 1e-9) and np.all(c <= 1.0 + 1e-9)

    def test_linear_in_roll_matches_paper_qualitative_result(self):
        # Liu et al. (2018) Fig. 2a: with yaw=pitch=declination=0,
        # theta(t) varies linearly with roll, and changing inclination
        # only shifts the curve (does not change its slope/shape).
        # Kept narrow enough to stay clear of the acute-angle wrap the
        # paper itself notes ("the sharp turns of the curve ... are
        # decided by the definition of Eq.7"), which is real, expected
        # behavior outside this range, not tested here.
        roll = np.linspace(-10.0, 10.0, 11)
        theta_1 = motion_coupling_angle(0.0, 0.0, roll, 30.0, 0.0)
        theta_2 = motion_coupling_angle(0.0, 0.0, roll, 50.0, 0.0)
        d_theta_1 = np.diff(theta_1)
        d_theta_2 = np.diff(theta_2)
        assert np.allclose(d_theta_1, d_theta_1[0], atol=1e-6)
        assert np.allclose(d_theta_2, d_theta_1, atol=1e-6)
        # different inclination only shifts theta(t), not its slope
        assert not np.allclose(theta_1, theta_2)

    def test_yaw_has_no_effect_at_zero_pitch_and_roll(self):
        # paper section 3: "the yaw angle has no effect ... for the
        # rotation axis is perpendicular to the coil plane" -- stated
        # for the pitch=roll=0 case (Fig. 2c/2f); with pitch/roll
        # nonzero, yaw generally does change theta(t), since the coil
        # normal is then no longer purely along the yaw rotation axis.
        yaw = np.linspace(-120.0, 120.0, 15)
        c = motion_coupling_cosine(yaw, 0.0, 0.0, 40.0, 20.0)
        assert np.allclose(c, c[0], atol=1e-10)


class TestSimulateMotionInducedVoltage:
    def test_zero_for_constant_cos_theta(self):
        cos_theta = np.full(20, 0.5)
        v = simulate_motion_induced_voltage(cos_theta, dt=0.01)
        assert np.allclose(v, 0.0)

    def test_rejects_bad_dt(self):
        with pytest.raises(ValueError):
            simulate_motion_induced_voltage(
                np.ones(5), dt=0.0,
            )


class TestCorrectMotionInducedNoise:
    def test_subtracts(self):
        measured = np.array([1.0, 2.0, 3.0])
        noise = np.array([0.1, 0.2, 0.3])
        out = correct_motion_induced_noise(measured, noise)
        assert np.allclose(out, [0.9, 1.8, 2.7])

    def test_shape_mismatch_raises(self):
        with pytest.raises(ValueError):
            correct_motion_induced_noise([1.0, 2.0], [1.0])


# ─────────────────────────────────────────────────────────────────────────
# Section 2/3 — Sites-based diagnostics
# ─────────────────────────────────────────────────────────────────────────


class TestAfmagTiltAngles:
    def test_table_shape_and_columns(self):
        sites = _profile(3)
        df = afmag_tilt_angles(sites)
        assert not df.empty
        expected_cols = {
            "station",
            "freq",
            "period",
            "tilt_real_deg",
            "tilt_real_azimuth_deg",
            "tilt_imag_deg",
            "tilt_imag_azimuth_deg",
            "tilt_resultant_deg",
        }
        assert expected_cols.issubset(df.columns)
        assert set(df["station"]) == {"S00", "S01", "S02"}

    def test_no_tipper_returns_empty(self):
        sites = [_site("S00", with_tipper=False)]
        df = afmag_tilt_angles(sites)
        assert df.empty


class TestMotionSusceptibilityTable:
    def test_one_row_per_station(self):
        sites = _profile(3)
        df = motion_susceptibility_table(
            sites,
            inclination=50.0,
            declination=0.0,
            roll_amplitude_deg=8.0,
            pitch_amplitude_deg=4.0,
        )
        assert len(df) == 3
        assert (df["susceptibility_score"] >= 0.0).all()

    def test_zero_amplitude_gives_zero_score(self):
        sites = _profile(1)
        df = motion_susceptibility_table(
            sites,
            inclination=50.0,
            declination=0.0,
            roll_amplitude_deg=0.0,
            pitch_amplitude_deg=0.0,
        )
        assert np.isclose(df["susceptibility_score"].iloc[0], 0.0, atol=1e-9)


class TestFlagMotionSusceptibleBand:
    def test_returns_sites_and_masks_band(self):
        sites = _profile(2)
        out = flag_motion_susceptible_band(
            sites,
            inclination=60.0,
            declination=5.0,
            roll_amplitude_deg=10.0,
            pitch_amplitude_deg=5.0,
            threshold=1e-6,
            band_hz=(1.0, 50.0),
        )
        assert out is not None
        df = afmag_tilt_angles(out)
        in_band = df[(df["freq"] >= 1.0) & (df["freq"] <= 50.0)]
        assert in_band["tilt_real_deg"].isna().all()

    def test_below_threshold_leaves_data_untouched(self):
        sites = _profile(1)
        before = afmag_tilt_angles(sites)
        out = flag_motion_susceptible_band(
            sites,
            inclination=60.0,
            declination=5.0,
            roll_amplitude_deg=10.0,
            pitch_amplitude_deg=5.0,
            threshold=10.0,  # unreachable score
            band_hz=(1.0, 50.0),
        )
        after = afmag_tilt_angles(out)
        assert np.allclose(
            before["tilt_real_deg"].to_numpy(),
            after["tilt_real_deg"].to_numpy(),
        )

    def test_invalid_action_raises(self):
        with pytest.raises(ValueError):
            flag_motion_susceptible_band(
                _profile(1),
                inclination=1.0,
                declination=1.0,
                roll_amplitude_deg=1.0,
                pitch_amplitude_deg=1.0,
                action="bogus",
            )

    def test_drop_action_removes_rows(self):
        sites = _profile(1)
        out = flag_motion_susceptible_band(
            sites,
            inclination=60.0,
            declination=5.0,
            roll_amplitude_deg=10.0,
            pitch_amplitude_deg=5.0,
            threshold=1e-6,
            band_hz=(1.0, 50.0),
            action="drop",
        )
        df = afmag_tilt_angles(out)
        assert not ((df["freq"] >= 1.0) & (df["freq"] <= 50.0)).any()


# ─────────────────────────────────────────────────────────────────────────
# Section 4 — plots
# ─────────────────────────────────────────────────────────────────────────


class TestPlots:
    def test_plot_afmag_tilt_profile(self):
        ax = plot_afmag_tilt_profile(_profile(4))
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_plot_afmag_tilt_profile_external_axes(self):
        _, ax_in = plt.subplots()
        ax = plot_afmag_tilt_profile(_profile(3), ax=ax_in)
        assert ax is ax_in
        plt.close("all")

    def test_plot_afmag_tilt_psection(self):
        ax = plot_afmag_tilt_psection(_profile(4))
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_plot_afmag_tilt_polar(self):
        ax = plot_afmag_tilt_polar(_profile(2))
        assert ax is not None
        plt.close("all")

    def test_plot_afmag_tilt_polar_no_tipper(self):
        sites = [_site("S00", with_tipper=False)]
        ax = plot_afmag_tilt_polar(sites)
        assert ax is not None
        plt.close("all")

    def test_plot_motion_coupling(self):
        t = np.linspace(0, 10, 40)
        ax = plot_motion_coupling(
            5 * np.sin(t), 3 * np.cos(t), 4 * np.sin(2 * t), 45.0, 10.0
        )
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_plot_motion_susceptibility_map(self):
        ax = plot_motion_susceptibility_map(
            _profile(3),
            inclination=55.0,
            declination=2.0,
            roll_amplitude_deg=8.0,
            pitch_amplitude_deg=4.0,
        )
        assert isinstance(ax, plt.Axes)
        plt.close("all")

    def test_plot_afmag_correction_comparison_auto_after(self):
        fig = plot_afmag_correction_comparison(
            _profile(3),
            inclination=60.0,
            declination=5.0,
            roll_amplitude_deg=10.0,
            pitch_amplitude_deg=5.0,
            threshold=1e-6,
            band_hz=(1.0, 50.0),
        )
        assert isinstance(fig, plt.Figure)
        plt.close("all")

    def test_plot_afmag_correction_comparison_requires_params(self):
        with pytest.raises(ValueError):
            plot_afmag_correction_comparison(_profile(2))

    def test_plot_afmag_correction_comparison_explicit_after(self):
        before = _profile(2)
        after = flag_motion_susceptible_band(
            before,
            inclination=60.0,
            declination=5.0,
            roll_amplitude_deg=10.0,
            pitch_amplitude_deg=5.0,
            threshold=1e-6,
            band_hz=(1.0, 50.0),
        )
        fig = plot_afmag_correction_comparison(before, after)
        assert isinstance(fig, plt.Figure)
        plt.close("all")
