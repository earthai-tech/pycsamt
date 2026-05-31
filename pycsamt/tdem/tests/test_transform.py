# -*- coding: utf-8 -*-
"""Tests for pycsamt.tdem transforms and data model."""

import numpy as np
import pytest

from pycsamt.tdem._base import TEMSounding
from pycsamt.tdem.transform import (
    LateTimeTransform,
    FourierTransform,
    _rho_a_late_time,
    _pseudo_freq,
    _phase_from_rho,
    _build_z_array,
    MU0,
)
from pycsamt.tdem.waveform import SquareWaveform, RampWaveform


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _synthetic_sounding(n=25, rho_true=100.0, current=8.0, loop_side=100.0):
    """
    Produce a synthetic TEMSounding for a uniform half-space with true
    apparent resistivity `rho_true` Ω·m.

    Forward formula (Ward & Hohmann 1988 / Christiansen 2009):
        |dBdt| = M × MU0^(5/2) / (10 × sqrt(π) × ρ^(3/2) × t^(5/2))
    """
    tx_area = loop_side ** 2
    M = current * tx_area
    t = np.logspace(-5, -2, n)

    dBdt = (M * MU0 ** 2.5) / (
        10.0 * np.sqrt(np.pi) * rho_true ** 1.5 * t ** 2.5
    )

    return TEMSounding(
        time_gates=t,
        data=dBdt,
        current=current,
        tx_area=tx_area,
        station_name="TEST",
        x=100.0,
        y=200.0,
        elevation=10.0,
    )


# ---------------------------------------------------------------------------
# TEMSounding tests
# ---------------------------------------------------------------------------

class TestTEMSounding:
    def test_basic_creation(self):
        snd = _synthetic_sounding()
        assert snd.n_gates == 25
        assert snd.moment == pytest.approx(8.0 * 100.0 ** 2)

    def test_dBdt_passthrough(self):
        snd = _synthetic_sounding()
        np.testing.assert_array_equal(snd.dBdt(), snd.data)

    def test_from_arrays_loop_side(self):
        t = np.logspace(-5, -2, 10)
        d = np.ones(10)
        snd = TEMSounding.from_arrays(t, d, current=5.0, loop_side=50.0)
        assert snd.tx_area == pytest.approx(2500.0)
        assert snd.moment == pytest.approx(5.0 * 2500.0)

    def test_from_arrays_loop_radius(self):
        t = np.logspace(-5, -2, 10)
        d = np.ones(10)
        snd = TEMSounding.from_arrays(t, d, current=1.0, loop_radius=10.0)
        assert snd.tx_area == pytest.approx(np.pi * 100.0)

    def test_shape_mismatch_raises(self):
        t = np.linspace(1e-5, 1e-2, 10)
        d = np.ones(8)
        with pytest.raises(ValueError, match="shape"):
            TEMSounding(t, d, current=1.0, tx_area=100.0)

    def test_invalid_data_type_raises(self):
        t = np.linspace(1e-5, 1e-2, 5)
        d = np.ones(5)
        with pytest.raises(ValueError, match="data_type"):
            TEMSounding(t, d, current=1.0, tx_area=100.0, data_type="bad")

    def test_voltage_to_dBdt(self):
        t = np.logspace(-5, -2, 10)
        rx_area, rx_turns = 0.5, 10
        dBdt_true = np.ones(10) * 1e-6
        voltage = dBdt_true * rx_area * rx_turns
        snd = TEMSounding(
            t, voltage, current=1.0, tx_area=100.0,
            data_type="voltage", rx_area=rx_area, rx_turns=rx_turns
        )
        np.testing.assert_allclose(snd.dBdt(), dBdt_true, rtol=1e-10)

    def test_repr(self):
        snd = _synthetic_sounding()
        assert "TEMSounding" in repr(snd)
        assert "TEST" in repr(snd)


# ---------------------------------------------------------------------------
# Late-time apparent resistivity
# ---------------------------------------------------------------------------

class TestRhoLatetime:
    def test_recovers_true_rho(self):
        rho_true = 100.0
        snd = _synthetic_sounding(n=30, rho_true=rho_true)
        rho_a = _rho_a_late_time(snd.dBdt(), snd.time_gates, snd.moment)
        np.testing.assert_allclose(rho_a, rho_true, rtol=1e-6)

    def test_half_space_50(self):
        rho_true = 50.0
        snd = _synthetic_sounding(n=20, rho_true=rho_true)
        rho_a = _rho_a_late_time(snd.dBdt(), snd.time_gates, snd.moment)
        np.testing.assert_allclose(rho_a, rho_true, rtol=1e-6)

    def test_nan_on_zero_data(self):
        t = np.array([1e-4, 1e-3, 1e-2])
        dBdt = np.array([0.0, 1e-8, 1e-9])
        rho = _rho_a_late_time(dBdt, t, 1e5)
        assert np.isnan(rho[0])
        assert np.isfinite(rho[1])


# ---------------------------------------------------------------------------
# Pseudo-frequency mapping
# ---------------------------------------------------------------------------

class TestPseudoFreq:
    def test_skin_depth_convention(self):
        t = np.array([1e-3])
        f = _pseudo_freq(t, "skin_depth")
        np.testing.assert_allclose(f, 1.0 / (2.0 * np.pi * 1e-3))

    def test_diffusion_convention(self):
        t = np.array([1e-3])
        f = _pseudo_freq(t, "diffusion")
        np.testing.assert_allclose(f, 1.0 / (4.0 * np.pi * 1e-3))

    def test_invalid_convention(self):
        with pytest.raises(ValueError):
            _pseudo_freq(np.array([1e-3]), "bad")


# ---------------------------------------------------------------------------
# Phase estimation
# ---------------------------------------------------------------------------

class TestPhaseEstimation:
    def test_homogeneous_is_45(self):
        rho = np.ones(10) * 100.0
        freq = np.logspace(1, 4, 10)
        phi = _phase_from_rho(rho, freq, mode="homogeneous")
        np.testing.assert_allclose(phi, 45.0)

    def test_weidelt_uniform_is_approx_45(self):
        # For a flat ρ_a curve, d(ln ρ)/d(ln ω) = 0, so φ ≈ 45°
        rho = np.ones(20) * 100.0
        freq = np.logspace(1, 4, 20)
        phi = _phase_from_rho(rho, freq, mode="weidelt")
        np.testing.assert_allclose(phi, 45.0, atol=1.0)

    def test_invalid_mode(self):
        with pytest.raises(ValueError):
            _phase_from_rho(np.ones(5), np.logspace(1, 3, 5), mode="xyz")


# ---------------------------------------------------------------------------
# Z array construction
# ---------------------------------------------------------------------------

class TestBuildZArray:
    def test_shape(self):
        n = 15
        rho = np.ones(n) * 50.0
        freq = np.logspace(0, 3, n)
        phi = np.full(n, 45.0)
        Z = _build_z_array(rho, phi, freq)
        assert Z.shape == (n, 2, 2)

    def test_1d_symmetry(self):
        n = 10
        rho = np.ones(n) * 100.0
        freq = np.logspace(0, 3, n)
        phi = np.full(n, 45.0)
        Z = _build_z_array(rho, phi, freq)
        # Zxx = Zyy = 0
        np.testing.assert_allclose(Z[:, 0, 0], 0.0)
        np.testing.assert_allclose(Z[:, 1, 1], 0.0)
        # Zyx = -Zxy
        np.testing.assert_allclose(Z[:, 1, 0], -Z[:, 0, 1])

    def test_rho_roundtrip(self):
        n = 20
        rho_true = 100.0
        freq = np.logspace(1, 4, n)
        rho = np.full(n, rho_true)
        phi = np.full(n, 45.0)
        Z = _build_z_array(rho, phi, freq)
        omega = 2.0 * np.pi * freq
        rho_back = np.abs(Z[:, 0, 1]) ** 2 / (omega * MU0)
        np.testing.assert_allclose(rho_back, rho_true, rtol=1e-10)


# ---------------------------------------------------------------------------
# LateTimeTransform integration
# ---------------------------------------------------------------------------

class TestLateTimeTransform:
    def test_result_keys(self):
        snd = _synthetic_sounding()
        tr = LateTimeTransform()
        res = tr.transform(snd)
        for key in ("freq", "Z", "Z_err", "rho_a", "phase_xy",
                    "station_name", "x", "y", "elevation"):
            assert key in res

    def test_freq_shape_matches_rho(self):
        snd = _synthetic_sounding(n=20)
        res = LateTimeTransform().transform(snd)
        assert res["freq"].shape == res["rho_a"].shape
        assert res["Z"].shape[0] == res["freq"].size

    def test_rho_roundtrip_100(self):
        snd = _synthetic_sounding(n=30, rho_true=100.0)
        res = LateTimeTransform().transform(snd)
        np.testing.assert_allclose(res["rho_a"], 100.0, rtol=1e-5)

    def test_rho_roundtrip_500(self):
        snd = _synthetic_sounding(n=30, rho_true=500.0)
        res = LateTimeTransform().transform(snd)
        np.testing.assert_allclose(res["rho_a"], 500.0, rtol=1e-5)

    def test_error_propagation(self):
        snd = _synthetic_sounding(n=20, rho_true=100.0)
        snd.error = snd.data * 0.05   # 5 % relative error
        res = LateTimeTransform().transform(snd)
        assert not np.all(np.isnan(res["Z_err"]))

    def test_station_metadata_preserved(self):
        snd = _synthetic_sounding()
        res = LateTimeTransform().transform(snd)
        assert res["station_name"] == "TEST"
        assert res["x"] == pytest.approx(100.0)
        assert res["y"] == pytest.approx(200.0)

    def test_freq_ascending(self):
        snd = _synthetic_sounding(n=20)
        res = LateTimeTransform().transform(snd)
        assert np.all(np.diff(res["freq"]) > 0)

    def test_weidelt_phase_mode(self):
        snd = _synthetic_sounding(n=30, rho_true=100.0)
        res = LateTimeTransform(phase_mode="weidelt").transform(snd)
        assert np.all(np.isfinite(res["phase_xy"]))

    def test_transform_many(self):
        snds = [_synthetic_sounding(rho_true=r) for r in [50, 100, 200]]
        results = LateTimeTransform().transform_many(snds)
        assert len(results) == 3

    def test_diffusion_convention(self):
        snd = _synthetic_sounding(n=20)
        res_sd = LateTimeTransform(freq_convention="skin_depth").transform(snd)
        res_df = LateTimeTransform(freq_convention="diffusion").transform(snd)
        # diffusion gives half the frequency of skin_depth
        np.testing.assert_allclose(
            res_df["freq"], res_sd["freq"] * 0.5, rtol=1e-6
        )


# ---------------------------------------------------------------------------
# FourierTransform stub
# ---------------------------------------------------------------------------

class TestFourierTransformStub:
    def test_raises_not_implemented(self):
        snd = _synthetic_sounding()
        with pytest.raises(NotImplementedError):
            FourierTransform().transform(snd)


# ---------------------------------------------------------------------------
# Waveform tests
# ---------------------------------------------------------------------------

class TestWaveforms:
    def test_square_off_time(self):
        wf = SquareWaveform(base_frequency=25.0)
        I = wf.current_at(np.array([0.011]))  # in the off-time of 25 Hz (hp=0.02)
        assert I[0] == pytest.approx(0.0)

    def test_square_on_time(self):
        wf = SquareWaveform(base_frequency=25.0)
        I = wf.current_at(np.array([0.005]))  # in the on-time (hp=0.02, duty=0.5)
        assert I[0] == pytest.approx(1.0)

    def test_ramp_midpoint(self):
        wf = RampWaveform(base_frequency=25.0, ramp_off=1e-4)
        I = wf.current_at(np.array([5e-5]))   # midpoint of ramp
        assert I[0] == pytest.approx(0.5)

    def test_ramp_post_ramp(self):
        wf = RampWaveform(base_frequency=25.0, ramp_off=1e-4)
        I = wf.current_at(np.array([2e-4]))
        assert I[0] == pytest.approx(0.0)
