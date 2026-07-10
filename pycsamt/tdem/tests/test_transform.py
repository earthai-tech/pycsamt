"""Tests for pycsamt.tdem transforms and data model."""

import numpy as np
import pytest

from pycsamt.tdem._base import TEMSounding
from pycsamt.tdem.transform import (
    MU0,
    FourierTransform,
    LateTimeTransform,
    _apply_waveform_correction,
    _biot_savart_circle_hz,
    _biot_savart_circle_hz_segments,
    _biot_savart_rect_hz,
    _biot_savart_rect_hz_segments,
    _build_z_array,
    _in_loop_geometry_factor,
    _in_loop_geometry_factor_td,
    _loop_inner_radius,
    _phase_from_rho,
    _pseudo_freq,
    _rho_a_in_loop,
    _rho_a_late_time,
    _rho_a_offset_loop,
    _waveform_moments,
)
from pycsamt.tdem.waveform import (
    CustomWaveform,
    HalfSineWaveform,
    RampWaveform,
    SquareWaveform,
)

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

# ---------------------------------------------------------------------------
# Biot-Savart helpers
# ---------------------------------------------------------------------------

class TestBiotSavartRect:
    def test_square_centre_known(self):
        """Hz at centre of square loop (side 2a) = I√2/(πa) per amp."""
        a = 50.0    # half-side
        hz = _biot_savart_rect_hz(0.0, 0.0, a, a)
        expected = np.sqrt(2.0) / (np.pi * a)
        np.testing.assert_allclose(hz, expected, rtol=1e-8)

    def test_symmetry_x(self):
        """Hz(r, 0) should equal Hz(-r, 0) for a square loop."""
        a = b = 40.0
        h1 = _biot_savart_rect_hz(10.0, 0.0, a, b)
        h2 = _biot_savart_rect_hz(-10.0, 0.0, a, b)
        np.testing.assert_allclose(h1, h2, rtol=1e-10)

    def test_symmetry_y(self):
        """Hz(0, r) should equal Hz(0, -r)."""
        a = b = 40.0
        h1 = _biot_savart_rect_hz(0.0, 15.0, a, b)
        h2 = _biot_savart_rect_hz(0.0, -15.0, a, b)
        np.testing.assert_allclose(h1, h2, rtol=1e-10)

    def test_off_centre_positive(self):
        """Hz must be positive for a receiver inside the loop."""
        a = b = 50.0
        hz = _biot_savart_rect_hz(20.0, 10.0, a, b)
        assert hz > 0.0

    def test_off_centre_larger_than_centre(self):
        """Hz at a point near the wire should be larger than at the centre."""
        a = b = 50.0
        hz_centre = _biot_savart_rect_hz(0.0, 0.0, a, b)
        hz_near   = _biot_savart_rect_hz(0.0, 45.0, a, b)  # close to top wire
        assert hz_near > hz_centre

    def test_rectangular_not_equal_square(self):
        """Hz at centre of a 2:1 rectangle ≠ Hz at centre of a square."""
        hz_sq   = _biot_savart_rect_hz(0.0, 0.0, 50.0, 50.0)
        hz_rect = _biot_savart_rect_hz(0.0, 0.0, 50.0, 25.0)
        assert not np.isclose(hz_sq, hz_rect)


class TestBiotSavartCircle:
    def test_circle_centre_known(self):
        """Hz at centre of circular loop = I/(2R) per amp."""
        R = 50.0
        hz = _biot_savart_circle_hz(0.0, 0.0, R)
        expected = 1.0 / (2.0 * R)
        np.testing.assert_allclose(hz, expected, rtol=1e-3)

    def test_off_centre_positive(self):
        R = 50.0
        hz = _biot_savart_circle_hz(20.0, 0.0, R)
        assert hz > 0.0

    def test_symmetry_radial(self):
        """Circular loop is azimuthally symmetric: Hz(r,0) = Hz(0,r)."""
        R = 50.0
        h1 = _biot_savart_circle_hz(15.0, 0.0, R)
        h2 = _biot_savart_circle_hz(0.0, 15.0, R)
        np.testing.assert_allclose(h1, h2, rtol=1e-3)


class TestInLoopGeometryFactor:
    def test_centre_is_one(self):
        eta = _in_loop_geometry_factor(0.0, 0.0, "square", (100.0,))
        assert eta == pytest.approx(1.0)

    def test_offset_increases_factor(self):
        """Moving toward the wire increases the local field → η > 1."""
        eta = _in_loop_geometry_factor(40.0, 0.0, "square", (100.0,))
        assert eta > 1.0

    def test_small_offset_near_unity(self):
        """For small offsets, η ≈ 1 (< 5 % deviation)."""
        eta = _in_loop_geometry_factor(5.0, 0.0, "square", (100.0,))
        assert abs(eta - 1.0) < 0.05

    def test_circle_shape(self):
        eta = _in_loop_geometry_factor(0.0, 0.0, "circle", (50.0,))
        assert eta == pytest.approx(1.0)

    def test_rectangular_shape(self):
        eta = _in_loop_geometry_factor(0.0, 0.0, "rectangular", (100.0, 50.0))
        assert eta == pytest.approx(1.0)


class TestLoopInnerRadius:
    def test_square(self):
        r = _loop_inner_radius("square", (100.0,))
        assert r == pytest.approx(50.0)

    def test_circle(self):
        r = _loop_inner_radius("circle", (75.0,))
        assert r == pytest.approx(75.0)

    def test_rectangular(self):
        r = _loop_inner_radius("rectangular", (200.0, 100.0))
        assert r == pytest.approx(50.0)  # half of the shorter side


# ---------------------------------------------------------------------------
# New apparent-resistivity formulas
# ---------------------------------------------------------------------------

class TestRhoAInLoop:
    def _make_central_snd(self, rho_true=100.0, n=25, loop_side=100.0):
        tx_area = loop_side ** 2
        M = 8.0 * tx_area
        t = np.logspace(-5, -2, n)
        dBdt = (M * MU0 ** 2.5) / (
            10.0 * np.sqrt(np.pi) * rho_true ** 1.5 * t ** 2.5
        )
        return dBdt, t, M

    def test_central_equals_late_time(self):
        """In-loop with (0,0) offset must give identical result to _rho_a_late_time."""
        dBdt, t, M = self._make_central_snd()
        rho_lt = _rho_a_late_time(dBdt, t, M)
        rho_il = _rho_a_in_loop(dBdt, t, M, "square", (100.0,), 0.0, 0.0)
        np.testing.assert_allclose(rho_il, rho_lt, rtol=1e-8)

    def test_offset_differs_from_central(self):
        """Off-centre receiver should produce a different ρ_a (higher near wire)."""
        dBdt, t, M = self._make_central_snd()
        rho_central = _rho_a_late_time(dBdt, t, M)
        rho_offset  = _rho_a_in_loop(dBdt, t, M, "square", (100.0,), 40.0, 0.0)
        # Near the wire the effective moment is larger → ρ_a is larger
        assert np.all(rho_offset > rho_central)

    def test_finite_results(self):
        dBdt, t, M = self._make_central_snd()
        rho = _rho_a_in_loop(dBdt, t, M, "circle", (50.0,), 10.0, 0.0)
        assert np.all(np.isfinite(rho))

    def test_rectangular_loop(self):
        dBdt, t, M = self._make_central_snd(loop_side=200.0)
        rho = _rho_a_in_loop(dBdt, t, M, "rectangular", (200.0, 100.0), 30.0, 0.0)
        assert np.all(np.isfinite(rho))


class TestRhoAOffsetLoop:
    def _make_offset_snd(self, rho_true=100.0, offset=500.0, n=25):
        M = 8.0 * 100.0 ** 2
        t = np.logspace(-5, -2, n)
        # Forward model: offset formula for σ = 1/rho_true
        dBdt = (M * MU0 ** 2.5) / (
            20.0 * np.sqrt(np.pi) * rho_true ** 1.5 * (offset ** 3) * t ** 2.5
        )
        return dBdt, t, M

    def test_recovers_true_rho(self):
        rho_true, offset = 100.0, 500.0
        dBdt, t, M = self._make_offset_snd(rho_true=rho_true, offset=offset)
        rho_a = _rho_a_offset_loop(dBdt, t, M, offset)
        np.testing.assert_allclose(rho_a, rho_true, rtol=1e-6)

    def test_different_offsets_same_model(self):
        """Changing the TX-RX separation changes ρ_a for fixed dBdt signal."""
        rho_true = 200.0
        dBdt, t, M = self._make_offset_snd(rho_true=rho_true, offset=500.0)
        rho500 = _rho_a_offset_loop(dBdt, t, M, 500.0)
        rho200 = _rho_a_offset_loop(dBdt, t, M, 200.0)
        assert not np.allclose(rho500, rho200)

    def test_invalid_offset_raises(self):
        dBdt = np.ones(10) * 1e-8
        t    = np.logspace(-5, -2, 10)
        with pytest.raises(ValueError, match="offset must be > 0"):
            _rho_a_offset_loop(dBdt, t, 1e4, 0.0)

    def test_finite_results(self):
        dBdt, t, M = self._make_offset_snd()
        rho = _rho_a_offset_loop(dBdt, t, M, 500.0)
        assert np.all(np.isfinite(rho))


# ---------------------------------------------------------------------------
# LateTimeTransform — config auto-detection
# ---------------------------------------------------------------------------

class TestLateTimeConfigDetection:
    def _make_snd(self, offset=0.0, loop_side=100.0):
        t = np.logspace(-5, -2, 20)
        M = 8.0 * loop_side ** 2
        dBdt = M * MU0 ** 2.5 / (10 * np.sqrt(np.pi) * 100.0 ** 1.5 * t ** 2.5)
        return TEMSounding(
            t, dBdt, current=8.0, tx_area=loop_side ** 2,
            offset=offset, loop_dims=(loop_side,), loop_shape="square",
        )

    def test_detects_central(self):
        snd = self._make_snd(offset=0.0)
        assert LateTimeTransform._detect_config(snd) == "central"

    def test_detects_in_loop(self):
        snd = self._make_snd(offset=30.0, loop_side=100.0)  # 30 < 50 (half-side)
        assert LateTimeTransform._detect_config(snd) == "in_loop"

    def test_detects_offset(self):
        snd = self._make_snd(offset=200.0, loop_side=100.0)  # 200 ≥ 50
        assert LateTimeTransform._detect_config(snd) == "offset"

    def test_result_has_loop_config_key(self):
        snd = self._make_snd(offset=0.0)
        res = LateTimeTransform().transform(snd)
        assert "loop_config" in res
        assert res["loop_config"] == "central"

    def test_in_loop_result_key(self):
        snd = self._make_snd(offset=30.0)
        res = LateTimeTransform().transform(snd)
        assert res["loop_config"] == "in_loop"

    def test_offset_loop_result_key(self):
        snd = self._make_snd(offset=200.0)
        res = LateTimeTransform().transform(snd)
        assert res["loop_config"] == "offset"

    def test_geometry_correction_false_uses_central(self):
        """With correction disabled, in-loop should give same ρ_a as central."""
        snd_central = self._make_snd(offset=0.0)
        snd_inloop  = self._make_snd(offset=30.0)
        tr = LateTimeTransform(loop_geometry_correction=False)
        res_c = tr.transform(snd_central)
        res_i = tr.transform(snd_inloop)
        np.testing.assert_allclose(res_c["rho_a"], res_i["rho_a"], rtol=1e-8)

    def test_in_loop_correction_increases_rho(self):
        """Off-centre Rx (near wire) → η > 1 → ρ_a_corrected > ρ_a_central."""
        snd_central = self._make_snd(offset=0.0)
        snd_inloop  = self._make_snd(offset=40.0)  # close to wire (half-side=50)
        tr = LateTimeTransform()
        res_c = tr.transform(snd_central)
        res_i = tr.transform(snd_inloop)
        # Both soundings have the same dBdt signal but different configurations;
        # the in-loop correction should produce a higher ρ_a near the wire.
        assert np.mean(res_i["rho_a"]) > np.mean(res_c["rho_a"])

    def test_rx_position_2d(self):
        """rx_position tuple (2-D) is honoured by LateTimeTransform."""
        t = np.logspace(-5, -2, 20)
        dBdt = 1e-8 * np.ones(20)
        snd = TEMSounding(
            t, dBdt, current=8.0, tx_area=100.0 ** 2,
            offset=30.0, loop_dims=(100.0,), loop_shape="square",
            rx_position=(20.0, 20.0),   # 2-D offset
        )
        res = LateTimeTransform().transform(snd)
        assert np.any(np.isfinite(res["rho_a"]))

    def test_offset_loop_recovers_rho(self):
        """Offset-loop transform should recover the true resistivity."""
        rho_true, offset = 150.0, 400.0
        M = 8.0 * 100.0 ** 2
        t = np.logspace(-5, -2, 25)
        dBdt = M * MU0 ** 2.5 / (
            20 * np.sqrt(np.pi) * rho_true ** 1.5 * offset ** 3 * t ** 2.5
        )
        snd = TEMSounding(
            t, dBdt, current=8.0, tx_area=100.0 ** 2,
            offset=offset, loop_dims=(100.0,), loop_shape="square",
        )
        res = LateTimeTransform().transform(snd)
        np.testing.assert_allclose(res["rho_a"], rho_true, rtol=1e-5)


# ---------------------------------------------------------------------------
# Waveform tests
# ---------------------------------------------------------------------------

def _analytical_halfsapce_dBdt(t, rho=100.0, a=50.0, current=1.0):
    """
    Exact central-loop step-off dBz/dt for a half-space (Ward & Hohmann 1988).
    Returns negative values (physically correct, decaying secondary field).
    """
    tx_area = np.pi * a ** 2
    moment = current * tx_area
    u = a * np.sqrt(MU0 / (4.0 * rho * t))
    dPhi_dt = -(4.0 * u ** 5 / np.sqrt(np.pi)) * np.exp(-(u ** 2)) / t
    return MU0 * moment / (2.0 * np.pi * a ** 3) * dPhi_dt


class TestFourierTransform:
    """Tests for FourierTransform: algorithmic correctness and output structure."""

    def _sounding_late_time(self, n=30, rho=100.0):
        """Late-time formula sounding (positive dBdt, standard helper)."""
        return _synthetic_sounding(n=n, rho_true=rho)

    def _sounding_analytical(self, n=30, rho=100.0, a=50.0):
        """Exact analytical half-space sounding (negative dBdt)."""
        t = np.logspace(-7, -2, n)
        dBdt = _analytical_halfsapce_dBdt(t, rho=rho, a=a)
        return TEMSounding(
            time_gates=t,
            data=dBdt,
            current=1.0,
            tx_area=np.pi * a ** 2,
        )

    # ── output structure ────────────────────────────────────────────────────

    def test_output_keys(self):
        snd = self._sounding_late_time()
        res = FourierTransform(n_freq=10).transform(snd)
        for key in ("freq", "Z", "Z_err", "rho_a", "phase_xy", "method"):
            assert key in res, f"missing key {key!r}"

    def test_method_key_is_fourier(self):
        snd = self._sounding_late_time()
        res = FourierTransform(n_freq=10).transform(snd)
        assert res["method"] == "fourier"

    def test_output_shapes(self):
        n_freq = 15
        snd = self._sounding_late_time()
        res = FourierTransform(n_freq=n_freq).transform(snd)
        assert res["freq"].shape[0] <= n_freq
        assert res["rho_a"].shape == res["freq"].shape
        assert res["phase_xy"].shape == res["freq"].shape
        assert res["Z"].ndim == 3
        assert res["Z"].shape[1:] == (2, 2)

    def test_freq_sorted(self):
        snd = self._sounding_late_time()
        res = FourierTransform(n_freq=15).transform(snd)
        assert np.all(np.diff(res["freq"]) > 0), "freqs must be ascending"

    # ── physical constraints ────────────────────────────────────────────────

    def test_rho_a_positive_and_finite(self):
        snd = self._sounding_late_time()
        res = FourierTransform(n_freq=20).transform(snd)
        assert np.all(res["rho_a"] > 0)
        assert np.all(np.isfinite(res["rho_a"]))

    def test_phase_in_range(self):
        """Phase must be in [0°, 90°] after clipping."""
        snd = self._sounding_analytical()
        res = FourierTransform(n_freq=20, waveform_correction=False).transform(snd)
        assert np.all(res["phase_xy"] >= 0.0)
        assert np.all(res["phase_xy"] <= 90.0)

    def test_analytical_negative_dBdt_gives_correct_im_sign(self):
        """Physically correct (negative) dBdt must yield Im[K] < 0."""
        from pycsamt.tdem.transform import (
            _cosine_transform_1d,
        )
        snd = self._sounding_analytical()
        t = snd.time_gates
        dBdt = snd.dBdt()
        omega_test = np.array([1e3])
        fc = _cosine_transform_1d(dBdt / snd.moment, t, omega_test, n_interp=0)
        # negative dBdt → negative cosine transform → im_k = fc/omega < 0
        assert fc[0] < 0

    def test_impedance_tensor_antisymmetric(self):
        """Z[i, 0, 1] == -Z[i, 1, 0] (MT anti-symmetry)."""
        snd = self._sounding_late_time()
        res = FourierTransform(n_freq=10).transform(snd)
        np.testing.assert_allclose(res["Z"][:, 0, 1], -res["Z"][:, 1, 0])

    # ── parameter variations ────────────────────────────────────────────────

    def test_custom_freq_range_honoured(self):
        snd = self._sounding_late_time()
        f_lo, f_hi = 1.0, 100.0
        res = FourierTransform(n_freq=20, freq_min=f_lo, freq_max=f_hi).transform(snd)
        assert res["freq"].min() >= f_lo * 0.99
        assert res["freq"].max() <= f_hi * 1.01

    def test_drop_nan_false_keeps_all_freqs(self):
        snd = self._sounding_late_time()
        res_drop = FourierTransform(n_freq=15, drop_nan=True).transform(snd)
        res_keep = FourierTransform(n_freq=15, drop_nan=False).transform(snd)
        assert len(res_keep["freq"]) >= len(res_drop["freq"])

    def test_n_aux_increased_changes_result(self):
        snd = self._sounding_late_time(n=40)
        res1 = FourierTransform(n_freq=10, n_aux=50).transform(snd)
        res2 = FourierTransform(n_freq=10, n_aux=200).transform(snd)
        # results differ when more KK quadrature points are used;
        # use relative comparison (values may be tiny, swamping allclose atol)
        rel_diff = np.abs(res1["rho_a"] - res2["rho_a"]) / (
            np.abs(res1["rho_a"]) + 1e-300
        )
        assert np.any(rel_diff > 1e-3)

    # ── ramp correction ─────────────────────────────────────────────────────

    def test_ramp_correction_changes_result(self):
        t = np.logspace(-4, -2, 25)
        dBdt = _analytical_halfsapce_dBdt(t)
        snd = TEMSounding(
            time_gates=t, data=dBdt, current=1.0, tx_area=np.pi * 50.0 ** 2,
            waveform=RampWaveform(base_frequency=25.0, ramp_off=5e-5),
        )
        res_no  = FourierTransform(n_freq=10, waveform_correction=False).transform(snd)
        res_yes = FourierTransform(n_freq=10, waveform_correction=True).transform(snd)
        # with ramp correction fewer gates may survive and/or freqs differ
        freqs_differ = not np.allclose(res_yes["freq"], res_no["freq"])
        shape_differs = res_yes["freq"].shape != res_no["freq"].shape
        assert shape_differs or freqs_differ

    def test_ramp_correction_reduces_gate_count(self):
        t = np.logspace(-4, -2, 30)
        dBdt = _analytical_halfsapce_dBdt(t)
        wf = RampWaveform(base_frequency=25.0, ramp_off=2e-4)
        TEMSounding(
            time_gates=t, data=dBdt, current=1.0, tx_area=np.pi * 50.0 ** 2,
            waveform=wf,
        )
        TEMSounding(time_gates=t, data=dBdt, current=1.0, tx_area=np.pi * 50.0 ** 2)
        # ramp correction drops early gates that fall within the ramp window
        corrected_d, corrected_t, _ = _apply_waveform_correction(dBdt, t, wf)
        assert len(corrected_t) < len(t)

    # ── error handling ──────────────────────────────────────────────────────

    def test_too_few_gates_raises(self):
        t = np.logspace(-4, -2, 3)
        snd = TEMSounding(time_gates=t, data=np.ones(3) * -1e-6,
                          current=1.0, tx_area=100.0)
        with pytest.raises(ValueError, match="4 time gates"):
            FourierTransform().transform(snd)

    # ── transform_many ──────────────────────────────────────────────────────

    def test_transform_many_returns_list(self):
        soundings = [self._sounding_late_time(n=20) for _ in range(3)]
        results = FourierTransform(n_freq=10).transform_many(soundings)
        assert isinstance(results, list)
        assert len(results) == 3

    def test_transform_many_consistent_with_single(self):
        snd = self._sounding_late_time(n=25)
        single = FourierTransform(n_freq=12).transform(snd)
        many = FourierTransform(n_freq=12).transform_many([snd])[0]
        np.testing.assert_array_equal(single["rho_a"], many["rho_a"])

    # ── station metadata propagation ────────────────────────────────────────

    def test_station_metadata_propagated(self):
        t = np.logspace(-5, -2, 20)
        snd = TEMSounding(
            time_gates=t, data=_synthetic_sounding(n=20).dBdt(),
            current=8.0, tx_area=1e4,
            station_name="ST01", x=10.0, y=20.0, elevation=5.0,
        )
        res = FourierTransform(n_freq=10).transform(snd)
        assert res["station_name"] == "ST01"
        assert res["x"] == pytest.approx(10.0)
        assert res["y"] == pytest.approx(20.0)
        assert res["elevation"] == pytest.approx(5.0)


# ---------------------------------------------------------------------------
# Waveform correction helpers
# ---------------------------------------------------------------------------

class TestWaveformMoments:
    """Tests for _waveform_moments across all waveform types."""

    def test_none_returns_zeros(self):
        assert _waveform_moments(None) == (0.0, 0.0)

    def test_square_returns_zeros(self):
        assert _waveform_moments(SquareWaveform(25.0)) == (0.0, 0.0)

    def test_ramp_analytic_tau_eff(self):
        wf = RampWaveform(base_frequency=25.0, ramp_off=2e-4)
        tau_eff, tau_window = _waveform_moments(wf)
        assert tau_eff == pytest.approx(1e-4)

    def test_ramp_analytic_tau_window(self):
        wf = RampWaveform(base_frequency=25.0, ramp_off=2e-4)
        tau_eff, tau_window = _waveform_moments(wf)
        assert tau_window == pytest.approx(2e-4)

    def test_halfsine_tau_eff_matches_analytic(self):
        # Analytic: tau_eff = hp*(pi-2)/(2*pi) for HalfSineWaveform
        hp = 0.5 / 25.0
        expected = hp * (np.pi - 2.0) / (2.0 * np.pi)
        tau_eff, tau_window = _waveform_moments(HalfSineWaveform(25.0))
        assert tau_eff == pytest.approx(expected, rel=1e-3)

    def test_halfsine_tau_window_is_zero(self):
        _, tau_window = _waveform_moments(HalfSineWaveform(25.0))
        assert tau_window == pytest.approx(0.0)

    def test_custom_ramp_matches_ramp_waveform(self):
        # A CustomWaveform with linear turn-off should give the same moments as RampWaveform
        tau_r = 3e-4
        t_wf = np.array([-1e-3, 0.0, tau_r, tau_r + 1e-5])
        I_wf = np.array([1.0,   1.0, 0.0,  0.0])
        cw   = CustomWaveform(t_wf, I_wf, base_frequency=25.0)
        tau_eff, tau_window = _waveform_moments(cw)
        assert tau_eff    == pytest.approx(tau_r / 2.0, rel=1e-2)
        # tau_window = 99th percentile ≈ 0.99*tau_r; allow 5% tolerance
        assert tau_window == pytest.approx(tau_r, rel=0.05)

    def test_tau_eff_positive(self):
        for wf in [RampWaveform(25.0, 1e-4), HalfSineWaveform(25.0)]:
            tau_eff, _ = _waveform_moments(wf)
            assert tau_eff >= 0.0


class TestApplyWaveformCorrection:
    """Tests for _apply_waveform_correction."""

    def _gates(self, n=30):
        return np.logspace(-4, -2, n)

    def _data(self, t):
        1.0 * (50.0 ** 2) * np.pi
        return _analytical_halfsapce_dBdt(t)

    def test_none_waveform_noop(self):
        t = self._gates()
        d = self._data(t)
        dc, tc, ec = _apply_waveform_correction(d, t, None)
        np.testing.assert_array_equal(d, dc)
        np.testing.assert_array_equal(t, tc)
        assert ec is None

    def test_square_waveform_noop(self):
        t = self._gates()
        d = self._data(t)
        dc, tc, _ = _apply_waveform_correction(d, t, SquareWaveform(25.0))
        np.testing.assert_array_equal(d, dc)
        np.testing.assert_array_equal(t, tc)

    def test_ramp_shifts_time(self):
        t = self._gates()
        d = self._data(t)
        tau_r = 1e-4
        dc, tc, _ = _apply_waveform_correction(d, t, RampWaveform(25.0, ramp_off=tau_r))
        # corrected times must be shifted by tau_r/2 for surviving gates
        surviving = t + tau_r / 2.0 > tau_r
        np.testing.assert_allclose(tc, t[surviving] + tau_r / 2.0)

    def test_ramp_amplitude_corrected(self):
        t = self._gates()
        d = self._data(t)
        tau_r = 1e-4
        dc, tc, _ = _apply_waveform_correction(d, t, RampWaveform(25.0, ramp_off=tau_r))
        surviving = t + tau_r / 2.0 > tau_r
        amp = 1.0 + (tau_r / 2.0) / (2.0 * tc)
        np.testing.assert_allclose(dc, d[surviving] * amp, rtol=1e-10)

    def test_ramp_rejects_early_gates(self):
        t = self._gates()
        d = self._data(t)
        tau_r = 5e-4  # large ramp to ensure some gates are rejected
        dc, tc, _ = _apply_waveform_correction(d, t, RampWaveform(25.0, ramp_off=tau_r))
        assert len(tc) < len(t)
        assert np.all(tc > tau_r / 2.0)

    def test_halfsine_shifts_time(self):
        t = self._gates()
        d = self._data(t)
        dc, tc, _ = _apply_waveform_correction(d, t, HalfSineWaveform(25.0))
        hp = 0.5 / 25.0
        expected_shift = hp * (np.pi - 2.0) / (2.0 * np.pi)
        np.testing.assert_allclose(tc, t + expected_shift, rtol=1e-3)

    def test_halfsine_no_gate_rejection(self):
        t = self._gates()
        d = self._data(t)
        dc, tc, _ = _apply_waveform_correction(d, t, HalfSineWaveform(25.0))
        assert len(tc) == len(t)

    def test_error_array_sliced_consistently(self):
        t = self._gates()
        d = self._data(t)
        err = np.abs(d) * 0.05
        tau_r = 5e-4
        dc, tc, ec = _apply_waveform_correction(d, t, RampWaveform(25.0, ramp_off=tau_r),
                                                 error=err)
        assert ec is not None
        assert len(ec) == len(tc)


class TestLateTimeTransformWaveform:
    """LateTimeTransform waveform_correction integration tests."""

    def _ramp_sounding(self, tau_r=1e-4, n=30):
        t = np.logspace(-5, -2, n)
        dBdt = _analytical_halfsapce_dBdt(t)
        return TEMSounding(
            time_gates=t, data=dBdt, current=1.0,
            tx_area=np.pi * 50.0 ** 2,
            waveform=RampWaveform(base_frequency=25.0, ramp_off=tau_r),
        )

    def test_correction_changes_result(self):
        snd = self._ramp_sounding(tau_r=5e-4)
        res_on  = LateTimeTransform(waveform_correction=True).transform(snd)
        res_off = LateTimeTransform(waveform_correction=False).transform(snd)
        # gate rejection changes shape or (if same shape) values differ
        assert (res_on["freq"].shape != res_off["freq"].shape
                or not np.allclose(res_on["freq"], res_off["freq"]))

    def test_correction_reduces_gate_count(self):
        snd = self._ramp_sounding(tau_r=5e-4)
        res_on  = LateTimeTransform(waveform_correction=True).transform(snd)
        res_off = LateTimeTransform(waveform_correction=False).transform(snd)
        assert len(res_on["freq"]) < len(res_off["freq"])

    def test_square_waveform_is_noop(self):
        t = np.logspace(-5, -2, 25)
        dBdt = _analytical_halfsapce_dBdt(t)
        snd_sq = TEMSounding(time_gates=t, data=dBdt, current=1.0,
                             tx_area=np.pi * 50.0 ** 2,
                             waveform=SquareWaveform(25.0))
        snd_no = TEMSounding(time_gates=t, data=dBdt, current=1.0,
                             tx_area=np.pi * 50.0 ** 2)
        res_sq = LateTimeTransform(waveform_correction=True).transform(snd_sq)
        res_no = LateTimeTransform(waveform_correction=True).transform(snd_no)
        np.testing.assert_array_equal(res_sq["rho_a"], res_no["rho_a"])

    def test_halfsine_waveform_shifts_freqs(self):
        t = np.logspace(-5, -2, 30)
        dBdt = _analytical_halfsapce_dBdt(t)
        snd = TEMSounding(time_gates=t, data=dBdt, current=1.0,
                          tx_area=np.pi * 50.0 ** 2,
                          waveform=HalfSineWaveform(25.0))
        res_on  = LateTimeTransform(waveform_correction=True).transform(snd)
        res_off = LateTimeTransform(waveform_correction=False).transform(snd)
        # time shift → all frequencies shift
        assert not np.allclose(res_on["freq"], res_off["freq"])

    def test_rho_a_positive_after_correction(self):
        snd = self._ramp_sounding()
        res = LateTimeTransform(waveform_correction=True).transform(snd)
        assert np.all(res["rho_a"] > 0)
        assert np.all(np.isfinite(res["rho_a"]))


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


# ---------------------------------------------------------------------------
# Per-segment Biot-Savart helpers
# ---------------------------------------------------------------------------

class TestBiotSavartRectHzSegments:
    """Tests for _biot_savart_rect_hz_segments."""

    def test_sum_matches_full_at_centre(self):
        # At the loop centre (0,0), segment sum must equal the static result
        a = b = 50.0
        hz_rx, hz_00, dist_rx, dist_00 = _biot_savart_rect_hz_segments(0.0, 0.0, a, b)
        assert pytest.approx(hz_rx.sum(), rel=1e-10) == hz_00.sum()

    def test_sum_matches_full_off_centre(self):
        # At (20,10) the segment sums must still equal the full Biot-Savart value
        a = b = 50.0
        hz_rx, hz_00, dist_rx, dist_00 = _biot_savart_rect_hz_segments(20.0, 10.0, a, b)
        full_rx  = _biot_savart_rect_hz(20.0, 10.0, a, b)
        full_00  = _biot_savart_rect_hz(0.0,  0.0,  a, b)
        assert pytest.approx(hz_rx.sum(), rel=1e-8) == full_rx
        assert pytest.approx(hz_00.sum(), rel=1e-8) == full_00

    def test_returns_four_segments(self):
        hz_rx, hz_00, dist_rx, dist_00 = _biot_savart_rect_hz_segments(5.0, 0.0, 50.0, 50.0)
        for arr in (hz_rx, hz_00, dist_rx, dist_00):
            assert len(arr) == 4

    def test_distances_bottom_wire(self):
        # Bottom wire of a 100×100 m loop is at y=-50; receiver at y=10 → dist = 60
        a = b = 50.0
        _, _, dist_rx, dist_00 = _biot_savart_rect_hz_segments(0.0, 10.0, a, b)
        assert dist_rx[0] == pytest.approx(60.0)   # |ry - (-b)|
        assert dist_00[0] == pytest.approx(50.0)    # b

    def test_distances_right_wire(self):
        # Right wire at x=+50; receiver at x=20 → dist = 30
        a = b = 50.0
        _, _, dist_rx, _ = _biot_savart_rect_hz_segments(20.0, 0.0, a, b)
        assert dist_rx[1] == pytest.approx(30.0)    # |a - rx|

    def test_rectangular_asymmetric(self):
        # 80×120 m loop; segment sums must match full Biot-Savart
        a, b = 40.0, 60.0
        hz_rx, hz_00, _, _ = _biot_savart_rect_hz_segments(10.0, -5.0, a, b)
        full_rx = _biot_savart_rect_hz(10.0, -5.0, a, b)
        full_00 = _biot_savart_rect_hz(0.0,   0.0, a, b)
        assert pytest.approx(hz_rx.sum(), rel=1e-8) == full_rx
        assert pytest.approx(hz_00.sum(), rel=1e-8) == full_00


class TestBiotSavartCircleHzSegments:
    """Tests for _biot_savart_circle_hz_segments."""

    def test_sum_matches_full_at_centre(self):
        r = 50.0
        hz_rx, hz_00, _, _ = _biot_savart_circle_hz_segments(0.0, 0.0, r)
        full = _biot_savart_circle_hz(0.0, 0.0, r)
        # centre: hz_rx ≈ hz_00 ≈ 1/(2r) per unit current in Gaussian units
        assert pytest.approx(hz_rx.sum(), rel=1e-3) == full

    def test_sum_matches_full_off_centre(self):
        r = 50.0
        rx, ry = 15.0, 10.0
        hz_rx, hz_00, _, _ = _biot_savart_circle_hz_segments(rx, ry, r, n_seg=720)
        full_rx = _biot_savart_circle_hz(rx, ry, r)
        full_00 = _biot_savart_circle_hz(0.0, 0.0, r)
        assert pytest.approx(hz_rx.sum(), rel=1e-3) == full_rx
        assert pytest.approx(hz_00.sum(), rel=1e-3) == full_00

    def test_n_seg_segments_returned(self):
        n = 180
        hz_rx, hz_00, dist_rx, dist_00 = _biot_savart_circle_hz_segments(0.0, 0.0, 30.0, n_seg=n)
        for arr in (hz_rx, hz_00, dist_rx, dist_00):
            assert len(arr) == n

    def test_dist_00_all_equal_radius(self):
        r = 40.0
        _, _, _, dist_00 = _biot_savart_circle_hz_segments(0.0, 0.0, r)
        np.testing.assert_allclose(dist_00, r)


# ---------------------------------------------------------------------------
# Time-dependent in-loop geometry factor
# ---------------------------------------------------------------------------

class TestInLoopGeometryFactorTD:
    """Tests for _in_loop_geometry_factor_td."""

    def _square_loop(self):
        return "square", [100.0]   # 100 m side

    def test_central_receiver_returns_ones(self):
        # rx=0, ry=0 must return η=1 for all time gates
        t = np.logspace(-5, -2, 20)
        rho_a = np.full(20, 100.0)
        eta = _in_loop_geometry_factor_td(0.0, 0.0, "square", [100.0], t, rho_a)
        np.testing.assert_allclose(eta, 1.0)

    def test_late_time_converges_to_static(self):
        # At very late time (large ρ_a) w_i → 1 for all segments → η_td → η_static
        t = np.logspace(-5, -2, 20)
        rho_a = np.full(20, 1e6)   # very high resistivity → large diffusion length
        eta_td = _in_loop_geometry_factor_td(20.0, 10.0, "square", [100.0], t, rho_a)
        eta_st = _in_loop_geometry_factor(20.0, 10.0, "square", [100.0])
        np.testing.assert_allclose(eta_td, eta_st, rtol=1e-4)

    def test_early_time_near_wire_exceeds_static(self):
        # Near-wire receiver: η_td(early) >> η_static because near segment dominates
        t = np.array([1e-6])       # very early gate
        rho_a = np.array([0.01])   # low resistivity → short diffusion length
        # rx=45 m, 5 m from right wire of 100 m loop (a=50)
        eta_td = _in_loop_geometry_factor_td(45.0, 0.0, "square", [100.0], t, rho_a)
        eta_st = _in_loop_geometry_factor(45.0, 0.0, "square", [100.0])
        assert float(eta_td[0]) > float(eta_st)

    def test_late_time_array_shape(self):
        n = 15
        t = np.logspace(-5, -2, n)
        rho_a = np.ones(n) * 50.0
        eta = _in_loop_geometry_factor_td(10.0, 0.0, "square", [100.0], t, rho_a)
        assert eta.shape == (n,)

    def test_rectangular_loop_shape(self):
        t = np.logspace(-5, -2, 10)
        rho_a = np.full(10, 1e6)
        eta_td = _in_loop_geometry_factor_td(10.0, 5.0, "rectangular", [80.0, 120.0], t, rho_a)
        eta_st = _in_loop_geometry_factor(10.0, 5.0, "rectangular", [80.0, 120.0])
        np.testing.assert_allclose(eta_td, eta_st, rtol=1e-4)

    def test_circle_loop_shape(self):
        t = np.logspace(-5, -2, 10)
        rho_a = np.full(10, 1e6)
        eta_td = _in_loop_geometry_factor_td(5.0, 0.0, "circle", [50.0], t, rho_a)
        eta_st = _in_loop_geometry_factor(5.0, 0.0, "circle", [50.0])
        np.testing.assert_allclose(eta_td, eta_st, rtol=1e-3)

    def test_eta_always_positive(self):
        t = np.logspace(-6, -2, 30)
        rho_a = np.logspace(-2, 4, 30)
        eta = _in_loop_geometry_factor_td(25.0, 15.0, "square", [100.0], t, rho_a)
        assert np.all(eta > 0.0)


# ---------------------------------------------------------------------------
# _rho_a_in_loop with time-dependent geometry (n_iter)
# ---------------------------------------------------------------------------

class TestRhoAInLoopTD:
    """Tests for the iterative geometry correction in _rho_a_in_loop."""

    def _make_data(self, n=25, rho_true=100.0):
        t = np.logspace(-5, -2, n)
        loop_side = 100.0
        M = 8.0 * loop_side ** 2
        dBdt = M * MU0 ** 2.5 / (10.0 * np.sqrt(np.pi) * rho_true ** 1.5 * t ** 2.5)
        return dBdt, t, M

    def test_central_receiver_unchanged_by_iter(self):
        # n_iter must not change result for on-centre receiver
        dBdt, t, M = self._make_data()
        r0 = _rho_a_in_loop(dBdt, t, M, "square", [100.0], 0.0, 0.0, n_iter=0)
        r3 = _rho_a_in_loop(dBdt, t, M, "square", [100.0], 0.0, 0.0, n_iter=3)
        np.testing.assert_allclose(r0, r3, rtol=1e-10)

    def test_n_iter_0_matches_static_only(self):
        # n_iter=0 should give the same result as using only η_static
        dBdt, t, M = self._make_data()
        r_iter0 = _rho_a_in_loop(dBdt, t, M, "square", [100.0], 20.0, 10.0, n_iter=0)
        eta_st   = _in_loop_geometry_factor(20.0, 10.0, "square", [100.0])
        r_static = _rho_a_late_time(dBdt, t, M * eta_st)
        np.testing.assert_allclose(r_iter0, r_static, rtol=1e-10)

    def test_n_iter_3_differs_from_static_at_early_time(self):
        # Near a wire, iterative correction should differ from static at early gates
        dBdt, t, M = self._make_data()
        r0 = _rho_a_in_loop(dBdt, t, M, "square", [100.0], 45.0, 0.0, n_iter=0)
        r3 = _rho_a_in_loop(dBdt, t, M, "square", [100.0], 45.0, 0.0, n_iter=3)
        # At least some early-time gates should differ
        assert not np.allclose(r0[:5], r3[:5], rtol=1e-4)

    def test_output_finite_and_positive(self):
        dBdt, t, M = self._make_data()
        r = _rho_a_in_loop(dBdt, t, M, "square", [100.0], 20.0, 10.0, n_iter=3)
        assert np.all(np.isfinite(r))
        assert np.all(r > 0.0)

    def test_n_iter_default_is_3(self):
        # Default n_iter=3 should match explicit n_iter=3
        dBdt, t, M = self._make_data()
        r_default = _rho_a_in_loop(dBdt, t, M, "square", [100.0], 20.0, 10.0)
        r3        = _rho_a_in_loop(dBdt, t, M, "square", [100.0], 20.0, 10.0, n_iter=3)
        np.testing.assert_array_equal(r_default, r3)


# ---------------------------------------------------------------------------
# LateTimeTransform with in_loop_n_iter parameter
# ---------------------------------------------------------------------------

class TestLateTimeTransformInLoopNIter:
    """Tests for LateTimeTransform.in_loop_n_iter propagation."""

    def _sounding_in_loop(self, rx=20.0, ry=10.0, n=25, rho_true=100.0):
        loop_side = 100.0
        M = 8.0 * loop_side ** 2
        t = np.logspace(-5, -2, n)
        dBdt = M * MU0 ** 2.5 / (10.0 * np.sqrt(np.pi) * rho_true ** 1.5 * t ** 2.5)
        return TEMSounding(
            time_gates=t, data=dBdt,
            current=8.0, tx_area=loop_side ** 2,
            loop_shape="square", loop_dims=[loop_side],
            rx_offset=(rx, ry),
        )

    def test_default_n_iter_is_3(self):
        ltt = LateTimeTransform()
        assert ltt.in_loop_n_iter == 3

    def test_n_iter_0_vs_3_differ_off_centre(self):
        snd = self._sounding_in_loop(rx=45.0, ry=0.0)
        res0 = LateTimeTransform(in_loop_n_iter=0).transform(snd)
        res3 = LateTimeTransform(in_loop_n_iter=3).transform(snd)
        # result dict is sorted by ascending frequency (late→early time)
        # the TD correction affects early-time (high-freq) gates → compare tail
        assert not np.allclose(res0["rho_a"][-5:], res3["rho_a"][-5:], rtol=1e-4)

    def test_n_iter_0_vs_3_same_on_centre(self):
        snd = self._sounding_in_loop(rx=0.0, ry=0.0)
        res0 = LateTimeTransform(in_loop_n_iter=0).transform(snd)
        res3 = LateTimeTransform(in_loop_n_iter=3).transform(snd)
        np.testing.assert_allclose(res0["rho_a"], res3["rho_a"], rtol=1e-8)

    def test_rho_a_positive_finite(self):
        snd = self._sounding_in_loop()
        res = LateTimeTransform(in_loop_n_iter=3).transform(snd)
        assert np.all(np.isfinite(res["rho_a"]))
        assert np.all(res["rho_a"] > 0.0)
