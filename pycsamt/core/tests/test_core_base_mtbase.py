import numpy as np
import pytest
from numpy.testing import assert_allclose

from pycsamt.core.base import MTBase


@pytest.fixture(scope="module")
def mt():
    return MTBase()


# ----------------------------- Constants --------------------------------- #


def test_constants_self_consistency(mt):
    # C0 == 1/sqrt(mu0*eps0) and ~ C
    c0 = 1.0 / np.sqrt(mt.MU0 * mt.EPS0)
    assert_allclose(mt.C0, c0, rtol=1e-12, atol=0.0)
    assert_allclose(mt.C, mt.C0, rtol=5e-11, atol=0.0)
    # or just compare relative error explicitly:
    rel = abs(mt.C - mt.C0) / mt.C
    assert rel < 5e-11

    # ETA0 == sqrt(mu0/eps0)
    eta0 = np.sqrt(mt.MU0 / mt.EPS0)
    assert_allclose(mt.ETA0, eta0, rtol=1e-12, atol=0.0)

    # RHO_FACTOR == 1/(mu0*2*pi)
    rho_fact = 1.0 / (mt.MU0 * 2.0 * np.pi)
    assert_allclose(mt.RHO_FACTOR, rho_fact, rtol=1e-15, atol=0.0)

    # Unit conversion helpers sanity
    assert_allclose(1e-6, mt.MICROVOLTS_TO_VOLTS)
    assert_allclose(1e-12, mt.PICOTESLA_TO_TESLA)
    assert_allclose(1e-9, mt.NANOTESLA_TO_TESLA)
    assert_allclose(1e-6, mt.MV_PER_KM_TO_V_PER_M)
    assert_allclose(1e-3, mt.METERS_TO_KILOMETERS)
    assert_allclose(1e3, mt.KILOMETERS_TO_METERS)
    assert mt.PERCENT_FACTOR == 100.0

    # (mV/km)/nT to (V/m)/T factor
    assert_allclose(mt.Z_UNIT_MVK_NT_TO_SI, 1e3)


# --------------------------- Basic utilities ------------------------------ #


def test_omega(mt):
    f = np.array([0.0, 0.25, 1.0])
    w = mt.omega(f)
    assert_allclose(w, 2.0 * np.pi * f, rtol=0, atol=0)


def test_freq_period_roundtrip(mt):
    f = np.array([0.0, 0.5, 1.0, 2.0])
    T = mt.freq_to_period(f)
    f2 = mt.period_to_freq(T)
    assert_allclose(f2, f, rtol=0, atol=0)

    # Edge cases
    assert np.isinf(mt.freq_to_period(0.0))
    assert_allclose(mt.period_to_freq(np.inf), 0.0)


# -------------------------- Z, rho, phase -------------------------------- #


def test_rho_phase_and_inverse(mt):
    # Choose rho, f, phi then rebuild Z and recover
    rho = np.array([10.0, 100.0])  # ohm·m
    f = np.array([1.0, 0.1])  # Hz
    phi = np.array([45.0, 30.0])  # deg

    Z = mt.z_from_rho_phase(rho, phi, f, phase_unit="deg")
    rho2, phi2 = mt.rho_phase_from_z(Z, f, phase_unit="deg")

    assert_allclose(rho2, rho, rtol=1e-12, atol=0)
    assert_allclose(phi2, phi, rtol=1e-12, atol=0)


def test_determinant_and_invariant(mt):
    # Build a simple 2D-like tensor: Zxy = A e^{iθ}, Zyx = -A e^{iθ}
    A = 5.0
    theta = np.deg2rad(30.0)
    z = np.zeros((2, 2), dtype=complex)
    z[0, 1] = A * np.exp(1j * theta)
    z[1, 0] = -A * np.exp(1j * theta)

    z_det = mt.determinant_z(z)
    # For this construction, Z_det = A e^{iθ}
    assert_allclose(z_det, A * np.exp(1j * theta), rtol=1e-12, atol=0)

    f = 0.5
    rho_det, phi_det = mt.rho_phase_from_det(z, f, phase_unit="deg")
    rho_expected = (A**2) / (mt.MU0 * mt.omega(f))
    assert_allclose(rho_det, rho_expected, rtol=1e-12, atol=0)
    assert_allclose(phi_det, np.degrees(theta), rtol=1e-12, atol=0)


def test_rotate_impedance_invariant_det(mt):
    # complex tensor; determinant invariant magnitude should hold
    A, theta = 2.0, np.deg2rad(15.0)
    z = np.array([[0, A * np.exp(1j * theta)], [-A * np.exp(1j * theta), 0]])
    det0 = mt.determinant_z(z)
    z_rot = mt.rotate_impedance(z, theta_deg=37.0)
    det1 = mt.determinant_z(z_rot)
    assert_allclose(det1, det0, rtol=1e-12, atol=1e-12)


# --------------------------- Diffusion scales ----------------------------- #


def test_skin_depth(mt):
    rho = 100.0  # ohm·m
    f = 1.0  # Hz
    delta = mt.skin_depth(f, rho)

    # Rule-of-thumb: δ ≈ 503 * sqrt(rho/f)  (meters)
    approx = 503.0 * np.sqrt(rho / f)
    assert_allclose(delta, approx, rtol=5e-3, atol=0.0)  # ~0.5%


# ------------------------------- Tipper ----------------------------------- #


def test_tipper_amp_phase(mt):
    # Tx = 3 + 4j, Ty = 0 -> amp = |Tx| = 5, phase = arg(Tx)
    T = np.array([[3.0 + 4.0j, 0.0 + 0.0j]])
    amp, phi = mt.tipper_amp_phase(T, phase_unit="deg")
    assert_allclose(amp, [5.0], rtol=0, atol=0)
    assert_allclose(phi, [np.degrees(np.arctan2(4.0, 3.0))], rtol=0, atol=0)


def test_tipper_rotate(mt):
    # Rotate [1,0] by +90° -> [0,-1] under chosen convention
    T = np.array([[1.0 + 0.0j, 0.0 + 0.0j]])
    Tr = MTBase.tipper_rotate(T, 90.0)
    assert_allclose(Tr, np.array([0.0, -1.0]), rtol=1e-12, atol=1e-12)
    # or
    assert_allclose(Tr.reshape(1, 2), np.array([[0.0, -1.0]]), rtol=1e-12, atol=1e-12)


def test_induction_arrows(mt):
    # Wiese arrows from T=(1+0j,0) -> (ax,ay) = (0, +1)
    T = np.array([[1.0 + 0.0j, 0.0 + 0.0j]])
    ax, ay = MTBase.induction_arrows(T, convention="wiese", use_imag=False)
    assert_allclose(ax, [0.0])
    assert_allclose(ay, [1.0])

    # Parkinson uses components directly -> (1, 0)
    ax2, ay2 = MTBase.induction_arrows(T, convention="parkinson", use_imag=False)
    assert_allclose(ax2, [1.0])
    assert_allclose(ay2, [0.0])


# --------------------- Apparent conductivity (sigma_a) -------------------- #


def test_apparent_conductivity_from_z(mt):
    Z = np.array([5.0, 10.0])  # ohms
    f = np.array([1.0, 0.1])
    sigma = mt.apparent_conductivity_from_z(Z, f)
    sigma_expected = (mt.MU0 * mt.omega(f)) / (np.abs(Z) ** 2)
    assert_allclose(sigma, sigma_expected, rtol=0, atol=0)


# --------------------------- Half-space impedance ------------------------- #


def test_halfspace_impedance(mt):
    f = np.logspace(-3, 3, 7)
    rho = 100.0
    Z = mt.halfspace_impedance(f, rho)

    # Magnitude should be sqrt(mu0*w*rho); phase ~ 45°
    mag = np.abs(Z)
    ang = np.degrees(np.angle(Z))
    assert_allclose(mag, np.sqrt(mt.MU0 * mt.omega(f) * rho), rtol=1e-12, atol=0.0)
    assert_allclose(ang, np.full_like(ang, 45.0), rtol=1e-12, atol=1e-12)


# ------------------------ Mixed-units conversion -------------------------- #


def test_z_mvk_nt_to_ohms(mt):
    z_field = np.array([50.0])  # (mV/km)/nT
    Z_ohm = mt.z_mvk_nt_to_ohms(z_field)
    # Expected: z * 1e3 * mu0  ( (mV/km)/nT → (V/m)/T, then ×μ0 to Ω )
    expected = z_field * mt.Z_UNIT_MVK_NT_TO_SI * mt.MU0
    assert_allclose(Z_ohm, expected, rtol=0, atol=0)


# ----------------------------- Rotations ---------------------------------- #


def test_rotate_fields(mt):
    e = np.array([[1.0, 0.0]])
    h = np.array([[0.0, 1.0]])
    e_r, h_r = MTBase.rotate_fields(e, h, 90.0)

    # With the chosen rotation matrix, +90° gives:
    # R = [[0, 1], [-1, 0]]
    assert_allclose(e_r, np.array([[0.0, -1.0]]), rtol=1e-12, atol=1e-12)
    assert_allclose(h_r, np.array([[1.0, 0.0]]), rtol=1e-12, atol=1e-12)


# ----------------------------- Phase tensor -------------------------------- #


def test_phase_tensor_basic(mt):
    # Choose X = I, Y = diag(a, a) -> Phi = Y
    a = 0.5
    Z = np.zeros((2, 2), dtype=complex)
    Z.real[...] = np.eye(2)
    Z.imag[...] = np.array([[a, 0.0], [0.0, a]])
    Phi = mt.phase_tensor(Z)
    assert_allclose(Phi, Z.imag, rtol=1e-12, atol=0.0)


def test_phase_tensor_params_shapes_and_values(mt):
    # Use the same well-behaved Phi (real symmetric, equal eigenvalues)
    a = 0.5
    Z = np.zeros((4, 2, 2), dtype=complex)
    Z.real[...] = np.eye(2)
    Z.imag[...] = np.array([[a, 0.0], [0.0, a]])
    phi_max, phi_min, alpha, beta, ellipt = mt.phase_tensor_params(Z, angle_unit="deg")

    # Shapes
    assert phi_max.shape == (4,)
    assert phi_min.shape == (4,)
    assert alpha.shape == (4,)
    assert beta.shape == (4,)
    assert ellipt.shape == (4,)

    # For this construction: alpha≈0, beta≈0, ellipt≈0,
    # principal phases equal arctan(a) in degrees.
    expect_phase = np.degrees(np.arctan(a))
    assert_allclose(phi_max, np.full(4, expect_phase), rtol=1e-12, atol=1e-12)
    assert_allclose(phi_min, np.full(4, expect_phase), rtol=1e-12, atol=1e-12)
    assert_allclose(alpha, np.zeros(4), rtol=1e-12, atol=1e-12)
    assert_allclose(beta, np.zeros(4), rtol=1e-12, atol=1e-12)
    assert_allclose(ellipt, np.zeros(4), rtol=1e-12, atol=1e-12)


def test_phase_tensor_azimuth(mt):
    # Construct a Phi with a clear azimuth (~45 deg)
    # Z = X + iY with X = I, Y = [[0, b], [c, 0]] with b>0, c<0 gives az ~ 45°
    b, c = 1.0, -1.0
    Z = np.zeros((2, 2), dtype=complex)
    Z.real[...] = np.eye(2)
    Z.imag[...] = np.array([[0.0, b], [c, 0.0]])
    az = mt.phase_tensor_azimuth(Z, unit="deg")

    # Single tensor -> scalar; wrap in array for uniformity
    assert np.ndim(az) == 0
    assert_allclose(az, 45.0, rtol=1e-12, atol=1e-12)


# ----------------------------- Swift skew --------------------------------- #


def test_swift_skew_ideal_2d(mt):
    # 2D ideal: zxx=zyy=0, zxy=-zyx=Ae^{iθ} -> skew s=0
    A, theta = 3.0, np.deg2rad(20.0)
    Z = np.zeros((2, 2), dtype=complex)
    Z[0, 1] = A * np.exp(1j * theta)
    Z[1, 0] = -A * np.exp(1j * theta)
    amp, ang = mt.swift_skew(Z)
    assert_allclose(amp, 0.0, rtol=0, atol=1e-14)
    # angle arbitrary when amp=0, but should be finite
    assert np.isfinite(ang)


# ----------------------------- Misc helpers -------------------------------- #


def test_z_field_unit_consistency_vs_si(mt):
    # Compare rho computed via SI vs legacy Zonge formula using matched units
    E_mVkm = 100.0
    B_nT = 50.0
    f = 1.0

    # # Convert to SI Z_si_like = (V/m)/T
    # Z_si_like = (E_mVkm * mt.MV_PER_KM_TO_V_PER_M) / (B_nT * mt.NANOTESLA_TO_TESLA)
    # rho_si = (Z_si_like ** 2) * mt.RHO_FACTOR / f

    # # Legacy form
    # rho_zonge = (mt.ZONGE_RHO_FACTOR / f) * (E_mVkm / B_nT) ** 2

    # # They are not exactly equal unless full unit bookkeeping is done,
    # # but they should scale similarly (same order of magnitude).
    # ratio = rho_si / rho_zonge

    Z_si = (
        (E_mVkm * mt.MV_PER_KM_TO_V_PER_M) / (B_nT * mt.NANOTESLA_TO_TESLA)
    ) * mt.H_TO_B  # multiply by μ0
    rho_si = (Z_si**2) * mt.RHO_FACTOR / f
    rho_zonge = (mt.ZONGE_RHO_FACTOR / f) * (E_mVkm / B_nT) ** 2
    assert_allclose(rho_si, rho_zonge, rtol=1e-2, atol=0)
