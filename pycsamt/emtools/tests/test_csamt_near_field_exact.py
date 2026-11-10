"""Regression tests: CSAMT near-field handling against the exact half-space.

Reference responses come from ``CSAMT1DForward`` with a short bipole, which
calls empymod's full quasi-static solution (independent of the tabulated
ratio used by the correction itself).
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.emtools.fieldzone import near_field_factor
from pycsamt.emtools.source_effects import correct_near_field
from pycsamt.forward import LayeredModel
from pycsamt.forward.em1d import CSAMT1DForward, MT1DForward

MU0 = 4e-7 * np.pi
RHO = 100.0
FREQ = 2.0 ** np.arange(0, 14)


class _FakeZ:
    def __init__(self, z, freq):
        self.z = np.asarray(z, dtype=complex)
        self.freq = np.asarray(freq, dtype=float)


class _FakeSite:
    def __init__(self, station, z, freq, offset):
        self.station = station
        self.Z = _FakeZ(z, freq)
        self.freq = np.asarray(freq, dtype=float)
        self.source_offset = float(offset)

    def get_section(self, *_, **__):
        return None


def _exact_site(offset):
    resp = CSAMT1DForward(FREQ, source_offset=offset, dipole_length=1.0).run(
        LayeredModel(resistivity=[RHO], thickness=[])
    )
    zf = resp.z * 1e-3 / MU0  # SI ohm -> EDI field units
    z = np.zeros((FREQ.size, 2, 2), complex)
    z[:, 0, 1], z[:, 1, 0] = zf, -zf
    return _FakeSite("S01", z, FREQ, offset), resp


@pytest.mark.parametrize("offset", [3000.0, 6000.0, 10000.0])
def test_forward_near_field_raises_rho_and_lowers_phase(offset):
    _, resp = _exact_site(offset)
    delta = 503.0 * np.sqrt(RHO / FREQ)
    near = offset / delta < 2.0
    far = offset / delta > 30.0
    assert np.all(resp.rho_a[near] > 1.2 * RHO)
    assert np.all(resp.phase[near] < 40.0)
    np.testing.assert_allclose(resp.rho_a[far], RHO, rtol=0.01)
    np.testing.assert_allclose(resp.phase[far], 45.0, atol=0.5)


def test_forward_far_field_matches_mt():
    m = LayeredModel(resistivity=[100.0, 10.0, 1000.0], thickness=[300.0, 500.0])
    f = np.array([256.0, 1024.0, 4096.0])
    cs = CSAMT1DForward(f, source_offset=50_000.0, dipole_length=1000.0).run(m)
    mt = MT1DForward(f).run(m)
    np.testing.assert_allclose(cs.rho_a, mt.rho_a, rtol=0.02)
    np.testing.assert_allclose(cs.phase, mt.phase, atol=0.5)


@pytest.mark.parametrize("offset", [3000.0, 6000.0, 10000.0])
def test_correct_near_field_recovers_halfspace(offset):
    site, _ = _exact_site(offset)
    out = list(correct_near_field([site], source_offset=offset))
    assert len(out) == 1
    ed = getattr(out[0], "edi", None) or out[0]
    z = np.asarray(ed.Z.z)[:, 0, 1]
    rho = 0.2 * np.abs(z) ** 2 / FREQ
    phase = np.degrees(np.angle(z)) % 180.0
    np.testing.assert_allclose(rho, RHO, rtol=0.02)
    np.testing.assert_allclose(phase, 45.0, atol=1.0)


def test_near_field_factor_above_one_in_near_field():
    site, _ = _exact_site(3000.0)
    df = near_field_factor([site], source_offset=3000.0)
    low = df[df["freq_hz"] <= 4.0]
    high = df[df["freq_hz"] >= 2048.0]
    assert (low["nf_factor"] > 1.2).all()
    np.testing.assert_allclose(high["nf_factor"], 1.0, atol=0.01)
