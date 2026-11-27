# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for :mod:`pycsamt.emtools.frames`."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.emtools.frames import (
    frame_azimuth,
    rotate_tensor,
    to_frame,
)
from pycsamt.utils.zmath import rotatematrix_incl_errors


def _rand(n=6, seed=0):
    rng = np.random.default_rng(seed)
    z = rng.normal(size=(n, 2, 2)) + 1j * rng.normal(size=(n, 2, 2))
    e = np.abs(rng.normal(size=(n, 2, 2))) * 0.1
    return z, e


@pytest.mark.parametrize("theta", [-90.0, 90.0, 180.0, 270.0, 360.0])
def test_quarter_turns_match_general_rotation(theta):
    z, e = _rand()
    zr, er = rotate_tensor(z, e, theta)
    for n in range(z.shape[0]):
        m, em = rotatematrix_incl_errors(z[n], theta, e[n])
        np.testing.assert_allclose(zr[n], m, atol=1e-12)
        np.testing.assert_allclose(er[n], em, atol=1e-12)


def test_quarter_turn_keeps_masks_per_element():
    z, e = _rand()
    z[2, 0, 1] = np.nan               # mask one Zxy
    zr, er = rotate_tensor(z, e, -90.0)
    # -90 deg: Z'yx = -Zxy, Z'xy = -Zyx, Z'xx = Zyy
    assert np.isnan(zr[2, 1, 0])
    assert np.isfinite(zr[2, 0, 1]) and np.isfinite(zr[2, 0, 0])
    np.testing.assert_allclose(zr[:, 0, 1], -z[:, 1, 0])
    np.testing.assert_allclose(zr[:, 0, 0], z[:, 1, 1])
    np.testing.assert_allclose(er[:, 1, 0], e[:, 0, 1])


def test_general_rotation_invariants_and_masking():
    z, e = _rand()
    z[1, 1, 1] = np.nan
    zr, _ = rotate_tensor(z, e, 37.0)
    assert np.all(np.isnan(zr[1]))               # whole freq masked
    ok = np.all(np.isfinite(z), axis=(1, 2))
    np.testing.assert_allclose(np.linalg.det(zr[ok]),
                               np.linalg.det(z[ok]), rtol=1e-10)
    back, _ = rotate_tensor(zr, None, -37.0)
    np.testing.assert_allclose(back[ok], z[ok], atol=1e-12)


def test_to_frame_roundtrip_on_real_edi(edi_imp_file: Path):
    from pycsamt.seg.edi import EDIFile

    ed = EDIFile(edi_imp_file, verbose=0)
    a0 = frame_azimuth(ed)
    z0 = ed.Z.z.copy()
    theta = to_frame(ed, 30.0)
    assert theta == pytest.approx(30.0 - a0)
    assert frame_azimuth(ed) == pytest.approx(30.0)
    hx = [h for h in ed.get_section("definemeas").hmeas
          if str(h.chtype).upper() == "HX"][0]
    assert float(hx.azm) == pytest.approx(0.0)
    to_frame(ed, a0)
    ok = np.all(np.isfinite(z0), axis=(1, 2))
    np.testing.assert_allclose(ed.Z.z[ok], z0[ok], rtol=1e-9,
                               atol=1e-12)
