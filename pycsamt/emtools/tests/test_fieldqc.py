# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for :mod:`pycsamt.emtools.fieldqc` (synthetic profiles)."""

from __future__ import annotations

from types import SimpleNamespace as NS

import numpy as np
import pytest

from pycsamt.emtools.fieldqc import (
    circular_mean_deg,
    correct_polarity,
    detect_band_break,
    detect_polarity,
    estimate_static_shift_neighbours,
    find_duplicate_stations,
    hampel_mask,
    point_qc_table,
)

F = np.logspace(np.log10(800.0), np.log10(1.0), 25)


def _z1d(rho=100.0, phase=45.0, shift=0.0):
    """1-D tensor with |Z| from rho (field units) and fixed phase."""
    mag = np.sqrt(rho * 10 ** shift * F / 0.2)
    z = np.zeros((F.size, 2, 2), complex)
    z[:, 0, 1] = mag * np.exp(1j * np.radians(phase))
    z[:, 1, 0] = -z[:, 0, 1]
    return z


def _site(name, z, rel_err=0.05):
    return NS(station=name, path=None,
              Z=NS(freq=F, z=z, z_err=np.abs(z) * rel_err))


def test_circular_mean():
    m, r = circular_mean_deg(np.array([10.0, 20.0, np.nan]))
    assert m == pytest.approx(15.0) and r == pytest.approx(
        np.cos(np.radians(5)))
    assert np.isnan(circular_mean_deg(np.array([np.nan]))[0])


@pytest.mark.parametrize("fxy, fyx", [(False, False), (True, False),
                                      (False, True), (True, True)])
def test_detect_and_correct_polarity(fxy, fyx):
    z = _z1d()
    bad = z.copy()
    if fxy:
        bad[:, 0, 1] *= -1
    if fyx:
        bad[:, 1, 0] *= -1
    d = detect_polarity(bad, F)
    assert (d["flip_xy"], d["flip_yx"]) == (fxy, fyx)
    assert d["trusted_xy"] and not d["anomalous_xy"]
    for assume in ("E", "H"):
        fixed = correct_polarity(bad, flip_xy=fxy, flip_yx=fyx,
                                 assume=assume)
        np.testing.assert_allclose(fixed[:, 0, 1], z[:, 0, 1])
        np.testing.assert_allclose(fixed[:, 1, 0], z[:, 1, 0])


def test_correct_polarity_row_vs_column():
    z = np.arange(1, 5, dtype=complex).reshape(1, 2, 2)
    e = correct_polarity(z, flip_xy=True, flip_yx=False, assume="E")
    h = correct_polarity(z, flip_xy=True, flip_yx=False, assume="H")
    np.testing.assert_array_equal(e[0], [[-1, -2], [3, 4]])
    np.testing.assert_array_equal(h[0], [[1, -2], [3, -4]])
    with pytest.raises(ValueError):
        correct_polarity(z, flip_xy=True, flip_yx=False, assume="X")


def test_detect_polarity_anomalous_and_empty_band():
    z = _z1d(phase=135.0)                # Q2: no sign flip fixes it
    assert detect_polarity(z, F)["anomalous_xy"]
    with pytest.raises(ValueError, match="band"):
        detect_polarity(z, F, band=(1e4, 1e5))


def test_find_duplicates():
    a = _site("A", _z1d())
    b = _site("B", _z1d())
    c = _site("C", _z1d(rho=10.0))
    assert find_duplicate_stations([a, b, c]) == [["A", "B"]]
    assert find_duplicate_stations([a, c]) == []


def test_hampel_mask():
    x = np.zeros(11)
    x[5] = 5.0
    m = hampel_mask(x, floor=0.5)
    assert m[5] and m.sum() == 1


def _profile(n=9, shifts=None, flip_every=3, noisy_band=None):
    rng = np.random.default_rng(0)
    sites = []
    for k in range(n):
        z = _z1d(shift=0.0 if shifts is None else shifts[k])
        z *= np.exp(1j * np.radians(rng.normal(0, 0.5, F.size)))[
            :, None, None]
        if flip_every and k % flip_every == 0:
            z[:, 0, :] *= -1                       # reversed Ex
        if noisy_band is not None:
            sel = (F >= noisy_band[0]) & (F <= noisy_band[1])
            z[sel] *= np.exp(1j * rng.uniform(-np.pi, np.pi,
                                              sel.sum()))[:, None, None]
        sites.append(_site(f"S{k}", z))
    return sites


def test_point_qc_table_clean_profile_has_no_flags():
    t = point_qc_table(_profile())
    assert not t.bad_xy.any() and not t.bad_yx.any()
    np.testing.assert_allclose(t.phase_xy, 45.0, atol=2.0)


def test_band_break_detected():
    t = point_qc_table(_profile(noisy_band=(100.0, 800.0)),
                       band=(1.0, 50.0))
    br = detect_band_break(t, min_jump=0.3)
    assert len(br) >= 1
    top = br.iloc[0]
    assert top.f_low < 100.0 <= top.f_high
    assert top.jump > 0.3


def test_static_shift_neighbours():
    shifts = np.zeros(9)
    shifts[4] = 0.5                       # station 4: +0.5 decade
    shifts[6] = 1.6                       # beyond max_abs
    t = point_qc_table(_profile(shifts=shifts, flip_every=0))
    s = estimate_static_shift_neighbours(t, band=(10.0, 800.0))
    s = s.set_index("station")
    assert s.loc["S4", "shift_xy"] == pytest.approx(0.5, abs=0.02)
    assert s.loc["S4", "apply_xy"]
    assert not s.loc["S6", "apply_xy"]    # 1.6 > max_abs = 1
    assert not s.loc["S0", "apply_xy"]    # ~0 < min_abs
