# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Regression guard for the ModEM apparent-resistivity unit convention.

ModEM impedance data are written in **field units** ``[mV/km]/[nT]``
(``ModEmConfig.units`` default), so the apparent resistivity is
``ρ_a = 0.2·|Z|²·T``. The response plot once used the SI form
``|Z|²/(μ₀·ω)`` directly on this field-unit Z, over-estimating ρ_a — and the
propagated error bars — by ``1/(0.2·2π·μ₀)`` ≈ 6.3·10⁵ (ρ_a shooting to
~10⁹ Ω·m on data whose true ρ_a is ~10²).
"""
from __future__ import annotations

import numpy as np

from pycsamt.models.modem.plot import _rho_phase_from_rows


def test_rho_uses_field_formula():
    # One synthetic field-unit impedance row: |Z| ≈ 3945 [mV/km]/[nT] at
    # T = 9.615e-5 s reproduces the real 23-18-002U ZXY (~299 Ω·m).
    T = 9.615384615384615e-05
    re, im = 2786.227, 2792.402
    err = 50.0
    p, rho, drho, phi, dphi = _rho_phase_from_rows([(T, re, im, err)])

    z2 = re ** 2 + im ** 2
    expected_rho = 0.2 * z2 * T
    np.testing.assert_allclose(rho[0], expected_rho, rtol=1e-9)

    # physical magnitude, NOT the 10⁹ produced by the SI-on-field-Z bug
    assert 50.0 < rho[0] < 1.0e4, rho[0]

    # error bar: 0.4·|Z|·σ_Z·T, a fraction of ρ_a for a small σ_Z
    expected_drho = 0.4 * np.sqrt(z2) * err * T
    np.testing.assert_allclose(drho[0], expected_drho, rtol=1e-9)
    assert drho[0] < rho[0]  # σ_Z=50 « |Z|≈3945 → sub-ρ error bar


def test_rho_not_si_formula():
    # Guard against silently reverting to |Z|²/(μ₀·ω).
    T = 1.0e-3
    re, im = 1000.0, 1000.0
    p, rho, drho, phi, dphi = _rho_phase_from_rows([(T, re, im, 10.0)])
    z2 = re ** 2 + im ** 2
    si_rho = z2 / ((4.0e-7 * np.pi) * (2.0 * np.pi / T))
    assert rho[0] < si_rho / 1.0e4  # field result is ~6.3e5× smaller
