# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Physics contract tests for :mod:`pycsamt.forward.em1d`.

The MT 1-D solver must honor the classical analytic limits:
a uniform half-space returns its own resistivity at every frequency
with a 45 degree phase, and layered models approach the top/bottom
layer resistivities at the high/low frequency asymptotes.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.forward.em1d import MT1DForward, TEM1DForward
from pycsamt.forward.synthetic import LayeredModel


def test_halfspace_recovers_resistivity_and_45deg_phase():
    freqs = np.logspace(-3, 4, 40)
    model = LayeredModel(resistivity=[100.0], thickness=[])
    resp = MT1DForward(freqs).run(model)

    assert resp.method == "MT1D"
    assert resp.rho_a.shape == freqs.shape
    assert np.allclose(resp.rho_a, 100.0, rtol=1e-8)
    assert np.allclose(resp.phase, 45.0, atol=1e-8)


def test_two_layer_asymptotic_limits():
    freqs = np.logspace(-5, 5, 60)
    model = LayeredModel(resistivity=[100.0, 10.0], thickness=[500.0])
    resp = MT1DForward(freqs).run(model)

    # high frequency -> shallow layer, low frequency -> halfspace
    assert resp.rho_a[-1] == pytest.approx(100.0, rel=0.05)
    assert resp.rho_a[0] == pytest.approx(10.0, rel=0.05)
    # phase stays within physical bounds for normal models
    assert np.all(resp.phase > 0.0)
    assert np.all(resp.phase < 90.0)


def test_rho_a_consistent_with_impedance():
    freqs = np.logspace(-2, 3, 20)
    model = LayeredModel(resistivity=[100.0, 10.0, 500.0], thickness=[300.0, 800.0])
    resp = MT1DForward(freqs).run(model)

    mu0 = 4e-7 * np.pi
    omega = 2 * np.pi * freqs
    assert np.allclose(resp.rho_a, np.abs(resp.z) ** 2 / (omega * mu0), rtol=1e-12)
    assert np.allclose(resp.phase, np.angle(resp.z, deg=True))


def test_conductor_lowers_apparent_resistivity():
    freqs = np.logspace(-3, 3, 30)
    background = LayeredModel(resistivity=[100.0], thickness=[])
    with_conductor = LayeredModel(
        resistivity=[100.0, 1.0, 100.0], thickness=[200.0, 300.0]
    )
    rho_bg = MT1DForward(freqs).run(background).rho_a
    rho_c = MT1DForward(freqs).run(with_conductor).rho_a
    # somewhere in the band the buried conductor must pull rho_a down
    assert rho_c.min() < 0.5 * rho_bg.min()


def test_forward_response_to_array_layout():
    freqs = np.logspace(-1, 2, 8)
    model = LayeredModel(resistivity=[50.0, 5.0], thickness=[250.0])
    resp = MT1DForward(freqs).run(model)

    vec = resp.to_array()
    assert vec.shape == (2 * freqs.size,)
    assert np.allclose(vec[: freqs.size], np.log10(resp.rho_a))
    assert np.allclose(vec[freqs.size :], resp.phase)

    vec_no_phase = resp.to_array(include_phase=False)
    assert vec_no_phase.shape == freqs.shape

    vec_linear = resp.to_array(log_rho=False)
    assert np.allclose(vec_linear[: freqs.size], resp.rho_a)


def test_layered_model_depth_and_validation():
    m = LayeredModel(resistivity=[100, 10, 500], thickness=[300, 800])
    assert m.n_layers == 3
    assert np.allclose(m.depth, [0.0, 300.0, 1100.0])
    with pytest.raises(Exception):
        LayeredModel(resistivity=[100, 10], thickness=[100, 200, 300])


def test_tem1d_matches_empymod_reference():
    """TEM1DForward's central-loop step-off response, cross-checked
    against a directly-constructed :func:`empymod.bipole` call using the
    same loop-as-scaled-dipole-segment trick empymod's own gallery uses
    to reproduce Ward & Hohmann (1988)'s central-loop figures (4.7-4.8).

    Exercises the moment -> strength conversion, the air-layer
    insertion, and the depth/resistivity array construction from a
    :class:`LayeredModel` -- the integration surface most likely to
    regress even though the transform itself is empymod's own.
    """
    import empymod

    radius = 25.0
    moment = 4.0
    times = np.logspace(-5, -3, 6)
    model = LayeredModel(resistivity=[60.0, 250.0, 900.0], thickness=[120.0, 700.0])

    resp = TEM1DForward(times=times, loop_radius=radius, moment=moment).run(model)

    strength = 2.0 * moment / radius
    ref = empymod.bipole(
        src=[radius, 0.0, 0.0, 90.0, 0.0],
        rec=[0.0, 0.0, 0.0, 0.0, 90.0],
        depth=[0.0, 120.0, 820.0],
        res=[2e14, 60.0, 250.0, 900.0],
        freqtime=times,
        signal=0,
        strength=strength,
        mrec=True,
        epermH=[0.0, 0.0, 0.0, 0.0],
        verb=0,
    ).real

    assert resp.method == "TEM1D"
    assert np.allclose(resp.dBz_dt, ref, rtol=1e-6)


def test_tem1d_halfspace_step_off_is_smooth_and_decaying():
    """A simple half-space's central-loop step-off response must decay
    smoothly and stay single-signed.

    An earlier, from-scratch Hankel/Fourier-transform implementation
    failed this even for its own default parameters: its frequency-domain
    kernel grew unboundedly at high frequency instead of decaying, and
    refining either quadrature's resolution changed the answer's *sign*
    at several times rather than converging. TEM1DForward now delegates
    to empymod's validated digital linear filters instead.
    """
    times = np.logspace(-6, -2, 20)
    model = LayeredModel(resistivity=[100.0], thickness=[])
    resp = TEM1DForward(times=times, loop_radius=50.0).run(model)

    assert np.all(np.isfinite(resp.dBz_dt))
    assert np.all(resp.dBz_dt > 0) or np.all(resp.dBz_dt < 0)
    assert np.all(np.diff(np.abs(resp.dBz_dt)) <= 0)
