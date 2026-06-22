"""Tests for package-wide plotting-view controls."""

from __future__ import annotations

import numpy as np

from pycsamt.api.control import (
    PYCSAMT_CONTROL,
    FrequencyAxisControl,
    PhaseViewControl,
    RhoViewControl,
    wrap_phase,
)


def test_wrap_phase_default_interval():
    """Phase wrapping should use the configured interval."""
    out = wrap_phase(np.array([-190.0, 0.0, 190.0]), (-180.0, 180.0))

    assert np.allclose(out, [170.0, 0.0, -170.0])


def test_control_context_restores_defaults():
    """Temporary control overrides should restore on exit."""
    old_range = PYCSAMT_CONTROL.phase.range

    with PYCSAMT_CONTROL.context(
        phase__range=(-90.0, 90.0),
        rho__view="linear",
        x__view="frequency",
    ):
        assert PYCSAMT_CONTROL.phase.range == (-90.0, 90.0)
        assert PYCSAMT_CONTROL.rho.view == "linear"
        assert PYCSAMT_CONTROL.x.view == "frequency"

    assert PYCSAMT_CONTROL.phase.range == old_range
    assert PYCSAMT_CONTROL.rho.view == "log10"
    assert PYCSAMT_CONTROL.x.view == "log10_period"


def test_view_transforms_and_labels():
    """View controls should transform rho and x consistently."""
    rho = RhoViewControl(view="log10")
    x = FrequencyAxisControl(view="log10_period")
    phase = PhaseViewControl(range=(-180.0, 180.0))

    assert np.allclose(rho.transform([10.0, 100.0]), [1.0, 2.0])
    assert rho.label() == r"$\log_{10}\rho_a$ ($\Omega\,\mathrm{m}$)"
    assert np.allclose(x.transform([1.0, 0.01]), [0.0, 2.0])
    assert x.label() == r"$\log_{10}T$ (s)"
    assert np.allclose(phase.transform([190.0]), [-170.0])
