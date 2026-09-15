"""Tests for package-wide plotting-view controls."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.api.control import (
    PYCSAMT_CONTROL,
    FrequencyAxisControl,
    PhaseViewControl,
    PyCSAMTControl,
    RhoViewControl,
    configure_control,
    reset_control,
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


# ─────────────────────────────────────────────────────────────────────────
# PhaseViewControl
# ─────────────────────────────────────────────────────────────────────────


def test_phase_view_control_rejects_bad_unit():
    with pytest.raises(ValueError, match="phase unit must be one of"):
        PhaseViewControl(unit="bogus").transform([0.0])


def test_phase_view_control_radian_round_trip():
    phase = PhaseViewControl(unit="radian", range=(-180.0, 180.0))
    out = phase.transform([np.pi + 0.1])  # slightly over pi -> wraps
    assert out[0] < np.pi
    assert phase.label() == r"Phase (rad)"


def test_phase_view_control_no_wrap_skips_wrapping():
    phase = PhaseViewControl(wrap=False)
    out = phase.transform([190.0])
    assert out[0] == 190.0


def test_phase_view_control_degree_label():
    assert PhaseViewControl(unit="degree").label() == r"Phase ($^\circ$)"


# ─────────────────────────────────────────────────────────────────────────
# RhoViewControl
# ─────────────────────────────────────────────────────────────────────────


def test_rho_view_control_rejects_bad_view():
    with pytest.raises(ValueError, match="rho view must be one of"):
        RhoViewControl(view="bogus").transform([1.0])


def test_rho_view_control_linear_view_passthrough():
    rho = RhoViewControl(view="linear")
    assert np.allclose(rho.transform([10.0, 100.0]), [10.0, 100.0])
    assert rho.label() == r"$\rho_a$ ($\Omega\,\mathrm{m}$)"


def test_rho_view_control_error_none_when_rho_err_none():
    rho = RhoViewControl(view="log10")
    assert rho.error([10.0], None) is None


def test_rho_view_control_error_log10_and_linear():
    rho_log = RhoViewControl(view="log10")
    err_log = rho_log.error([10.0], [1.0])
    assert err_log[0] > 0

    rho_linear = RhoViewControl(view="linear")
    err_linear = rho_linear.error([10.0], [1.0])
    assert err_linear[0] == 1.0


# ─────────────────────────────────────────────────────────────────────────
# FrequencyAxisControl
# ─────────────────────────────────────────────────────────────────────────


def test_frequency_axis_control_rejects_bad_view():
    with pytest.raises(ValueError, match="x view must be one of"):
        FrequencyAxisControl(view="bogus").transform([1.0])


def test_frequency_axis_control_every_view_and_label():
    freq = [1.0, 10.0]
    period_view = FrequencyAxisControl(view="period")
    assert np.allclose(period_view.transform(freq), [1.0, 0.1])
    assert period_view.label() == "Period (s)"
    assert period_view.use_log_scale() is True

    logfreq_view = FrequencyAxisControl(view="log10_frequency")
    assert np.allclose(logfreq_view.transform(freq), np.log10(freq))
    assert logfreq_view.label() == r"$\log_{10}f$ (Hz)"

    freq_view = FrequencyAxisControl(view="frequency")
    assert np.allclose(freq_view.transform(freq), freq)
    assert freq_view.label() == "Freq (Hz)"
    assert freq_view.use_log_scale() is True

    logperiod_view = FrequencyAxisControl(view="log10_period")
    assert logperiod_view.use_log_scale() is False


# ─────────────────────────────────────────────────────────────────────────
# PyCSAMTControl container
# ─────────────────────────────────────────────────────────────────────────


def test_pycsamt_control_reset_restores_defaults():
    control = PyCSAMTControl()
    control.configure(rho__view="linear")
    control.reset()
    assert control.rho.view == "log10"


def test_pycsamt_control_summary_and_repr():
    control = PyCSAMTControl()
    summary = control.summary()
    assert "PyCSAMTControl" in summary
    assert "phase.range" in summary
    assert repr(control) == summary


def test_pycsamt_control_context_without_kwargs_is_a_noop_block():
    control = PyCSAMTControl()
    original = control.rho.view
    with control.context():
        assert control.rho.view == original
    assert control.rho.view == original


# ─────────────────────────────────────────────────────────────────────────
# wrap_phase / module-level convenience functions
# ─────────────────────────────────────────────────────────────────────────


def test_wrap_phase_rejects_zero_width_range():
    with pytest.raises(ValueError, match="upper bound must be greater"):
        wrap_phase([0.0], (10.0, 10.0))


def test_module_level_configure_and_reset_control():
    configure_control(rho__view="linear")
    assert PYCSAMT_CONTROL.rho.view == "linear"
    reset_control()
    assert PYCSAMT_CONTROL.rho.view == "log10"
