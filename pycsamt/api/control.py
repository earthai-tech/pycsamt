"""Global plotting-view controls for pyCSAMT."""

from __future__ import annotations

import copy
from collections.abc import Generator
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Any

import numpy as np

_PHASE_UNITS = {"degree", "radian"}
_RHO_VIEWS = {"log10", "linear", "rho"}
_X_VIEWS = {
    "log10_period",
    "period",
    "log10_frequency",
    "frequency",
}


@dataclass
class PhaseViewControl:
    """Control how phase values are wrapped and displayed."""

    range: tuple[float, float] = (-180.0, 180.0)
    unit: str = "degree"
    wrap: bool = True

    def transform(self, phase: Any) -> np.ndarray:
        """Return phase values in the configured display interval."""
        values = np.asarray(phase, dtype=float)
        unit = self.unit.lower()
        if unit not in _PHASE_UNITS:
            msg = f"phase unit must be one of {sorted(_PHASE_UNITS)}."
            raise ValueError(msg)
        if unit == "radian":
            values = np.rad2deg(values)
        if self.wrap and self.range is not None:
            values = wrap_phase(values, self.range)
        if unit == "radian":
            values = np.deg2rad(values)
        return values

    def label(self) -> str:
        """Return the phase-axis label."""
        if self.unit.lower() == "radian":
            return r"Phase (rad)"
        return r"Phase ($^\circ$)"


@dataclass
class RhoViewControl:
    """Control apparent-resistivity display scale."""

    view: str = "log10"

    def transform(self, rho: Any) -> np.ndarray:
        """Return apparent resistivity in the configured display scale."""
        values = np.asarray(rho, dtype=float)
        view = self.view.lower()
        if view not in _RHO_VIEWS:
            msg = f"rho view must be one of {sorted(_RHO_VIEWS)}."
            raise ValueError(msg)
        if view == "log10":
            return np.log10(values)
        return values

    def error(self, rho: Any, rho_err: Any | None) -> np.ndarray | None:
        """Return display-scale rho errors."""
        if rho_err is None:
            return None
        rho = np.asarray(rho, dtype=float)
        err = np.asarray(rho_err, dtype=float)
        if self.view.lower() == "log10":
            return err / (np.maximum(np.abs(rho), 1e-24) * np.log(10.0))
        return err

    def label(self) -> str:
        """Return the apparent-resistivity axis label."""
        if self.view.lower() == "log10":
            return r"$\log_{10}\rho_a$ ($\Omega\,\mathrm{m}$)"
        return r"$\rho_a$ ($\Omega\,\mathrm{m}$)"


@dataclass
class FrequencyAxisControl:
    """Control the 1-D frequency/period axis."""

    view: str = "log10_period"

    def transform(self, freq: Any) -> np.ndarray:
        """Return x-axis values from frequency in hertz."""
        freq = np.asarray(freq, dtype=float)
        view = self.view.lower()
        if view not in _X_VIEWS:
            msg = f"x view must be one of {sorted(_X_VIEWS)}."
            raise ValueError(msg)
        period = 1.0 / np.maximum(freq, 1e-24)
        if view == "log10_period":
            return np.log10(period)
        if view == "period":
            return period
        if view == "log10_frequency":
            return np.log10(freq)
        return freq

    def label(self) -> str:
        """Return the x-axis label."""
        labels = {
            "log10_period": r"$\log_{10}T$ (s)",
            "period": "Period (s)",
            "log10_frequency": r"$\log_{10}f$ (Hz)",
            "frequency": "Freq (Hz)",
        }
        return labels[self.view.lower()]

    def use_log_scale(self) -> bool:
        """Return whether Matplotlib should use log scaling."""
        return self.view.lower() in {"period", "frequency"}


class PyCSAMTControl:
    """Package-wide plotting-view control container."""

    def __init__(self) -> None:
        self._init_defaults()

    def _init_defaults(self) -> None:
        self.phase = PhaseViewControl()
        self.rho = RhoViewControl()
        self.x = FrequencyAxisControl()

    def configure(self, **kw: Any) -> None:
        """Configure controls using ``section__attribute`` paths."""
        for path, value in kw.items():
            parts = path.split("__")
            obj = self
            for part in parts[:-1]:
                obj = getattr(obj, part)
            setattr(obj, parts[-1], value)

    @contextmanager
    def context(self, **kw: Any) -> Generator[PyCSAMTControl, None, None]:
        """Temporarily override controls, then restore them."""
        snapshot = self._snapshot()
        try:
            if kw:
                self.configure(**kw)
            yield self
        finally:
            self._restore(snapshot)

    def reset(self) -> None:
        """Reset controls to package defaults."""
        self._init_defaults()

    def summary(self) -> str:
        """Return a human-readable summary."""
        lines = [
            "PyCSAMTControl",
            f"  phase.range = {self.phase.range!r}",
            f"  phase.unit  = {self.phase.unit!r}",
            f"  phase.wrap  = {self.phase.wrap}",
            f"  rho.view    = {self.rho.view!r}",
            f"  x.view      = {self.x.view!r}",
        ]
        return "\n".join(lines)

    def __repr__(self) -> str:  # noqa: D105
        return self.summary()

    def _snapshot(self) -> dict[str, Any]:
        return {
            "phase": copy.deepcopy(self.phase),
            "rho": copy.deepcopy(self.rho),
            "x": copy.deepcopy(self.x),
        }

    def _restore(self, snapshot: dict[str, Any]) -> None:
        self.phase = snapshot["phase"]
        self.rho = snapshot["rho"]
        self.x = snapshot["x"]


def wrap_phase(
    phase: Any,
    phase_range: tuple[float, float] = (-180.0, 180.0),
) -> np.ndarray:
    """Wrap phase values into ``phase_range``."""
    values = np.asarray(phase, dtype=float)
    lo, hi = float(phase_range[0]), float(phase_range[1])
    width = hi - lo
    if width <= 0:
        msg = "phase_range upper bound must be greater than lower bound."
        raise ValueError(msg)
    return (values - lo) % width + lo


PYCSAMT_CONTROL = PyCSAMTControl()


def configure_control(**kw: Any) -> None:
    """Configure :data:`PYCSAMT_CONTROL` with dotted-path arguments."""
    PYCSAMT_CONTROL.configure(**kw)


def reset_control() -> None:
    """Reset :data:`PYCSAMT_CONTROL` to package defaults."""
    PYCSAMT_CONTROL.reset()


__all__ = [
    "FrequencyAxisControl",
    "PYCSAMT_CONTROL",
    "PhaseViewControl",
    "PyCSAMTControl",
    "RhoViewControl",
    "configure_control",
    "reset_control",
    "wrap_phase",
]
