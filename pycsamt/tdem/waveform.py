# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Transmitter waveform models for TEM systems."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

__all__ = ["SquareWaveform", "RampWaveform", "HalfSineWaveform", "CustomWaveform"]


class _WaveformBase:
    """Shared interface for all waveform types."""

    def current_at(self, t: np.ndarray) -> np.ndarray:
        r"""
        Return the normalised transmitter current (0–1) at times *t*.

        Parameters
        ----------
        t : array-like
            Times in seconds relative to the switch-off event at t = 0.
            Negative values are on-time; positive values are off-time.

        Returns
        -------
        np.ndarray
            Current envelope (dimensionless, in [0, 1]).
        """
        raise NotImplementedError

    @property
    def base_frequency(self) -> float:
        raise NotImplementedError

    @property
    def period(self) -> float:
        return 1.0 / self.base_frequency

    @property
    def half_period(self) -> float:
        return 0.5 / self.base_frequency


@dataclass
class SquareWaveform(_WaveformBase):
    r"""
    Ideal square-wave transmitter current (zero ramp time).

    The waveform alternates between +I and −I at ``base_frequency``
    Hz.  For late-time TEM processing the waveform is treated as an
    ideal step turn-off, so this is the appropriate choice when ramp
    effects are negligible.

    Parameters
    ----------
    base_frequency : float
        Fundamental repetition rate in Hz (e.g. 25 Hz for ZongeGDP).
    duty_cycle : float, optional
        Fraction of the half-period spent at full current.  Default 0.5
        (50 % duty cycle → symmetric square wave).

    Examples
    --------
    >>> from pycsamt.tdem.waveform import SquareWaveform
    >>> wf = SquareWaveform(base_frequency=25.0)
    >>> wf.half_period
    0.02
    """

    base_frequency: float = 25.0
    duty_cycle: float = 0.5

    def current_at(self, t: np.ndarray) -> np.ndarray:
        t = np.asarray(t, float)
        hp = self.half_period
        # normalised position within the positive half-period
        t_mod = np.mod(t, hp)
        on_time = hp * self.duty_cycle
        return np.where(t_mod < on_time, 1.0, 0.0)


@dataclass
class RampWaveform(_WaveformBase):
    r"""
    Transmitter current with a finite linear ramp on switch-off.

    The current decays linearly from 1 to 0 over ``ramp_off`` seconds
    starting at t = 0 (the nominal switch-off instant).  This is the
    waveform required for accurate early-time TEM processing where the
    ramp duration is comparable to the first measurement gates.

    Parameters
    ----------
    base_frequency : float
        Fundamental repetition rate in Hz.
    ramp_off : float
        Duration of the current turn-off ramp in seconds.
    ramp_on : float, optional
        Duration of the current turn-on ramp in seconds.  Rarely
        needed; default 0.0 (ideal step turn-on).
    duty_cycle : float, optional
        Fraction of the half-period at full current (before ramp).
        Default 0.5.

    Examples
    --------
    >>> from pycsamt.tdem.waveform import RampWaveform
    >>> wf = RampWaveform(base_frequency=25.0, ramp_off=1e-4)
    >>> wf.ramp_off
    0.0001
    """

    base_frequency: float = 25.0
    ramp_off: float = 1e-4
    ramp_on: float = 0.0
    duty_cycle: float = 0.5

    def current_at(self, t: np.ndarray) -> np.ndarray:
        t = np.asarray(t, float)
        I = np.ones_like(t)
        # ramp-off region: 0 ≤ t < ramp_off
        mask_ramp = (t >= 0.0) & (t < self.ramp_off)
        I[mask_ramp] = 1.0 - t[mask_ramp] / self.ramp_off
        # post-ramp: fully off
        I[t >= self.ramp_off] = 0.0
        # negative t (on-time): full current
        I[t < 0.0] = 1.0
        return I


@dataclass
class HalfSineWaveform(_WaveformBase):
    r"""
    Half-sine transmitter current (used by some CSEM / airborne systems).

    Parameters
    ----------
    base_frequency : float
        Fundamental repetition rate in Hz.

    Examples
    --------
    >>> from pycsamt.tdem.waveform import HalfSineWaveform
    >>> wf = HalfSineWaveform(base_frequency=30.0)
    >>> round(wf.half_period, 6)
    0.016667
    """

    base_frequency: float = 30.0

    def current_at(self, t: np.ndarray) -> np.ndarray:
        t = np.asarray(t, float)
        hp = self.half_period
        t_mod = np.mod(t, hp)
        return np.maximum(0.0, np.sin(np.pi * t_mod / hp))


class CustomWaveform(_WaveformBase):
    r"""
    User-defined waveform supplied as paired time/current arrays.

    Parameters
    ----------
    t_waveform : array-like of float
        Time samples in seconds (relative to switch-off at t = 0).
    i_waveform : array-like of float
        Normalised current values (0–1) at each ``t_waveform`` sample.
    base_frequency : float
        Fundamental repetition rate in Hz.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.tdem.waveform import CustomWaveform
    >>> t = np.linspace(-0.02, 0.02, 200)
    >>> I = np.where(t < 0, 1.0, 0.0)
    >>> wf = CustomWaveform(t, I, base_frequency=25.0)
    """

    def __init__(
        self,
        t_waveform,
        i_waveform,
        *,
        base_frequency: float = 25.0,
    ) -> None:
        self._t = np.asarray(t_waveform, float)
        self._I = np.asarray(i_waveform, float)
        self._base_frequency = float(base_frequency)
        if self._t.shape != self._I.shape or self._t.ndim != 1:
            raise ValueError("t_waveform and i_waveform must be 1-D arrays of equal length")

    @property
    def base_frequency(self) -> float:
        return self._base_frequency

    def current_at(self, t: np.ndarray) -> np.ndarray:
        return np.interp(np.asarray(t, float), self._t, self._I,
                         left=self._I[0], right=self._I[-1])
