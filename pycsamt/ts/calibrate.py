# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
Sensor (coil) transfer-function correction for time series.

Induction coils and some electric-field preamplifiers record a
voltage :math:`v(t)` that is a filtered version of the physical
field :math:`b(t)`. In the frequency domain

.. math::

    V(f) = S(f)\,B(f), \qquad S(f) = A(f)\,e^{i\varphi(f)},

where :math:`A` is the sensitivity (output units per input unit,
e.g. mV/nT) and :math:`\varphi` the phase lead of the output (deg in
calibration tables). :func:`apply_sensor_response` estimates
:math:`b(t)` by spectral division

.. math::

    \hat B(f) = \frac{V(f)}{S_w(f)}, \qquad
    S_w = \begin{cases} S & |S| \ge w\,\max|S| \\
          w\,\max|S|\, S/|S| & \text{otherwise,} \end{cases}

with a water level :math:`w` that bounds the gain where the sensor
barely responds (near DC for coils). The zero-frequency term is set
to zero (coils do not pass DC). The convention matches
:func:`numpy.fft.rfft` (forward kernel :math:`e^{-2\pi i f t}`):
for :math:`b = \cos(2\pi f t)`, the output is
:math:`A\cos(2\pi f t + \varphi)`.

Calibration tables are interpolated with :math:`\log A` linear in
:math:`\log f` and :math:`\varphi` linear in :math:`\log f`; outside
the table the edge values are held (``extrapolate="edge"``) or an
error is raised (``"raise"``).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ._base import TSData

__all__ = ["SensorResponse", "read_sensor_response",
           "apply_sensor_response"]


@dataclass
class SensorResponse:
    """Tabulated sensor transfer function.

    Parameters
    ----------
    freq : array_like
        Frequencies (Hz), strictly positive.
    amplitude : array_like
        Sensitivity :math:`A(f)` (output / input units, > 0).
    phase_deg : array_like
        Phase lead :math:`\\varphi(f)` of the output (degrees).
    input_units, output_units : str
        Units of the physical field and of the recorded signal,
        e.g. ``"nT"`` and ``"mV"``.
    name : str, optional
        Sensor identifier (serial number).

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.ts.calibrate import SensorResponse
    >>> s = SensorResponse([1, 10, 100], [0.1, 1.0, 1.0],
    ...                    [90.0, 45.0, 0.0])
    >>> np.round(abs(s.at(np.array([10.0]))), 3)
    array([1.])
    """

    freq: np.ndarray
    amplitude: np.ndarray
    phase_deg: np.ndarray
    input_units: str = "nT"
    output_units: str = "mV"
    name: str = ""

    def __post_init__(self) -> None:
        f = np.asarray(self.freq, float).ravel()
        a = np.asarray(self.amplitude, float).ravel()
        ph = np.asarray(self.phase_deg, float).ravel()
        if not (f.size == a.size == ph.size) or f.size < 2:
            raise ValueError(
                "freq, amplitude and phase_deg need the same length "
                "(>= 2).")
        if np.any(~np.isfinite(f) | (f <= 0)):
            raise ValueError("frequencies must be finite and > 0.")
        if np.any(~np.isfinite(a) | (a <= 0)):
            raise ValueError("amplitudes must be finite and > 0.")
        o = np.argsort(f)
        if np.any(np.diff(f[o]) == 0):
            raise ValueError("duplicate calibration frequencies.")
        self.freq, self.amplitude, self.phase_deg = f[o], a[o], ph[o]

    def at(self, f: np.ndarray, extrapolate: str = "edge"
           ) -> np.ndarray:
        """Complex response :math:`S(f)` at frequencies *f* (> 0).

        Parameters
        ----------
        f : ndarray
            Frequencies (Hz), > 0.
        extrapolate : {"edge", "raise"}, default "edge"
            Behaviour outside the tabulated band.

        Returns
        -------
        ndarray of complex
        """
        f = np.asarray(f, float)
        if np.any(f <= 0):
            raise ValueError("frequencies must be > 0.")
        lo, hi = self.freq[0], self.freq[-1]
        if extrapolate == "raise" and (f.min() < lo or f.max() > hi):
            raise ValueError(
                f"frequencies outside the calibration band "
                f"[{lo:g}, {hi:g}] Hz.")
        if extrapolate not in ("edge", "raise"):
            raise ValueError("extrapolate must be 'edge' or 'raise'.")
        lf = np.log10(np.clip(f, lo, hi))
        lfreq = np.log10(self.freq)
        amp = 10 ** np.interp(lf, lfreq, np.log10(self.amplitude))
        ph = np.interp(lf, lfreq, self.phase_deg)
        return amp * np.exp(1j * np.radians(ph))


def read_sensor_response(
    path: str | Path,
    *,
    columns: tuple[int, int, int] = (0, 1, 2),
    delimiter: str | None = None,
    comments: str = "#",
    input_units: str = "nT",
    output_units: str = "mV",
) -> SensorResponse:
    """Read a ``freq amplitude phase_deg`` calibration table.

    Lines that do not parse as numbers (headers, comments) are
    skipped.

    Parameters
    ----------
    path : str or Path
        Text table, whitespace- or *delimiter*-separated.
    columns : tuple of int, default (0, 1, 2)
        Column indices of frequency (Hz), amplitude and phase (deg).
    delimiter : str, optional
        Column separator (default: any whitespace; use ``","`` for
        CSV).
    comments : str, default "#"
        Comment marker.
    input_units, output_units : str
        Stored on the returned :class:`SensorResponse`.

    Returns
    -------
    SensorResponse

    Raises
    ------
    ValueError
        If fewer than two numeric rows are found.
    """
    rows = []
    for ln in Path(path).read_text(errors="replace").splitlines():
        ln = ln.split(comments, 1)[0].strip()
        if not ln:
            continue
        tok = ln.split(delimiter) if delimiter else ln.split()
        try:
            rows.append([float(tok[c]) for c in columns])
        except (ValueError, IndexError):
            continue
    if len(rows) < 2:
        raise ValueError(f"{path}: fewer than two calibration rows.")
    a = np.asarray(rows)
    return SensorResponse(a[:, 0], a[:, 1], a[:, 2],
                          input_units=input_units,
                          output_units=output_units,
                          name=Path(path).stem)


def apply_sensor_response(
    ts: TSData,
    responses: dict[str, SensorResponse],
    *,
    water_level: float = 1e-3,
    extrapolate: str = "edge",
    inplace: bool = False,
) -> TSData:
    """Remove sensor transfer functions from time-series channels.

    Parameters
    ----------
    ts : TSData
        Time series with a sampling interval ``dt``.
    responses : dict of str to SensorResponse
        Channel id (``"HX"`` ...) to its sensor response.
    water_level : float, default 1e-3
        Relative water level :math:`w` (see module notes); must be
        in (0, 1).
    extrapolate : {"edge", "raise"}, default "edge"
        Passed to :meth:`SensorResponse.at`.
    inplace : bool, default False
        Modify *ts* instead of a copy.

    Returns
    -------
    TSData
        Corrected series. ``units[ch]`` becomes the response's
        ``input_units``; ``meta["calibration"][ch]`` records the
        sensor name and water level.

    Raises
    ------
    ValueError
        If ``dt`` is missing, a channel is absent, contains NaN
        (fill gaps first, e.g. :func:`pycsamt.ts.preprocess`), or
        ``water_level`` is out of range.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.ts import TSData
    >>> from pycsamt.ts.calibrate import (SensorResponse,
    ...                                   apply_sensor_response)
    >>> t = np.arange(4096) / 256.0
    >>> ts = TSData({"HX": 2.0 * np.cos(2 * np.pi * 8 * t)},
    ...             dt=1 / 256.0)
    >>> flat = SensorResponse([0.1, 1000.0], [2.0, 2.0], [0.0, 0.0])
    >>> out = apply_sensor_response(ts, {"HX": flat})
    >>> round(float(np.abs(out.get("HX")).max()), 3)
    1.0
    """
    if ts.dt is None or ts.dt <= 0:
        raise ValueError("ts.dt (sampling interval) is required.")
    if not 0 < water_level < 1:
        raise ValueError("water_level must be in (0, 1).")
    out = ts if inplace else _copy_ts(ts)
    cal = dict(out.meta.get("calibration", {}) or {})
    for cid, resp in responses.items():
        key = cid.upper()
        if key not in out.data:
            raise ValueError(f"channel {key!r} not in time series.")
        x = out.data[key]
        if np.isnan(x).any():
            raise ValueError(
                f"channel {key} contains NaN; fill gaps first.")
        n = x.size
        spec = np.fft.rfft(x)
        f = np.fft.rfftfreq(n, d=out.dt)
        s = np.zeros(f.size, dtype=complex)
        s[1:] = resp.at(f[1:], extrapolate=extrapolate)
        mag = np.abs(s[1:])
        floor = water_level * mag.max()
        small = mag < floor
        s[1:][small] = floor * s[1:][small] / mag[small]
        spec[0] = 0.0
        spec[1:] /= s[1:]
        out.data[key] = np.fft.irfft(spec, n=n)
        out.units[key] = resp.input_units
        cal[key] = {"sensor": resp.name, "water_level": water_level}
    out.meta["calibration"] = cal
    out.meta["calibrated"] = True
    return out


def _copy_ts(ts: TSData) -> TSData:
    new = ts.copy_meta()
    for cid in ts.channels():
        new.add_channel(cid, ts.data[cid].copy())
    return new
