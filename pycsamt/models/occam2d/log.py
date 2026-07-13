# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""OccamLog — parse Occam2D log files.

The log file (typically ``LogFile.logfile``) records per-iteration
convergence statistics written by the Fortran binary.

Entry point
-----------
``OccamLog.read(path)``
    Parse the log and populate per-iteration arrays.

Log-file structure (one block per iteration)
--------------------------------------------
Each iteration block looks like::

    ** ITERATION  N  **
    TOFMU: MISFIT = ...  AMU = ...   ← internal search steps (ignored)
    MINIMUM TOL FROM fminocc IS AT MU = <lagrange>
    AND IS = <achieved_rms>
    USING N EVALUATIONS OF FUNCTION
    [DIVERGENCE PROBLEMS, CUTTING STEP SIZE   ← may repeat]
    [more TOFMU / MINIMUM TOL / AND IS lines  ← keep the *last* values]
    STEPSIZE IS = <stepsize>
    ROUGHNESS IS = <roughness>

The "STARTING R.M.S." line before each "** ITERATION **" is the ending
RMS from the previous iteration and is used only for cross-validation.
Some iterations end with "INTERCEPT IS AT MU" instead of
"MINIMUM TOL FROM fminocc".  The last iteration may have neither
STEPSIZE nor ROUGHNESS if the run stopped with "Convergence problems".
"""

from __future__ import annotations

from pathlib import Path
from typing import Union

import numpy as np

from .base import OccamBase
from .doc import _occam_param_docs as _params

PathLike = Union[str, Path]

__all__ = ["OccamLog"]

_NAN = float("nan")


def _parse_float(line: str) -> float:
    """Return the float after the last equals sign."""
    try:
        return float(line.rsplit("=", 1)[1].strip())
    except (ValueError, IndexError):
        return _NAN


class OccamLog(OccamBase):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.iterations: np.ndarray = np.array([], dtype=int)
        self.rms: np.ndarray = np.array([])
        self.roughness: np.ndarray = np.array([])
        self.lagrange: np.ndarray = np.array([])
        self.stepsize: np.ndarray = np.array([])

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> OccamLog:
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"Occam log file not found: {p}")

        obj = cls(**kwargs)

        # Per-iteration accumulators
        iter_nums: list[int] = []
        rms_vals: list[float] = []
        rough_vals: list[float] = []
        step_vals: list[float] = []
        lagr_vals: list[float] = []

        # Mutable state for the iteration currently being parsed.
        # Using a one-element list so the inner helper can mutate them.
        cur = {
            "iter": None,
            "rms": _NAN,
            "roughness": _NAN,
            "stepsize": _NAN,
            "lagrange": _NAN,
        }

        def _save() -> None:
            """Flush current iteration state into the accumulators."""
            if cur["iter"] is None:
                return
            iter_nums.append(cur["iter"])
            rms_vals.append(cur["rms"])
            rough_vals.append(cur["roughness"])
            step_vals.append(cur["stepsize"])
            lagr_vals.append(cur["lagrange"])

        def _reset() -> None:
            cur["iter"] = None
            cur["rms"] = _NAN
            cur["roughness"] = _NAN
            cur["stepsize"] = _NAN
            cur["lagrange"] = _NAN

        with p.open("r", errors="replace") as fh:
            for raw in fh:
                line = raw.strip()
                if not line:
                    continue
                up = line.upper()

                if "** ITERATION" in up:
                    # Save completed iteration, start a new one
                    _save()
                    _reset()
                    # Extract the integer after "ITERATION"
                    tokens = line.replace("*", "").split()
                    for i, tok in enumerate(tokens):
                        if tok.upper() == "ITERATION" and i + 1 < len(tokens):
                            try:
                                cur["iter"] = int(tokens[i + 1])
                            except ValueError:
                                pass
                            break

                elif "MINIMUM TOL FROM" in up or "INTERCEPT IS AT MU" in up:
                    # Lagrange multiplier at the accepted step.
                    # May appear multiple times (divergence sub-rounds);
                    # keep the last value for this iteration.
                    v = _parse_float(line)
                    if not np.isnan(v):
                        cur["lagrange"] = v

                elif up.startswith("AND IS"):
                    # Achieved misfit.  Also appears multiple times;
                    # keep the last value (closest to ROUGHNESS IS).
                    v = _parse_float(line)
                    if not np.isnan(v):
                        cur["rms"] = v

                elif "STEPSIZE IS" in up:
                    v = _parse_float(line)
                    if not np.isnan(v):
                        cur["stepsize"] = v

                elif "ROUGHNESS IS" in up:
                    v = _parse_float(line)
                    if not np.isnan(v):
                        cur["roughness"] = v

        # Save the last iteration (no following "** ITERATION **" to trigger it)
        _save()

        obj.iterations = np.array(iter_nums, dtype=int)
        obj.rms = np.array(rms_vals)
        obj.roughness = np.array(rough_vals)
        obj.stepsize = np.array(step_vals)
        obj.lagrange = np.array(lagr_vals)

        if obj.verbose:
            obj.logger.info(
                "OccamLog.read: %d iterations parsed from %s", obj.n_iter, p
            )

        return obj

    # ------------------------------------------------------------------
    # Derived
    # ------------------------------------------------------------------
    @property
    def n_iter(self) -> int:
        """Number of parsed iterations."""
        return len(self.iterations)

    @property
    def converged(self) -> bool:
        """True if any iteration achieved RMS ≤ 1.0 (Occam target)."""
        if not self.rms.size:
            return False
        finite = self.rms[np.isfinite(self.rms)]
        return bool(finite.size and finite.min() <= 1.0)

    @property
    def best_iteration(self) -> int:
        """1-based iteration number with the lowest finite RMS misfit."""
        if not self.rms.size:
            return 0
        finite_mask = np.isfinite(self.rms)
        if not finite_mask.any():
            return 0
        idx = int(np.argmin(np.where(finite_mask, self.rms, np.inf)))
        return int(self.iterations[idx])

    def summary(self) -> str:
        """Return a one-paragraph convergence summary."""
        if not self.n_iter:
            return "OccamLog: no iterations parsed."
        return (
            f"OccamLog: {self.n_iter} iterations | "
            f"initial RMS {self.rms[0]:.4f} → "
            f"final RMS {self.rms[-1]:.4f} | "
            f"best iter {self.best_iteration} "
            f"(RMS {self.rms[self.best_iteration - 1]:.4f}) | "
            f"converged: {self.converged}"
        )


OccamLog.__doc__ = rf"""
Represent an Occam2D convergence log.

``OccamLog`` parses the text log written by the Occam2D
Fortran executable. The file records one block per inversion
iteration, including accepted misfit, model roughness,
Lagrange multiplier, and line-search step size. Parsed arrays
align by index, so ``iterations[i]``, ``rms[i]``,
``roughness[i]``, ``lagrange[i]``, and ``stepsize[i]``
describe the same iteration.

The main convergence statistic is the normalized RMS data
misfit. If :math:`r_i` are weighted residuals for :math:`N`
data, the reported quantity is commonly interpreted as

.. math::

    \phi_d =
    \sqrt{{\frac{{1}}{{N}}\sum_{{i=1}}^N r_i^2}} .

An inversion is usually considered well weighted when
:math:`\phi_d \approx 1`. The practical target still depends
on error estimates and modeling assumptions [1]_.

Parameters
----------
{_params.common.verbose}
{_params.common.logger}

Attributes
----------
iterations : ndarray of int, shape (n_iter,)
    One-based iteration numbers parsed from ``** ITERATION``
    blocks. The values are preserved as written by Occam.
rms : ndarray of float, shape (n_iter,)
    Accepted normalized RMS misfit for each iteration. When
    an iteration contains repeated search steps, the parser
    keeps the last ``AND IS =`` value in that block.
roughness : ndarray of float, shape (n_iter,)
    Model roughness reported by Occam. The final entry may be
    ``nan`` when a run stops before writing ``ROUGHNESS IS``.
lagrange : ndarray of float, shape (n_iter,)
    Accepted Lagrange multiplier, :math:`\mu`, for each
    iteration. Values are read from ``MINIMUM TOL FROM`` or
    ``INTERCEPT IS AT MU`` lines.
stepsize : ndarray of float, shape (n_iter,)
    Accepted step size for each iteration. The final entry may
    be ``nan`` if convergence problems stop the run early.

Notes
-----
The parser is intentionally tolerant of Occam2D log variants.
It ignores intermediate ``TOFMU`` search lines and keeps the
last accepted values in each iteration block. This behavior
matches logs where divergence problems trigger repeated
Lagrange searches before a step is accepted.

See Also
--------
OccamRunner
    Produces the log file by launching the executable.
InversionResult
    Loads logs with model, iteration, and response files.
Plot2D.misfit
    Visualizes RMS convergence from an ``OccamLog`` object.

Examples
--------
Read a log and inspect the best iteration:

>>> from pycsamt.models.occam2d import OccamLog
>>> log = OccamLog.read("occam_run/LogFile.logfile")
>>> log.best_iteration

Check whether the inversion reached the common RMS target:

>>> log.converged

Print a compact report for scripts:

>>> log.summary()

References
----------
.. [1] Constable, S. C., Parker, R. L., and Constable,
   C. G., "Occam's inversion: A practical algorithm for
   generating smooth models from electromagnetic sounding
   data", Geophysics, 52(3), 289-300, 1987.
.. [2] deGroot-Hedlin, C., and Constable, S.,
   "Occam's inversion to generate smooth, two-dimensional
   models from magnetotelluric data", Geophysics, 55(12),
   1613-1624, 1990.
"""

OccamLog.read.__func__.__doc__ = rf"""
Read an Occam2D log file.

The reader scans the file line by line with a small state
machine. A ``** ITERATION`` line starts a new block and saves
the previous block. Within each block, the last accepted RMS,
roughness, Lagrange multiplier, and step size are retained.

This is useful for logs where Occam cuts the step size or
repeats the Lagrange search. Intermediate trial values are not
stored because they do not describe the accepted model.

Parameters
----------
{_params.common.path}
**kwargs
    Additional keyword arguments forwarded to the ``OccamLog``
    constructor. Use this for ``verbose`` or ``logger`` when
    integrating the parser into a larger workflow.

Returns
-------
OccamLog
    Parsed convergence-log container with one array entry per
    completed iteration block.

Raises
------
FileNotFoundError
    Raised when ``path`` does not exist.

See Also
--------
OccamLog.summary
    Returns a short text summary of convergence.
OccamLog.best_iteration
    Reports the iteration with the lowest finite RMS misfit.

Examples
--------
>>> from pycsamt.models.occam2d import OccamLog
>>> log = OccamLog.read("occam_run/LogFile.logfile")
>>> log.n_iter
>>> log.rms[-1]
"""
