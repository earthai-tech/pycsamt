# -*- coding: utf-8 -*-
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

PathLike = Union[str, Path]

__all__ = ["OccamLog"]

_NAN = float("nan")


def _parse_float(line: str) -> float:
    """Return the float after the last '=' on *line*, or NaN on failure."""
    try:
        return float(line.rsplit("=", 1)[1].strip())
    except (ValueError, IndexError):
        return _NAN


class OccamLog(OccamBase):
    """Occam2D log-file container.

    One entry per iteration in every array; arrays are aligned so that
    ``iterations[i]``, ``rms[i]``, ``roughness[i]``, etc. all describe
    the same iteration.

    Attributes
    ----------
    iterations : np.ndarray[int], shape (n_iter,)
        Iteration numbers (1-based, as written by Occam).
    rms : np.ndarray[float], shape (n_iter,)
        Achieved RMS misfit at the end of each iteration.
        (The last "AND IS = …" value before ROUGHNESS IS for that block.)
    roughness : np.ndarray[float], shape (n_iter,)
        Roughness value (deGroot-Hedlin measure).  ``NaN`` for the last
        iteration when the run ended without writing ROUGHNESS IS.
    lagrange : np.ndarray[float], shape (n_iter,)
        Lagrange multiplier (µ) at the accepted step of each iteration.
        Sourced from "MINIMUM TOL FROM fminocc IS AT MU" or
        "INTERCEPT IS AT MU".
    stepsize : np.ndarray[float], shape (n_iter,)
        Stepsize value written at the end of each iteration.  ``NaN``
        if the run ended without writing STEPSIZE IS.
    n_iter : int
        Total number of parsed iterations.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.iterations: np.ndarray = np.array([], dtype=int)
        self.rms:        np.ndarray = np.array([])
        self.roughness:  np.ndarray = np.array([])
        self.lagrange:   np.ndarray = np.array([])
        self.stepsize:   np.ndarray = np.array([])

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> "OccamLog":
        """Parse an Occam log file at *path*.

        Reads the file line-by-line using a state machine:  a new
        "** ITERATION **" line triggers saving the previous iteration's
        accumulated values.  The final iteration is saved after the loop.

        Parameters
        ----------
        path : path-like
            Path to the Occam log file (e.g. ``LogFile.logfile``).

        Returns
        -------
        OccamLog
            Populated instance with one array entry per iteration.

        Raises
        ------
        FileNotFoundError
            If *path* does not exist.
        """
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"Occam log file not found: {p}")

        obj = cls(**kwargs)

        # Per-iteration accumulators
        iter_nums:   list[int]   = []
        rms_vals:    list[float] = []
        rough_vals:  list[float] = []
        step_vals:   list[float] = []
        lagr_vals:   list[float] = []

        # Mutable state for the iteration currently being parsed.
        # Using a one-element list so the inner helper can mutate them.
        cur = {
            "iter":      None,
            "rms":       _NAN,
            "roughness": _NAN,
            "stepsize":  _NAN,
            "lagrange":  _NAN,
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
            cur["iter"]      = None
            cur["rms"]       = _NAN
            cur["roughness"] = _NAN
            cur["stepsize"]  = _NAN
            cur["lagrange"]  = _NAN

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
        obj.rms        = np.array(rms_vals)
        obj.roughness  = np.array(rough_vals)
        obj.stepsize   = np.array(step_vals)
        obj.lagrange   = np.array(lagr_vals)

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
