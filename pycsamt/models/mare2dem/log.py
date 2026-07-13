# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""MARE2DEM log file parsers.

Provides:

* :class:`Mare2DEMLog` — parse the main ``OccamLog.2012.0`` per-iteration
  log written by MARE2DEM during inversion.
* :class:`GroupRMSLog` (re-exported) — parse the per-group RMS CSV log.
* :class:`IterationRecord` — one row of the per-iteration log.

MARE2DEM log format (``OccamLog.2012.0``)::

    Format:     OccamLog.2012.0
    ...
    ** Iteration     1 **
    ...
                  Model Misfit:        6.875
                     Roughness:        1.407
                    Optimal Mu:        4.000
            Convergence Status:        0
    ...
    ** Iteration     2 **
    ...

Convergence is detected when ``Convergence Status: 1`` or when the
ending line contains ``"target misfit achieved"`` / ``"misfit reached"``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from .iotools.group_rms import GroupRMSLog, read_group_rms_log

__all__ = [
    "IterationRecord",
    "Mare2DEMLog",
    "GroupRMSLog",
    "read_group_rms_log",
]


@dataclass
class IterationRecord:
    """One completed iteration from a MARE2DEM log file.

    Attributes
    ----------
    iteration : int
        Iteration number.
    rms : float
        Normalized RMS misfit at this iteration.
    roughness : float
        Model roughness.
    lambda_ : float
        Log10 of the Lagrange (regularization) multiplier.
    """

    iteration: int
    rms: float
    roughness: float
    lambda_: float


class Mare2DEMLog:
    """Parse the MARE2DEM per-iteration ``OccamLog.2012.0`` log file.

    MARE2DEM writes one block per completed iteration containing
    ``Model Misfit``, ``Roughness``, and ``Optimal Mu`` lines.  This
    parser extracts those values and exposes them as
    :class:`IterationRecord` objects.

    Parameters
    ----------
    path : path-like
        Path to the log file (usually ``*.logfile`` or ``*.log``).

    Attributes
    ----------
    path : pathlib.Path
    iterations : list of IterationRecord
    converged : bool
    """

    # Iteration header: "** Iteration     N **"
    _ITER_HDR = re.compile(r"\*\*\s+Iteration\s+(\d+)\s+\*\*")
    # Summary value lines
    _MISFIT = re.compile(r"Model Misfit\s*:\s*([\d.eE+\-]+)")
    _ROUGH = re.compile(r"Roughness\s*:\s*([\d.eE+\-]+)")
    _MU = re.compile(r"Optimal Mu\s*:\s*([\d.eE+\-]+)")
    _CONV = re.compile(r"Convergence Status\s*:\s*(\d+)")

    def __init__(self, path: str | Path):
        self.path = Path(path).resolve()
        self.iterations: list[IterationRecord] = []
        self.converged: bool = False
        if self.path.exists():
            self._parse()

    def _parse(self) -> None:
        text = self.path.read_text(errors="replace")
        lines = text.splitlines()

        current_iter: int | None = None
        misfit: float | None = None
        roughness: float | None = None
        mu: float | None = None

        for line in lines:
            # check for iteration header
            m = self._ITER_HDR.search(line)
            if m:
                # flush previous iteration if complete
                if current_iter is not None and misfit is not None:
                    self.iterations.append(
                        IterationRecord(
                            iteration=current_iter,
                            rms=misfit,
                            roughness=roughness or 0.0,
                            lambda_=mu or 0.0,
                        )
                    )
                current_iter = int(m.group(1))
                misfit = roughness = mu = None
                continue

            m = self._MISFIT.search(line)
            if m:
                misfit = float(m.group(1))
                continue

            m = self._ROUGH.search(line)
            if m:
                roughness = float(m.group(1))
                continue

            m = self._MU.search(line)
            if m:
                mu = float(m.group(1))
                continue

            m = self._CONV.search(line)
            if m:
                if int(m.group(1)) == 1:
                    self.converged = True
                continue

            low = line.lower()
            if (
                "target misfit achieved" in low
                or "misfit reached" in low
                or "converged" in low
            ):
                self.converged = True

        # flush last iteration
        if current_iter is not None and misfit is not None:
            self.iterations.append(
                IterationRecord(
                    iteration=current_iter,
                    rms=misfit,
                    roughness=roughness or 0.0,
                    lambda_=mu or 0.0,
                )
            )

    @property
    def final_rms(self) -> float | None:
        """RMS at the last logged iteration."""
        return self.iterations[-1].rms if self.iterations else None

    @property
    def n_iterations(self) -> int:
        """Number of completed iterations in the log."""
        return len(self.iterations)

    def rms_history(self) -> list[float]:
        """Return per-iteration RMS values in order."""
        return [r.rms for r in self.iterations]

    def __repr__(self) -> str:
        return (
            f"Mare2DEMLog(n_iterations={self.n_iterations}, "
            f"final_rms={self.final_rms}, converged={self.converged})"
        )
