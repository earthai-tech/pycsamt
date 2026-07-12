# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Parse ModEM nonlinear conjugate-gradient log files.

The module defines :class:`ModEmLog`, a small container for the
iteration history written by ModEM during NLCG inversion.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Union

import numpy as np

from .base import ModEmBase
from .doc import _modem_param_docs as _params

PathLike = Union[str, Path]

__all__ = ["ModEmLog"]

# Patterns
_RE_COMPLETED = re.compile(
    r"Completed\s+NLCG\s+iteration\s+(\d+)", re.IGNORECASE
)
_RE_WITH = re.compile(
    r"with:\s+f=\s*([\d.E+\-]+)\s+m2=\s*([\d.E+\-]+)\s+rms=\s*([\d.E+\-]+)"
    r"\s+lambda=\s*([\d.E+\-]+)\s+alpha=\s*([\d.E+\-]+)",
    re.IGNORECASE,
)
_RE_START = re.compile(
    r"START:\s+f=\s*([\d.E+\-]+)\s+m2=\s*([\d.E+\-]+)\s+rms=\s*([\d.E+\-]+)"
    r"\s+lambda=\s*([\d.E+\-]+)\s+alpha=\s*([\d.E+\-]+)",
    re.IGNORECASE,
)


def _parse_log(path: Path) -> dict:
    """Parse iteration records from a ModEM NLCG log."""
    records: list[dict] = []

    with path.open("r", errors="replace") as fh:
        lines = fh.readlines()

    pending_iter: int | None = None

    for ln in lines:
        # START line -> iteration 0
        m = _RE_START.search(ln)
        if m and not records:
            records.append(
                {
                    "iteration": 0,
                    "f": float(m.group(1)),
                    "m2": float(m.group(2)),
                    "rms": float(m.group(3)),
                    "lambda": float(m.group(4)),
                    "alpha": float(m.group(5)),
                }
            )
            continue

        m = _RE_COMPLETED.search(ln)
        if m:
            pending_iter = int(m.group(1))
            continue

        if pending_iter is not None:
            m = _RE_WITH.search(ln)
            if m:
                records.append(
                    {
                        "iteration": pending_iter,
                        "f": float(m.group(1)),
                        "m2": float(m.group(2)),
                        "rms": float(m.group(3)),
                        "lambda": float(m.group(4)),
                        "alpha": float(m.group(5)),
                    }
                )
                pending_iter = None

    return {"records": records}


class ModEmLog(ModEmBase):
    def __init__(self, **kwargs):
        """Initialize an empty ModEM log container."""
        super().__init__(**kwargs)
        self.iterations: np.ndarray = np.array([], dtype=int)
        self.rms: np.ndarray = np.array([])
        self.objective: np.ndarray = np.array([])
        self.model_norm: np.ndarray = np.array([])
        self.lagrange: np.ndarray = np.array([])
        self.alpha: np.ndarray = np.array([])

    @property
    def n_iter(self) -> int:
        """Number of parsed log records."""
        return len(self.iterations)

    @property
    def final_rms(self) -> float:
        """Final parsed RMS value, or ``nan`` for an empty log."""
        return float(self.rms[-1]) if self.rms.size else float("nan")

    @property
    def best_iter(self) -> int:
        """Iteration number with the lowest parsed RMS value."""
        if self.rms.size == 0:
            return 0
        return int(self.iterations[int(np.argmin(self.rms))])

    @classmethod
    def read(cls, path: PathLike, **kwargs) -> ModEmLog:
        """Parse a ModEM NLCG log file.

        Parameters
        ----------
        path : path-like
            Path to a ModEM log file, commonly
            ``"Modular_NLCG.log"`` for bundled examples or
            ``"inverse.log"`` for some 3-D runs.
        **kwargs : dict
            Additional keyword arguments forwarded to
            :class:`ModEmLog`, commonly ``verbose`` or
            ``logger``.

        Returns
        -------
        ModEmLog
            Parsed log object containing iteration numbers,
            RMS values, objective values, model norms, lambda
            values, and line-search step lengths.

        Raises
        ------
        FileNotFoundError
            If ``path`` does not exist.

        Examples
        --------
        >>> from pycsamt.models.modem.log import ModEmLog
        >>> log = ModEmLog.read("Modular_NLCG.log")
        >>> log.n_iter > 0
        True
        >>> log.final_rms == log.rms[-1]
        True
        """
        p = Path(path)
        if not p.exists():
            msg = f"ModEM log file not found: {p}"
            raise FileNotFoundError(msg)

        parsed = _parse_log(p)
        obj = cls(**kwargs)
        recs = parsed["records"]

        if recs:
            obj.iterations = np.array(
                [r["iteration"] for r in recs],
                dtype=int,
            )
            obj.rms = np.array([r["rms"] for r in recs])
            obj.objective = np.array([r["f"] for r in recs])
            obj.model_norm = np.array([r["m2"] for r in recs])
            obj.lagrange = np.array([r["lambda"] for r in recs])
            obj.alpha = np.array([r["alpha"] for r in recs])

        if obj.verbose:
            obj.logger.info(
                "ModEmLog.read: %d iterations, final RMS=%.4f from %s",
                obj.n_iter, obj.final_rms, p,
            )
        return obj


ModEmLog.__doc__ = rf"""
Represent a parsed ModEM NLCG iteration log.

``ModEmLog`` stores the convergence history written by ModEM
during nonlinear conjugate-gradient inversion. It extracts the
initial ``START`` record and completed iteration records, then
exposes each tracked quantity as a NumPy array. The object is
used directly by result loaders and plotting helpers to inspect
misfit reduction and inversion progress.

Each parsed record corresponds to a line with values similar to
``f``, ``m2``, ``rms``, ``lambda``, and ``alpha``. These values
summarize the trade-off between data fit and regularization:

.. math::

   \Phi(m) =
   \Phi_d(m) + \lambda \Phi_m(m),

where :math:`\Phi_d` measures data misfit, :math:`\Phi_m`
measures model roughness or size, and :math:`\lambda` is the
damping parameter reported by the log.

Parameters
----------
{_params.common.verbose}
{_params.common.logger}

Attributes
----------
iterations : numpy.ndarray, shape (n_iter,)
    Parsed iteration numbers. The first record is often ``0``
    from the ``START`` line. Some ModEM logs can contain
    restarts, so iteration numbers are not guaranteed to be
    unique or strictly increasing.
rms : numpy.ndarray, shape (n_iter,)
    Normalized root-mean-square misfit for each parsed record.
    Values near one indicate data fit comparable to the
    assigned errors when uncertainties are realistic.
objective : numpy.ndarray, shape (n_iter,)
    Total objective-function value ``f`` reported by ModEM.
model_norm : numpy.ndarray, shape (n_iter,)
    Model roughness or model norm ``m2`` reported by ModEM.
lagrange : numpy.ndarray, shape (n_iter,)
    Lambda damping values used during inversion.
alpha : numpy.ndarray, shape (n_iter,)
    Line-search step length values reported by ModEM.

Notes
-----
``final_rms`` returns the last parsed RMS value, while
``best_iter`` returns the iteration number associated with the
lowest parsed RMS value. If no records are parsed, ``final_rms``
is ``nan`` and ``best_iter`` is ``0``.

The parser is intentionally focused on the standard lines used
by the ModEM examples and common NLCG outputs. Additional log
messages are ignored.

See Also
--------
InversionResult
    Loads :class:`ModEmLog` while scanning a run directory.
PlotMisfit
    Plots RMS history from a parsed log.
ModEmControl
    Defines target RMS and lambda controls used by the run.
ModEmRunner
    Launches the ModEM executable that writes the log.

Examples
--------
Read a ModEM log and inspect convergence:

>>> from pycsamt.models.modem.log import ModEmLog
>>> log = ModEmLog.read("Modular_NLCG.log")
>>> log.n_iter > 0
True
>>> log.best_iter in set(log.iterations)
True

Use the RMS history for plotting:

>>> rms = log.rms
>>> rms.shape == log.iterations.shape
True

References
----------
.. [1] Egbert, G. D., and Kelbert, A., "Computational
   recipes for electromagnetic inverse problems", Geophysical
   Journal International, 189(1), 251-267, 2012,
   doi:10.1111/j.1365-246X.2011.05347.x.
.. [2] Kelbert, A., Meqbel, N., Egbert, G. D., and Tandon,
   K., "ModEM: A modular system for inversion of
   electromagnetic geophysical data", Computers and
   Geosciences, 66, 40-53, 2014,
   doi:10.1016/j.cageo.2014.01.010.
"""
