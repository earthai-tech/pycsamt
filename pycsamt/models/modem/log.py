# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""ModEmLog — parse the ModEM NLCG log file.

The log file (``Modular_NLCG.log`` or ``inverse.log``) records one or more
lines per iteration.  The key lines we extract are::

    Completed NLCG iteration     N
         with: f=<f>  m2=<m2>  rms=<rms>  lambda=<λ>  alpha=<α>

as well as the ``START:`` line for iteration 0.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Union

import numpy as np

from .base import ModEmBase

PathLike = Union[str, Path]

__all__ = ["ModEmLog"]

# Patterns
_RE_COMPLETED = re.compile(
    r"Completed\s+NLCG\s+iteration\s+(\d+)", re.IGNORECASE
)
_RE_WITH = re.compile(
    r"with:\s+f=([\d.E+\-]+)\s+m2=([\d.E+\-]+)\s+rms=\s*([\d.E+\-]+)"
    r"\s+lambda=\s*([\d.E+\-]+)\s+alpha=\s*([\d.E+\-]+)",
    re.IGNORECASE,
)
_RE_START = re.compile(
    r"START:\s+f=([\d.E+\-]+)\s+m2=([\d.E+\-]+)\s+rms=\s*([\d.E+\-]+)"
    r"\s+lambda=\s*([\d.E+\-]+)\s+alpha=\s*([\d.E+\-]+)",
    re.IGNORECASE,
)


def _parse_log(path: Path) -> dict:
    records: list[dict] = []

    with path.open("r", errors="replace") as fh:
        lines = fh.readlines()

    pending_iter: int | None = None

    for ln in lines:
        # START line → iteration 0
        m = _RE_START.search(ln)
        if m and not records:
            records.append({
                "iteration": 0,
                "f":         float(m.group(1)),
                "m2":        float(m.group(2)),
                "rms":       float(m.group(3)),
                "lambda":    float(m.group(4)),
                "alpha":     float(m.group(5)),
            })
            continue

        m = _RE_COMPLETED.search(ln)
        if m:
            pending_iter = int(m.group(1))
            continue

        if pending_iter is not None:
            m = _RE_WITH.search(ln)
            if m:
                records.append({
                    "iteration": pending_iter,
                    "f":         float(m.group(1)),
                    "m2":        float(m.group(2)),
                    "rms":       float(m.group(3)),
                    "lambda":    float(m.group(4)),
                    "alpha":     float(m.group(5)),
                })
                pending_iter = None

    return {"records": records}


class ModEmLog(ModEmBase):
    """ModEM NLCG log container.

    Attributes
    ----------
    iterations : np.ndarray, shape (n_iter,)
        Iteration numbers (0-based; 0 = initial START).
    rms : np.ndarray
        RMS misfit per iteration.
    objective : np.ndarray
        Total objective function value f per iteration.
    model_norm : np.ndarray
        Model roughness / norm m² per iteration.
    lagrange : np.ndarray
        Damping parameter λ per iteration.
    alpha : np.ndarray
        Step-size α per iteration.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.iterations:  np.ndarray = np.array([], dtype=int)
        self.rms:         np.ndarray = np.array([])
        self.objective:   np.ndarray = np.array([])
        self.model_norm:  np.ndarray = np.array([])
        self.lagrange:    np.ndarray = np.array([])
        self.alpha:       np.ndarray = np.array([])

    @property
    def n_iter(self) -> int:
        return len(self.iterations)

    @property
    def final_rms(self) -> float:
        return float(self.rms[-1]) if self.rms.size else float("nan")

    @property
    def best_iter(self) -> int:
        if self.rms.size == 0:
            return 0
        return int(self.iterations[int(np.argmin(self.rms))])

    @classmethod
    def read(cls, path: PathLike, **kwargs) -> "ModEmLog":
        """Parse a ModEM NLCG log file."""
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"ModEM log file not found: {p}")

        parsed = _parse_log(p)
        obj    = cls(**kwargs)
        recs   = parsed["records"]

        if recs:
            obj.iterations = np.array([r["iteration"] for r in recs], dtype=int)
            obj.rms        = np.array([r["rms"]       for r in recs])
            obj.objective  = np.array([r["f"]         for r in recs])
            obj.model_norm = np.array([r["m2"]        for r in recs])
            obj.lagrange   = np.array([r["lambda"]    for r in recs])
            obj.alpha      = np.array([r["alpha"]     for r in recs])

        if obj.verbose:
            obj.logger.info(
                "ModEmLog.read: %d iterations, final RMS=%.4f from %s",
                obj.n_iter, obj.final_rms, p,
            )
        return obj
