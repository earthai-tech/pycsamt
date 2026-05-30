# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""OccamStartup and OccamIter — Startup / .iter file handling.

Both file types share the ``OCCAMITER_FLEX`` format header.  The
distinction is that a Startup file has ``Iteration: 0`` and a uniform
initial model vector, while an ``.iter`` file carries a real iteration
number and the model parameters from that iteration.

OCCAMITER_FLEX format
---------------------
The file is split into two sections separated by ``Param Count: N``:

Header (keyword: value pairs, one per line)::

    Format:             OCCAMITER_FLEX
    Description:        <text>
    Model File:         <filename>
    Data File:          <filename>
    Date/Time:          <text>
    Iterations to run:  <int>
    Target Misfit:      <float>
    Roughness Type:     <int>
    Diagonal Penalties: <int>
    Stepsize Cut Count: <int>
    !Model Limits:      <text>       ← comment lines (!) are skipped
    !Model Value Steps: <text>
    Debug Level:        <int>
    Iteration:          <int>        ← 0 = Startup, N ≥ 1 = iter file
    Lagrange Value:     <float>
    Roughness Value:    <float>
    Misfit Value:       <float>
    Misfit Reached:     <int>        ← 0 or 1
    Param Count:        <int>        ← triggers switch to param block

Parameter block (space-separated floats, 4 per line)::

    log10(ρ₀)  log10(ρ₁)  log10(ρ₂)  log10(ρ₃)
    ...

Classes
-------
OccamStartup
    Reads / writes / builds the initial Startup file (Iteration: 0).
OccamIter
    Reads any ``.iter`` file (Iteration: N ≥ 1) produced by Occam.
"""

from __future__ import annotations

import datetime
from pathlib import Path
from typing import Dict, Optional, Union

import numpy as np

from .base   import OccamBase
from .config import OccamConfig

PathLike = Union[str, Path]

__all__ = ["OccamStartup", "OccamIter"]

_FORMAT_TAG = "OCCAMITER_FLEX"

# -----------------------------------------------------------------------
# Shared low-level parser
# -----------------------------------------------------------------------

# Canonical key map: upper-cased header keyword → attribute name
_KEY_MAP: Dict[str, str] = {
    "FORMAT":             "format_str",
    "DESCRIPTION":        "description",
    "MODEL FILE":         "model_file",
    "DATA FILE":          "data_file",
    "DATE/TIME":          "datetime_str",
    "ITERATIONS TO RUN":  "max_iterations",
    "TARGET MISFIT":      "target_misfit",
    "ROUGHNESS TYPE":     "roughness_type",
    "DIAGONAL PENALTIES": "diagonal_penalties",
    "STEPSIZE CUT COUNT": "stepsize_cut_count",
    "DEBUG LEVEL":        "debug_level",
    "ITERATION":          "iteration",
    "LAGRANGE VALUE":     "lagrange_value",
    "ROUGHNESS VALUE":    "roughness_value",
    "MISFIT VALUE":       "misfit_value",
    "MISFIT REACHED":     "misfit_reached",
    "PARAM COUNT":        "n_params",
}

# Fields that need integer coercion
_INT_FIELDS = {
    "max_iterations", "roughness_type", "diagonal_penalties",
    "stepsize_cut_count", "debug_level", "iteration", "n_params",
}


def _parse_iter_flex(path: Path) -> dict:
    """Parse an OCCAMITER_FLEX file and return a plain dict.

    Parameters
    ----------
    path : Path
        Resolved path to the file.

    Returns
    -------
    dict
        Keys matching ``_KEY_MAP`` values plus ``"param_values"``
        (``np.ndarray`` of log₁₀ resistivity, length ``n_params``).

    Raises
    ------
    FileNotFoundError
        If *path* does not exist.
    ValueError
        If the file does not contain a valid OCCAMITER_FLEX header.
    """
    if not path.exists():
        raise FileNotFoundError(f"OCCAMITER_FLEX file not found: {path}")

    result: dict = {v: None for v in _KEY_MAP.values()}
    param_tokens: list[str] = []
    reading_params = False
    n_params = 0

    with path.open("r", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("!"):
                continue

            if reading_params:
                param_tokens.extend(line.split())
                continue

            # Header key: value
            if ":" not in line:
                continue
            raw_key, _, raw_val = line.partition(":")
            key_up = raw_key.strip().upper()
            val    = raw_val.strip()

            attr = _KEY_MAP.get(key_up)
            if attr is None:
                continue

            if attr in _INT_FIELDS:
                try:
                    result[attr] = int(float(val))
                except ValueError:
                    result[attr] = 0
            elif attr == "misfit_reached":
                try:
                    result[attr] = bool(int(float(val)))
                except ValueError:
                    result[attr] = False
            elif attr == "format_str":
                result[attr] = val
                if val.upper() != _FORMAT_TAG:
                    raise ValueError(
                        f"Expected format '{_FORMAT_TAG}', got '{val}' in {path}"
                    )
            else:
                # float fields and string fields
                try:
                    result[attr] = float(val)
                except ValueError:
                    result[attr] = val  # store as string (description, filenames…)

            if attr == "n_params":
                n_params = result["n_params"] or 0
                reading_params = True

    if result["format_str"] is None:
        raise ValueError(f"File does not contain a valid OCCAMITER_FLEX header: {path}")

    # Parse parameter vector
    try:
        arr = np.array([float(t) for t in param_tokens[:n_params]], dtype=float)
    except ValueError:
        arr = np.array([], dtype=float)

    if arr.size < n_params and n_params > 0:
        # Pad with NaN if fewer tokens than declared (truncated file)
        pad = np.full(n_params - arr.size, np.nan)
        arr = np.concatenate([arr, pad])

    result["param_values"] = arr
    return result


# -----------------------------------------------------------------------
# OccamStartup
# -----------------------------------------------------------------------

class OccamStartup(OccamBase):
    """Occam2D startup (iteration-0) file.

    Attributes
    ----------
    description : str
    model_file : str
    data_file : str
    datetime_str : str
    max_iterations : int
    target_misfit : float
    roughness_type : int
    diagonal_penalties : int
    stepsize_cut_count : int
    debug_level : int
    iteration : int
        Always 0 for a valid Startup file.
    lagrange_value : float
        Initial Lagrange multiplier (written as ``Lagrange Value``).
    roughness_value : float
    misfit_value : float
    misfit_reached : bool
    n_params : int
    param_values : np.ndarray, shape (n_params,)
        Log₁₀ resistivity values; uniform ``log10(initial_rho)`` at startup.
    """

    def __init__(
        self,
        config: Optional[OccamConfig] = None,
        description: str = "startup created by pycsamt",
        **kwargs,
    ):
        super().__init__(**kwargs)
        cfg = config or OccamConfig()
        self.config             = cfg
        self.format_str         = _FORMAT_TAG
        self.description        = description
        self.model_file         = cfg.model_file
        self.data_file          = cfg.data_file
        self.datetime_str       = str(datetime.datetime.now())
        self.max_iterations     = cfg.max_iterations
        self.target_misfit      = cfg.target_misfit
        self.roughness_type     = cfg.roughness_type
        self.diagonal_penalties = cfg.diagonal_penalties
        self.stepsize_cut_count = cfg.stepsize_cut_count
        self.debug_level        = cfg.debug_level
        self.iteration          = 0
        self.lagrange_value     = cfg.lagrange_start
        self.roughness_value    = 1e10
        self.misfit_value       = 1000.0
        self.misfit_reached     = False
        self.n_params           = 0
        self.param_values: np.ndarray = np.array([])

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------
    @classmethod
    def from_model(
        cls,
        model: "OccamModel",  # noqa: F821
        config: Optional[OccamConfig] = None,
        **kwargs,
    ) -> "OccamStartup":
        """Build a startup from a populated model object.

        Sets the parameter vector to ``log10(config.initial_rho)``
        repeated ``model.n_params`` times.

        Parameters
        ----------
        model : OccamModel
            Populated model (must have ``n_params > 0``).
        config : OccamConfig, optional
            Run configuration.  Provides ``initial_rho``, file names, …

        Returns
        -------
        OccamStartup
        """
        cfg = config or OccamConfig()
        obj = cls(config=cfg, **kwargs)

        n = model.n_params
        if n <= 0:
            raise ValueError(
                "OccamStartup.from_model: model has no parameters"
            )

        obj.n_params     = n
        obj.param_values = np.full(n, np.log10(cfg.initial_rho), dtype=float)
        obj.model_file   = cfg.model_file
        obj.data_file    = cfg.data_file

        if obj.verbose:
            obj.logger.info(
                "OccamStartup.from_model: %d params, rho0=%.1f Ω·m",
                n, cfg.initial_rho,
            )
        return obj

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> "OccamStartup":
        """Parse a Startup file (``Iteration: 0``).

        Parameters
        ----------
        path : path-like
            Path to the Startup file.

        Returns
        -------
        OccamStartup

        Raises
        ------
        FileNotFoundError
            If *path* does not exist.
        ValueError
            If the file is not OCCAMITER_FLEX or ``Iteration`` ≠ 0.
        """
        p = Path(path)
        parsed = _parse_iter_flex(p)

        if parsed.get("iteration", 0) != 0:
            raise ValueError(
                f"Expected Iteration: 0 (Startup) but got "
                f"Iteration: {parsed['iteration']} in {p}"
            )

        obj = cls(**kwargs)
        _apply_parsed(obj, parsed)

        if obj.verbose:
            obj.logger.info(
                "OccamStartup.read: %d params loaded from %s", obj.n_params, p
            )
        return obj

    def write(self, path: PathLike) -> Path:
        """Write startup file to *path* in OCCAMITER_FLEX format.

        Returns the resolved path.
        """
        return _write_iter_flex(self, Path(path))


# -----------------------------------------------------------------------
# OccamIter
# -----------------------------------------------------------------------

class OccamIter(OccamBase):
    """Reader for Occam ``.iter`` files (``Iteration: N ≥ 1``).

    Attributes
    ----------
    format_str : str
    description : str
    model_file : str
    data_file : str
    datetime_str : str
    max_iterations : int
    target_misfit : float
    roughness_type : int
    diagonal_penalties : int
    stepsize_cut_count : int
    debug_level : int
    iteration : int
        Iteration number (≥ 1).
    lagrange_value : float
        Lagrange multiplier accepted at this iteration.
    roughness_value : float
    misfit_value : float
        Normalised RMS misfit at this iteration.
    misfit_reached : bool
        ``True`` if the target misfit was achieved.
    n_params : int
    param_values : np.ndarray, shape (n_params,)
        Log₁₀ resistivity values at this iteration.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.format_str:      str        = _FORMAT_TAG
        self.description:     str        = ""
        self.model_file:      str        = ""
        self.data_file:       str        = ""
        self.datetime_str:    str        = ""
        self.max_iterations:  int        = 0
        self.target_misfit:   float      = 1.0
        self.roughness_type:  int        = 1
        self.diagonal_penalties: int     = 0
        self.stepsize_cut_count: int     = 8
        self.debug_level:     int        = 1
        self.iteration:       int        = 0
        self.lagrange_value:  float      = 0.0
        self.roughness_value: float      = 0.0
        self.misfit_value:    float      = 0.0
        self.misfit_reached:  bool       = False
        self.n_params:        int        = 0
        self.param_values:    np.ndarray = np.array([])

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> "OccamIter":
        """Parse an ``.iter`` file at *path*.

        Parameters
        ----------
        path : path-like
            Path to the ``.iter`` file.

        Returns
        -------
        OccamIter

        Raises
        ------
        FileNotFoundError
            If *path* does not exist.
        ValueError
            If the file is not OCCAMITER_FLEX or ``Iteration`` == 0
            (use ``OccamStartup.read()`` for the Startup file).
        """
        p = Path(path)
        parsed = _parse_iter_flex(p)

        if parsed.get("iteration", 0) == 0:
            raise ValueError(
                f"File has Iteration: 0 — this is a Startup file. "
                f"Use OccamStartup.read() instead: {p}"
            )

        obj = cls(**kwargs)
        _apply_parsed(obj, parsed)

        if obj.verbose:
            obj.logger.info(
                "OccamIter.read: iter %d, misfit=%.4f, %d params from %s",
                obj.iteration, obj.misfit_value, obj.n_params, p,
            )
        return obj

    def write(self, path: PathLike) -> Path:
        """Write ``.iter`` file to *path* in OCCAMITER_FLEX format.

        Returns the resolved path.
        """
        return _write_iter_flex(self, Path(path))

    # ------------------------------------------------------------------
    # Derived
    # ------------------------------------------------------------------
    def to_resistivity(self) -> np.ndarray:
        """Return resistivity model (Ω·m) from log₁₀ param_values."""
        if not self.param_values.size:
            return np.array([])
        return 10.0 ** self.param_values

    @property
    def log10_rho_stats(self) -> dict:
        """Basic statistics of the log₁₀ resistivity parameter vector."""
        if not self.param_values.size:
            return {}
        v = self.param_values
        return {
            "min":    float(np.nanmin(v)),
            "max":    float(np.nanmax(v)),
            "mean":   float(np.nanmean(v)),
            "median": float(np.nanmedian(v)),
        }


# -----------------------------------------------------------------------
# Shared writer (used by both OccamStartup.write and OccamIter.write)
# -----------------------------------------------------------------------

# Ordered write sequence: (file keyword, attribute name, value type)
_WRITE_ORDER = [
    ("Format",             "format_str",         "str"),
    ("Description",        "description",        "str"),
    ("Model File",         "model_file",         "str"),
    ("Data File",          "data_file",          "str"),
    ("Date/Time",          "datetime_str",       "str"),
    ("Iterations to run",  "max_iterations",     "int"),
    ("Target Misfit",      "target_misfit",      "float"),
    ("Roughness Type",     "roughness_type",     "int"),
    ("Diagonal Penalties", "diagonal_penalties", "int"),
    ("Stepsize Cut Count", "stepsize_cut_count", "int"),
    ("Debug Level",        "debug_level",        "int"),
    ("Iteration",          "iteration",          "int"),
    ("Lagrange Value",     "lagrange_value",     "float"),
    ("Roughness Value",    "roughness_value",    "float"),
    ("Misfit Value",       "misfit_value",       "float"),
    ("Misfit Reached",     "misfit_reached",     "int"),
    ("Param Count",        "n_params",           "int"),
]

_KEY_WIDTH = 20   # columns reserved for "Key:" + padding


def _write_iter_flex(obj: OccamBase, path: Path) -> Path:
    """Serialise an OCCAMITER_FLEX object to *path*."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []

    for kw, attr, vtype in _WRITE_ORDER:
        val = getattr(obj, attr, None)
        if val is None:
            val = 0
        key_str = f"{kw}:"
        pad      = max(1, _KEY_WIDTH - len(key_str))
        if vtype == "int":
            val_str = str(int(val))
        elif vtype == "float":
            val_str = str(float(val))
        else:
            val_str = str(val)
        lines.append(f"{key_str}{' ' * pad}{val_str}\n")

        if attr == "stepsize_cut_count":
            # Insert comment lines between Stepsize Cut Count and Debug Level
            lines.append("!Model Limits:      none\n")
            lines.append("!Model Value Steps: none\n")

    # Param block: 4 values per line
    pv = getattr(obj, "param_values", np.array([]))
    for i in range(0, len(pv), 4):
        chunk = pv[i: i + 4]
        lines.append("  " + "    ".join(f"{v:.4f}" for v in chunk) + "\n")

    with path.open("w") as fh:
        fh.writelines(lines)

    return path


# -----------------------------------------------------------------------
# Helper: apply parsed dict to any object that has the right attributes
# -----------------------------------------------------------------------

def _apply_parsed(obj: OccamBase, parsed: dict) -> None:
    """Copy all non-None values from *parsed* onto *obj*."""
    for attr, val in parsed.items():
        if val is not None and hasattr(obj, attr):
            setattr(obj, attr, val)
    # Ensure n_params is consistent with param_values length
    pv = parsed.get("param_values")
    if pv is not None and isinstance(pv, np.ndarray):
        obj.param_values = pv
        if not parsed.get("n_params"):
            obj.n_params = len(pv)
