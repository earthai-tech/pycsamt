# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""ModEmControl — read and write the ModEM inversion control file.

The control file is a plain key–value text file.  The same format is used
for both 2-D and 3-D runs.  Example::

    Model and data output file name    : ModEM_out
    Initial damping factor lambda      : 10.
    To update lambda divide by         : 100.
    Initial search step in model units : 10.
    Restart when rms diff is less than : 5.0e-4
    Exit search when rms is less than  : 1.05
    Exit when lambda is less than      : 1.0e-4
    Maximum number of iterations       : 100
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

from .base   import ModEmBase
from .config import ModEmConfig

PathLike = Union[str, Path]

__all__ = ["ModEmControl"]

_KEYS = [
    ("output_stem",       "Model and data output file name"),
    ("initial_lambda",    "Initial damping factor lambda"),
    ("lambda_divisor",    "To update lambda divide by"),
    ("initial_alpha",     "Initial search step in model units"),
    ("rms_diff_tol",      "Restart when rms diff is less than"),
    ("target_rms",        "Exit search when rms is less than"),
    ("lambda_exit",       "Exit when lambda is less than"),
    ("max_iterations",    "Maximum number of iterations"),
]
_FLOAT_ATTRS = {"initial_lambda", "lambda_divisor", "initial_alpha",
                "rms_diff_tol",   "target_rms",     "lambda_exit"}
_INT_ATTRS   = {"max_iterations"}


class ModEmControl(ModEmBase):
    """ModEM inversion control file container.

    Attributes
    ----------
    output_stem : str
        Base name for all ModEM output files.
    initial_lambda : float
    lambda_divisor : float
    initial_alpha : float
    rms_diff_tol : float
    target_rms : float
    lambda_exit : float
    max_iterations : int
    """

    def __init__(self, config: Optional[ModEmConfig] = None, **kwargs):
        super().__init__(**kwargs)
        cfg = config or ModEmConfig()
        self.config:         ModEmConfig = cfg
        self.output_stem:    str         = cfg.output_stem
        self.initial_lambda: float       = cfg.initial_lambda
        self.lambda_divisor: float       = cfg.lambda_divisor
        self.initial_alpha:  float       = cfg.initial_alpha
        self.rms_diff_tol:   float       = cfg.rms_diff_tol
        self.target_rms:     float       = cfg.target_rms
        self.lambda_exit:    float       = cfg.lambda_exit
        self.max_iterations: int         = cfg.max_iterations

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------
    @classmethod
    def from_config(
        cls,
        config: Optional[ModEmConfig] = None,
        **kwargs,
    ) -> "ModEmControl":
        """Build a control object from a ``ModEmConfig``."""
        return cls(config=config, **kwargs)

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> "ModEmControl":
        """Parse an existing ModEM control file."""
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"ModEM control file not found: {p}")

        obj = cls(**kwargs)
        with p.open("r", errors="replace") as fh:
            for raw in fh:
                if ":" not in raw:
                    continue
                key_raw, _, val_raw = raw.partition(":")
                key = key_raw.strip().upper()
                val = val_raw.strip()
                for attr, label in _KEYS:
                    if label.upper() == key:
                        if attr in _FLOAT_ATTRS:
                            try:
                                setattr(obj, attr, float(val))
                            except ValueError:
                                pass
                        elif attr in _INT_ATTRS:
                            try:
                                setattr(obj, attr, int(float(val)))
                            except ValueError:
                                pass
                        else:
                            setattr(obj, attr, val)
                        break
        if obj.verbose:
            obj.logger.info("ModEmControl.read: loaded from %s", p)
        return obj

    def write(self, path: PathLike) -> Path:
        """Write control file to *path*.

        Returns the resolved path.
        """
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)

        _W = 44   # key field width
        lines: list[str] = []
        for attr, label in _KEYS:
            val = getattr(self, attr)
            key = f"{label:<{_W}}"
            if attr in _FLOAT_ATTRS:
                lines.append(f"{key}: {val:.4g}\n")
            elif attr in _INT_ATTRS:
                lines.append(f"{key}: {int(val)}\n")
            else:
                lines.append(f"{key}: {val}\n")

        with p.open("w") as fh:
            fh.writelines(lines)
        return p
