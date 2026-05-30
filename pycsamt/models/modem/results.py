# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""InversionResult — aggregated view of a completed ModEM run.

After an inversion finishes the working directory contains a collection of
files:

* log file   (``*.log``)         — RMS history
* control    (``*.inv``)         — run parameters
* models     (``m*.rho`` / ``m*.ws``) — initial, iteration, final
* data       (``d*.dat``)        — observed and predicted responses

``InversionResult`` scans a workdir, auto-detects 2D vs 3D, and loads
all available artefacts into a single container.

Quick start
-----------
>>> from pycsamt.models.modem.results import InversionResult
>>> r = InversionResult("./modem_run")
>>> print(r.final_rms)
>>> r.model_final.rho_linear   # 2D: (nz, nx)  /  3D: (nz, ny, nx)
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Optional, Union

from .base       import ModEmBase
from .config     import ModEmConfig
from .data       import ModEmData
from .log        import ModEmLog
from .control    import ModEmControl
from .model2d    import ModEmModel2D
from .model3d    import ModEmModel3D
from .covariance import ModEmCovariance

PathLike = Union[str, Path]

__all__ = ["InversionResult"]

# Stem patterns ModEM uses for model files
_MODEL_STEM_RE = re.compile(r"^m(\d+|i)$", re.IGNORECASE)
# Data stems
_DATA_STEM_RE  = re.compile(r"^d(\d+|i|e|[a-z]+_\w+)$", re.IGNORECASE)


def _stem(p: Path) -> str:
    return p.stem.split("_")[0]          # strip "_JT" or similar suffixes


class InversionResult(ModEmBase):
    """Container for all artefacts produced by a ModEM run.

    Parameters
    ----------
    workdir : path-like
        Directory containing ModEM output files.
    config : ModEmConfig, optional

    Attributes
    ----------
    workdir : Path
    mode : str
        ``'2d'`` or ``'3d'``, detected from model file extensions.
    log : ModEmLog or None
    control : ModEmControl or None
    model_initial : ModEmModel2D or ModEmModel3D or None
        The ``m0`` starting model.
    model_final : ModEmModel2D or ModEmModel3D or None
        The ``mi`` (or highest-numbered) model.
    models : dict[str, ModEmModel2D | ModEmModel3D]
        All model files keyed by stem (``'m0'``, ``'m1'``, ``'mi'``, …).
    data_obs : ModEmData or None
        The ``d0`` observed-data file.
    data_pred : ModEmData or None
        The ``di`` (or ``pr``) predicted-response file.
    covariance : ModEmCovariance or None
    """

    def __init__(
        self,
        workdir: PathLike,
        config: Optional[ModEmConfig] = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.workdir: Path = Path(workdir)
        if not self.workdir.is_dir():
            raise FileNotFoundError(f"ModEM workdir not found: {self.workdir}")

        self.config: ModEmConfig = config or ModEmConfig()

        self.mode:          str     = "unknown"
        self.log:           Optional[ModEmLog]                        = None
        self.control:       Optional[ModEmControl]                    = None
        self.model_initial: Optional[Union[ModEmModel2D, ModEmModel3D]] = None
        self.model_final:   Optional[Union[ModEmModel2D, ModEmModel3D]] = None
        self.models:        Dict[str, Union[ModEmModel2D, ModEmModel3D]] = {}
        self.data_obs:      Optional[ModEmData]                       = None
        self.data_pred:     Optional[ModEmData]                       = None
        self.covariance:    Optional[ModEmCovariance]                 = None

        self._scan()

    # ------------------------------------------------------------------
    # Scan
    # ------------------------------------------------------------------

    def _scan(self) -> None:
        wd = self.workdir

        # -- detect mode -----------------------------------------------
        if list(wd.glob("*.ws")):
            self.mode = "3d"
        elif list(wd.glob("*.rho")):
            self.mode = "2d"

        # -- log -------------------------------------------------------
        log_files = sorted(wd.glob("*.log"))
        if log_files:
            try:
                self.log = ModEmLog.read(log_files[0])
            except Exception:
                pass

        # -- control ---------------------------------------------------
        inv_files = sorted(wd.glob("*.inv"))
        if inv_files:
            try:
                self.control = ModEmControl.read(inv_files[0])
            except Exception:
                pass

        # -- covariance ------------------------------------------------
        cov_files = sorted(wd.glob("*.cov"))
        if cov_files:
            try:
                self.covariance = ModEmCovariance.read(cov_files[0])
            except Exception:
                pass

        # -- models ----------------------------------------------------
        ext = "*.ws" if self.mode == "3d" else "*.rho"
        model_cls = ModEmModel3D if self.mode == "3d" else ModEmModel2D
        for mf in sorted(wd.glob(ext)):
            st = mf.stem
            # accept m0, m1, …, mi; reject data-model hybrids
            if not _MODEL_STEM_RE.match(_stem(mf)):
                continue
            try:
                m = model_cls.read(mf)
                self.models[st] = m
            except Exception:
                pass

        self.model_initial = self.models.get("m0")
        # prefer "mi" (final), otherwise the highest numbered model
        if "mi" in self.models:
            self.model_final = self.models["mi"]
        elif self.models:
            def _sort_key(k: str) -> int:
                m = re.match(r"m(\d+)", k)
                return int(m.group(1)) if m else -1
            last = max(self.models, key=_sort_key)
            self.model_final = self.models[last]

        # -- data files ------------------------------------------------
        for df in sorted(wd.glob("*.dat")):
            st = _stem(df)
            if st == "d0":
                try:
                    self.data_obs = ModEmData.read(df)
                except Exception:
                    pass
            elif st in ("di", "pr"):
                try:
                    self.data_pred = ModEmData.read(df)
                except Exception:
                    pass

        if self.verbose:
            self.logger.info(
                "InversionResult: mode=%s, %d models, log=%s",
                self.mode, len(self.models), "yes" if self.log else "no",
            )

    # ------------------------------------------------------------------
    # Convenience properties
    # ------------------------------------------------------------------

    @property
    def n_iter(self) -> int:
        """Total number of inversion iterations recorded in the log."""
        return self.log.n_iter if self.log else 0

    @property
    def final_rms(self) -> float:
        """Final data RMS misfit (from log)."""
        return self.log.final_rms if self.log else float("nan")

    @property
    def rms_history(self):
        """RMS misfit per iteration as a numpy array."""
        import numpy as np
        return self.log.rms if self.log else np.array([])

    @property
    def best_rms(self) -> float:
        """Minimum RMS achieved during the inversion."""
        import numpy as np
        h = self.rms_history
        return float(np.min(h)) if h.size else float("nan")

    @property
    def iteration_numbers(self):
        """Iteration index array from the log."""
        import numpy as np
        return self.log.iterations if self.log else np.array([], dtype=int)

    def __repr__(self) -> str:
        return (
            f"InversionResult(mode={self.mode!r}, "
            f"n_iter={self.n_iter}, "
            f"final_rms={self.final_rms:.4f}, "
            f"models={list(self.models)})"
        )
