# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Result containers for completed ModEM inversion runs."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from .base import ModEmBase
from .config import ModEmConfig
from .control import ModEmControl
from .covariance import ModEmCovariance
from .data import ModEmData
from .doc import _modem_param_docs as _params
from .log import ModEmLog
from .model2d import ModEmModel2D
from .model3d import ModEmModel3D

PathLike = str | Path
ModelObject = ModEmModel2D | ModEmModel3D

__all__ = ["InversionResult"]

# Stem patterns ModEM uses for model files
_MODEL_STEM_RE = re.compile(r"^m(\d+|i)$", re.IGNORECASE)


def _stem(p: Path) -> str:
    """Return a ModEM stem stripped of known generated suffixes."""
    return p.stem.split("_")[0]


class InversionResult(ModEmBase):
    def __init__(
        self,
        workdir: PathLike,
        config: ModEmConfig | None = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.workdir: Path = Path(workdir)
        if not self.workdir.is_dir():
            msg = f"ModEM workdir not found: {self.workdir}"
            raise FileNotFoundError(msg)

        self.config: ModEmConfig = config or ModEmConfig()

        self.mode: str = "unknown"
        self.log: ModEmLog | None = None
        self.control: ModEmControl | None = None
        self.model_initial: ModelObject | None = None
        self.model_final: ModelObject | None = None
        self.models: dict[str, ModelObject] = {}
        self.data_obs: ModEmData | None = None
        self.data_pred: ModEmData | None = None
        self.covariance: ModEmCovariance | None = None

        self._scan()

    # ------------------------------------------------------------------
    # Scan
    # ------------------------------------------------------------------

    def _scan(self) -> None:
        """Discover and load available ModEM output artefacts."""
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
            # Accept m0, m1, and mi; reject data-model hybrids.
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
                self.mode,
                len(self.models),
                "yes" if self.log else "no",
            )

    # ------------------------------------------------------------------
    # Convenience properties
    # ------------------------------------------------------------------

    @property
    def n_iter(self) -> int:
        """Total number of inversion iterations parsed from the log."""
        return self.log.n_iter if self.log else 0

    @property
    def final_rms(self) -> float:
        """Final root-mean-square data misfit reported by ModEM.

        Returns
        -------
        float
            Last parsed RMS value. If no log was loaded, the value is
            ``nan`` so callers can use :func:`numpy.isfinite` to check
            availability.
        """
        return self.log.final_rms if self.log else float("nan")

    @property
    def rms_history(self) -> np.ndarray:
        """RMS misfit values for each parsed inversion iteration."""
        return self.log.rms if self.log else np.array([])

    @property
    def best_rms(self) -> float:
        """Minimum RMS value reached by the inversion history."""
        h = self.rms_history
        return float(np.min(h)) if h.size else float("nan")

    @property
    def iteration_numbers(self) -> np.ndarray:
        """Iteration numbers associated with ``rms_history``."""
        if self.log:
            return self.log.iterations
        return np.array([], dtype=int)

    def __repr__(self) -> str:
        return (
            f"InversionResult(mode={self.mode!r}, "
            f"n_iter={self.n_iter}, "
            f"final_rms={self.final_rms:.4f}, "
            f"models={list(self.models)})"
        )


InversionResult.__doc__ = rf"""
Aggregate the files produced by a ModEM inversion run.

``InversionResult`` scans one ModEM working directory, detects
whether it contains two-dimensional or three-dimensional model
files, and loads every recognized artefact into a single
Python object. It is intended for post-processing, plotting,
quality control, and workflow scripts that need a consistent
view of logs, control files, models, data, and covariance
settings after a run has finished.

The main scalar quality measure exposed by the class is the
root-mean-square data misfit. For :math:`N` weighted data
residuals, the value reported by ModEM is commonly read as

.. math::
   \mathrm{{RMS}} =
   \sqrt{{\frac{{1}}{{N}}
   \sum_{{i=1}}^N
   \left(\frac{{d_i^{{obs}} - d_i^{{pred}}}}
   {{\sigma_i}}\right)^2}}.

Here :math:`d_i^{{obs}}` and :math:`d_i^{{pred}}` are observed
and predicted data values, and :math:`\sigma_i` is the data
uncertainty used for weighting. The exact composition of the
sum depends on the component family selected in the ModEM data
file [1]_.

Parameters
----------
{_params.common.workdir}
{_params.common.config}
**kwargs : dict
    Additional keyword arguments passed to
    :class:`~pycsamt.models.modem.base.ModEmBase`. Typical
    values include logging and verbosity options shared by
    the ModEM wrapper classes.

Attributes
----------
workdir : pathlib.Path
    Resolved directory scanned for ModEM artefacts. The
    directory must exist before the result object is created.
mode : {{"2d", "3d", "unknown"}}
    Detected inversion dimensionality. A directory containing
    ``*.ws`` files is treated as 3-D, a directory containing
    ``*.rho`` files is treated as 2-D, and an empty or
    unrecognized directory remains ``"unknown"``.
log : ModEmLog or None
    Parsed run log containing iteration numbers, objective
    values, model norms, trade-off parameters, and RMS
    history. The attribute is ``None`` when no readable log
    file is present.
control : ModEmControl or None
    Parsed inversion-control file, usually read from the first
    ``*.inv`` file found in ``workdir``. It records nonlinear
    solver limits, target misfit, lambda controls, and output
    naming settings.
{_params.result.model_initial}
{_params.result.model_final}
{_params.result.models}
{_params.result.data_obs}
{_params.result.data_pred}
covariance : ModEmCovariance or None
    Parsed covariance file for 3-D runs. The object contains
    smoothing weights, smoothing iteration count, and active
    masks when a readable ``*.cov`` file is available.

Notes
-----
The scanner is deliberately tolerant. If a recognized file is
missing or cannot be parsed, the corresponding attribute is
left as ``None`` and scanning continues. This makes the class
useful for inspecting incomplete, interrupted, or partially
copied run directories.

Model files are keyed by their file stem. The initial model is
selected from ``m0`` when present. The final model is selected
from ``mi`` when present; otherwise the highest numbered
``m<N>`` model is used.

Examples
--------
Load a finished inversion directory and inspect the misfit:

>>> from pycsamt.models.modem.results import InversionResult
>>> result = InversionResult("modem_run")
>>> result.mode in {{"2d", "3d", "unknown"}}
True
>>> float(result.final_rms) == result.final_rms
True

Access model and data products when they were loaded:

>>> if result.model_final is not None:
...     rho = result.model_final.rho_linear
...     rho.ndim in (2, 3)
... else:
...     rho = None

Compare the best and final RMS values:

>>> if result.rms_history.size:
...     result.best_rms <= result.final_rms
... else:
...     result.n_iter == 0
True

See Also
--------
ModEmLog
    Parser for ModEM iteration logs and RMS histories.
ModEmData
    Reader and writer for observed and predicted data files.
ModEmModel2D
    Two-dimensional ModEM resistivity model container.
ModEmModel3D
    Three-dimensional ModEM resistivity model container.
ModEmControl
    Reader and writer for ModEM inversion-control files.
ModEmCovariance
    Reader and writer for 3-D covariance smoothing files.

References
----------
.. [1] Kelbert, A., Meqbel, N., Egbert, G. D., and
   Tandon, K. (2014). ModEM: A modular system for inversion
   of electromagnetic geophysical data. *Computers &
   Geosciences*, 66, 40-53.
   doi:10.1016/j.cageo.2014.01.010.
"""
