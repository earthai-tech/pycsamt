# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Loader for MARE2DEM inversion output files."""

from __future__ import annotations

from pathlib import Path

from .base import Mare2DEMBase
from .config import Mare2DEMConfig
from .data import EMData, read_emdata
from .log import Mare2DEMLog
from .mesh import ResistivityModel, read_resistivity
from .validation import (
    is_emdata_file, is_log_file, is_resistivity_file, is_response_file,
    Mare2DEMFileType, detect_file_type,
)

__all__ = ["InversionResult"]


class InversionResult(Mare2DEMBase):
    """Load and expose MARE2DEM inversion output files.

    ``InversionResult`` scans a MARE2DEM run directory after the
    binary has finished and loads the iteration log, final
    resistivity model, observed-data file, and predicted-response
    file.

    Parameters
    ----------
    workdir : path-like
        Directory to scan for MARE2DEM output files.
    config : Mare2DEMConfig, optional
        Configuration providing default file stems. When omitted,
        the scanner looks for any ``.log``, ``.resistivity``,
        ``.emdata``, and ``*_MARE2DEM.emdata`` files.
    **kwargs :
        Forwarded to :class:`Mare2DEMBase`.

    Attributes
    ----------
    workdir : pathlib.Path
        Absolute path of the scanned run directory.
    config : Mare2DEMConfig
        Configuration used for file-name hints.
    log : Mare2DEMLog or None
        Parsed iteration log.
    model : ResistivityModel or None
        Final inverted resistivity model.
    data : EMData or None
        Observed data file.
    response : EMData or None
        Predicted-response file (``*_MARE2DEM.emdata``).

    Examples
    --------
    >>> from pycsamt.models.mare2dem import InversionResult
    >>> result = InversionResult("./mare2dem_run")
    >>> result.log.final_rms
    0.98
    >>> result.log.converged
    True
    """

    def __init__(
        self,
        workdir: str | Path,
        config: Mare2DEMConfig | None = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.workdir: Path = Path(workdir).resolve()
        self.config: Mare2DEMConfig = config or Mare2DEMConfig()

        self.log: Mare2DEMLog | None = None
        self.model: ResistivityModel | None = None
        self.data: EMData | None = None
        self.response: EMData | None = None

        self._scan()

    # ------------------------------------------------------------------
    # Internal scanning
    # ------------------------------------------------------------------

    def _scan(self) -> None:
        """Scan *workdir* and populate output attributes."""
        if not self.workdir.exists():
            if self.verbose:
                self.logger.warning(
                    "InversionResult: workdir does not exist: %s",
                    self.workdir,
                )
            return

        all_files = list(self.workdir.iterdir())

        for f in all_files:
            if not f.is_file():
                continue
            if is_log_file(f) and self.log is None:
                self.log = Mare2DEMLog(f)
            elif is_response_file(f) and self.response is None:
                self.response = read_emdata(f)
            elif is_emdata_file(f) and self.data is None:
                self.data = read_emdata(f)
            elif is_resistivity_file(f) and self.model is None:
                self.model = read_resistivity(f)

        if self.verbose:
            self.logger.info(
                "InversionResult: scanned %s — log=%s, model=%s, "
                "data=%s, response=%s",
                self.workdir,
                self.log is not None,
                self.model is not None,
                self.data is not None,
                self.response is not None,
            )

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------

    @property
    def converged(self) -> bool:
        """``True`` when the log reports successful convergence."""
        return bool(self.log and self.log.converged)

    @property
    def final_rms(self) -> float | None:
        """Final normalized RMS from the log, or ``None``."""
        return self.log.final_rms if self.log else None

    @property
    def n_iterations(self) -> int:
        """Number of completed inversion iterations."""
        return self.log.n_iterations if self.log else 0

    def summary(self) -> str:
        """Return a human-readable summary string."""
        lines = [
            f"InversionResult({self.workdir})",
            f"  converged   : {self.converged}",
            f"  final RMS   : {self.final_rms}",
            f"  n_iterations: {self.n_iterations}",
            f"  model       : {self.model}",
            f"  data        : {self.data}",
            f"  response    : {self.response}",
        ]
        return "\n".join(lines)

    def print_summary(self) -> None:
        """Print :meth:`summary` to stdout."""
        print(self.summary())

    def __repr__(self) -> str:
        return (
            f"InversionResult(workdir={self.workdir}, "
            f"converged={self.converged}, "
            f"final_rms={self.final_rms})"
        )
