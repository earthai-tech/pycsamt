# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""OccamData — build and parse OccamDataFile.dat from EDI sources.

The data file is the primary Occam2D input.  It encodes apparent
resistivity, phase (and optionally tipper) for selected TE/TM modes
at each station and frequency.

Entry points
------------
``OccamData.from_edi(edi_collection_or_sites, ...)``
    Build a new data file object directly from EDI data.
``OccamData.read(path)``
    Parse an existing ``OccamDataFile.dat``.
``OccamData.write(path)``
    Serialise to disk in the ``OCCAM2MTDATA_1.0`` format.

File format (OCCAM2MTDATA_1.0)
-------------------------------
FORMAT:           OCCAM2MTDATA_1.0
TITLE:            <title>
SITES:            <n>
   S00 ... S<n-1>
OFFSETS (M):
   <x0> ... <x_{n-1}>
FREQUENCIES:      <nf>
   <f0> ... <f_{nf-1}>
DATA BLOCKS:      <nd>
   <site_idx> <freq_idx> <type_code> <value> <error>
   ...
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Union

import numpy as np

from .base   import OccamBase
from .config import OccamConfig

PathLike = Union[str, Path]

__all__ = ["OccamData"]

# Data-type codes used in the OCCAM2MTDATA format
DATA_TYPE_CODES = {
    "RhoTE":   1,
    "PhsTE":   2,
    "RhoTM":   5,
    "PhsTM":   6,
    "ReZXX":   11,
    "ImZXX":   12,
    "ReZXY":   13,
    "ImZXY":   14,
    "ReZYX":   15,
    "ImZYX":   16,
    "ReZYY":   17,
    "ImZYY":   18,
    "ReTZX":   3,
    "ImTZX":   4,
    "ReTZY":   7,
    "ImTZY":   8,
}


class OccamData(OccamBase):
    """Occam2DMT data file container.

    Attributes
    ----------
    title : str
        Free-text title written into the file header.
    sites : list[str]
        Ordered station names.
    offsets : np.ndarray, shape (n_sites,)
        Horizontal profile offsets in metres.
    frequencies : np.ndarray, shape (n_freq,)
        Frequencies in Hz (descending order, as Occam expects).
    data_blocks : np.ndarray, shape (n_data, 5)
        Columns: site_index (1-based), freq_index (1-based), type_code,
        value, error.
    """

    def __init__(
        self,
        title: str = "pycsamt Occam2D data file",
        config: Optional[OccamConfig] = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.title      = title
        self.config     = config or OccamConfig()
        self.sites:       List[str]    = []
        self.offsets:     np.ndarray   = np.array([])
        self.frequencies: np.ndarray   = np.array([])
        self.data_blocks: np.ndarray   = np.empty((0, 5))

    # ------------------------------------------------------------------
    # Construction from EDI
    # ------------------------------------------------------------------
    @classmethod
    def from_edi(
        cls,
        source,
        modes: Optional[List[str]] = None,
        config: Optional[OccamConfig] = None,
        title: str = "pycsamt Occam2D data file",
        **kwargs,
    ) -> "OccamData":
        """Build an ``OccamData`` from an ``EDICollection`` or ``Sites``.

        Parameters
        ----------
        source : EDICollection | Sites
            Loaded EDI data (single source of truth).
        modes : list[str], optional
            Subset of ``["TE", "TM", "Tipper"]``.  Default from config.
        config : OccamConfig, optional
            Run configuration (error floors, frequency band, …).
        """
        # TODO: implement
        # 1. Extract TFBundle / z arrays from source
        # 2. Apply frequency band and error-floor filters from config
        # 3. Compute profile offsets from station coordinates
        # 4. Populate self.sites, self.offsets, self.frequencies, self.data_blocks
        raise NotImplementedError("OccamData.from_edi — not yet implemented")

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @classmethod
    def read(cls, path: PathLike, **kwargs) -> "OccamData":
        """Parse an existing ``OccamDataFile.dat``."""
        # TODO: implement parsing of OCCAM2MTDATA_1.0 format
        raise NotImplementedError("OccamData.read — not yet implemented")

    def write(self, path: PathLike) -> Path:
        """Write the data file to *path* in OCCAM2MTDATA_1.0 format.

        Returns the resolved path.
        """
        # TODO: implement
        raise NotImplementedError("OccamData.write — not yet implemented")

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------
    @property
    def n_sites(self) -> int:
        return len(self.sites)

    @property
    def n_frequencies(self) -> int:
        return len(self.frequencies)

    @property
    def n_data(self) -> int:
        return len(self.data_blocks)
