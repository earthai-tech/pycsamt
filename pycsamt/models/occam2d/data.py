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
::

    FORMAT:           OCCAM2MTDATA_1.0
    TITLE:            <title>
    SITES:            <n>
       S00
       ...
       S<n-1>
    OFFSETS (M):
       <x0>
       ...
       <x_{n-1}>
    FREQUENCIES:      <nf>
       <f0>
       ...
       <f_{nf-1}>
    DATA BLOCKS:      <nd>
    SITE  FREQ  TYPE   DATUM    ERROR
       <site_idx>  <freq_idx>  <type_code>  <value>  <error>
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

_FORMAT_TAG = "OCCAM2MTDATA_1.0"

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


# -----------------------------------------------------------------------
# Low-level parser
# -----------------------------------------------------------------------

def _parse_data(path: Path) -> dict:
    """Parse an OCCAM2MTDATA_1.0 file.

    Returns
    -------
    dict with keys:
        format_str, title, sites (list[str]), offsets (list[float]),
        frequencies (list[float]), data_rows (list[list]) where each
        row is [site_idx, freq_idx, type_code, datum, error].

    Raises
    ------
    FileNotFoundError
    ValueError
        If the Format tag is absent or wrong.
    """
    if not path.exists():
        raise FileNotFoundError(f"OccamDataFile not found: {path}")

    with path.open("r", errors="replace") as fh:
        lines = [l.rstrip("\n") for l in fh]

    result: dict = {
        "format_str":  None,
        "title":       "",
        "sites":       [],
        "offsets":     [],
        "frequencies": [],
        "data_rows":   [],
    }

    i = 0
    N = len(lines)
    n_sites = 0
    n_freq  = 0

    while i < N:
        raw     = lines[i]
        stripped = raw.strip()
        i       += 1

        if not stripped:
            continue

        if ":" not in stripped:
            continue

        raw_key, _, raw_val = stripped.partition(":")
        key_up = raw_key.strip().upper()
        val    = raw_val.strip()

        if key_up == "FORMAT":
            result["format_str"] = val
            if val.upper() != _FORMAT_TAG:
                raise ValueError(
                    f"Expected format '{_FORMAT_TAG}', got '{val}' in {path}"
                )

        elif key_up == "TITLE":
            result["title"] = val

        elif key_up == "SITES":
            n_sites = int(val)
            sites: list[str] = []
            while len(sites) < n_sites and i < N:
                s = lines[i].strip()
                i += 1
                if s:
                    sites.append(s)
            result["sites"] = sites

        elif key_up.startswith("OFFSETS"):
            offsets: list[float] = []
            while len(offsets) < n_sites and i < N:
                s = lines[i].strip()
                i += 1
                if not s:
                    continue
                try:
                    offsets.append(float(s))
                except ValueError:
                    i -= 1   # not a float → put line back
                    break
            result["offsets"] = offsets

        elif key_up == "FREQUENCIES":
            n_freq = int(val)
            freqs: list[float] = []
            while len(freqs) < n_freq and i < N:
                s = lines[i].strip()
                i += 1
                if not s:
                    continue
                try:
                    freqs.append(float(s))
                except ValueError:
                    i -= 1
                    break
            result["frequencies"] = freqs

        elif key_up.startswith("DATA BLOCKS"):
            # Skip the column-header line ("SITE  FREQ  TYPE   DATUM    ERROR")
            while i < N:
                s = lines[i].strip()
                i += 1
                if s and s.upper().startswith("SITE"):
                    break
            # Read data rows until end of file
            data_rows: list[list] = []
            while i < N:
                s = lines[i].strip()
                i += 1
                if not s:
                    continue
                tokens = s.split()
                if len(tokens) == 5:
                    try:
                        data_rows.append([
                            int(tokens[0]), int(tokens[1]), int(tokens[2]),
                            float(tokens[3]), float(tokens[4]),
                        ])
                    except ValueError:
                        pass
            result["data_rows"] = data_rows

    if result["format_str"] is None:
        raise ValueError(
            f"File does not contain a valid OCCAM2MTDATA_1.0 header: {path}"
        )

    return result


# -----------------------------------------------------------------------
# OccamData
# -----------------------------------------------------------------------

class OccamData(OccamBase):
    """Occam2DMT data file container.

    Attributes
    ----------
    format_str : str
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
        self.format_str:  str          = _FORMAT_TAG
        self.title:       str          = title
        self.config:      OccamConfig  = config or OccamConfig()
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
        """Parse an existing ``OccamDataFile.dat``.

        Parameters
        ----------
        path : path-like
            Path to the data file.

        Returns
        -------
        OccamData

        Raises
        ------
        FileNotFoundError
            If *path* does not exist.
        ValueError
            If the Format tag is absent or wrong.
        """
        p      = Path(path)
        parsed = _parse_data(p)
        obj    = cls(**kwargs)

        obj.format_str  = parsed["format_str"]
        obj.title       = parsed["title"]
        obj.sites       = parsed["sites"]
        obj.offsets     = np.array(parsed["offsets"], dtype=float)
        obj.frequencies = np.array(parsed["frequencies"], dtype=float)

        rows = parsed["data_rows"]
        if rows:
            obj.data_blocks = np.array(rows, dtype=float)
        else:
            obj.data_blocks = np.empty((0, 5))

        if obj.verbose:
            obj.logger.info(
                "OccamData.read: %d sites, %d freqs, %d data blocks from %s",
                obj.n_sites, obj.n_frequencies, obj.n_data, p,
            )
        return obj

    def write(self, path: PathLike) -> Path:
        """Write the data file to *path* in OCCAM2MTDATA_1.0 format.

        Returns the resolved path.
        """
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)

        lines: list[str] = [
            f"FORMAT:           {self.format_str}\n",
            f"TITLE:            {self.title}\n",
            f"SITES:            {self.n_sites}\n",
        ]
        for s in self.sites:
            lines.append(f"   {s}\n")

        lines.append("OFFSETS (M):\n")
        for off in self.offsets:
            lines.append(f"   {off:.1f}\n")

        lines.append(f"FREQUENCIES:      {self.n_frequencies}\n")
        for f in self.frequencies:
            lines.append(f"   {f:.7e}\n")

        lines.append(f"DATA BLOCKS:      {self.n_data}\n")
        lines.append("SITE  FREQ  TYPE   DATUM    ERROR\n")
        for row in self.data_blocks:
            s_i = int(row[0])
            f_i = int(row[1])
            t_c = int(row[2])
            dat = row[3]
            err = row[4]
            lines.append(f"{s_i:5d} {f_i:5d} {t_c:5d} {dat:12.4f} {err:12.4f}\n")

        with p.open("w") as fh:
            fh.writelines(lines)
        return p

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

    @property
    def type_codes(self) -> np.ndarray:
        """Unique data-type codes present in this dataset."""
        if not self.data_blocks.size:
            return np.array([], dtype=int)
        return np.unique(self.data_blocks[:, 2].astype(int))
