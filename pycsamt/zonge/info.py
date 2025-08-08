# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
Container that fan-outs an AVG *data block* into the
domain-specific helper classes implemented in the surrounding
sub-packages (``meas``, ``var``, ``resphase``, …).

Typical usage
-------------
>>> from pycsamt.zonge.utils import load_avg
>>> from pycsamt.zonge.info  import DataInfo
>>>
>>> df, meta = load_avg("LCS01.avg")
>>> info     = DataInfo(df)
>>> info.station.value          # → ndarray
>>> info.phase["S07"]           # → per-station Phase slice
>>> info.z.mag.shape            # → (n_freq, n_station, 2, 2)
"""

from __future__ import annotations

from pathlib import Path
from typing  import Any, Dict #, Sequence, List

import numpy as np
import pandas as pd

from ..exceptions  import AvgDataError
from ..log.logger import get_logger

from .meas  import CompMeas, Amps, Frequency
from .survey import Station
from .var import (
    PcEmag, 
    PcHmag, 
    PcRho, 
    SPhz, 
    SHphz, 
    SEphz
)
from .properties import FieldAliases 
from .resphase   import Phase, Resistivity
from .utils import load_avg
from .z import Z

__all__ = ["DataInfo"]

logger = get_logger(__name__)


class DataInfo:
    """Parse a *kind-1/2* AVG **DataFrame** into typed sub-components."""

    def __init__(self,
        frame: pd.DataFrame | None = None,
        *,
        to_degree: bool = False
        ) -> None:

        self.station     : Station      | None = None
        self.frequency   : Frequency    | None = None
        self.amps        : Amps         | None = None
        self.comp        : CompMeas     | None = None
        self.phase       : Phase        | None = None
        self.rho         : Resistivity  | None = None
        self.z           : Z            | None = None

        # quality metrics
        self.pc_emag : PcEmag | None = None
        self.pc_hmag : PcHmag | None = None
        self.pc_rho  : PcRho  | None = None
        self.s_ephz  : SEphz  | None = None
        self.s_hphz  : SHphz  | None = None
        self.s_phz   : SPhz   | None = None

        if frame is not None:
            self._populate_from_frame(frame.copy(), to_degree)


    def _populate_from_frame(
            self, df: pd.DataFrame, to_degree: bool
        ) -> None:
        """Dispatch every recognised column to its helper class."""
        #  mandatory columns 
        try:
            stn = df[self._col("station")].to_numpy()
            frq = df[self._col("freq")   ].to_numpy()
        except KeyError as exc:
            raise AvgDataError(f"Missing essential column: {exc}") from exc

        # grid dimensions
        unique_freq, n_per_stn = np.unique(frq, return_counts=True)
        n_freq     = unique_freq.size
        n_station  = n_per_stn[0] if n_per_stn.size else 0

        # ---- core helpers 
        self.station   = Station(stn, normalize=True)
        self.frequency = Frequency(frq)
        self.comp      = CompMeas(df.get(self._col("component")))

        # amps, E/H field magnitudes & phases 
        if (col := self._maybe(df, "amps")) is not None:
            self.amps = Amps(df[col], n_freq=n_freq, n_station=n_station)

        if (col := self._maybe(df, "phase")) is not None:
            self.phase = Phase(df[col], n_freq=n_freq,
                               n_stations=n_station, to_degree=to_degree)

        if (col := self._maybe(df, "rho")) is not None:
            self.rho = Resistivity(df[col], n_freq=n_freq,
                                   n_stations=n_station)
            # SRes column?
            if (sc := self._maybe(df, "sres")) is not None:
                self.rho.set_sres(df[sc], n_freq=n_freq,
                                  n_stations=n_station)

        # quality / variation metrics 
        self._set_metric(df, "pcemag", PcEmag,  "pc_emag",
                         n_freq, n_station)
        self._set_metric(df, "pchmag", PcHmag,  "pc_hmag",
                         n_freq, n_station)
        self._set_metric(df, "pcrho",  PcRho,   "pc_rho",
                         n_freq, n_station)
        self._set_metric(df, "sphz",   SPhz,    "s_phz",
                         n_freq, n_station, to_degree=to_degree)
        self._set_metric(df, "sephz",  SEphz,   "s_ephz",
                         n_freq, n_station, to_degree=to_degree)
        self._set_metric(df, "shphz",  SHphz,   "s_hphz",
                         n_freq, n_station, to_degree=to_degree)

        # Z tensor (requires the four sub-components) 
        z_cols = {k: self._maybe(df, f"z{k}") for k in ("xx", "xy", "yx", "yy")}
        if all(v is not None for v in z_cols.values()):
            z_stack = {k: df[col] for k, col in z_cols.items()}  # type: ignore[arg-type]
            # reshape to (n_freq, n_station)
            z_arr   = np.stack(
                [v.to_numpy().reshape(n_station, n_freq).T for v in z_stack.values()],
                axis=-1
            ).reshape(n_freq, n_station, 2, 2)
            self.z = Z(z_arr, freq=unique_freq,
                       station_ids=list(self.station.names))

    @classmethod
    def from_avg(
        cls,
        path:        str | Path,
        *,
        to_degree:   bool                    = False,
        load_kwargs: dict[str, Any] | None   = None,
    ) -> "DataInfo":
        """
        Quick shortcut::

            info = DataInfo.from_avg("LCS01.avg", to_degree=True)

        Parameters
        ----------
        path
            Filesystem location of the ``.avg`` file.
        to_degree
            Convert internal phase caches to **degrees**.
        load_kwargs
            Extra keyword arguments forwarded verbatim to
            :func:`pycsamt.zonge.utils.load_avg`
            (e.g. ``ll_columns=('lat', 'lon')`` or ``utm_zone=30``).

        Returns
        -------
        DataInfo
            Fully populated instance.
        """
        
        df, _meta = load_avg(path, **(load_kwargs or {}))
        logger.info("AVG loaded (%d rows) – building DataInfo …", len(df))
        return cls(df, to_degree=to_degree)

    @staticmethod
    def _col(df: pd.DataFrame, alias: str) -> str:
        """
        Return the first concrete column name that exists in *df* and
        matches the FieldAliases entry for *alias*.
        """
        for cand in getattr(FieldAliases, alias, ()):
            if cand in df.columns:            # ← check presence
                return cand
        raise KeyError(f"No column found for alias '{alias}'")
    
    @staticmethod
    def _maybe(df: pd.DataFrame, alias: str) -> str | None:
        for cand in getattr(FieldAliases, alias, ()):
            if cand in df.columns:
                return cand
        return None


    def _set_metric(
        self,
        df:        pd.DataFrame,
        alias:     str,
        cls,                    # type: ignore[valid-type]
        attr_name: str,
        n_freq:    int,
        n_station: int,
        **kws) -> None:
        """Instantiate *cls* if *alias* is available and attach it."""
        if (col := self._maybe(df, alias)) is None:
            return
        obj = cls(df[col], n_freq=n_freq, n_station=n_station, **kws)
        setattr(self, attr_name, obj)

    def __repr__(self) -> str:
        n_stn = len(self.station.names) if self.station else 0
        n_frq = len(self.frequency.value) if self.frequency else 0
        return f"DataInfo(stations={n_stn}, freqs={n_frq})"

    # quick attribute proxy (e.g. info['S03']['rho'])
    def __getitem__(self, station: str) -> Dict[str, Any]:
        return {
            "rho"   : self.rho.loc[station]   if self.rho   else None,
            "phase" : self.phase.loc[station] if self.phase else None,
            "pcE"   : self.pc_emag.loc[station] if self.pc_emag else None,
            "pcH"   : self.pc_hmag.loc[station] if self.pc_hmag else None,
            "sϕZ"   : self.s_phz.loc[station]  if self.s_phz  else None,
            "z"     : self.z[station]          if self.z      else None,
    }

