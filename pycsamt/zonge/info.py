# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
DataInfo - High-level AVG data aggregator.

This module provides the DataInfo class, which serves as a
primary facade for interacting with a complete Zonge AVG dataset.
It composes all other components (Header, Z, Resistivity, Phase,
and various QC metrics) into a single, convenient container.
"""
from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Mapping, Optional, Tuple, Union

import pandas as pd

from ..exceptions import AvgDataError
from .base import AVGFrame
from .heads import Header
from .meas import Amps, CompMeas, Frequency
from .resphase import Phase, Resistivity
from .survey import Station
from .utils import load_avg
from .var_pc import PcEmag, PcHmag, PcRho
from .var_std import SEphz, SHphz, SPhz
from .z import Z


__all__ = ["DataInfo"]


class DataInfo:
    """
    High-level aggregator for a complete Zonge AVG dataset.

    This class provides a unified interface to all parsed data
    components from an AVG file, including header metadata,
    impedance, resistivity, phase, and quality control metrics.
    """

    def __init__(self) -> None:
        # Core data holders
        self._frame: Optional[AVGFrame] = None
        self.df: Optional[pd.DataFrame] = None
        self.meta: Optional[Mapping[str, Any]] = None

        # Component containers
        self.header: Header = Header()
        self.station: Station = Station()
        self.z: Z = Z()
        self.resistivity: Resistivity = Resistivity()
        self.phase: Phase = Phase()
        self.frequency: Frequency = Frequency()
        self.amps: Amps = Amps()
        self.comp: CompMeas = CompMeas()

        # QC components
        self.pc_emag: PcEmag = PcEmag()
        self.pc_hmag: PcHmag = PcHmag()
        self.pc_rho: PcRho = PcRho()
        self.s_ephz: SEphz = SEphz()
        self.s_hphz: SHphz = SHphz()
        self.s_phz: SPhz = SPhz()

    @classmethod
    def from_avg(
        cls,
        avg: Union[
            str, Path, AVGFrame, pd.DataFrame,
            Tuple[pd.DataFrame, Mapping[str, Any]]
        ],
        *,
        meta: Optional[Mapping[str, Any]] = None,
    ) -> "DataInfo":
        """
        Build a DataInfo object from a path, AVGFrame, or DataFrame.
        """
        if isinstance(avg, (str, Path)):
            df, m = load_avg(Path(avg))
            frame = AVGFrame(df, m, Path(avg))
        elif isinstance(avg, AVGFrame):
            frame = avg
        elif isinstance(avg, tuple) and len(avg) == 2:
            df, m = avg
            frame = AVGFrame(df, dict(m))
        elif isinstance(avg, pd.DataFrame):
            frame = AVGFrame(avg, dict(meta or {}))
        else:
            raise TypeError(
                "from_avg expects Path|AVGFrame|DataFrame|"
                "(DataFrame, meta) tuple."
            )

        obj = cls()
        obj.read(frame.data, frame.meta)
        obj._frame = frame
        return obj

    def read(
        self,
        source: pd.DataFrame,
        meta: Optional[Mapping[str, Any]] = None,
    ) -> None:
        """
        Orchestrate reading data into all sub-components.
        """
        self.df = source
        self.meta = meta or {}

        # Populate header from metadata
        self.header.read(meta=self.meta)

        # Populate data components from the DataFrame
        components = [
            self.station, self.z, self.resistivity, self.phase,
            self.frequency, self.amps, self.comp, self.pc_emag,
            self.pc_hmag, self.pc_rho, self.s_ephz,
            self.s_hphz, self.s_phz,
        ]

        for comp in components:
            try:
                comp.read(self.df, self.meta)
            except AvgDataError as e:
                # Gracefully skip components if their data is missing
                # For example, a file might not have %Hmag
                warnings.warn(
                    f"Notice: Could not load component "
                    f"'{comp.__class__.__name__}': {e}"
                )
            except Exception as e:
                warnings.warn(
                    f"Warning: Unexpected error loading "
                    f"'{comp.__class__.__name__}': {e}"
                )

    def __str__(self) -> str:
        if self.df is None:
            return "DataInfo(empty)"

        n_st = (
            self.df["station"].nunique()
            if "station" in self.df.columns
            else 0
        )
        n_f = (
            self.df["freq"].nunique()
            if "freq" in self.df.columns
            else 0
        )
        return (
            f"DataInfo(stations={n_st}, freqs={n_f}, "
            f"rows={len(self.df)})"
        )

    __repr__ = __str__
