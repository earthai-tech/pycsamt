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
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import pandas as pd

from ..exceptions import AvgDataError
from .base import AVGFrame
from .config import Zonge
from .heads import Header
from .meas import Amps, CompMeas, Frequency
from .resphase import Phase, Resistivity
from .survey import Station
from .utils import _standardise_columns, load_avg
from .var_pc import PcEmag, PcHmag, PcRho
from .var_std import SEphz, SHphz, SPhz
from .z import Z

__all__ = ["DataInfo"]


class DataInfo(Zonge):
    r"""High-level aggregator for a complete Zonge AVG dataset.

    This class acts as the primary container and orchestrator for
    all data and metadata parsed from a Zonge AVG file. It
    composes all other data components (e.g., Header, Z,
    Resistivity, Phase, and QC metrics) into a single,
    convenient object.

    Its main role is to provide a unified interface, holding all
    the structured information in one place after the initial
    parsing is complete.

    Attributes
    ----------
    df : pandas.DataFrame or None
        The core tidy DataFrame containing all available data
        columns after parsing and standardization.
    meta : mapping or None
        The raw metadata dictionary extracted from the file's
        header section.
    header : :class:`~.heads.Header`
        A component that aggregates all header-level metadata.
    station : :class:`~.survey.Station`
        A component that manages survey line geometry.
    z : :class:`~.z.Z`
        The component for computing the complex impedance tensor.
    resistivity, phase : :class:`~.resphase.Resistivity`, :class:`~.resphase.Phase`
        Components for apparent resistivity and phase data.
    frequency, amps, comp : :class:`~.meas.Frequency`, etc.
        Components for core measurement quantities.
    pc_emag, pc_hmag, pc_rho : :class:`~.var_pc.PcEmag`, etc.
        Components for percent-error quality control metrics.
    s_ephz, s_hphz, s_phz : :class:`~.var_std.SEphz`, etc.
        Components for phase standard deviation QC metrics.

    Methods
    -------
    from_avg(avg, meta=None)
        A classmethod factory to build a `DataInfo` object from
        various sources, including a file path or DataFrame.
    read(source, meta=None)
        Orchestrates the population of all sub-components from a
        standardized DataFrame and metadata dictionary.

    Notes
    -----
    The `read` method is the core of this class. It iterates
    through all its component attributes and calls their
    respective `read` methods. It includes a robust error-
    handling mechanism that will issue a warning and skip any
    component that fails to load (e.g., due to missing data in
    the source file), making the loading process resilient.

    Examples
    --------
    While typically used internally by the `AVG` class, you could
    use `DataInfo` directly:

    >>> from pycsamt.zonge.info import DataInfo
    >>> from pycsamt.zonge.utils import load_avg
    >>> df, meta = load_avg("data/avg/K2.avg")
    >>> data_info = DataInfo()
    >>> data_info.read(df, meta)
    >>> print(data_info.station)
    Station(n=28, span=25.0–1375.0 m, inc=50.0)

    See Also
    --------
    pycsamt.zonge.avg.AVG : The main user-facing class that uses
        `DataInfo`.
    pycsamt.zonge.base.AVGComponentBase : The base class for all
        components held by `DataInfo`.
    """

    def __init__(self, verbose: bool = False) -> None:
        super().__init__(verbose=verbose)

        # Core data holders
        self._frame: AVGFrame | None = None
        self.df: pd.DataFrame | None = None
        self.meta: Mapping[str, Any] | None = None

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
        avg: str
        | Path
        | AVGFrame
        | pd.DataFrame
        | tuple[pd.DataFrame, Mapping[str, Any]],
        *,
        meta: Mapping[str, Any] | None = None,
    ) -> DataInfo:
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
        meta: Mapping[str, Any] | None = None,
    ) -> None:
        """
        Orchestrate reading data into all sub-components.
        """
        self.df = _standardise_columns(source.copy())
        if "comp" not in self.df.columns:
            self.df["comp"] = "ExHy"
        self.meta = meta or {}

        # Populate header from metadata
        self.header.read(meta=self.meta)

        # Populate data components from the DataFrame
        components = [
            self.station,
            self.z,
            self.resistivity,
            self.phase,
            self.frequency,
            self.amps,
            self.comp,
            self.pc_emag,
            self.pc_hmag,
            self.pc_rho,
            self.s_ephz,
            self.s_hphz,
            self.s_phz,
        ]

        for comp in components:
            try:
                comp.read(self.df, self.meta)
            except AvgDataError as e:
                # Gracefully skip components if their data is missing
                msg = (
                    f"Notice: Could not load component "
                    f"'{comp.__class__.__name__}': {e}"
                )
                if self.verbose:
                    self._logger.warning(msg)
                else:
                    warnings.warn(msg, stacklevel=2)

            except Exception as e:
                msg = (
                    f"Warning: Unexpected error loading "
                    f"'{comp.__class__.__name__}': {e}"
                )
                if self.verbose:
                    self._logger.error(msg)
                else:
                    warnings.warn(msg, stacklevel=2)

        return self

    def __str__(self) -> str:
        if self.df is None:
            return "DataInfo(empty)"

        n_st = (
            self.df["station"].nunique() if "station" in self.df.columns else 0
        )
        n_f = self.df["freq"].nunique() if "freq" in self.df.columns else 0
        return f"DataInfo(stations={n_st}, freqs={n_f}, rows={len(self.df)})"

    __repr__ = __str__
