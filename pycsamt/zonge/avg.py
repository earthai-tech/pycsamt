# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later
"""
AVG - Main user-facing data container for Zonge datasets.

This module provides the `AVG` class, which serves as the
primary entry point and high-level facade for interacting with a
complete Zonge AVG dataset. It composes all other components
(Header, Z, Resistivity, Phase, and QC metrics) into a single,
convenient container.
"""
from __future__ import annotations

import warnings
from pathlib import Path
from typing import (
    Any,
    Mapping,
    Optional,
    Tuple,
    Union,
    Literal, 
    Dict
)

import numpy as np
import pandas as pd

from ..api.bunch import Bunch 
from ..constants import MU_0, PI
from ..decorators import has_fit
from ..log.logger import get_logger
from ..utils.deps import ensure_pkg 
from ..utils.validation import has_read 
from ._transfer import LegacyAVGBase
from .base import AVGFrame
from .info import DataInfo
from .utils import ( 
    classify_avg_format, 
    load_avg, 
    write_avg
)


__all__ = ["BaseAVG","AVG", "AMTAVG"]

logger = get_logger(__name__)

class BaseAVG:
    r"""Base class for AVG data handling and file writing.

    This class serves as the core engine for reading, parsing,
    and writing Zonge AVG data files. It handles the logic for
    identifying file types (legacy vs. modern), transforming
    legacy data, and populating a structured data model.

    Parameters
    ----------
    verbose : bool, default False
        If ``True``, log messages will be printed to the console
        during file reading and writing operations, providing
        insight into the process.

    Attributes
    ----------
    info : DataInfo
        An aggregator object that holds all the individual data
        components (e.g., Header, Station, Z, Resistivity).
    verbose : bool
        Controls the verbosity of logging output.
    _kind : {1, 2} or None
        Indicates the detected file format: 1 for legacy, 2 for
        modern. ``None`` if no file has been read.
    _source_path : pathlib.Path or None
        The path to the source file that was read.

    Notes
    -----
    This class is typically not instantiated directly by users.
    Instead, the :class:`~pycsamt.zonge.AVG` class, which
    inherits from `BaseAVG`, is the primary user-facing entry
    point.

    The `read` method is the main entry point for data ingestion
    and is designed to be flexible, accepting file paths or
    in-memory data structures.

    Examples
    --------
    While direct use is uncommon, one could use `BaseAVG` like so:

    >>> from pycsamt.zonge.avg import BaseAVG
    >>> avg_processor = BaseAVG(verbose=True)
    >>> avg_processor.read('path/to/your/data.avg')
    >>> avg_processor.write('path/to/output.avg', fmt='modern')

    See Also
    --------
    AVG : The primary user-facing class for AVG data.
    DataInfo : The main component aggregator used by this class.
    """
    def __init__(self, verbose: bool = False):
        self.info: DataInfo = DataInfo()
        self.verbose: bool = verbose
        self._kind: Optional[int] = None
        self._source_path: Optional[Path] = None

    def read(
        self,
        source: Union[str, Path, pd.DataFrame, AVGFrame],
        meta: Optional[Mapping[str, Any]] = None,
    ):
        """
        Read, parse, and populate components from a data source.
        """
        df: pd.DataFrame
        final_meta: Mapping[str, Any]

        if isinstance(source, (str, Path)):
            self._source_path = Path(source)
            if self.verbose:
                logger.info(
                    f"Reading from file: {self._source_path}"
                )
            lines = self._source_path.read_text(
                errors="replace"
            ).splitlines()
            self._kind = classify_avg_format(lines)
            df, final_meta = load_avg(self._source_path)

            if self._kind == 1:
                if self.verbose:
                    logger.info("Transforming legacy data.")
                transformer = LegacyAVGBase()
                ds = transformer.from_dataframe(df, meta=final_meta)
                df = ds.to_dataframe().reset_index()
                final_meta = ds.attrs

        elif isinstance(source, pd.DataFrame):
            if self.verbose:
                logger.info("Reading from pandas DataFrame.")
            df = source
            final_meta = meta or {}
            self._kind = 2  # Assume modern if from DataFrame

        elif isinstance(source, AVGFrame):
            if self.verbose:
                logger.info("Reading from AVGFrame object.")
            df = source.data
            final_meta = source.meta
            self._source_path = source.source
            self._kind = 2  # Assume modern

        else:
            raise TypeError(
                "Unsupported source type. Expected Path, str, "
                "DataFrame, or AVGFrame."
            )

        if self.verbose:
            logger.info("Populating data components...")
        self.info.read(df, final_meta)
        if self.verbose:
            logger.info("AVG data successfully loaded.")
        
        return self 
    
    def to_modern(
        self,
        path: Optional[Union[str, Path]] = None,
        *,
        stamp: bool = True,
        float_fmt: str = "%.6g",
        na_rep: str = "",
        header_spaces: bool = False,
    ):
        """
        Write data to a modern (kind-2) AVG file.
        """
        if self.info.df is None:
            raise ValueError("Data frame is not loaded.")
        
        has_read (self, attributes="info") 
        
        if path is None:
            base = (
                self._source_path.stem
                if self._source_path
                else "export"
            )
            path = Path.cwd() / f"{base}_modern.avg"

        if self.verbose:
            logger.info(f"Writing modern AVG file to: {path}")

        meta = self.info.header.to_keywords()

        write_avg(
            core=self.info.df,
            extra=None,
            meta=meta,
            path=path,
            stamp=stamp,
            float_fmt=float_fmt,
            na_rep=na_rep,
            header_spaces=header_spaces,
        )
        
    def to_legacy(
        self,
        path: Optional[Union[str, Path]] = None,
        *,
        precision: int = 4,
        na_rep: str = "*",
    ):
        """
        Write data to a legacy (kind-1) AVG file.
        """
        if self.info.df is None:
            raise ValueError("Data frame is not loaded.")
        
        has_read (self, attributes="info") 
        
        if path is None:
            base = (
                self._source_path.stem
                if self._source_path
                else "export"
            )
            path = Path.cwd() / f"{base}_legacy.avg"

        if self.verbose:
            logger.info(f"Writing legacy AVG file to: {path}")

        lines = []
        h = self.info.header
        lines.append(
            f"\\ AMTAVG 7.76: "
            f"\"{h.hardware.source_file or 'pycsamt.export'}\", "
            f"Dated {h.hardware.dated or '...'}, "
            f"Processed {h.hardware.processed or '...'}"
        )
        if h.rx.length_m is not None:
            lines.append(f"$ ASPACE= {h.rx.length_m}m")
        if h.tx.gdp_station is not None:
            lines.append(f"$ XMTR  = {h.tx.gdp_station:>5.0f}.")

        header_str = (
            " skp Station Freq  Comp Amps     Emag     Ephz      "
            "Hmag     Hphz  Resistivity   Phase  %Emag  sEphz  "
            "%Hmag  sHphz   %Rho   sPhz"
        )
        separator_str = (
            "\\-++------++----++---++----++---------++------++"
            "---------++------++---------++------++-----++-----++"
            "-----++-----++-----++-----+"
        )
        lines.append(header_str)
        lines.append(separator_str)

        df = self.info.df.sort_values(by=["station", "freq"])

        # Define format specifier to avoid linter confusion
        sci_fspec = f"%9.{precision}e"

        def _format_val(val, fspec, na_char):
            width = len(fspec % 0)
            if pd.isna(val):
                return na_char.center(width)
            return fspec % val
        
        for _, row in df.iterrows():
            skp = 2 if row.get("use", True) else 1

            line = (
                f" {skp:1d} "
                f"{_format_val(row.get('station'), '%7.1f', na_rep)} "
                f"{_format_val(row.get('freq'), '%5.0f', na_rep)} "
                f"{str(row.get('comp', 'ExHy')):<4s} "
                f"{_format_val(row.get('amps'), '%4.1f', na_rep)}   "
                f"{_format_val(row.get('emag'), sci_fspec, na_rep)} "
                f"{_format_val(row.get('ephz'), '%7.1f', na_rep)}  "
                f"{_format_val(row.get('hmag'), sci_fspec, na_rep)} "
                f"{_format_val(row.get('hphz'), '%7.1f', na_rep)}  "
                f"{_format_val(row.get('rho'), sci_fspec, na_rep)} "
                f"{_format_val(row.get('phase'), '%7.1f', na_rep)}   "
                f"{_format_val(row.get('pc_emag'), '%4.1f', na_rep)}  "
                f"{_format_val(row.get('s_ephz'), '%5.1f', na_rep)}  "
                f"{_format_val(row.get('pc_hmag'), '%4.1f', na_rep)}  "
                f"{_format_val(row.get('s_hphz'), '%5.1f', na_rep)}   "
                f"{_format_val(row.get('pc_rho'), '%4.1f', na_rep)}  "
                f"{_format_val(row.get('s_phz'), '%5.1f', na_rep)}"
            )
            lines.append(line)

        Path(path).write_text("\n".join(lines))

    def write(
        self,
        path: Optional[Union[str, Path]] = None,
        *,
        fmt: Literal[
            "kind1", "kind2", "legacy", "modern", "auto"
        ] = "auto",
        **kwargs,
    ):
        """
        Write AVG data to a file in the specified format.
        """

        fmt_lower = fmt.lower()
        if fmt_lower in ("legacy", "kind1"):
            self.to_legacy(path, **kwargs)
        elif fmt_lower in ("modern", "kind2", "auto"):
            self.to_modern(path, **kwargs)
        else:
            raise ValueError(f"Unknown format '{fmt}' specified.")
            
    @property
    def summary(self):
        """
        Provide a Bunch object with a summary of key attributes.
        """
        if self.info.df is None:
            return Bunch(status="Data not loaded")
        
        has_read (self, attributes="info") 
        
        # Safely get values, providing defaults if not loaded
        hdr = self.info.header
        st = self.info.station
        frq = self.info.frequency

        # Smartly select useful information
        info_dict = {
            "source_file": (
                self._source_path.name if self._source_path else "N/A"
            ),
            "data_kind": f"Kind-{self._kind}" if self._kind else "N/A",
            "project": hdr.annotation.project_name or "N/A",
            "survey_type": hdr.config.survey_type or "N/A",
            "line_name": hdr.config.line_name or "N/A",
            "num_stations": st.n_unique if st else 0,
            "num_frequencies": frq.n_unique if frq else 0,
            "station_range": (
                f"{st.span[0]} - {st.span[1]} {st.unit}"
                if st and st.span
                else "N/A"
            ),
            "frequency_range": (
                f"{frq.unique().min():.4g} - {frq.unique().max():.4g} Hz"
                if frq and frq.n_unique > 0
                else "N/A"
            ),
            "total_rows": len(self.info.df),
        }

        return Bunch(**info_dict)
    
    def asdict(self) -> Dict[str, Any]:
        r"""Return a shallow, JSON-serializable representation.

        This method aggregates the keyword dictionaries from all
        header components into a single dictionary, providing a
        complete, serializable view of the survey's metadata.

        Returns
        -------
        dict
            A dictionary containing all header information.
        """
        if self.info.df is None:
            return {"status": "Data not loaded"}

        return self.header.to_keywords()   
    
    def __str__(self) -> str:
        """Provide a concise, human-readable representation."""
        if self.info.df is None or self.info.df.empty:
            status = "empty"
        else:
            n_st = self.info.station.n_unique
            n_f = self.info.frequency.n_unique
            status = (
                f"stations={n_st}, freqs={n_f}, "
                f"rows={len(self.info.df)}"
            )

        src = (
            f", source='{self._source_path.name}'"
            if self._source_path
            else ""
        )
        return f"{self.__class__.__name__}({status}{src})"

    def __repr__(self) -> str:
        """Provide an unambiguous developer representation."""
        path_repr = (
            f"Path('{self._source_path}')"
            if self._source_path
            else "None"
        )
        return (
            f"{self.__class__.__name__}("
            f"source_path={path_repr}, "
            f"verbose={self.verbose}, "
            f"loaded={'is not None' if self.info.df is not None else 'False'}"
            ")"
        )


@has_fit("raise")
class AVG(BaseAVG):
    r"""High-level façade for a Zonge AVG/AMTAVG dataset.

    An :class:`AVG` instance loads a raw text file, normalizes
    its content, and hydrates a full suite of data components,
    including stations, frequencies, impedance, quality metrics,
    and header blocks. This makes the entire survey accessible
    through strongly-typed attributes [1]_.

    Parameters
    ----------
    verbose : bool, default False
        If ``True``, log messages will be printed to the console
        during file reading and writing operations.

    Attributes
    ----------
    header : :class:`pycsamt.zonge.heads.Header`
        Aggregates *Hardware*, *SurveyAnnotation*,
        *SurveyConfiguration*, and *Rx/Tx* property blocks.
    station : :class:`pycsamt.zonge.survey.Station`
        Manages survey line geometry and station coordinates.
    frequency : :class:`pycsamt.zonge.meas.Frequency`
        Manages the frequency axis for the measurements.
    resistivity, phase : :class:`~.resphase.Resistivity`, :class:`~.resphase.Phase`
        Containers for apparent resistivity and phase data.
    z : :class:`pycsamt.zonge.z.Z`
        Component for computing the complex impedance tensor.
    df : pd.DataFrame or None
        The core tidy DataFrame containing all available data
        columns after parsing and normalization.

    Methods
    -------
    from_file(path, verbose=False)
        Classmethod to load and parse an AVG file. This is the
        primary entry point for creating an `AVG` instance.
    read(source, meta=None)
        Populates the object from a source (file, DataFrame, etc.).
    write(path, fmt='auto', **kwargs)
        Writes the data to a file, dispatching to `to_modern` or
        `to_legacy` based on the `fmt` argument.
    to_xarray()
        Exports the entire dataset into a single, comprehensive
        :class:`xarray.Dataset`.
    to_tensor(var='z', **kwargs)
        Exports a specific variable as a 2x2 NumPy tensor.
    asdict()
        Returns a JSON-serializable dictionary of all header
        metadata.

    Notes
    -----
    Zonge Engineering produced two AVG file formats:

    * **Kind-1 (Legacy)**: Fixed-width text tables with minimal
        metadata.
    * **Kind-2 (Modern)**: CSV-based format with rich
        ``$key=value`` headers, used by modern AMTAVG and ASTATIC.

    The loader automatically detects the format and transforms it
    into a canonical, tidy DataFrame. All subsequent processing is
    format-agnostic. Zonge placeholders (``*``) are converted to
    ``NaN``.

    Examples
    --------
    Load a file, access a component, and write back:

    >>> from pycsamt.zonge import AVG
    >>> avg = AVG.from_file('LCS01.avg', verbose=True)
    >>> # Access resistivity for a specific component
    >>> rho_xy = avg.resistivity.frame[
    ...     avg.resistivity.frame.comp == 'ExHy'
    ... ]
    >>> avg.write('LCS01_clean.avg')

    Build from an in-memory DataFrame:

    >>> from pycsamt.zonge.utils import load_avg
    >>> df, meta = load_avg('raw.avg')
    >>> avg = AVG()
    >>> avg.read(df, meta=meta)
    >>> print(avg.station.span)

    References
    ----------
    .. [1] Zonge International, Inc. (2014). *ASTATIC v3.70
           User Manual*.
    """
    def __init__(self, verbose: bool = False):
        super().__init__(verbose=verbose)

    @classmethod
    def from_file(
        cls,
        path: Union[str, Path],
        *,
        verbose: bool = False,
    ) -> "AVG":
        """
        Load and parse an AVG file (legacy or modern).

        This factory method handles file classification, parsing,
        and optional transformation of legacy data into a modern,
        consistent structure.
        """
        obj = cls(verbose=verbose)
        obj.read(path)
        return obj

    # --- Properties to expose components ---
    @property
    def header(self):
        """Access the Header component."""
        return self.info.header

    @property
    def station(self):
        """Access the Station component."""
        return self.info.station

    @property
    def z(self):
        """Access the impedance (Z) component."""
        return self.info.z

    @property
    def resistivity(self):
        """Access the Resistivity component."""
        return self.info.resistivity

    @property
    def phase(self):
        """Access the Phase component."""
        return self.info.phase

    @property
    def frequency(self):
        """Access the Frequency component."""
        return self.info.frequency
        
    @property
    def df(self) -> Optional[pd.DataFrame]:
        """
        Access the core tidy DataFrame containing all available
        data columns after parsing and normalization.
        """
        return self.info.df

    @ensure_pkg ('xarray', extra="xarray is required.")
    def to_xarray(
        self,
        *,
        include_qc: bool = True,
        include_z: bool = True,
    ):
        """
        Export the complete dataset to a single xarray.Dataset.

        This method aggregates primary data (resistivity, phase),
        computed impedance (Z), and optional quality control (QC)
        metrics into a unified, multi-dimensional dataset.
        """
        if self.info.df is None:
            warnings.warn(
                "Cannot create xarray.Dataset from empty data.")
            return None
        has_read (self, attributes="info") 
        
        # Start with the primary data variables
        ds_rho = self.resistivity.to_xarray()
        ds_phase = self.phase.to_xarray()

        # Merge primary datasets
        ds = ds_rho.merge(ds_phase)

        # Optionally compute and merge complex impedance
        if include_z:
            try:
                z_real = self.z.to_xarray(var="z_real")
                z_imag = self.z.to_xarray(var="z_imag")
                z_err = self.z.to_xarray(var="z_err")
                ds = ds.merge(z_real)
                ds = ds.merge(z_imag)
                ds = ds.merge(z_err)
            except Exception as e:
                if self.verbose:
                    logger.warning(
                        f"Could not compute impedance Z: {e}")

        # Optionally merge all available QC metrics
        if include_qc:
            qc_components = {
                "pc_emag": self.info.pc_emag,
                "pc_hmag": self.info.pc_hmag,
                "pc_rho": self.info.pc_rho,
                "s_ephz": self.info.s_ephz,
                "s_hphz": self.info.s_hphz,
                "s_phz": self.info.s_phz,
            }
            for name, comp in qc_components.items():
                try:
                    ds_qc = comp.to_xarray()
                    # Rename data var to avoid conflicts
                    var_name = list(ds_qc.data_vars)[0]
                    ds = ds.merge(ds_qc.rename({var_name: name}))
                except Exception:
                    if self.verbose:
                        logger.info(
                            f"Skipping QC component '{name}': "
                            "data not available."
                        )

        # Attach comprehensive header metadata
        ds.attrs.update(self.header.to_keywords())
        ds.attrs["source_file"] = (
            str(self._source_path) 
            if self._source_path else "Unknown"
        )
        return ds

    def to_tensor(
        self,
        var: str = "z",
        *,
        station: Optional[Union[int, float]] = None,
        agg: str | None = "mean",
        fill_value: float = np.nan,
        sort_freq: bool = True,
        align: str = "union",
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Export a specific variable as a 2x2 tensor.

        This is a convenience method that delegates to the
        appropriate component's `to_tensor` method.
        """
        has_read (self, attributes="info") 
        component_map = {
            "z": self.z,
            "z_real": self.z,
            "z_imag": self.z,
            "z_err": self.z,
            "rho": self.resistivity,
            "phase": self.phase,
        }
        if var not in component_map:
            raise ValueError(
                f"Variable '{var}' not supported for tensor export. "
                f"Choose from: {list(component_map.keys())}"
            )

        component = component_map[var]
        return component.to_tensor(
            var=var,
            station=station,
            agg=agg,
            fill_value=fill_value,
            sort_freq=sort_freq,
            align=align,
        )

    def __str__(self) -> str:
        """Provide a concise string representation."""
        if self.info.df is None:
            return "AVG(empty)"
        return f"AVG(source='{self._source_path.name}')"

    __repr__ = __str__


        
class AMTAVG(AVG):
    r"""Extends AVG with tensor component properties and computations.

    This class inherits from :class:`~pycsamt.zonge.avg.AVG` and
    provides a more specialized interface for accessing individual
    components of the impedance, resistivity, and phase tensors as
    pandas Series. It is designed for users who need to work
    directly with specific tensor elements (e.g., z_xy, rho_yx).

    In addition to providing direct access to tensor components,
    this class includes methods for advanced data manipulation,
    such as recalculating resistivity and phase from the impedance
    tensor.

    Attributes
    ----------
    z_xx, z_xy, z_yx, z_yy : pd.Series
        The complex impedance for the respective tensor component.
    z_xx_err, z_xy_err, z_yx_err, z_yy_err : pd.Series
        The propagated error for the respective impedance component.
    res_xx, res_xy, res_yx, res_yy : pd.Series
        The apparent resistivity for the respective tensor component.
    res_xx_err, res_xy_err, res_yx_err, res_yy_err : pd.Series
        The percent error for the respective resistivity component.
    phase_xx, phase_xy, phase_yx, phase_yy : pd.Series
        The impedance phase for the respective tensor component.
    phase_xx_err, phase_xy_err, phase_yx_err, phase_yy_err : pd.Series
        The standard deviation of the respective phase component.

    Methods
    -------
    compute_resistivity_phase()
        Calculates apparent resistivity and phase from the complex
        impedance tensor Z.
    set_resistivity_phase(rho, phi, rho_err=None, phi_err=None)
        Updates the dataset with new resistivity and phase values,
        triggering a recalculation of the impedance tensor.
    get_tensor_by_station(station_id, var='z')
        Fetches a 3D tensor for a single station using xarray.

    Examples
    --------
    >>> from pycsamt.zonge import AMTAVG
    >>> amt_avg = AMTAVG.from_file('path/to/your/data.avg')
    >>> # Access a specific component of the resistivity tensor
    >>> rho_xy_values = amt_avg.res_xy
    >>> print(rho_xy_values.head())

    See Also
    --------
    AVG : The parent class providing the main data loading and
          exporting functionality.
    Z : The component that handles impedance calculations.
    Resistivity : The component that manages resistivity data.
    Phase : The component that manages phase data.
    """
    
    def compute_resistivity_phase(self):
        """
        Compute rho and phi from the complex impedance Z.
        """
        if self.info.df is None:
            raise ValueError("Data frame is not loaded.")
        
        has_read (self, attributes="info") 
        
        z_complex = self.z.z
        freq = self.info.df["freq"]
        omega = 2 * PI * freq

        rho = (np.abs(z_complex)**2) / (omega * MU_0)
        phi_rad = np.arctan2(z_complex.imag, z_complex.real)
        phi_mrad = phi_rad * 1000.0

        # Error propagation
        z_err = self.z.z_err
        # Avoid division by zero if z_complex is zero
        abs_z = np.abs(z_complex)
        safe_abs_z = np.where(abs_z == 0, np.nan, abs_z)

        rho_err = 200 * (z_err / safe_abs_z)
        phi_err = (
            1000
            * (z_err / safe_abs_z)
            * np.abs(np.sin(phi_rad))
        )

        return rho, phi_mrad, rho_err, phi_err

    def set_resistivity_phase(
        self,
        rho: pd.Series,
        phi: pd.Series,
        rho_err: Optional[pd.Series] = None,
        phi_err: Optional[pd.Series] = None,
    ):
        """
        Attach new rho/phi data and reconstruct Z.
        """
        if self.info.df is None:
            raise ValueError("Data frame is not loaded.")
        
        has_read (self, attributes="info") 
        
        # Work on a copy to avoid modifying the original df
        df = self.info.df.copy()

        # Use .loc to ensure alignment and avoid warnings
        # df.loc[rho.index, "rho"] = rho
        # df.loc[phi.index, "phase"] = phi
        # if rho_err is not None:
        #     df.loc[rho_err.index, "pc_rho"] = rho_err
        # if phi_err is not None:
        #     df.loc[phi_err.index, "sphz"] = phi_err
        
        # Direct assignment is safer and respects index alignment
        df["rho"] = rho
        df["phase"] = phi
        if rho_err is not None:
            df["pc_rho"] = rho_err
        if phi_err is not None:
            df["sphz"] = phi_err
            
        # Re-read components to update their internal state
        self.info.read(df, self.info.meta)
        
    def get_tensor_by_station(
        self,
        station_id: Union[int, float],
        *,
        var: Literal["z", "rho", "phase"] = "z",
    ):
        """
        Fetch a 3D tensor for a single station using xarray.

        This method returns an xarray.DataArray for the specified
        variable, indexed by frequency and the 2x2 tensor axes.
        """
        has_read (self, attributes="info") 
        
        component_map = {
            "z": self.z,
            "rho": self.resistivity,
            "phase": self.phase,
        }
        if var not in component_map:
            raise ValueError(
                f"Variable '{var}' not supported. "
                f"Choose from: {list(component_map.keys())}"
            )

        component = component_map[var]
        # The component's to_xarray method handles tensor shape
        ds = component.to_xarray()

        try:
            # Use .sel() for label-based selection
            station_tensor = ds.sel(station=station_id)
            return station_tensor
        except KeyError:
            logger.error(f"Station ID '{station_id}' not found.")
            raise


    # --- Z Tensor Components ---
    @property
    def z_xx(self): 
        return self.z.z_xx
    @property
    def z_xy(self): 
        return self.z.z_xy
    @property
    def z_yx(self):
        return self.z.z_yx
    @property
    def z_yy(self): 
        return self.z.z_yy
    @property
    def z_xx_err(self):
        return self.z.z_xx_err
    @property
    def z_xy_err(self): 
        return self.z.z_xy_err
    @property
    def z_yx_err(self): 
        return self.z.z_yx_err
    @property
    def z_yy_err(self): 
        return self.z.z_yy_err

    # --- Resistivity Tensor Components ---
    @property
    def res_xx(self):
        df = self.resistivity.frame
        return df[df.comp == "ExHx"]["rho"]
    @property
    def res_xy(self):
        df = self.resistivity.frame
        return df[df.comp == "ExHy"]["rho"]
    @property
    def res_yx(self):
        df = self.resistivity.frame
        return df[df.comp == "EyHx"]["rho"]
    @property
    def res_yy(self):
        df = self.resistivity.frame
        return df[df.comp == "EyHy"]["rho"]

    # --- Resistivity Error Components ---
    @property
    def res_xx_err(self):
        df = self.info.pc_rho.frame
        return df[df.comp == "ExHx"]["pc_rho"]
    @property
    def res_xy_err(self):
        df = self.info.pc_rho.frame
        return df[df.comp == "ExHy"]["pc_rho"]
    @property
    def res_yx_err(self):
        df = self.info.pc_rho.frame
        return df[df.comp == "EyHx"]["pc_rho"]
    @property
    def res_yy_err(self):
        df = self.info.pc_rho.frame
        return df[df.comp == "EyHy"]["pc_rho"]

    # --- Phase Tensor Components ---
    @property
    def phase_xx(self):
        df = self.phase.frame
        return df[df.comp == "ExHx"]["phase"]
    @property
    def phase_xy(self):
        df = self.phase.frame
        return df[df.comp == "ExHy"]["phase"]
    @property
    def phase_yx(self):
        df = self.phase.frame
        return df[df.comp == "EyHx"]["phase"]
    @property
    def phase_yy(self):
        df = self.phase.frame
        return df[df.comp == "EyHy"]["phase"]

    # --- Phase Error Components ---
    @property
    def phase_xx_err(self):
        df = self.info.s_phz.frame
        return df[df.comp == "ExHx"]["sphz"]
    @property
    def phase_xy_err(self):
        df = self.info.s_phz.frame
        return df[df.comp == "ExHy"]["sphz"]
    @property
    def phase_yx_err(self):
        df = self.info.s_phz.frame
        return df[df.comp == "EyHx"]["sphz"]
    @property
    def phase_yy_err(self):
        df = self.info.s_phz.frame
        return df[df.comp == "EyHy"]["sphz"]

    def __str__(self) -> str:
        """Provide a concise, human-readable representation."""
        if self.info.df is None or self.info.df.empty:
            status = "empty"
        else:
            n_st = self.info.station.n_unique
            n_f = self.info.frequency.n_unique
            status = (
                f"stations={n_st}, freqs={n_f}, "
                f"rows={len(self.info.df)}"
            )

        src = (
            f", source='{self._source_path.name}'"
            if self._source_path
            else ""
        )
        return f"AMTAVG({status}{src})"

    def __repr__(self) -> str:
        """Provide an unambiguous developer representation."""
        if self._source_path:
            return (
                f"AMTAVG.from_file("
                f"'{self._source_path!s}', "
                f"verbose={self.verbose})"
            )
        return f"AMTAVG(verbose={self.verbose}, loaded=False)"       
    