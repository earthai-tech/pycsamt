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
from collections.abc import Mapping
from datetime import datetime
from pathlib import Path
from typing import (
    Any,
    Literal,
)

import numpy as np
import pandas as pd

from ..api.bunch import Bunch
from ..constants import MU_0, PI
from ..decorators import has_fit
from ..exceptions import AvgDataError
from ..utils._dependency import import_optional_dependency
from ..utils.validation import has_read
from ._transfer import LegacyAVGBase
from .base import AVGFrame, guess_kind_from_df
from .config import Zonge
from .heads import _emit_hardware_banner
from .info import DataInfo
from .survey import Topography
from .utils import classify_avg_format, load_avg, write_avg

__all__ = ["BaseAVG","AVG", "AMTAVG"]



class BaseAVG(Zonge):
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
        super().__init__(verbose=verbose)

        self.info: DataInfo = DataInfo(verbose=verbose)
        self._kind: int | None = None
        self._source_path: Path | None = None
        self.topo: Topography | None = None

    def read(
        self,
        source: str | Path | AVGFrame | pd.DataFrame,
        meta: Mapping[str, Any] | None = None,
    ) -> BaseAVG:
        r"""Read, parse, and populate components from a data source.

        This is the primary data ingestion method for the class. It
        is designed to be flexible, accepting various data sources
        and orchestrating the entire process of parsing, transforming,
        and populating the internal data components.

        Parameters
        ----------
        source : str, Path, AVGFrame, or pd.DataFrame
            The data source to load. This can be:

            - A string or `pathlib.Path` pointing to a Zonge AVG
              file (both legacy and modern formats are supported).
            - A `pandas.DataFrame` containing AVG data. Note that
              if the DataFrame is not already standardized, it will
              be automatically standardized by `AVGFrame`.
            - An existing :class:`~.base.AVGFrame` object.

        meta : mapping, optional
            An optional dictionary of metadata. This is primarily
            used when `source` is a `pd.DataFrame` to provide
            header-like information.

        Returns
        -------
        self : BaseAVG
            The method returns the instance of the class, allowing
            for convenient method chaining.

        Raises
        ------
        TypeError
            If the `source` is of an unsupported type.
        FileNotFoundError
            If `source` is a path that does not exist.
        AvgFileError
            If a file source cannot be reliably classified as either
            legacy (kind-1) or modern (kind-2).
        AvgDataError
            If parsing fails due to malformed or missing data.

        Notes
        -----
        The `read` method encapsulates a complete data processing
        pipeline:

        1.  **File Handling**: If a path is provided, it reads the
            file and uses :func:`~.utils.classify_avg_format` to
            determine the file kind.
        2.  **Transformation**: If a legacy (kind-1) file is
            detected, it automatically uses the
            :class:`~._transfer.LegacyAVGBase` transformer to
            convert the data to a modern, structured format.
        3.  **Standardization**: If a raw `pd.DataFrame` is
            provided, it is wrapped in an `AVGFrame`, which
            triggers the automatic standardization of column names
            to the internal canonical schema.
        4.  **Component Population**: The final, clean DataFrame and
            metadata are passed to the `DataInfo` component, which
            in turn populates all the individual data components
            (e.g., `Station`, `Z`, `Resistivity`).

        Examples
        --------
        >>> from pycsamt.zonge import AVG
        >>> # Read from a modern file
        >>> avg = AVG().read('data/avg/K2.avg')
        >>>
        >>> # Read from a legacy file (will be auto-transformed)
        >>> avg_legacy = AVG().read('data/avg/K1.avg')
        >>>
        >>> # Read from an in-memory DataFrame
        >>> import pandas as pd
        >>> df = pd.DataFrame({'Freq': [1024], 'Resistivity': [100]})
        >>> avg_from_df = AVG().read(df)
        """
        df: pd.DataFrame
        final_meta: Mapping[str, Any]

        if isinstance(source, (str, Path)):
            self._source_path = Path(source)
            if self.verbose:
                self._logger.info(
                    f"Reading from file: {self._source_path}"
                )
            lines = self._source_path.read_text(
                errors="replace"
            ).splitlines()
            self._kind = classify_avg_format(lines)
            df, final_meta = load_avg(self._source_path)

            if self._kind == 1:
                if self.verbose:
                    self._logger.info("Transforming legacy data.")
                transformer = LegacyAVGBase()
                ds = transformer.from_dataframe(df, meta=final_meta)
                df = ds.to_dataframe().reset_index()
                final_meta = ds.attrs

        elif isinstance(source, pd.DataFrame):
            if self.verbose:
                self._logger.info("Reading from pandas DataFrame.")
            frame = AVGFrame (data=source, meta=meta or {})

            # Guess the kind from the now-standardized frame
            # The default 'soft' mode will handle this correctly
            df, final_meta, self._kind = guess_kind_from_df(
                frame, transform=True,
                verbose = self.verbose
            )
            self._source_path = frame.source

        elif isinstance(source, AVGFrame):
            if self.verbose:
                self._logger.info("Reading from AVGFrame object.")
            df = source.data
            final_meta = source.meta
            self._source_path = source.source
            self._kind = guess_kind_from_df(df)

        else:
            raise TypeError(
                "Unsupported source type. Expected Path, str, "
                "DataFrame, or AVGFrame."
            )

        if self.verbose:
            self._logger.info("Populating data components...")
        self.info.read(df, final_meta)
        if self.verbose:
            self._logger.info("AVG data successfully loaded.")

        return self

    def add_topography(
        self,
        stn_file: str | Path | pd.DataFrame,
        utm_zone: str= None,
        epsg: int =None,
    ) -> BaseAVG:
        r"""Read and attach station topography data.

        This method reads a Zonge ``.stn`` file (or a DataFrame
        with the same structure) and attaches a populated
        :class:`~.survey.Topography` object to the instance.

        Parameters
        ----------
        stn_file : str, Path, or DataFrame
            The path to the ``.stn`` file or a pre-loaded
            DataFrame containing the station location data.

        Returns
        -------
        self : BaseAVG
            The method returns the instance, allowing for
            method chaining.
        """

        has_read(self)

        if self.verbose:
            self._logger.info(f"Reading topography from: {stn_file}")

        # Create and read the Topography component
        self.topo = Topography(
            verbose=self.verbose,
            utm_zone= utm_zone,
            epsg= epsg
            ).read(stn_file)

        # XXX TODO: Optional: add logic here to merge elevation
        # into the main df if needed for specific calculations,
        # but keeping it separate seems generally better.
        # For example:
        # topo_map = self.topo.frame.set_index('station')['elevation']
        # self.info.df['elevation'] = self.info.df['station'].map(topo_map)
        # let keep safe for this version now.

        return self

    def to_modern(
        self,
        path: str | Path | None = None,
        *,
        stamp: bool = True,
        float_fmt: str = "%.6g",
        na_rep: str = "*",
        header_spaces: bool = False,
    ):
        r"""Write data to a modern (kind-2) AVG file.

        This method serializes the object's data into a modern,
        CSAVGW-style ``.avg`` file. It delegates the core writing
        and formatting logic to the :func:`~.utils.write_avg`
        function.

        Parameters
        ----------
        path : str or pathlib.Path, optional
            The output file path. If ``None``, a default filename
            (e.g., ``K2_modern.avg``) is created in the current
            working directory.
        stamp : bool, default True
            If ``True``, a ``$Written=<timestamp>`` line is added to
            the header for provenance.
        float_fmt : str, default "%.6g"
            The format specifier for floating-point numbers.
        na_rep : str, default "*"
            The string representation for missing (NaN) values.
        header_spaces : bool, default False
            If ``True``, adds spaces around the equals sign in header
            keywords (e.g., ``$Key = Value``).

        Notes
        -----
        This method intelligently prepares the data for writing by:
        1.  Assembling a complete header by calling the
            :meth:`~.heads.Header.to_keywords` method on the `header`
            component.
        2.  Generating a file banner from the `hardware` component.
        3.  Filtering out any artificial data rows (e.g., those
            with ``NaN`` in the 'rho' column) that may have been
            created during a legacy-to-modern transformation.

        See Also
        --------
        to_legacy : Write data to a legacy (kind-1) file.
        pycsamt.zonge.utils.write_avg : The underlying file writing
            function.
        """
        has_read (self)

        if path is None:
            base = (
                self._source_path.stem
                if self._source_path
                else "export"
            )
            path = Path.cwd() / f"{base}_modern.avg"

        if self.verbose:
            self._logger.info(f"Writing modern AVG file to: {path}")

        # Filter out artificial NaN rows before writing ---
        df_to_write = self.info.df.dropna(subset=['rho']).copy()

        meta = self.info.header.to_keywords()
        banner = _emit_hardware_banner(self.info.header.hardware)

        write_avg(
            core=df_to_write,
            extra=None,
            meta=meta,
            path=path,
            stamp=stamp,
            float_fmt=float_fmt,
            na_rep=na_rep,
            header_spaces=header_spaces,
            banner_lines=banner,
        )

    def to_legacy(
        self,
        path: str | Path | None = None,
        *,
        precision: int = 4,
        na_rep: str = "*",
    ):
        r"""Write data to a legacy (kind-1) AVG file.

        This method reconstructs the classic, fixed-width format of
        older Zonge AVG files. It uses the populated data components
        to generate the header and formatted data rows.

        Parameters
        ----------
        path : str or pathlib.Path, optional
            The output file path. If ``None``, a default filename
            (e.g., ``K2_legacy.avg``) is created in the current
            working directory.
        precision : int, default 4
            The number of decimal places for scientific notation
            fields (e.g., 'Emag', 'Resistivity').
        na_rep : str, default "*"
            The string representation for missing (NaN) values.

        Notes
        -----
        The conversion from a modern, comprehensive data structure to
        the legacy format involves several key steps:

        - **Header Synthesis**: It creates a minimal legacy header
          using ``$ASPACE`` and ``$XMTR`` keywords derived from the
          `header` component.
        - **Column Translation**: It translates the internal canonical
          column names (e.g., 'pc_emag') back to their legacy
          equivalents (e.g., '%Emag') for correct formatting.
        - **Data Filtering**: It first filters out any artificial
          data rows (e.g., those with ``NaN`` in the 'rho' column) to
          ensure the output is clean.

        Be aware that some information from a modern dataset may be
        lost in this conversion, as the legacy format does not
        support all the fields of the modern format.

        See Also
        --------
        to_modern : Write data to a modern (kind-2) file.
        """
        has_read (self)

        if path is None:
            base = (
                self._source_path.stem
                if self._source_path
                else "export"
            )
            path = Path.cwd() / f"{base}_legacy.avg"

        lines = []
        h = self.info.header
        # Today's date for 'Processed' if not available
        today_str = datetime.now().strftime("%d %b %y")
        lines.append(
            f"\\ AMTAVG 7.76: "
            f"\"{h.hardware.source_file or 'pycsamt.export'}\", "
            f"Dated {h.hardware.dated or '...'}, "
            f"Processed {h.hardware.processed or today_str}"
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

        # Map from the canonical lowercase names in self.info.df
        # to the names expected by the row.get() calls below.
        rename_map = {
            "e.%err": "pc_emag", "e.perr": "s_ephz",
            "h.%err": "pc_hmag", "h.perr": "s_hphz",
            "rho.%err": "pc_rho", "z.perr": "s_phz",
        }
        df = self.info.df.dropna(subset=['rho']).copy()
        df = df.rename(
            columns=rename_map
        ).sort_values(by=["station", "freq"])

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

        if self.verbose:
            self._logger.info(f"Writing legacy AVG file to: {path}")


    def write(
        self,
        path: str | Path | None = None,
        *,
        fmt: Literal[
            "kind1", "kind2", "legacy", "modern", "auto"
        ] = "auto",
        **kwargs,
    ):
        r"""Write AVG data to a file in the specified format.

        This method acts as a high-level dispatcher, calling either
        the :meth:`to_legacy` or :meth:`to_modern` method based on the
        `fmt` parameter. It provides a single, convenient entry point
        for file serialization.

        Parameters
        ----------
        path : str or pathlib.Path, optional
            The output file path. If ``None``, a default filename will
            be generated based on the source file's name.
        fmt : {'legacy', 'kind1', 'modern', 'kind2', 'auto'}, default 'auto'
            The desired output format.
            - 'legacy' or 'kind1': Writes a fixed-width legacy file.
            - 'modern' or 'kind2': Writes a CSV-based modern file.
            - 'auto': Defaults to the modern format.
        **kwargs
            Additional keyword arguments to be passed to the specific
            writing method (`to_legacy` or `to_modern`). For example,
            `precision` for `to_legacy` or `float_fmt` for
            `to_modern`.

        Raises
        ------
        ValueError
            If an unsupported `fmt` is specified.
        NotReadError
            If this method is called before data has been loaded.

        See Also
        --------
        to_legacy : Writes data to a legacy (kind-1) file.
        to_modern : Writes data to a modern (kind-2) file.
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
        r"""Provide a `Bunch` object with a summary of key attributes.

        This property aggregates the most important descriptive
        metadata and statistics from the loaded dataset into a single,
        easy-to-read object. It is useful for quickly inspecting the
        contents and overall scope of the survey data.

        Returns
        -------
        pycsamt.api.bunch.Bunch
            A `Bunch` object (which behaves like a dictionary with
            attribute-style access) containing key information such
            as the source filename, data format kind, project name,
            number of stations and frequencies, and their respective
            ranges.

        Notes
        -----
        The `summary` property is read-only and is generated on-demand
        by introspecting the various data components (e.g., `header`,
        `station`, `frequency`). If no data has been loaded, it will
        return a `Bunch` with a 'status' key indicating this.

        The `Bunch` object has a user-friendly string representation,
        making it ideal for printing in interactive sessions.

        Examples
        --------
        >>> from pycsamt.zonge import AVG
        >>> avg = AVG.from_file('data/avg/K2.avg')
        >>> print(avg.summary)
        <Bunch object with keys: source_file, data_kind, project, ...>
        >>> avg.summary.num_stations
        28
        """
        if self.info.df is None:
            return Bunch(status="Data not loaded")

        has_read (self)

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

    def asdict(self) -> dict[str, Any]:
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

    def __has_read__(self) -> bool:
        """
        Checks if the object has been successfully populated with
        data.
        """
        # The most reliable check is whether the core DataFrame
        # exists and is not empty.
        return (
            self.info is not None
            and self.info.df is not None
            and not self.info.df.empty
        )


    def __repr__(self) -> str:
        """Provide an unambiguous developer representation."""
        path_repr = (
            f"Path('{self._source_path}')"
            if self._source_path
            else "None"
        )
        return (
            f"{self.__class__.__name__}.from_file("
            f"path={path_repr}, "
            f"verbose={self.verbose})"
            if self.__has_read__() else
            f"{self.__class__.__name__}(verbose={self.verbose}, loaded=False)"
        )

    def _summary_stats(self) -> str:
        """Creates a statistical summary block for the AVG data."""
        if not self.__has_read__():
            return "  (No statistics available: data not loaded)\n"

        summary_bunch = self.summary
        n_stations = summary_bunch.get("num_stations", 0)
        n_freqs = summary_bunch.get("num_frequencies", 0)

        # Initialize coordinate and elevation data as None
        lats, lons, elevs = None, None, None

        # Get location data ONLY if a topography file is loaded
        if self.topo is not None and not self.topo.frame.empty:
            # Ensure lat/lon are available for the summary
            if 'latitude' not in self.topo.frame.columns or (
                    'longitude' not in self.topo.frame.columns):
                self.topo.convert_coords(to='ll', inplace=True)

            lats = self.topo.latitude
            lons = self.topo.longitude
            elevs = self.topo.elevation

        freqs = self.info.frequency.unique()

        # Helper to safely format data ranges
        def format_range(
            arr: np.ndarray | None, fmt: str
        ) -> str:
            """Safely formats the min/max of an array."""
            # Check if array is valid and has finite numbers
            if (
                arr is not None and arr.size > 0
                and np.any(np.isfinite(arr))
            ):
                return (
                    f"{np.nanmin(arr):{fmt}} to "
                    f"{np.nanmax(arr):{fmt}}"
                )
            return "N/A"

        # Generate range strings for all variables
        lat_range = format_range(lats, ".4f")
        lon_range = format_range(lons, ".4f")
        elev_range = format_range(elevs, ".2f") # Fixed
        freq_range = format_range(freqs, ".2E")

        lines = [
            "  " + "-" * 68,
            "  Statistical Summary:",
            ("    Survey Type: "
             f"{summary_bunch.get('survey_type', 'N/A')}"),
            ("    Line Name:   "
             f"{summary_bunch.get('line_name', 'N/A')}"),
            (f"    Stations:    {n_stations} "
             f"({summary_bunch.get('station_range', 'N/A')})"),
            f"    Frequencies: {n_freqs}",
            f"    Freq Range (Hz): {freq_range}",
            f"    Latitude Range:  {lat_range}",
            f"    Longitude Range: {lon_range}",
            f"    Elevation Range (m): {elev_range}",
            "  " + "-" * 68,
        ]
        return "\n".join(lines)


    def __str__(self) -> str:  # pragma: no cover
        """Provides a detailed, robust summary of the AVG data."""
        if not self.__has_read__():
            src = (
                f"'{self._source_path.name}'"
                if self._source_path else "Not Loaded"
            )
            return f"{self.__class__.__name__}(source={src}, status=empty)"

        # --- Header ---
        title = (
            f" {self.__class__.__name__} Summary"
            f" (Source: {self._source_path.name}) "
        )
        width = 95  # Adjusted width for new table format
        header = [
            "=" * width,
            title.center(width),
            "=" * width
        ]

        # --- Statistical Summary ---
        stats_str = self._summary_stats()

        # --- Per-Site Details (from topography component) ---
        details = ["\nStation Details:"]

        if self.topo is not None and not self.topo.frame.empty:
            topo = self.topo

            # 1. Calculate all the summary statistics from topo data
            n_stations = len(topo.stations)
            steps = topo.get_step()
            valid_steps = steps[steps > 0]
            total_length = valid_steps.sum()
            mean_azimuth = topo.get_azimuth(mode='mean')
            elevations = topo.elevation

            # 2. Format the calculated values for clean display
            length_str = (
                f"{total_length / 1000:.2f} km" if total_length > 1000
                else f"{total_length:.1f} m"
            )

            def to_cardinal(deg: float) -> str:
                """Helper to convert azimuth degree to cardinal direction."""
                if np.isnan(deg): return ""
                dirs = [
                    "N", "N-NE", "NE", "E-NE", "E", "E-SE", "SE", "S-SE",
                    "S", "S-SW", "SW", "W-SW", "W", "W-NW", "NW", "N-NW"
                ]
                return dirs[int(round(deg / 22.5)) % 16]

            azimuth_str = (
                f"{mean_azimuth:.1f}° ({to_cardinal(mean_azimuth)})"
            )

            spacing_str = (
                f"Min: {valid_steps.min():.1f}, "
                f"Mean: {valid_steps.mean():.1f}, "
                f"Max: {valid_steps.max():.1f}"
            ) if not valid_steps.empty else "N/A"

            elev_range_str = (
                f"{elevations.min():.1f} - {elevations.max():.1f} "
                f"(Δ: {elevations.max() - elevations.min():.1f})"
            ) if elevations.size > 1 else "N/A"

            start_pt_str = f"({topo.easting[0]:.1f}, {topo.northing[0]:.1f})"
            end_pt_str = f"({topo.easting[-1]:.1f}, {topo.northing[-1]:.1f})"

            # 3. Structure the data for a formatted key-value table
            table_data = [
                ("Number of Stations", f"{n_stations}"),
                ("Total Line Length", length_str),
                ("Line Azimuth (Mean)", azimuth_str),
                ("Station Spacing (m)", spacing_str),
                ("Elevation Range (m)", elev_range_str),
                ("Start Point (E, N)", start_pt_str),
                ("End Point (E, N)", end_pt_str),
            ]

            # 4. Build the formatted table string for display
            header1, header2 = "Geometry", "Value"

            # max_key_len = max(len(key) for key, _ in table_data)
            # for key, val in table_data:
            #     details.append(f"  {key.ljust(max_key_len)}: {val}")
            # 4. Calculate dynamic column widths based on content
            col1_width = max(len(header1), max(len(k) for k, _ in table_data))
            col2_width = max(len(header2), max(len(v) for _, v in table_data))

            # 5. Build the formatted table strings
            separator = (f"  {'-' * (col1_width + col2_width + 7)}")

            details.append(separator)
            details.append(f"  | {header1.ljust(col1_width)} | "
                           f"{header2.ljust(col2_width)} |")
            details.append(separator)

            for key, val in table_data:
                details.append(f"  | {key.ljust(col1_width)} | "
                               f"{val.ljust(col2_width)} |")
            details.append(separator)


        else:
            # Message when no topography data is loaded
            details.append(
                "  (Not available: topography file not loaded)"
            )
        # --- New Section: Dataset Details ---
        dataset_details = ["\nDataset Details (Per-Station Summary):"]
        df = self.info.df.copy()
        agg_cols = ['rho', 'phase', 'pc_rho', 's_phz']

        # 1. Aggregate main data to get per-station stats
        summary_df = df.groupby('station')[agg_cols].agg(
            ['min', 'max'])

        # 2. FIX: Robustly flatten the MultiIndex columns into strings
        # This is the most important change. It turns ('rho','min') into 'rho_min'
        summary_df.columns = [
            '_'.join(filter(None, col)).strip()
            for col in summary_df.columns.values
        ]
        summary_df.reset_index(inplace=True)

        # 3. Create the combined string columns (now using simple string keys)
        summary_df['Station'] = summary_df['station'].apply('{:,.1f}'.format)
        summary_df['ρ Range'] = summary_df.apply(
            lambda r: f"{r['rho_min']:.2e} - {r['rho_max']:.2e}", axis=1
        )
        summary_df['Φ Range'] = summary_df.apply(
            lambda r: f"{r['phase_min']:.2f} - {r['phase_max']:.2f}", axis=1
        )
        summary_df['%err Range'] = summary_df.apply(
            lambda r: f"{r['pc_rho_min']:.2f} - {r['pc_rho_max']:.2f}", axis=1
        )
        summary_df['σΦ Range'] = summary_df.apply(
            lambda r: f"{r['s_phz_min']:.2f} - {r['s_phz_max']:.2f}", axis=1
        )

        # Select only the formatted string columns for display
        display_df = summary_df[[
            'Station', 'ρ Range', 'Φ Range', '%err Range', 'σΦ Range'
        ]]

        # 4. Calculate dynamic column widths (this code now works correctly)
        headers = {
            'Station': 'Station', 'ρ Range': 'ρ Min - ρ Max',
            'Φ Range': 'Φ Min - Φ Max', '%err Range': '%err Min - %err Max',
            'σΦ Range': 'σΦ Min - σΦ Max'
        }
        col_widths = {
            col: max(len(headers[col]), display_df[col].str.len().max())
            for col in display_df.columns
        }

        # 4. Build the table strings
        header_row = " | ".join([
            headers[col].center(col_widths[col])
            for col in display_df.columns
        ])
        separator_row = "-+-".join(['-' * col_widths[col]
                                    for col in display_df.columns])

        dataset_details.append(header_row)
        dataset_details.append(separator_row)

        # Function to format a single data row
        def format_row(row):
            return " | ".join([
                row['Station'].ljust(col_widths['Station']),
                row['ρ Range'].center(col_widths['ρ Range']),
                row['Φ Range'].center(col_widths['Φ Range']),
                row['%err Range'].center(col_widths['%err Range']),
                row['σΦ Range'].center(col_widths['σΦ Range'])
            ])

        # 5. Handle Truncation and build data rows
        n_stations = len(display_df)
        if n_stations > 30:
            for _, row in display_df.head(5).iterrows():
                dataset_details.append(format_row(row))

            ellipsis_row = " | ".join(['...'.center(col_widths[col])
                                      for col in display_df.columns])
            dataset_details.append(ellipsis_row)

            for _, row in display_df.tail(5).iterrows():
                dataset_details.append(format_row(row))
        else:
            for _, row in display_df.iterrows():
                dataset_details.append(format_row(row))

        dataset_details.append(separator_row)

        # --- Combine and Return ---
        return "\n".join(
            header + [stats_str] + details + dataset_details
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
        path: str | Path,
        *,
        verbose: bool = False,
    ) -> AVG:
        r"""Load and parse an AVG file from a path.

        This classmethod is the primary factory for creating a
        fully populated :class:`AVG` instance directly from a file
        on disk. It handles the entire data loading pipeline.

        Parameters
        ----------
        path : str or pathlib.Path
            The filesystem path to the ``.avg`` file. Both legacy
            (kind-1) and modern (kind-2) formats are supported.
        verbose : bool, default False
            If ``True``, log messages will be printed to the console
            during the file loading and parsing process.

        Returns
        -------
        AVG
            A new, fully initialized instance of the :class:`AVG`
            class containing all the data and metadata from the
            file.

        Notes
        -----
        This method is a convenient wrapper around the more general
        :meth:`read` method. It creates an instance of the class
        and then immediately calls ``obj.read(path)`` to perform
        the actual data loading and component population.

        Examples
        --------
        >>> from pycsamt.zonge import AVG
        >>> # Load a modern AVG file
        >>> avg = AVG.from_file('data/avg/K2.avg', verbose=True)
        >>> print(avg.summary)

        See Also
        --------
        read : The underlying method that performs the data
            ingestion.
        """

        obj = cls(verbose=verbose)
        obj.read(path)
        return obj

    def to_xarray(
        self,
        *,
        include_qc: bool = True,
        include_z: bool = True,
    ):
        r"""Export the complete dataset to a single xarray.Dataset.

        This is the primary export method for creating a comprehensive,
        multi-dimensional, and self-describing dataset. It aggregates
        the primary data variables (resistivity, phase), the computed
        impedance tensor (Z), and all available quality control (QC)
        metrics into a single object.

        Parameters
        ----------
        include_qc : bool, default True
            If ``True``, all available QC components (e.g., `pc_emag`,
            `s_phz`) will be included as data variables in the
            dataset.
        include_z : bool, default True
            If ``True``, the real and imaginary parts of the complex
            impedance tensor (Z) and its propagated error will be
            included as data variables.

        Returns
        -------
        xarray.Dataset
            A dataset with coordinates ``(station, freq, comp)`` and
            data variables for all included measurements and QC
            metrics. The dataset's attributes are populated from the
            AVG file's header.

        Notes
        -----
        This method provides a high-level, convenient way to get all
        the processed data into a standard format for analysis and
        plotting. It gracefully handles cases where optional
        components (like certain QC metrics) are not available in the
        source data by skipping them.

        Examples
        --------
        >>> from pycsamt.zonge import AVG
        >>> avg = AVG.from_file('data/avg/K2.avg')
        >>> # Export all available data
        >>> ds = avg.to_xarray()
        >>> print(ds.data_vars)
        Data variables:
            rho        (station, freq, comp) float64 ...
            phase      (station, freq, comp) float64 ...
            z_real     (station, freq, comp) float64 ...
            ...

        See Also
        --------
        to_tensor : Export a single variable to a NumPy tensor.
        """
        import_optional_dependency(
            "xarray",
            extra="xarray is required for to_xarray()",
            errors="raise",
        )
        has_read (self)

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
                    self._logger.warning(
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
                        self._logger.info(
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
        station: int | float | None = None,
        agg: str | None = "mean",
        fill_value: float = np.nan,
        sort_freq: bool = True,
        align: str = "union",
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        r"""Export a specific variable as a 2x2 tensor.

        This is a convenience method that provides direct access to the
        underlying tensor-shaping capabilities of the data components.
        It delegates the call to the appropriate component (e.g., `Z`,
        `Resistivity`) to reshape the selected variable into a NumPy
        array.

        Parameters
        ----------
        var : {'z', 'z_real', 'z_imag', 'z_err', 'rho', 'phase'}, default 'z'
            The specific data variable to export as a tensor.
        station : int or float, optional
            If provided, the method returns a 3D tensor for only this
            station. If ``None``, it returns a 4D tensor containing
            data for all stations.
        agg : str or None, default 'mean'
            The aggregation method to use if duplicate measurements
            exist for the same station, frequency, and component.
        fill_value : float, default np.nan
            The value to use for missing elements in the tensor grid.
        sort_freq : bool, default True
            If ``True``, the frequency axis of the tensor will be
            sorted.
        align : {'union', 'intersection'}, default 'union'
            For multi-station tensors, determines how to handle
            different frequency sets between stations.

        Returns
        -------
        tensor : numpy.ndarray
            The data reshaped into a tensor. Shape is
            ``(n_freq, 2, 2)`` for a single station or
            ``(n_station, n_freq, 2, 2)`` for multiple stations.
        freqs : numpy.ndarray
            The frequency axis for the tensor.
        stations : numpy.ndarray
            The station axis for the tensor (empty for a single-
            station request).

        See Also
        --------
        to_xarray : Export the entire dataset to a labeled xarray object.
        pycsamt.zonge.tensor.TensorBase.to_tensor : The underlying
            implementation.
        """
        has_read (self)
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
    def df(self) -> pd.DataFrame | None:
        """
        Access the core tidy DataFrame containing all available
        data columns after parsing and normalization.
        """
        return self.info.df



class AMTAVG(AVG):
    r"""Extends AVG with tensor components and analytical methods.

    This class inherits from :class:`~pycsamt.zonge.avg.AVG` and
    provides a more specialized, analysis-focused interface. It is
    designed for users who need to perform advanced geophysical
    data processing and manipulation. It is designed for working
    directly with specific tensor elements (e.g., z_xy, rho_yx).

    The primary features of this class include:

    - Direct, property-based access to individual components of
      the impedance, resistivity, and phase tensors as pandas
      Series (e.g., ``z_xy``, ``res_yx_err``).
    - A suite of methods for advanced data processing, such as
      tensor rotation, phase unwrapping, and statistical analysis.

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
    calculate_statistics()
        Computes statistical QC metrics from repeated measurements.
    unwrap_phase()
        Corrects for 2π phase wrapping in the impedance phase data.
    rotate(angle_deg)
        Performs a 2D rotation of the impedance tensor.
    calculate_tipper(hz_data)
        Computes the Tipper transfer function from vertical field data.

    Examples
    --------
    >>> from pycsamt.zonge import AMTAVG
    >>> amt_avg = AMTAVG.from_file('data/avg/K2.avg')
    >>> # Rotate the data by 30 degrees clockwise
    >>> amt_avg.rotate(30)
    >>> # Access the newly rotated Z_xy component
    >>> z_xy_rotated = amt_avg.z_xy

    See Also
    --------
    AVG : The parent class providing the main data loading and
          exporting functionality.
    Z : The component that handles impedance calculations.
    Resistivity : The component that manages resistivity data.
    Phase : The component that manages phase data.
    """

    def compute_resistivity_phase(
        self, todeg: bool = False
    ) -> tuple[pd.Series, pd.Series, pd.Series, pd.Series]:
        r"""Compute rho and phi from the complex impedance Z.

        This method calculates the apparent resistivity (:math:`\rho_a`)
        and impedance phase (:math:`\phi`) from the complex
        impedance tensor (:math:`\mathbf{Z}`). It also propagates
        the errors from :math:`\mathbf{Z}` to the new values.

        Parameters
        ----------
        todeg : bool, default False
            If ``True``, the output impedance phase and its error
            will be converted from milliradians to degrees.

        Returns
        -------
        rho : pandas.Series
            The calculated apparent resistivity in ohm-meters.
        phi : pandas.Series
            The calculated impedance phase in milliradians or
            degrees.
        rho_err : pandas.Series
            The propagated relative error in apparent resistivity
            (in percent).
        phi_err : pandas.Series
            The propagated error in impedance phase (in the same
            units as `phi`).

        Notes
        -----
        The calculations are based on the fundamental MT
        relationships [1]_:

        .. math::
            \rho_a = \frac{1}{\omega \mu_0} |Z|^2
            \quad \text{and} \quad
            \phi = \arctan\left(\frac{\text{Im}(Z)}{\text{Re}(Z)}\right)

        Error propagation for :math:`\rho_a` is given by:

        .. math::
            \frac{\delta\rho_a}{\rho_a} = 2 \frac{\delta|Z|}{|Z|}

        References
        ----------
        .. [1] Chave, A. D., and Jones, A. G. (2012). *The
               Magnetotelluric Method: Theory and Practice*.
               Cambridge University Press.
        """

        has_read(self)

        z_complex = self.z.z
        freq = self.info.df["freq"]
        omega = 2 * PI * freq

        rho = (np.abs(z_complex)**2) / (omega * MU_0)
        phi_rad = np.arctan2(np.imag(z_complex), np.real(z_complex))

        # Error propagation
        z_err = self.z.z_err
        abs_z = np.abs(z_complex)
        safe_abs_z = np.where(abs_z == 0, np.nan, abs_z)

        rho_err = 200 * (z_err / safe_abs_z)
        phi_err_mrad = (
            1000
            * (z_err / safe_abs_z)
            * np.abs(np.sin(phi_rad))
        )

        if todeg:
            phi = np.rad2deg(phi_rad)
            phi_err = phi_err_mrad * (180.0 / (PI * 1000.0))
        else:
            phi = phi_rad * 1000.0
            phi_err = phi_err_mrad

        return rho, phi, rho_err, phi_err

    def set_resistivity_phase(
        self,
        rho: pd.Series,
        phi: pd.Series,
        rho_err: pd.Series | None = None,
        phi_err: pd.Series | None = None,
    ):
        r"""Attach new rho/phi data and reconstruct Z.

        This method updates the core DataFrame with new apparent
        resistivity and impedance phase values. After updating, it
        automatically re-initializes all data components, ensuring
        that dependent values (like the impedance tensor Z) are
        recalculated and the entire object state remains consistent.

        Parameters
        ----------
        rho : pandas.Series
            A Series containing the new apparent resistivity values.
            The index must align with the internal DataFrame.
        phi : pandas.Series
            A Series containing the new impedance phase values,
            assumed to be in milliradians.
        rho_err : pandas.Series, optional
            A Series containing the new relative error values for
            resistivity (in percent).
        phi_err : pandas.Series, optional
            A Series containing the new error values for phase (in
            milliradians).

        Notes
        -----
        This is a powerful method for data manipulation, such as
        applying static shifts or corrections. By calling
        `self.info.read` at the end, it ensures that all changes
        are propagated correctly throughout the entire data model.
        """

        has_read(self)
        df = self.info.df.copy()

        # Direct assignment is safe and respects index alignment
        df["rho"] = rho
        df["phase"] = phi
        if rho_err is not None:
            df["pc_rho"] = rho_err
        if phi_err is not None:
            df["s_phz"] = phi_err

        # Re-read components to update their internal state
        self.info.read(df, self.info.meta)

    def get_tensor_by_station(
        self,
        station_id: int | float,
        *,
        var: Literal["z", "rho", "phase"] = "z",
    ):
        r"""Fetch a 3D tensor for a single station using xarray.

        This method extracts the data for a single station from the
        full multi-station dataset and returns it as a 3D tensor
        in an :class:`xarray.DataArray`.

        Parameters
        ----------
        station_id : int or float
            The unique identifier for the station to retrieve. This
            must match one of the values in the 'station'
            coordinate of the dataset.
        var : {'z', 'rho', 'phase'}, default 'z'
            The specific data variable to retrieve.
            - 'z': Complex impedance tensor.
            - 'rho': Apparent resistivity tensor.
            - 'phase': Impedance phase tensor.

        Returns
        -------
        xarray.DataArray
            A 3D DataArray with dimensions ``(frequency, e, h)``
            containing the tensor data for the requested station.

        Raises
        ------
        KeyError
            If the specified `station_id` is not found in the
            dataset's station coordinates.
        ValueError
            If an unsupported `var` is requested.
        NotReadError
            If this method is called before data has been loaded.

        Notes
        -----
        This method leverages the powerful label-based indexing of
        xarray's ``.sel()`` method. It first constructs the full
        multi-station tensor in memory and then selects the slice
        corresponding to the given `station_id`.

        Examples
        --------
        >>> from pycsamt.zonge import AMTAVG
        >>> avg = AMTAVG.from_file('data/avg/K2.avg')
        >>> # Get the resistivity tensor for station 25.0
        >>> rho_tensor_25 = avg.get_tensor_by_station(25.0, var='rho')
        >>> print(rho_tensor_25.shape)
        (28, 2, 2)
        """
        has_read(self)

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
        # Create the full multi-station tensor first
        full_tensor = component.to_xarray_tensor(var=var)

        try:
            # Now, select the station. This will raise KeyError
            # if the station_id is not found in the coordinate.
            station_tensor = full_tensor.sel(station=station_id)
            return station_tensor
        except KeyError:
            # Use the instance's logger for error reporting
            if self.verbose:
                self._logger.error(
                    f"Station ID '{station_id}' not found in the "
                    "dataset's coordinates."
                )
            raise

    def calculate_statistics(
        self,
        *,
        drop_on_failure: bool = True,
        update_components: bool = True,
    ) -> pd.DataFrame:
        r"""Calculate statistics from grouped data.

        This method computes statistical metrics for repeated
        measurements, mirroring the process of the original Zonge
        AMTAVG program. It groups the data by station and frequency
        and calculates the mean, standard deviation (σ), and
        coefficient of variation (C-var) for key measurement
        columns.

        Parameters
        ----------
        drop_on_failure : bool, default True
            If ``True``, rows where statistical calculations
            cannot be performed (e.g., only one measurement in a
            group) will be excluded from the output DataFrame.
        update_components : bool, default True
            If ``True``, after calculating the statistics, this
            method will automatically update the corresponding QC
            components (e.g., `self.pc_rho`, `self.s_phz`) with the
            newly computed values.

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the calculated statistical
            metrics, including mean, standard deviation, and
            coefficient of variation for ``rho``, ``phase``,
            ``emag``, and ``hmag``.

        Notes
        -----
        The calculations are based on the formulas described in the
        AMTAVG manual [1]_. The coefficient of variation (C-var),
        which corresponds to the legacy ``%`` columns, is
        calculated as:

        .. math::
            C_{var} = 100 \cdot \frac{\sigma}{\bar{x}}

        where :math:`\sigma` is the standard deviation and
        :math:`\bar{x}` is the mean of the measurements.

        References
        ----------
        .. [1] Zonge Engineering (1996). *AMTAVG v7.2x User
               Manual*, Appendix B.
        """
        has_read(self)
        df = self.info.df

        # Define the columns for which to calculate statistics
        agg_cols = ["rho", "phase", "emag", "hmag"]
        # Ensure all required columns are present
        if not all(col in df.columns for col in agg_cols):
            raise ValueError(
                "DataFrame must contain 'rho', 'phase', "
                "'emag', and 'hmag' columns."
            )

        # Group by station and frequency to find repeated measurements
        grouped = df.groupby(["station", "freq"])

        # Calculate mean and standard deviation for each group
        stats = grouped[agg_cols].agg(['mean', 'std'])

        # Calculate Coefficient of Variation where applicable
        stats[('rho', 'cvar')] = (
            100 * stats[('rho', 'std')] / stats[('rho', 'mean')]
        )
        stats[('emag', 'cvar')] = (
            100 * stats[('emag', 'std')] / stats[('emag', 'mean')]
        )
        stats[('hmag', 'cvar')] = (
            100 * stats[('hmag', 'std')] / stats[('hmag', 'mean')]
        )

        # Flatten the multi-level column index
        stats.columns = ['_'.join(col) for col in stats.columns]
        stats = stats.reset_index()

        if drop_on_failure:
            stats = stats.dropna()

        if update_components:
            # Create a map from the new stats columns to the
            # canonical component names
            update_map = {
                'rho_cvar': 'pc_rho',
                'phase_std': 's_phz',
                'emag_cvar': 'pc_emag',
                'hmag_cvar': 'pc_hmag',
            }
            # Update the main DataFrame with the new QC values
            merged_df = pd.merge(
                self.info.df,
                stats,
                on=['station', 'freq'],
                how='left',
                suffixes=('', '_stat')
            )
            for stat_col, canon_col in update_map.items():
                # Check if column exists before updating
                # In principle, this is always true.
                if canon_col in self.info.df.columns:
                    if stat_col in merged_df.columns:
                        # Use combine_first to fill NaNs in the original
                        # column with new values from the stats column
                        self.info.df[canon_col] = merged_df[
                            canon_col
                        ].combine_first(merged_df[stat_col])
                elif self.verbose:
                    self._logger.info(
                        f"Skipping update for '{canon_col}': "
                        "column not found in source DataFrame."
                    )

            # Re-read the info object to update all components
            self.info.read(self.info.df, self.info.meta)

        return stats

    def unwrap_phase(
        self,
        *,
        todeg: bool = False,
        update_components: bool = True,
    ) -> pd.DataFrame:
        r"""Unwraps the impedance phase to correct for 2π jumps.

        This method implements a phase referencing/unwrapping
        algorithm, similar to the `PHASEREF` mode in the original
        Zonge AMTAVG program. It corrects for phase wrapping by
        ensuring that the absolute difference between consecutive
        phase values (sorted by frequency) does not exceed π
        radians (1000π milliradians).

        Parameters
        ----------
        todeg : bool, default False
            If ``True``, the final unwrapped phase values in the
            DataFrame will be converted from milliradians to
            degrees.
        update_components : bool, default True
            If ``True``, the main DataFrame (`self.info.df`) is
            updated in place with the unwrapped phase values, and
            all data components are re-read to reflect this change.

        Returns
        -------
        pandas.DataFrame
            A new DataFrame with the unwrapped phase values. The
            main `self.info.df` is also updated if
            `update_components` is ``True``.

        Notes
        -----
        The unwrapping is performed independently for each unique
        station and component pair. The data is first sorted by
        frequency to ensure the correct sequence for unwrapping.

        The core of the algorithm uses :func:`numpy.unwrap`, which
        operates in radians. The phase values, assumed to be in
        milliradians, are converted to radians before unwrapping
        and then converted back to the desired output unit.

        Examples
        --------
        >>> from pycsamt.zonge import AMTAVG
        >>> avg = AMTAVG.from_file('data/avg/K2.avg')
        >>> # Unwrap phase and update the object in place
        >>> avg.unwrap_phase()
        >>> # Get the unwrapped phase for the xy component
        >>> unwrapped_phi_xy = avg.phase_xy

        See Also
        --------
        numpy.unwrap : The underlying function used for unwrapping.
        """

        has_read(self)
        df = self.info.df.copy()

        # # Convert phase from milliradians to radians for unwrapping
        # phase_rad = df["phase"] * 1e-3

        # Group by station and component, then apply unwrap to each
        unwrapped_rad = df.groupby(
            ["station", "comp"]
        )["phase"].transform(
            lambda s: np.unwrap(s * 1e-3) # Convert to rad before unwrap
        )

        if todeg:
            # Convert unwrapped radians to degrees
            df["phase"] = unwrapped_rad * (180.0 / PI)
            unit = "deg"
        else:
            # Convert unwrapped radians back to milliradians
            df["phase"] = unwrapped_rad * 1000.0
            unit = "mrad"

        if update_components:
            # Update the main DataFrame and all components
            self.info.df["phase"] = df["phase"]
            self.info.meta["Unit.Phase"] = unit
            self.info.read(self.info.df, self.info.meta)

        return df

    def rotate(
        self,
        angle: float,
        *,
        update_components: bool = True,
    ):
        r"""Rotate the impedance tensor by a given angle.

        This method performs a 2D rotation of the impedance tensor
        :math:`\mathbf{Z}` for all measurements. The rotation is
        applied in a clockwise-positive direction, which is standard
        in many geophysical applications.

        Parameters
        ----------
        angle : float
            The rotation angle in degrees, where a positive angle
            corresponds to a clockwise rotation of the coordinate
            system.
        update_components : bool, default True
            If ``True``, the main DataFrame (`self.info.df`) is
            updated in place with the new, rotated impedance values,
            and all dependent components (like resistivity and
            phase) are recalculated and updated.

        Returns
        -------
        pandas.DataFrame
            A new DataFrame containing the four rotated complex
            impedance components (``z_xx_rot``, ``z_xy_rot``,
            ``z_yx_rot``, ``z_yy_rot``). The main `self.info.df` is
            also updated if `update_components` is ``True``.

        Notes
        -----
        The rotation is performed using the standard 2D tensor
        rotation matrix [1]_:

        .. math::
            \mathbf{Z'} = \mathbf{R}(\theta) \mathbf{Z}
            \mathbf{R}^T(\theta)

        where :math:`\mathbf{Z}` is the original impedance tensor,
        :math:`\mathbf{Z'}` is the rotated tensor, :math:`\theta`
        is the rotation angle, and :math:`\mathbf{R}(\theta)` is
        the rotation matrix:

        .. math::
            \mathbf{R}(\theta) = \begin{pmatrix}
            \cos\theta & \sin\theta \\
            -\sin\theta & \cos\theta
            \end{pmatrix}

        This operation updates the ``z_xx``, ``z_xy``, ``z_yx``,
        and ``z_yy`` components. If `update_components` is
        ``True``, the apparent resistivity and phase are then
        recalculated from this new rotated impedance tensor.

        References
        ----------
        .. [1] Simpson, F., & Bahr, K. (2005). *Practical
               Magnetotellurics*. Cambridge University Press.

        See Also
        --------
        compute_resistivity_phase : Recalculates rho and phi from Z.
        """
        has_read(self)
        df = self.info.df.copy()

        # Convert angle to radians for trigonometric functions
        angle_rad = np.deg2rad(angle)
        c = np.cos(angle_rad)
        s = np.sin(angle_rad)
        c2 = c * c
        s2 = s * s
        cs = c * s

        # The component accessors are row-subsets of ``df`` (they keep
        # the original row index), so the four series carry *disjoint*
        # indices. Re-key everything on (station, freq) before mixing
        # components; a component that was not measured for a sounding
        # contributes zero to the rotation (scalar-CSAMT convention).
        keys = pd.MultiIndex.from_frame(df[["station", "freq"]])

        def _keyed(series: pd.Series) -> pd.Series:
            sub = df.loc[series.index, ["station", "freq"]]
            out = pd.Series(
                series.to_numpy(),
                index=pd.MultiIndex.from_frame(sub),
                dtype=complex,
            )
            return out[~out.index.duplicated(keep="first")]

        pair_index = keys.unique()
        zxx = _keyed(self.z.z_xx).reindex(pair_index, fill_value=0j)
        zxy = _keyed(self.z.z_xy).reindex(pair_index, fill_value=0j)
        zyx = _keyed(self.z.z_yx).reindex(pair_index, fill_value=0j)
        zyy = _keyed(self.z.z_yy).reindex(pair_index, fill_value=0j)

        # Perform the 2D tensor rotation (aligned on station/freq)
        zxx_rot = c2 * zxx + cs * (zxy + zyx) + s2 * zyy
        zxy_rot = -cs * zxx + c2 * zxy - s2 * zyx + cs * zyy
        zyx_rot = -cs * zxx - s2 * zxy + c2 * zyx + cs * zyy
        zyy_rot = s2 * zxx - cs * (zxy + zyx) + c2 * zyy

        # Create a new DataFrame for the rotated values
        rotated_df = pd.DataFrame({
            "z_xx_rot": zxx_rot,
            "z_xy_rot": zxy_rot,
            "z_yx_rot": zyx_rot,
            "z_yy_rot": zyy_rot,
        })

        if update_components:
            # Pick each row's rotated tensor element by its component
            # label (labels are stored upper-cased) and its
            # (station, freq) key.
            comp_map = {
                "EXHX": zxx_rot, "ZXX": zxx_rot,
                "EXHY": zxy_rot, "ZXY": zxy_rot,
                "EYHX": zyx_rot, "ZYX": zyx_rot,
                "EYHY": zyy_rot, "ZYY": zyy_rot,
            }
            z_new = np.full(len(df), np.nan + 0j, dtype=complex)
            comp_norm = df["comp"].astype(str).str.upper().to_numpy()
            for label, rot in comp_map.items():
                mask = comp_norm == label
                if mask.any():
                    z_new[mask] = rot.reindex(keys[mask]).to_numpy()

            # Recalculate rho and phi from the new complex Z
            omega = 2 * PI * df["freq"].to_numpy(dtype=float)
            abs_z_sq = np.abs(z_new) ** 2

            self.info.df["rho"] = abs_z_sq / (omega * MU_0)
            self.info.df["phase"] = (
                np.arctan2(z_new.imag, z_new.real) * 1000.0  # to mrad
            )

            # Re-read the info object to update all components
            self.info.read(self.info.df, self.info.meta)

        return rotated_df

    def calculate_tipper(
        self,
        hz_data: pd.Series,
        *,
        update_components: bool = True,
    ) -> pd.DataFrame:
        r"""Calculate the Tipper transfer function.

        This method computes the complex Tipper components (Tx, Ty)
        by solving the linear relationship:
        :math:`H_z = T_x H_x + T_y H_y`

        Parameters
        ----------
        hz_data : pandas.Series
            A Series containing the complex vertical magnetic field
            (:math:`H_z`) data. The index of this Series must match
            the index of the main DataFrame (`self.info.df`).
        update_components : bool, default True
            If ``True``, a new `Tipper` component will be created and
            attached to the `self.info` object, populated with the
            calculated Tx and Ty values.

        Returns
        -------
        pandas.DataFrame
            A new DataFrame containing the calculated complex Tipper
            components, `tx` and `ty`.

        Notes
        -----
        The calculation is performed for each unique station and
        frequency using a least-squares approach to solve for
        :math:`T_x` and :math:`T_y`. This requires at least two
        linearly independent horizontal field measurements for a
        stable solution.

        This method requires the horizontal magnetic field data
        (`hmag`, `hphz`) to be present in the dataset.
        """
        from .tipper import Tipper

        has_read(self)
        df = self.info.df.copy()

        if "hmag" not in df.columns or "hphz" not in df.columns:
            raise AvgDataError(
                "Horizontal magnetic field data (hmag, hphz) is "
                "required to calculate the Tipper."
            )
        # Reconstruct complex Hx and Hy from magnitude and phase
        h_mag = df["hmag"]
        h_phase_rad = df["hphz"] * 1e-3
        h_complex = h_mag * np.exp(1j * h_phase_rad)

        # Create Hx and Hy based on component type
        df["hx"] = np.where(
            df["comp"].isin(["ExHx", "EyHx"]), h_complex, 0
        )
        df["hy"] = np.where(
            df["comp"].isin(["ExHy", "EyHy"]), h_complex, 0
        )
        df["hz"] = hz_data

        tipper_results = []
        # Group by each measurement point
        for (station, freq), group in df.groupby(
            ["station", "freq"]
        ):
            # We need at least two equations (two components)
            if len(group) < 2:
                continue

            # Set up the linear system A*x = b
            # A = [[Hx1, Hy1], [Hx2, Hy2], ...]
            # x = [[Tx], [Ty]]
            # b = [[Hz1], [Hz2], ...]
            A = group[["hx", "hy"]].values
            b = group["hz"].values

            try:
                # Solve for x using least squares
                x, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
                tipper_results.append({
                    "station": station,
                    "freq": freq,
                    "tx": x[0],
                    "ty": x[1],
                })
            except np.linalg.LinAlgError:
                if self.verbose:
                    self._logger.warning(
                        f"Could not solve for Tipper at station "
                        f"{station}, freq {freq}."
                    )

        if not tipper_results:
            warnings.warn("Tipper calculation resulted in no data.", stacklevel=2)
            return pd.DataFrame()

        tipper_df = pd.DataFrame(tipper_results)

        if update_components:
            # Create and populate the new Tipper component
            self.info.tipper = Tipper(verbose=self.verbose)
            self.info.tipper.read(tipper_df, self.info.meta)

        return tipper_df

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
        return df[df.comp == "ExHx"]["s_phz"]
    @property
    def phase_xy_err(self):
        df = self.info.s_phz.frame
        return df[df.comp == "ExHy"]["s_phz"]
    @property
    def phase_yx_err(self):
        df = self.info.s_phz.frame
        return df[df.comp == "EyHx"]["s_phz"]
    @property
    def phase_yy_err(self):
        df = self.info.s_phz.frame
        return df[df.comp == "EyHy"]["s_phz"]
