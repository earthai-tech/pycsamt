# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
stratagem.survey
================

End-to-end Stratagem AMT survey processing pipeline.

:class:`StratagemSurvey` replaces the legacy ``watex``-based
``stratagem_edi_process_script.py`` with a single composable class that
covers the full post-acquisition workflow:

.. code-block:: text

   WinGLink EDI directory
     ↓  EDIBatch              load & natural-sort
     ↓  CoordinateInjector    inject GPS coords from CSV / XLS / XLSX
     ↓  QualityController     per-station QC report (optional)
     ↓  StaticShiftCorrector  AMA static-shift removal
     ↓  FrequencyFilter       band select + incoherent-freq masking
     ↓  NoiseRemover          powerline notch + Hampel + smoothing
     ↓  export()              write corrected EDIs
     ↓  rename()              standardise filenames (T2.000.edi …)

Every processing step is optional and chainable.  The pipeline stores
intermediate objects as attributes so individual results are always
accessible.

Usage example
-------------
Replicate the old ``stratagem_edi_process_script.py`` in four lines:

>>> from pycsamt.stratagem import StratagemSurvey
>>> sv = StratagemSurvey(
...     edi_dir="2/2EDI",
...     coord_file="2.csv",
...     raw_dir="原始数据/2HX",
...     epsg=32649,
... ).fit()
>>> sv.remove_static_shift().remove_noises().drop_frequencies(fmin=10.0)
>>> sv.export("2/2EDIP").rename(basename="T2.", dst_path="2/renamedEDIs")
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Union

from ..api.property import MetadataMixin, PyCSAMTObject
from ..emtools._core import ensure_sites
from ..exceptions import NotFittedError
from ..seg.edi import EDIFile
from ..site.base import Sites
from .gis_correct import CoordinateInjector
from .io import EDIBatch, StratagemRawReader
from .process import NoiseRemover, StaticShiftCorrector
from .qc import FrequencyFilter, QualityController
from .rename import EDIRenamer, EDIWriter

__all__ = ["StratagemSurvey"]

#: Accepted forms for the EDI source: a WinGLink export directory (uses
#: Stratagem's own natural 3-digit sort), or anything the conventional
#: ``pycsamt.emtools``/``pycsamt.site`` stack already understands (a
#: :class:`~pycsamt.site.base.Sites`, a list of
#: :class:`~pycsamt.seg.edi.EDIFile`, or an ``EDICollection``).
EdiSource = Union[str, Path, Sites, Sequence[EDIFile]]


class StratagemSurvey(PyCSAMTObject, MetadataMixin):
    """End-to-end Stratagem AMT survey processing pipeline.

    Parameters
    ----------
    edi_dir : path-like, Sites, or sequence of EDIFile
        Directory of WinGLink-exported EDI files (uses Stratagem's own
        natural 3-digit sort via :class:`~pycsamt.stratagem.io.EDIBatch`),
        or an already-loaded :class:`~pycsamt.site.base.Sites` /
        ``EDICollection`` / list of :class:`~pycsamt.seg.edi.EDIFile` —
        normalised through :func:`~pycsamt.emtools._core.ensure_sites`,
        the same entry point the rest of ``pycsamt.emtools`` uses. In the
        latter case :attr:`batch_` stays ``None`` (there's no directory
        to report) and ordering is whatever the source already has.
    coord_file : path-like
        GPS coordinate table (CSV / XLS / XLSX).
    raw_dir : path-like, optional
        Directory of raw Stratagem hardware files (``X*.NNN``, …).
        When supplied, hardware SNR masks are used in QC and frequency
        filtering.
    epsg : int, default 32649
        EPSG code of the projected CRS of *coord_file*.  Use ``32649``
        (UTM Zone 49N WGS84) for the standard south-China survey area.
    utm_zone : str, default ``'49N'``
        UTM zone string for ``project_point_utm2ll`` when *epsg* is not
        sufficient.
    coordinate_system : str, default ``'utm'``
    order : str, default ``'auto'``
        Station-to-GPS row ordering for
        :class:`~pycsamt.stratagem.gis_correct.StationLocator`.
    drop_stations : list of int, optional
        0-based indices into the loaded ``EDIBatch`` to exclude *before*
        coordinate injection — e.g. a calibration/test shot that isn't a
        real profile position and has no matching row in *coord_file*.
    easting_col, northing_col : str, optional
        Column names in *coord_file* for the projected E-W / N-S
        coordinates.  Forwarded to
        :meth:`~pycsamt.stratagem.gis_correct.CoordinateInjector.fit`.
        Required whenever *coord_file* has more than two numeric columns
        besides *elev_col* — auto-detection raises rather than guessing
        in that case (see :mod:`pycsamt.stratagem.gis_correct`).
    elev_col : str, default ``'elev'``
    station_col : str, default ``'stations'``
    read_kwargs : dict, optional
        Extra keyword arguments forwarded to the *coord_file* reader.
    verbose : int, default 0

    Attributes
    ----------
    batch_ : EDIBatch or None
        Loaded EDI collection when *edi_dir* was a directory path;
        ``None`` when *edi_dir* was already a ``Sites``/list of EDIFile.
    raw_reader_ : StratagemRawReader or None
        Hardware file reader (None when *raw_dir* not supplied).
    injector_ : CoordinateInjector
        Coordinate-injected EDI wrapper.
    qc_ : QualityController or None
        QC report (populated after :meth:`run_qc`).
    edi_objects_ : list of EDIFile
        Current working set of EDI objects.  Modified in-place by each
        processing step.

    Examples
    --------
    Full pipeline, one fluent expression:

    >>> sv = (
    ...     StratagemSurvey(
    ...         edi_dir="2/2EDI",
    ...         coord_file="2.csv",
    ...         raw_dir="原始数据/2HX",
    ...         epsg=32649,
    ...     )
    ...     .fit()
    ...     .run_qc()
    ...     .remove_static_shift()
    ...     .drop_frequencies(fmin=10.0)
    ...     .remove_noises()
    ...     .export("2/2EDIP")
    ...     .rename(basename="T2.", dst_path="2/renamedEDIs")
    ... )
    >>> print(sv.qc_.summary())
    """

    __repr_fields__ = (
        "edi_dir",
        "epsg",
        "utm_zone",
        "n_stations_",
    )

    def __init__(
        self,
        edi_dir: EdiSource,
        coord_file: str | Path,
        *,
        raw_dir: str | Path | None = None,
        epsg: int = 32649,
        utm_zone: str = "49N",
        coordinate_system: str = "utm",
        order: str = "auto",
        drop_stations: list[int] | None = None,
        easting_col: str | None = None,
        northing_col: str | None = None,
        elev_col: str = "elev",
        station_col: str = "stations",
        read_kwargs: dict | None = None,
        verbose: int = 0,
    ) -> None:
        self.edi_dir = edi_dir
        self.coord_file = coord_file
        self.raw_dir = raw_dir
        self.epsg = epsg
        self.utm_zone = utm_zone
        self.coordinate_system = coordinate_system
        self.order = order
        self.drop_stations = drop_stations
        self.easting_col = easting_col
        self.northing_col = northing_col
        self.elev_col = elev_col
        self.station_col = station_col
        self.read_kwargs = read_kwargs
        self.verbose = verbose

        # populated by fit()
        self.batch_: EDIBatch | None = None
        self.raw_reader_: StratagemRawReader | None = None
        self.injector_: CoordinateInjector | None = None
        self.qc_: QualityController | None = None
        self.edi_objects_: list | None = None

    # ------------------------------------------------------------------
    # mandatory step
    # ------------------------------------------------------------------

    def fit(self) -> StratagemSurvey:
        """Load EDIs, optional raw files, and inject GPS coordinates.

        This is the only mandatory step.  All processing methods
        (:meth:`remove_static_shift`, :meth:`drop_frequencies`, etc.)
        must be called after :meth:`fit`.

        Returns
        -------
        self
        """
        # ── 1. load EDIs ─────────────────────────────────────────────
        if isinstance(self.edi_dir, (str, Path)):
            self.batch_ = EDIBatch(
                self.edi_dir,
                verbose=self.verbose,
            ).fit()
            n_loaded = len(self.batch_)
            edi_list = list(self.batch_.edi_objects_)
        else:
            # already a Sites / EDICollection / list of EDIFile — hand it
            # to the same normaliser the rest of pycsamt.emtools uses.
            self.batch_ = None
            sites = ensure_sites(self.edi_dir, verbose=self.verbose)
            n_loaded = len(sites)
            edi_list = sites.as_list()

        # drop non-survey stations (e.g. a calibration/test shot) by
        # 0-based index *before* coordinate injection, so the coordinate
        # table only needs one row per real station.
        if self.drop_stations:
            drop = set(self.drop_stations)
            edi_list = [e for i, e in enumerate(edi_list) if i not in drop]

        self.n_stations_ = len(edi_list)

        if self.verbose:
            print(
                f"[StratagemSurvey] loaded {n_loaded} EDI files"
                + (
                    f", using {self.n_stations_} after dropping "
                    f"{n_loaded - self.n_stations_} station(s)"
                    if self.drop_stations
                    else ""
                )
            )

        # ── 2. load raw hardware files (optional) ─────────────────────
        if self.raw_dir is not None:
            try:
                self.raw_reader_ = StratagemRawReader(
                    self.raw_dir,
                    verbose=self.verbose,
                ).fit()
                if self.verbose:
                    cov = self.raw_reader_.station_coverage()
                    print(
                        f"[StratagemSurvey] raw reader: "
                        f"{self.raw_reader_.n_stations_} stations, "
                        f"coverage={cov:.1%}"
                    )
            except Exception as exc:
                if self.verbose:
                    print(f"[StratagemSurvey] raw reader failed: {exc}")
                self.raw_reader_ = None

        # ── 3. inject GPS coordinates ─────────────────────────────────
        self.injector_ = CoordinateInjector(
            coordinate_system=self.coordinate_system,
            epsg=self.epsg,
            utm_zone=self.utm_zone,
            order=self.order,
            verbose=self.verbose,
        ).fit(
            edi_list,
            self.coord_file,
            easting_col=self.easting_col,
            northing_col=self.northing_col,
            elev_col=self.elev_col,
            station_col=self.station_col,
            read_kwargs=self.read_kwargs,
        )

        self.edi_objects_ = self.injector_.edi_objects_

        if self.verbose:
            print(
                f"[StratagemSurvey] coordinates injected, "
                f"reversed={self.injector_.reversed_}"
            )

        return self

    # ------------------------------------------------------------------
    # optional processing steps (all chainable)
    # ------------------------------------------------------------------

    def run_qc(
        self,
        *,
        min_frac_ok: float = 0.6,
        min_snr_med: float = 2.0,
        max_skew_med: float = 6.0,
        include_skew: bool = True,
    ) -> StratagemSurvey:
        """Run the station-level QC report.

        Results stored in :attr:`qc_`.  Does not modify ``Z`` data.

        Returns
        -------
        self
        """
        self._require_fit()
        self.qc_ = QualityController(
            min_frac_ok=min_frac_ok,
            min_snr_med=min_snr_med,
            max_skew_med=max_skew_med,
            include_skew=include_skew,
            verbose=self.verbose,
        ).fit(self.edi_objects_, raw_reader=self.raw_reader_)

        if self.verbose:
            n_flagged = len(self.qc_.flagged_stations())
            print(
                f"[StratagemSurvey] QC: {len(self.qc_.report_)} stations, "
                f"{n_flagged} flagged"
            )
        return self

    def remove_static_shift(
        self,
        *,
        sort_by: str = "lon",
        half_window: int = 3,
        weights: str = "tri",
        pband: tuple | None = None,
        max_skew: float | None = 6.0,
    ) -> StratagemSurvey:
        """Apply AMA static-shift correction.

        .. important::
           Call this **before** :meth:`drop_frequencies` to ensure the
           full frequency range is available for spatial averaging.

        Returns
        -------
        self
        """
        self._require_fit()
        sc = StaticShiftCorrector(
            sort_by=sort_by,
            half_window=half_window,
            weights=weights,
            pband=pband,
            max_skew=max_skew,
            verbose=self.verbose,
        ).fit(self.edi_objects_)
        self._ss_corrector_ = sc

        if self.verbose:
            med = sc.factors_["fac_z"].median()
            print(
                f"[StratagemSurvey] static shift: "
                f"{len(sc.factors_)} stations, median fac_z={med:.3f}"
            )
        return self

    def drop_frequencies(
        self,
        *,
        fmin: float | None = None,
        fmax: float | None = None,
        snr_thresh: float = 2.5,
        min_frac: float = 0.4,
        use_hardware_mask: bool = True,
    ) -> StratagemSurvey:
        """Filter frequency bands and mask incoherent bins.

        Parameters
        ----------
        fmin, fmax : float, optional
            Frequency band limits in Hz.
        snr_thresh : float
            Per-station SNR threshold for incoherence masking.
        min_frac : float
            Minimum fraction of stations that must pass *snr_thresh*.
        use_hardware_mask : bool
            Apply hardware SNR mask when raw files were loaded.

        Returns
        -------
        self
        """
        self._require_fit()
        ff = FrequencyFilter(
            fmin=fmin,
            fmax=fmax,
            snr_thresh=snr_thresh,
            min_frac=min_frac,
            use_hardware_mask=use_hardware_mask,
            verbose=self.verbose,
        ).fit(self.edi_objects_, raw_reader=self.raw_reader_)
        self._freq_filter_ = ff

        if self.verbose:
            print(
                f"[StratagemSurvey] freq filter: "
                f"hw={ff.n_masked_hw_}, band={ff.n_dropped_band_}, "
                f"incoherent={ff.n_masked_stat_}"
            )
        return self

    def remove_noises(
        self,
        *,
        mains_hz: float = 50.0,
        n_harm: int = 30,
        tol_hz: float = 0.08,
        notch_mode: str = "interp",
        hampel_win: int = 3,
        hampel_nsig: float = 3.0,
        smooth: bool = False,
        smooth_win: int = 3,
    ) -> StratagemSurvey:
        """Apply powerline notch + Hampel outlier + optional smoothing.

        Returns
        -------
        self
        """
        self._require_fit()
        nr = NoiseRemover(
            mains_hz=mains_hz,
            n_harm=n_harm,
            tol_hz=tol_hz,
            notch_mode=notch_mode,
            hampel_win=hampel_win,
            hampel_nsig=hampel_nsig,
            smooth=smooth,
            smooth_win=smooth_win,
            verbose=self.verbose,
        ).fit(self.edi_objects_)
        self._noise_remover_ = nr

        if self.verbose:
            print(
                f"[StratagemSurvey] noise removal: "
                f"{nr.n_stations_} stations processed"
            )
        return self

    # ------------------------------------------------------------------
    # output steps
    # ------------------------------------------------------------------

    def export(
        self,
        savepath: str | Path,
        *,
        dataid_prefix: str | None = None,
        overwrite: bool = False,
    ) -> StratagemSurvey:
        """Write the current ``edi_objects_`` to *savepath*.

        Parameters
        ----------
        savepath : path-like
            Output directory (created if absent).
        dataid_prefix : str, optional
            When given, ``>HEAD.DATAID`` is standardised to
            ``{dataid_prefix}{i:03d}`` before writing.
        overwrite : bool, default False

        Returns
        -------
        self
        """
        self._require_fit()
        writer = EDIWriter(
            dataid_prefix=dataid_prefix,
            overwrite=overwrite,
            verbose=self.verbose,
        ).fit(self.edi_objects_, savepath)
        self._last_export_dir_ = Path(savepath).expanduser().resolve()
        self._writer_ = writer

        if self.verbose:
            print(
                f"[StratagemSurvey] export: "
                f"{writer.n_written_} files → {self._last_export_dir_}"
            )
        return self

    def rename(
        self,
        basename: str,
        dst_path: str | Path,
        *,
        zero_pad: int = 3,
        trailer: str = "",
        overwrite: bool = False,
        source: str | Path | None = None,
    ) -> StratagemSurvey:
        """Rename EDI files with a standardised basename.

        Parameters
        ----------
        basename : str
            Filename prefix, e.g. ``'T2.'`` → ``T2.000.edi``.
        dst_path : path-like
            Output directory for renamed files.
        zero_pad : int, default 3
        trailer : str, default ``''``
        overwrite : bool, default False
        source : path-like, optional
            Source directory or list.  Defaults to the directory written
            by the most recent :meth:`export` call; falls back to the
            current ``edi_objects_``.

        Returns
        -------
        self
        """
        self._require_fit()

        src: str | Path | list
        if source is not None:
            src = source
        elif hasattr(self, "_last_export_dir_"):
            src = self._last_export_dir_
        else:
            src = self.edi_objects_

        renamer = EDIRenamer(
            basename=basename,
            zero_pad=zero_pad,
            trailer=trailer,
            overwrite=overwrite,
            verbose=self.verbose,
        ).fit(src, dst_path)
        self._renamer_ = renamer

        if self.verbose:
            print(
                f"[StratagemSurvey] renamed {renamer.n_renamed_} files → "
                f"{Path(dst_path).expanduser().resolve()}"
            )
        return self

    # ------------------------------------------------------------------
    # convenience helpers
    # ------------------------------------------------------------------

    def summary(self) -> str:
        """Return a human-readable pipeline status summary."""
        self._require_fit()

        lines = [
            "StratagemSurvey",
            f"  edi_dir    : {self.edi_dir}",
            f"  coord_file : {self.coord_file}",
            f"  n_stations : {self.n_stations_}",
            f"  epsg / zone: {self.epsg} / {self.utm_zone}",
            f"  raw_reader : {'yes (' + str(self.raw_reader_.n_stations_) + ' stations)' if self.raw_reader_ else 'not loaded'}",
        ]
        if self.injector_ is not None:
            lines.append(
                f"  coord order: {'reversed' if self.injector_.reversed_ else 'forward'}"
            )
        if self.qc_ is not None:
            n_flag = len(self.qc_.flagged_stations())
            lines.append(
                f"  QC flags   : {n_flag} / {len(self.qc_.report_)} stations"
            )
        if hasattr(self, "_ss_corrector_"):
            med = self._ss_corrector_.factors_["fac_z"].median()
            lines.append(f"  SS fac_z   : median={med:.3f}")
        if hasattr(self, "_freq_filter_"):
            ff = self._freq_filter_
            lines.append(
                f"  freq filter: hw={ff.n_masked_hw_}, "
                f"band={ff.n_dropped_band_}, incoh={ff.n_masked_stat_}"
            )
        if hasattr(self, "_noise_remover_"):
            lines.append(
                f"  noise rm   : {self._noise_remover_.n_stations_} stations"
            )
        if hasattr(self, "_writer_"):
            lines.append(
                f"  export     : {self._writer_.n_written_} files → {self._last_export_dir_}"
            )
        if hasattr(self, "_renamer_"):
            lines.append(f"  rename     : {self._renamer_.n_renamed_} files")
        return "\n".join(lines)

    @property
    def coordinate_frame(self):
        """DataFrame of WGS84 coordinates (requires :meth:`fit`)."""
        self._require_fit()
        return self.injector_.coordinate_frame()

    @property
    def sites_(self) -> Sites:
        """Current :attr:`edi_objects_` wrapped as a :class:`~pycsamt.site.base.Sites`.

        A fresh view built on every access, so it always reflects the
        current pipeline state (post-QC, post-static-shift, etc.). This is
        the interop point with the conventional ``pycsamt.emtools`` /
        ``pycsamt.site`` stack — e.g. use ``sv.sites_.write(outdir)`` for
        the generic ``{station}.edi`` writer instead of Stratagem's own
        :meth:`export`/:meth:`rename` (which additionally handle DATAID
        prefixing and Stratagem's zero-padded naming convention).

        Examples
        --------
        >>> sv.sites_.write("out_dir", exist_ok=True)
        """
        self._require_fit()
        return Sites(self.edi_objects_)

    # ------------------------------------------------------------------
    # internal
    # ------------------------------------------------------------------

    def _require_fit(self) -> None:
        if self.edi_objects_ is None:
            raise NotFittedError(
                "Call StratagemSurvey.fit() before running processing steps."
            )
