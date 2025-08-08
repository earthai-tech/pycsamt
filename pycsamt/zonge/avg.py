# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""
Level façade around :pymod:`pycsamt.zonge.utils`.

* Loads any 'kind-1 / kind-2' Zonge AVG file into a tidy
  :class:`pandas.DataFrame`.
* Keeps header metadata in a plain ``dict``.
* Provides a symmetrical :py:meth:`write` helper to re-emit a
  *kind-2* CSV AVG.
"""
from __future__ import annotations

from pathlib import Path
from typing  import ( 
    Iterable, 
    Literal, 
    Sequence, 
    Any, 
    Callable
)
from datetime import datetime, timezone
import numpy as np 
import pandas as pd
from scipy.ndimage import uniform_filter1d

from ..decorators import has_fit
from ..exceptions import ( 
    AvgFileError, 
    AvgDataError, 
    ResistivityError, 
    PhaseError, 
    ZError, 
    
)
from ..log.logger  import get_logger
from ..utils.validation import check_is_fitted 

from .base import BaseAVG, OpsBase
from .heads import Head
from .properties import (
    Hardware,
    SurveyAnnotation,
    SurveyConfiguration,
    Receiver,
    Transmitter,
    SkipFlag,
)
from .meas import (
    CompMeas,
    Amps,
    Frequency,
)
from .resphase import Resistivity, Phase
from .survey   import Station
from .var      import (
    PcEmag,
    PcHmag,
    PcRho,
    SPhz,
    SHphz,
    SEphz,
)
from .info import DataInfo 
from .utils        import (
    load_avg, 
    write_avg, 
    extract_core_columns,
)
from .z import Z

logger = get_logger(__name__)

__all__ = ["AVG", "AVGOps"]


@has_fit("raise")
class AVG(BaseAVG):
    """
    High-level façade around a Zonge **AVG / AMTAVG** file.

    An :class:`AVG` instance loads the raw text file (or an already
    parsed :class:`~pandas.DataFrame`), normalises the column names,
    *and* hydrates every helper object – stations, frequencies,
    impedance, quality metrics, header blocks, and so on.  This makes
    the whole survey instantly accessible through strongly-typed
    attributes instead of ad-hoc column look-ups [1]_.

    Parameters
    ----------
    obj_or_path : str or pathlib.Path or pandas.DataFrame, optional
        * Path to a ``.avg`` file *or*
        * A prepared DataFrame that already follows the canonical
          schema.  When this form is used the caller becomes
          responsible for column consistency.
        * If *None* the instance is created empty; call
          :py:meth:`read` or :py:meth:`fit` later.
    meta : dict, optional
        Header key–value pairs that overwrite entries read from the
        file.  Ignored when *obj_or_path* is a DataFrame unless you
        explicitly supply it.
    verbose : bool, default ``False``
        Emit INFO messages during load / write operations.


    Attributes
    ----------
    Head : :class:`pycsamt.zonge.heads.Head`
        Aggregates *Hardware*, *SurveyAnnotation*, *SurveyConfiguration*
        and the *Rx / Tx* property blocks plus the skip flag.
    Station, Frequency : :class:`pycsamt.zonge.survey.Station`,
        :class:`pycsamt.zonge.meas.Frequency`
        One-dimensional vectors keyed by station id.
    Amps, Emag, Ephz, Hmag, Hphz : measurement vectors.
    Resistivity, Phase : 1-D apparent-resistivity and phase containers.
    PcEmag, PcHmag, PcRho, SEphz, SHphz, SPhz : quality / variation
        metrics.
    Z : :class:`pycsamt.zonge.z.Z`
        Full impedance tensor stack.
    DataInfo : :class:`pycsamt.zonge.info.DataInfo`
        Convenience bundle mirroring all of the above.

    The *values* and *loc* maps exposed by each component follow the
    same API: ``obj.loc['S05']`` returns the per-frequency 1-D array
    for station *S05*.


    Methods
    ----------
    read / fit
        Parse a file or DataFrame and populate every component.
    write
        Emit a *kind-2* CSV AVG (see ``to_kind2`` and ``to_kind1`` for
        direct control over the flavour).
    to_kind1, to_kind2
        Low-level helpers for serialisation without touching *self*.
    asdict
        Shallow JSON-serialisable representation of the whole object
        (delegated to :class:`BaseAVG`).


    Notes
    -------
    Quick background on AVG files:

    Zonge Engineering introduced two slightly different flavours:

    * *kind-1* – legacy, fixed-width white-space tables.  No metadata
      apart from a simple hardware banner.  Each row holds a single
      `(station, freq, component)` triplet [2]_.

    * *kind-2* – modern CSV with a rich ``$key = value`` header.  All
      data for a given receiver station are grouped together.  This
      is the format emitted by the current **AMTAVG** utility and by
      *ASTATIC* after static correction.

    The loader auto-detects the flavour and returns a tidy DataFrame
    with canonical lower-case column names (``station, freq, emag,
    rho, phase, ...``).  All subsequent processing is flavour-agnostic.
    
    * The class is decorated with :func:`pycsamt.decorators.has_fit`.
      Therefore ``avg.fit(src)`` and ``avg.read(src)`` are synonyms
      unless a real ``fit`` method is later implemented.
    * Numeric columns are *never* kept as strings.  The loader replaces
      Zonge placeholders (``*`` or blanks) with ``NaN`` and converts
      the rest to ``float64``.
    * Phases are stored internally in **radians**; the *to_degree*
      flag propagated down the component tree controls the public
      presentation.

    Examples
    ----------
    Load, tweak a component, and write back::

        >>> from pycsamt.zonge.avg import AVG
        >>> avg = AVG('LCS01.avg', verbose=True).fit()
        >>> avg.PcEmag.values *= 0.75       # re-scale E-error
        >>> avg.write('LCS01_clean.avg')

    Build from an in-memory DataFrame::

        >>> df, meta = load_avg('raw.avg')
        >>> avg = AVG(df, meta=meta).read()
        >>> print(avg.DataInfo['S10']['rho'])


    References
    ----------
    .. [1] Zonge Engineering (2016). *AMTAVG 7.76 User Manual*.

    .. [2] Zonge Engineering (1991). *Electromagnetic Methods in Applied
           Geophysics*, Vol. 2, pp. 713-809.

    """

    def __init__(
        self,
        obj_or_path: str | Path | pd.DataFrame | None = None,
        *,
        meta: dict[str, str] | None = None,
        verbose: bool = False
    ) -> None:
  
        # 1) create *all* component shells – they stay empty until
        #    _populate_components() fills them
        self.Head   = Head()

        # header-property blocks
        self.Hardware            = Hardware()
        self.SurveyAnnotation    = SurveyAnnotation()
        self.SurveyConfiguration = SurveyConfiguration()
        self.Receiver            = Receiver()
        self.Transmitter         = Transmitter()
        self.SkipFlag            = SkipFlag()

        # per-station / per-frequency measurements
        self.Station     = Station()
        self.Frequency   = Frequency()
        self.Amps        = Amps()
        self.Emag        = CompMeas("ExHy")
        self.Ephz        = CompMeas("ExHy")
        self.Hmag        = CompMeas("ExHy")
        self.Hphz        = CompMeas("ExHy")

        # variations & QC
        self.PcEmag      = PcEmag()
        self.PcHmag      = PcHmag()
        self.PcRho       = PcRho()
        self.SEphz       = SEphz()
        self.SHphz       = SHphz()
        self.SPhz        = SPhz()

        # resistivity / phase
        self.Resistivity = Resistivity()
        self.Phase       = Phase()

        # impedance tensor – created empty; maybe filled later
        self.Z           = Z()

        # misc
        self.DataInfo    = DataInfo()

        # 2) base initialisation (empty frame for now)
        super().__init__(pd.DataFrame(), meta or {}, None)
        self.verbose = verbose
        self.obj_or_path = obj_or_path 
        self.meta = meta 

        
    def read(
        self,
        obj_or_path: str | Path | pd.DataFrame |None,
        *,
        meta: dict[str, str] | None = None,
        inplace: bool = False
    ) -> "AVG":
        """
        Load an AVG file **or** a ready-made :class:`pandas.DataFrame`
        *and hydrate every component object*.

        • When *obj_or_path* is a path the file is parsed by
          :func:`pycsamt.zonge.utils.load_avg`.
        • When it is a DataFrame the caller is assumed to have done the
          parsing already and *meta* becomes the header dictionary.
        """
        # ---- 1. acquire DataFrame + header ------------------------
        self.obj_or_path = obj_or_path or self.obj_or_path 
        self.meta = meta or self.meta 
        
        if isinstance(self.obj_or_path, (str, Path)):
            df, hdr   = load_avg(Path(self.obj_or_path), inplace=inplace)
            _source   = Path(self.obj_or_path)
            if meta:                                   # user overrides
                hdr.update(self.meta)
        elif isinstance(self.obj_or_path, pd.DataFrame):
            df, hdr, _source = self.obj_or_path.copy(), (meta or {}), None
        else:
            raise TypeError("read() expects AVG-path or DataFrame")

        # 2. keep pristine copy for round-trips 
        self._df_core = df.copy()

        # 3. hydrate header / property blocks 
        self._populate_header(hdr)

        # ---- 4. create numeric component views 
        self._populate_components(df)

        # 5. update BaseAVG bookkeeping 
        super().__init__(data=df, meta=hdr, source=_source)

        if self.verbose:
            logger.info("AVG loaded – %s", self)
        return self

    def _populate_header(self, hdr: dict[str, str]) -> None:
        """
        Map the *hdr* dict (key=value pairs returned by ``load_avg``)
        onto the lightweight property containers.
        """
        
        # -- hardware 
        self.Hardware.set(
            version      = hdr.get("AMTAVG",    self.Hardware.version),
            source_file  = self.source,
            dated        = hdr.get("Dated"),
            processed    = hdr.get("Processed"),
            astatic_ver  = hdr.get("ASTATIC",   self.Hardware.astatic_ver),
            updated      = hdr.get("Updated"),
        )

        # -- survey configuration 
        self.SurveyConfiguration.set(
            survey_type   = hdr.get("Survey.Type",   "CSAMT"),
            array_type    = hdr.get("Survey.Array"),
            line_name     = hdr.get("Line.Name"),
            line_number   = _try_float(hdr.get("Line.Number")),
            line_azim_deg = _try_float(hdr.get("Line.Azimuth")),
            unit_length   = hdr.get("Unit.Length",   "m"),
            unit_emag     = hdr.get("Unit.E",        "nV/m"),
            unit_hfield   = hdr.get("Unit.B",        "pT"),
            unit_phase    = hdr.get("Unit.Phase",    "mrad"),
            utm_zone      = _try_int(hdr.get("UTM_Zone")),
        )

        # -- survey annotation 
        self.SurveyAnnotation.set(
            project_name  = hdr.get("Job.Name",  
                                    self.SurveyAnnotation.project_name),
            project_area  = hdr.get("Area"),
            customer_name = hdr.get("Customer"),
            contractor_name = hdr.get("Contractor"),
            project_label = hdr.get("Label"),
        )

        # -- transmitter & receiver 
        self.Transmitter.set(
            station   = _try_int(hdr.get("Tx.Stn",       0)),
            gdp_station = _try_int(hdr.get("Tx.GdpStn")),
            tx_type   = hdr.get("Tx.Type"),
        )
        self.Receiver.set(
            station   = _try_int(hdr.get("Rx.Stn",       0)),
            gdp_station = _try_int(hdr.get("Rx.GdpStn")),
            length_m  = _try_float(hdr.get("Rx.Length")),
            comps     = hdr.get("Rx.Cmp", "ExHy"),
        )

        # -- skip-flag (fall back to “2 = good”) 
        self.SkipFlag.set(hdr.get("Skip", "2"))

    def _populate_components(self, df: pd.DataFrame) -> None:
        """
        Slice the canonical DataFrame columns and feed every component
        object.  The logic is intentionally strict: **all** mandatory
        columns must exist; if something is missing a clear exception is
        raised so the caller knows the AVG is incomplete or wrongly
        parsed.
        """
        try:
            stn   = df["station"].values
            freq  = df["freq"].values
        except KeyError as exc:
            raise AvgDataError(f"mandatory column missing: {exc}") from None

        n_freq     = len(np.unique(freq))
        n_station  = len(np.unique(stn))

        # ── station & frequency first 
        self.Station.read(  stn,  unit="m")
        self.Frequency.read(freq, n_freq=n_freq)

        # ── core scalar / vector components
        _req = ("amps", "emag", "ephz", "hmag", "hphz", "rho", "phase")
        for col in _req:
            if col not in df.columns:
                raise AvgDataError(f"AVG file lacks <{col}> column")
                
        self.Amps.read( df["amps"].values,    n_freq, n_station)
        self.Emag.read( df["emag"].values,    n_freq, n_station)
        self.Ephz.read( df["ephz"].values,    n_freq, n_station)
        self.Hmag.read( df["hmag"].values,    n_freq, n_station)
        self.Hphz.read( df["hphz"].values,    n_freq, n_station)
        self.Resistivity.read(df["rho"].values,   n_freq, n_station)
        self.Phase.read( df["phase"].values, n_freq, n_station)

        # ── optional variations / QC columns 
        self.PcEmag.read(df.get("e.%err", np.zeros_like(freq)),
                         n_freq, n_station)
        self.PcHmag.read(df.get("h.%err", np.zeros_like(freq)),
                         n_freq, n_station)
        self.PcRho .read(df.get("rho.%err", np.zeros_like(freq)),
                         n_freq, n_station)
        self.SEphz.read(df.get("e.perr", np.zeros_like(freq)),
                        n_freq, n_station)
        self.SHphz.read(df.get("h.perr", np.zeros_like(freq)),
                        n_freq, n_station)
        self.SPhz .read(df.get("phase.%err", np.zeros_like(freq)),
                        n_freq, n_station)

        # ── impedance tensor (kind-2 only) ─────────────────────────────
        if {"z.mag", "z.phz"}.issubset(df.columns):
            z_abs   = df["z.mag"].values
            z_phase = df["z.phz"].values * 1e-3         # mrad → rad   (✓)
        
            # complex |Z|·e^{jφ}
            z_complex = z_abs * (np.cos(z_phase) + 1j * np.sin(z_phase))
        
            # assume ExHy only → xx = yy = 0,  xy = Z,  yx = 0
            z_stack = np.zeros((n_freq * n_station, 2, 2), dtype=complex)
            z_stack[:, 0, 1] = z_complex
        
            # (n_freq, n_station, 2, 2)
            z_stack = (
                z_stack.reshape(n_station, n_freq, 2, 2)
                       .swapaxes(0, 1)
            )
        
            self.Z.read(
                z_stack,
                freq        = np.unique(freq),
                station_ids = self.Station.names,
            )

    def _native_writer(
        self,
        path: str | Path | None,
        core_only: bool = False,
        stamp: bool = True,
        keep: Iterable[str] | None = None, 
        **kw
        ) -> Path:
        """
        Serialise the current frame back to disk (*kind-2 CSV*).

        Parameters
        ----------
        path       : str | Path
            Destination filename (``.avg`` will **not** be added
            automatically).
        core_only  : bool, default *False*
            When *True* only the *core* columns (see
            :func:`~pycsamt.zonge.utils.extract_core_columns`) are
            written; ancillary fields go to the *extra* block.
        stamp      : bool, default *True*
            Append a `$Written = …` UTC time-stamp.
        keep       : Iterable[str] | None
            Custom column subset passed to *extract_core_columns*
            (ignored if *core_only* is *False*).
        """
        if core_only:
            core, extra = extract_core_columns(self.data, keep=keep)
        else:
            core, extra = self.data, pd.DataFrame()

        path = write_avg(
            core, extra, 
            self.meta, path=path, 
            stamp=stamp, 
            **kw
            )
        if self.verbose:
        
            logger.info("AVG saved -> %s", Path(path).resolve())
            
        return path 
    
    def write(
        self,
        path: str | Path | None = None,
        *,
        fmt: Literal["kind1", "kind2", "auto", "native"] = "auto",
        when_unspecified: Literal["file", "string"] = "file",
        core_only: bool = False,
        stamp: bool = True,
        keep: Iterable[str] | None = None,
        **export_kw,
    ) -> Path | str:
        """
        Serialise the current AVG **components** back to disk or return
        the generated text.
    
        Parameters
        ----------
        path : str | Path | None, default *None*
            Destination.  If *None* the behaviour depends on
            *when_unspecified*.
        fmt : {"kind1", "kind2", "auto"}, default ``"auto"``
            Which flavour to emit.
    
            ``"auto"`` → decide from *path* suffix  
            ``"kind1"`` → legacy whitespace table  
            ``"kind2"`` → modern CSV with metadata
        when_unspecified : {"file", "string"}, default ``"file"``
            Fallback when *path* is *None*:
    
            * ``"file"``  → create a file called
              ``<basename>_<fmt>.avg`` in *cwd*.
            * ``"string"`` → return the text only.
        core_only, stamp, keep
            Forwarded to the underlying exporter **when applicable**.
        export_kw
            Extra keyword arguments propagated to *to_kind1* / *to_kind2*
            (e.g. *float_fmt*, *na_rep* …).
    
        Returns
        -------
        Path | str
            The written file **or** the AVG text.
        """

        # 0) Sanity
        if self._df_core.empty:
            raise AvgFileError(
                "Nothing to write – load or build data first"
            )
        
        fmt = str(fmt).lower() 
        
        if fmt=='native': 
            return self._native_write(
                path= path, 
                core_only=core_only, 
                stamp=stamp, 
                keep=keep, 
                **export_kw
             )
        
        # 1) Resolve format
        if fmt == "auto":
            if isinstance(path, (str, Path)):
                fmt = "kind1" if str(path).lower().endswith(
                    "_kind1.avg") else "kind2"
            else:
                fmt = "kind2"                       # sensible default
    
        if fmt not in {"kind1", "kind2"}:
            raise ValueError(
                "fmt must be 'kind1', 'kind2' or 'auto'")
    
        # 2) Delegate to the right exporter
        if fmt == "kind1":
            result = self.to_kind1(
                savefile=path,
                when_unspecified=when_unspecified,
                stamp=stamp,
                core_only=core_only,
                keep=keep,
                **export_kw,
            )
        else:  # kind-2
            result = self.to_kind2(
                savefile=path,
                when_unspecified=when_unspecified,
                stamp=stamp,
                core_keep=keep if core_only else None,
                **export_kw,
            )
    
        return result

    @property
    def core(self) -> pd.DataFrame:
        """Return *kind-2* core columns (always a **copy**)."""
        core, _ = extract_core_columns(self.data)
        return core

    @property
    def extra(self) -> pd.DataFrame:
        """Return ancillary columns (copy)."""
        _, extra = extract_core_columns(self.data)
        return extra


    @classmethod
    def from_avg(cls, path: str | Path, **kw) -> "AVG":
        """Shorthand for ``AVG().read(path, **kw)``."""
        return cls().read(path, **kw)


    def __getitem__(self, key: str) -> pd.Series:
        """Column selection: ``avg['rho']``."""
        return self.data[key]

    def __len__(self) -> int:        # row count
        return self.nrows

    def _build_frame_from_components(self) -> pd.DataFrame:
        """
        Recompose a tidy DataFrame from the current component objects.
    
        The routine assumes that *all “vector” components* share the same
        grid size ``(n_freq × n_station)``.  Station & frequency provide
        the master axes.
        """
        if not (self.Station and self.Frequency):
            raise AvgDataError("Station and/or Frequency not populated")
    
        stn = self.Station.value.ravel()
        frq = np.repeat(self.Frequency.value, len(
            self.Station.value) // len(self.Frequency.value))
    
        # helper: flatten component → 1-D array or np.nan
        def _vec(obj, attr="value"):
            return getattr(obj, attr).ravel(
                ) if obj is not None else np.full_like(frq, np.nan, dtype=float)
    
        frame_dict = {
            "Station"    : stn,
            "Freq"       : frq,
            "Comp"       : self.CompMeas.name if self.CompMeas else "ExHy",
            "Amps"       : _vec(self.Amps),
            "Emag"       : _vec(self.Emag),
            "Ephz"       : _vec(self.Ephz),
            "Hmag"       : _vec(self.Hmag),
            "Hphz"       : _vec(self.Hphz),
            "Resistivity": _vec(self.Resistivity),
            "Phase"      : _vec(self.Phase),
            "%Emag"      : _vec(self.PcEmag),
            "sEphz"      : _vec(self.SEphz),
            "%Hmag"      : _vec(self.PcHmag),
            "sHphz"      : _vec(self.SHphz),
            "%Rho"       : _vec(self.PcRho),
            "sPhz"       : _vec(self.SPhz),
        }
    
        df = pd.DataFrame(frame_dict)
        # leading skip flag (single value broadcasted)
        df.insert(0, "skp", self.Head.skip.code)
        return df


    def to_kind2(
        self,
        savefile: str | Path | None = None,
        *,
        when_unspecified: Literal["file", "string"] = "file",
        core_keep: Iterable[str] | None = None,
        stamp: bool = True,
        na_rep: str = "*",
        float_fmt: str = "%.6g",
    ) -> str | Path:
        """
        Emit a **kind-2** AVG (metadata header + CSV data block).
    
        Parameters
        ----------
        savefile : str | Path | None, default *None*
            Where to write.  See *when_unspecified* for the fallback rule.
        when_unspecified : {"file", "string"}, default ``"file"``
            Behaviour when *savefile* is *None*.
    
            ``"file"``   → create ``<basename>_kind2.avg`` in *cwd*.  
            ``"string"`` → return the text only.
        core_keep : Iterable[str] | None, optional
            Column names to keep in the **core** section
            (defaults to the helper’s internal minimal set).
        stamp : bool, default *True*
            Append a ``$Written=<UTC-timestamp>`` record.
        na_rep, float_fmt : str
            Place-holder for NaNs and numeric formatting (passed to
            :pandas:`DataFrame.to_csv`).
    
        Returns
        -------
        Path | str
            Written file path **or** the raw text.
        """

        # 1) Data     – built from live components, *not* the raw DF
        core, extra = extract_core_columns(
            self._build_frame_from_components(), keep=core_keep
        )
        data_block  = pd.concat([core, extra], axis=1)
    
        csv_txt = data_block.to_csv(
            index=False,
            na_rep=na_rep,
            float_format=float_fmt,
        )
    
        # 2) Metadata – combine Head + self.meta
        header_lines: list[str] = []
    
        # 2-a  pycsamt “Head” container   → $key = value
        if hasattr(self, "Head") and isinstance(self.Head, Head):
            for ln in self.Head.write():
                header_lines.append(f"${ln}")
    
        # 2-b  original AVG metadata      → $key = value
        for k, v in (self.meta or {}).items():
            header_lines.append(f"${k} = {v}")
    
        # 2-c  optional time-stamp
        if stamp:
            
            ts = datetime.now(tz=timezone.utc).isoformat(timespec="seconds")
            header_lines.append(f"$Written = {ts}")
    
        header_lines.append("")                      # blank line before CSV
        full_text = "\n".join(header_lines) + csv_txt
    

        # 3) Decide output
        if savefile is not None:
            out = Path(savefile).with_suffix(".avg")
            out.write_text(full_text)
            logger.info("kind-2 AVG written → %s", out)
            return out
    
        if when_unspecified == "string":
            logger.info("kind-2 AVG emitted as string (no file written)")
            return full_text
    
        # default: write to current directory
        base = (self.source.stem if getattr(self, "source", None)
                else "exported")
        out  = Path.cwd() / f"{base}_kind2.avg"
        out.write_text(full_text)
        logger.info("kind-2 AVG written (implicit) → %s", out)
        return out


    def to_kind1(
        self,
        savefile: str | Path | None = None,
        *,
        when_unspecified: Literal["file", "string"] = "file",
        precision: int = 3,
        na_rep: str = "*",
    ) -> str | Path:
        """
        Emit a *legacy* kind-1 AVG table built from **live components**.
    
        Parameters
        ----------
        savefile : str | Path | None, default *None*
            • Path → write exactly there.  
            • *None* → behaviour controlled by *when_unspecified*.
        when_unspecified : {"file", "string"}, default ``"file"``
            Action when *savefile* is *None*:
    
            ``"file"``   → create ``<basename>_kind1.avg`` in *cwd*.  
            ``"string"`` → return the text (do **not** write to disk).
        precision : int, default *3*
            Significant digits for numeric fields.
        na_rep : str, default ``"*"``
            Placeholder for missing values.
    
        Returns
        -------
        Path | str
            • :class:`pathlib.Path` when a file is produced.  
            • ``str`` when the text is returned.
        """

        # 1) rebuild DataFrame from the *current* component state
        df = self._build_frame_from_components()  # internal helper
    
        k1_order = [
            "skp", "Station", "Freq", "Comp", "Amps",
            "Emag", "Ephz", "Hmag", "Hphz",
            "Resistivity", "Phase",
            "%Emag", "sEphz", "%Hmag", "sHphz",
            "%Rho", "sPhz",
        ]
        df = df[k1_order]
    
        num_fmt = f"{{:.{precision}g}}".format
        body = df.applymap(
            lambda v: na_rep if pd.isna(v)
            else num_fmt(v) if isinstance(v, (int, float, np.floating))
            else str(v)
        )
    
        header = " ".join(k1_order)
        text   = header + "\n" + "\n".join(
            " ".join(r) for r in body.to_numpy()) + "\n"
    

        # 2) Decide what to do with *text*
        # a) user gave an explicit path  → honour it
        if savefile is not None:
            path = Path(savefile)
            if path.suffix.lower() != ".avg":
                path = path.with_suffix(".avg")
            path.write_text(text)
            logger.info("kind-1 AVG written -> %s", path)
            return path
    
        # b) *savefile* omitted
        if when_unspecified == "string":
            logger.info("kind-1 AVG emitted as string (no file written)")
            return text
    
        # default branch → produce a file in cwd
        base   = (self.source.stem if getattr(self, "source", None)
                  else "exported")
        path   = Path.cwd() / f"{base}_kind1.avg"
        path.write_text(text)
        logger.info("kind-1 AVG written (implicit) → %s", path)
        return path


@has_fit("warn")                                      # a real `.fit()` is present
class AVGOps(OpsBase):
    """
    Light **post-processing façade** around an :class:`~pycsamt.zonge.avg.AVG`
    instance.

    The class extracts the *Z-tensor*, apparent resistivity (ρ) and phase (ϕ)
    stacks into fast NumPy arrays and exposes a handful of convenience
    operations:

    * :meth:`get_reference_frequency` – pick the *best* (or closest) frequency
      according to a user criterion.
    * :meth:`apply_filter` – run a generic 1-D smoothing / custom callable on
      ρ and/or ϕ along the *frequency* axis.
    * :meth:`update` – push modified **z / rho / phase** back into the cached
      arrays **and** into the underlying :class:`AVG` (if any).

    Parameters
    ----------
    obj : str | pathlib.Path | pandas.DataFrame | AVG | None
        Any payload accepted by :pyclass:`AVG`.  When *None* the object must
        be provided later to :meth:`fit`.
    to_degree : bool, default *False*
        Convert phase to **degrees** upon caching.
    verbose : bool, default *False*
        INFO-level progress logging.
    **avg_kwargs
        Extra arguments forwarded to :class:`AVG` when the loader is invoked
        internally (e.g. ``verbose=False``).

    Notes
    -----
    The class is **stateless** with respect to the original source – you may
    freely modify the cached ``_z``, ``_rho``, ``_phase`` NumPy arrays and
    then call :meth:`update` to propagate changes back to :class:`AVG`.
    """

    def __init__(
        self,
        obj: str | Path | pd.DataFrame | AVG | None = None,
        *,
        to_degree: bool = False,
        verbose:   bool = False,
        **avg_kwargs: Any,
    ) -> None:

        super().__init__(obj, verbose=verbose)
        self._to_degree = to_degree

        # mutable working copies (shape: n_freq × n_station)
        self._z:     np.ndarray | None = None        # complex (…,2,2)
        self._rho:   np.ndarray | None = None
        self._phase: np.ndarray | None = None

        # AVG loader may need extra parameters later
        self._avg_kwargs = avg_kwargs

        # user handed something → try an immediate fit
        if obj is not None:
            self.fit(obj, **avg_kwargs)


    def fit(self, obj: Any | None = None, **avg_kwargs) -> "AVGOps":  # type: ignore[override]
        """
        Load *obj* (unless an :class:`AVG` is already held) and cache the
        essential arrays.  **Idempotent** – repeated calls are ignored once
        fitted.
        """
        if hasattr(self, "_is_fitted") and self._is_fitted:        # noqa: B023
            return self                                            # already done


        obj = obj or self.obj
        if isinstance(obj, AVG):
            self.avg: AVG = obj
            if self.avg.data.empty:
                raise AvgDataError("AVG instance has no data – call `.fit()` "
                                   "on it beforehand.")
        elif isinstance(obj, (pd.DataFrame, str, Path)):
            avg_kwargs |= self._avg_kwargs
            self.avg = AVG(obj, verbose=self.verbose, **avg_kwargs)
            self.avg.fit(obj)                                      # guarantee data
        else:
            raise TypeError("Unsupported payload for AVGOps.fit")


        # 2) cache NumPy views – makes later operations blazingly fast
        self._cache_raw()
        self._is_fitted = True
        if self.verbose:
            logger.info("AVGOps fitted: %d freq × %d station",
                        self.n_freq, self.n_station)
        return self


    @property
    def n_freq(self) -> int:
        """Number of discrete transmit frequencies."""
        return self._z.shape[0]                                

    @property
    def n_station(self) -> int:
        """Number of stations in the line / grid."""
        return self._z.shape[1]                                   

    def get_reference_frequency(
        self,
        target: float | None = None,
        *,
        strategy: Literal["nearest", "median", "max_rho"] = "nearest",
    ) -> int:
        """
        Return the *row-index* (0-based) of a “reference” frequency.

        Parameters
        ----------
        target : float | None, default *None*
            Desired frequency in **Hz** when *strategy* = ``'nearest'``.
        strategy : {"nearest", "median", "max_rho"}, default "nearest"
            Heuristic to pick the row:

            * ``'nearest'`` – closest to *target*.
            * ``'median'``  – 50-th percentile of the survey band.
            * ``'max_rho'`` – frequency at which median(ρ) is maximal.

        Raises
        ------
        NotFittedError
            When :meth:`fit` has not been called yet.
        """
        check_is_fitted(self, "_is_fitted")

        freq = self.avg.Frequency.value                      # 1-D unique array
        if strategy == "nearest":
            if target is None:
                raise ValueError("`target` must be supplied for 'nearest'")
            idx = int(np.abs(freq - target).argmin())
            return idx

        if strategy == "median":
            return int(len(freq) // 2)

        if strategy == "max_rho":
            med_rho = np.nanmedian(self._rho, axis=1)        # shape n_freq
            return int(np.nanargmax(med_rho))

        raise ValueError(f"Unknown strategy: {strategy}")

    def apply_filter(
        self,
        method: str | Callable[[np.ndarray], np.ndarray] = "moving",
        *,
        size:   int = 3
    ) -> None:
        """
        Run a *very* light 1-D filter along the **frequency axis**.

        Parameters
        ----------
        method : {"moving"} | callable
            * ``"moving"`` – centred moving average with *size* samples.
            * *callable*  – must accept / return a NumPy array of same
              shape *(n_freq, n_station)*.
        size : int, default *3*
            Window length for the moving average.

        Notes
        -----
        The operation **modifies** the cached ``_rho`` and ``_phase`` in-place.
        Call :meth:`update` afterwards if you want the underlying
        :class:`AVG` object to reflect the changes.
        """
        check_is_fitted(self, "_is_fitted")

        if callable(method):
            self._rho   = method(self._rho)
            self._phase = method(self._phase)
            return

        if method == "moving":
            if size < 2:
                raise ValueError("`size` must be ≥ 2 for moving-average")
            self._rho   = uniform_filter1d(self._rho,   size, axis=0,
                                           mode="nearest")
            self._phase = uniform_filter1d(self._phase, size, axis=0,
                                           mode="nearest")
            return

        raise ValueError(f"Unsupported filter: {method}")

    def update(
        self,
        *,
        z:     np.ndarray | None = None,
        rho:   np.ndarray | None = None,
        phase: np.ndarray | None = None
    ) -> None:
        """
        Replace **one or more** cached arrays and propagate the changes
        back to the underlying :pyattr:`avg` object.

        Notes
        -----
        Shapes **must** match the originals::

            (n_freq, n_station, 2, 2)  for *z*
            (n_freq, n_station)        for *rho* / *phase*
        """
        check_is_fitted(self, "_is_fitted")

        if z is not None:
            if z.shape != self._z.shape:
                raise ZError("Shape mismatch for `z`")
            self._z = z
            self.avg.Z.read(
                z, freq=self.avg.Frequency.value,
                station_ids=list(self.avg.Station.names)
            )

        if rho is not None:
            if rho.shape != self._rho.shape:
                raise ResistivityError("Shape mismatch for `rho`")
            self._rho = rho
            self.avg.Resistivity.read(
                rho.ravel(), n_freq=self.n_freq, n_stations=self.n_station)

        if phase is not None:
            if phase.shape != self._phase.shape:
                raise PhaseError("Shape mismatch for `phase`")
            self._phase = phase
            ph = np.deg2rad(phase) if self._to_degree else phase
            self.avg.Phase.read(
                ph.ravel(),
                n_freq=self.n_freq,
                n_stations=self.n_station,
                to_degree=self._to_degree
        )

    def _cache_raw(self) -> None:
        """
        Extract *Z*, ρ and ϕ from :pyattr:`avg` and store them in
        contiguous NumPy arrays – ready for fast numerical work.
        """
        z_ten = self.avg.Z.as_tensor()               # 3-D ImpedanceTensor
        self._z     = z_ten.z.copy()                 # (n_freq,n_stn,2,2)
        self._rho   = self.avg.Resistivity.values.reshape(  # type: ignore[attr-defined]
                          self.n_freq, self.n_station)
        raw_phase   = self.avg.Phase.values.reshape( # type: ignore[attr-defined]
                          self.n_freq, self.n_station)
        self._phase = (np.rad2deg(raw_phase) if self._to_degree
                       else raw_phase)


    def __str__(self) -> str:                        # noqa: D401
        if not hasattr(self, "_is_fitted"):
            return "AVGOps(<unfitted>)"
        return (f"AVGOps({self.n_freq} freq × {self.n_station} station, "
                f"degree={self._to_degree})")

    __repr__ = __str__


def _validate_shape(arr: np.ndarray, expected: Sequence[int], label: str) -> None:
    if arr.shape != tuple(expected):
        raise ValueError(f"{label} has shape {arr.shape}, expected {expected}")

# helper function 
def _try_float(val: str | int | float | None) -> float | None:
    try:
        return None if val in (None, "") else float(val)
    except (TypeError, ValueError):
        return None

def _try_int(val: str | int | None) -> int | None:
    try:
        return None if val in (None, "") else int(float(val))
    except (TypeError, ValueError):
        return None

