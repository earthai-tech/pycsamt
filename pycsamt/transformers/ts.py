# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.transformers.ts
=======================

Transform raw MT field time series to MT-impedance EDI files.

The core operation is:

    field time series  →  cross-spectra  →  Z (+ tipper)  →  EDI

built on :mod:`pycsamt.ts` (readers + spectral estimation) and
exposed here with the same transformer API as
:class:`~pycsamt.transformers.SpectraToEDI`,
:class:`~pycsamt.transformers.AVGtoEDI` and
:class:`~pycsamt.transformers.JtoEDI`: flexible input (single
file, directory, :class:`~pycsamt.ts.TSData`, or a list), an
:class:`~pycsamt.seg.collection.EDICollection` result, and a
:class:`~pycsamt.transformers.TransformResult` batch report.

Supported inputs are the :func:`pycsamt.ts.readers.read_ts`
formats: SAMTEX/LiMS ``.ts``, EMSLAB ``.asc``, EDI-embedded
``>TSERIES`` and generic ASCII tables.

Quick start
-----------
::

    from pycsamt.transformers import TStoEDI

    # --- single field file ---
    col = TStoEDI(nfft=10240).transform(
        "data/MT/TS/kap103as.ts/kap103as.ts",
        output_dir="imp_edis/",
    )
    print(col[0].station, col[0].Z.n_freq)

    # --- whole directory (recursive) ---
    col = TStoEDI().transform("data/MT/TS/", output_dir="imp_edis/")

    # --- in-memory record, remote reference ---
    from pycsamt.ts import read_ts

    local = read_ts("kap103as.ts")
    remote = read_ts("kap112as.ts")
    col = TStoEDI(estimator="rr", remote=remote).transform(local)

    # --- inspect failures ---
    result = TStoEDI().transform_batch("data/MT/TS/")
    for r in result.failures:
        print(r.source, r.error)

See Also
--------
pycsamt.ts.ts_to_edifile
    The single-record assembly core reused here.
pycsamt.transformers.SpectraToEDI
    Analogous transformer for spectra-EDI inputs.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from ..seg.collection import EDICollection
from ..seg.edi import EDIFile
from ..ts._base import TSData
from ._base import TransformerMixin
from .spectra import TransformResult, _FailRecord

__all__ = ["TStoEDI"]

_log = logging.getLogger(__name__)

#: file patterns collected when *source* is a directory
_TS_PATTERNS = ("*.ts", "*.asc", "*.dat", "*.txt")


class TStoEDI(TransformerMixin):
    r"""Transform raw MT field time series to impedance EDI files.

    Each source record is processed with
    :func:`pycsamt.ts.ts_to_edifile`: band-averaged cross-power
    spectra (Hann-tapered overlapping segments, Huber-weighted
    robust stacking), impedance recovery ``Z = S_EH · S_HH⁻¹``
    (or remote reference ``Z = S_ER · S_HR⁻¹``), and EDI header
    synthesis from the recorded site metadata (coordinates,
    channel azimuths, dipole lengths).  All successfully
    converted records are gathered into a single
    :class:`~pycsamt.seg.collection.EDICollection`.

    Parameters
    ----------
    estimator : {'ls', 'rr'}, default 'ls'
        Single-site least squares, or remote reference (needs
        ``remote``).
    remote : TSData, optional
        Simultaneous remote-reference record whose horizontal
        magnetic channels are appended as ``RHX``/``RHY``.
    nfft : int, optional
        FFT segment length.  Defaults to an automatic choice
        targeting ~60+ segments at 50 % overlap.
    overlap : float, default 0.5
        Fractional segment overlap.
    per_decade : int, default 7
        Estimation frequencies per decade.
    fmin, fmax : float, optional
        Target frequency range (Hz).
    robust : {'huber', None}, default 'huber'
        Segment stacking scheme.
    ridge : float, optional
        Tikhonov regularisation of the inverted spectral block.
    estimate_error : bool, default False
        Propagate 1-σ impedance uncertainties into
        ``>ZXX.VAR`` / … blocks (DoF from the band averaging).
    include_spectra : bool, default False
        Embed the estimated cross-spectra as ``>SPECTRA``
        blocks in each output EDI.
    include_tseries : bool, default False
        Embed a decimated head of the raw series as
        ``>TSERIES`` blocks.
    reader_kws : dict, optional
        Extra options for :func:`pycsamt.ts.readers.read_ts`
        (e.g. ``{"declination": -19.5}`` for EMSLAB, or
        ``{"format": "ascii", "dt": 1.0}``).
    station_suffix : str, default ""
        Appended to the station name of every converted EDI
        (e.g. ``"_TS"``).
    skip_errors : bool, default True
        Catch and record per-record failures instead of
        aborting the batch.
    verbose : int, default 0
        Verbosity level.
    **spectra_kws
        Any further option forwarded to
        :func:`pycsamt.ts.process.ts_to_spectra`
        (``freqs``, ``n_iter``, ``huber_k``,
        ``min_segments``, ``fill_gap``, ...).

    Examples
    --------
    Single LiMS file, write to disk::

        col = TStoEDI(nfft=10240).transform(
            "data/MT/TS/kap103as.ts/kap103as.ts",
            output_dir="imp_edis/",
        )

    EMSLAB file with reader options and errors::

        col = TStoEDI(
            estimate_error=True,
            reader_kws={"declination": -19.5},
        ).transform("data/MT/TS/emsl01.asc/emsl01.asc")

    Batch a directory and inspect failures::

        result = TStoEDI().transform_batch("data/MT/TS/")
        print(result.n_ok, "ok,", result.n_fail, "failed")

    See Also
    --------
    pycsamt.ts.ts_to_edi
        Functional one-record counterpart.
    pycsamt.transformers.SpectraToEDI
        Same API for spectra-EDI inputs.
    """

    def __init__(
        self,
        *,
        estimator: str = "ls",
        remote: TSData | None = None,
        nfft: int | None = None,
        overlap: float = 0.5,
        per_decade: int = 7,
        fmin: float | None = None,
        fmax: float | None = None,
        robust: str | None = "huber",
        ridge: float | None = None,
        estimate_error: bool = False,
        include_spectra: bool = False,
        include_tseries: bool = False,
        reader_kws: dict | None = None,
        station_suffix: str = "",
        skip_errors: bool = True,
        verbose: int = 0,
        **spectra_kws: Any,
    ) -> None:
        super().__init__()
        self.estimator = estimator
        self.remote = remote
        self.nfft = nfft
        self.overlap = overlap
        self.per_decade = per_decade
        self.fmin = fmin
        self.fmax = fmax
        self.robust = robust
        self.ridge = ridge
        self.estimate_error = estimate_error
        self.include_spectra = include_spectra
        self.include_tseries = include_tseries
        self.reader_kws = dict(reader_kws or {})
        self.station_suffix = station_suffix
        self.skip_errors = skip_errors
        self.verbose = verbose
        self.spectra_kws = dict(spectra_kws)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def transform(
        self,
        source: Any,
        *,
        output_dir: str | Path | None = None,
        station_name: str | None = None,
    ) -> EDICollection:
        """Transform one or multiple time-series records to EDIs.

        Parameters
        ----------
        source : path-like, TSData, or list
            Accepts:

            * A single field file path (LiMS ``.ts``, EMSLAB
              ``.asc``, generic ASCII, or an EDI carrying
              ``>TSERIES`` blocks).
            * A directory — files matching ``*.ts``, ``*.asc``,
              ``*.dat``, ``*.txt`` are collected recursively.
            * An in-memory :class:`~pycsamt.ts.TSData`.
            * A list mixing any of the above.
        output_dir : path-like, optional
            When given, each converted EDI is written there as
            ``<station>_ts.edi`` (directory created if needed).
        station_name : str, optional
            Station override for **single-record** input only.
            Ignored when *source* resolves to multiple records.

        Returns
        -------
        EDICollection
            Successfully converted impedance EDI files.

        Raises
        ------
        RuntimeError
            When ``skip_errors=False`` and any record fails.
        ValueError
            When *source* resolves to zero records.
        """
        result = self.transform_batch(
            source,
            output_dir=output_dir,
            station_name=station_name,
        )
        if result.n_fail and not self.skip_errors:
            msgs = "\n  ".join(f"{r.source}: {r.error}" for r in result.failures)
            raise RuntimeError(
                f"{result.n_fail} record(s) failed to convert:\n  {msgs}"
            )
        return result.collection

    def transform_batch(
        self,
        source: Any,
        *,
        output_dir: str | Path | None = None,
        station_name: str | None = None,
    ) -> TransformResult:
        """Like :meth:`transform` but always returns a
        :class:`~pycsamt.transformers.TransformResult`.

        Parameters
        ----------
        source, output_dir, station_name
            Same as :meth:`transform`.

        Returns
        -------
        TransformResult
        """
        records = self._resolve_sources(source)
        if not records:
            raise ValueError(
                f"No time-series records found in {source!r}. "
                "Provide a field file (.ts/.asc/.dat/.txt), a "
                "directory, or a TSData object."
            )

        out_dir = Path(output_dir) if output_dir is not None else None
        if out_dir is not None:
            out_dir.mkdir(parents=True, exist_ok=True)

        edis: list[EDIFile] = []
        failures: list[_FailRecord] = []
        single = len(records) == 1

        for rec in records:
            sname = station_name if single else None
            label = self._label(rec)
            try:
                ed = self._transform_one(rec, station_name=sname)
                if out_dir is not None:
                    self._write_edi(ed, out_dir, verbose=self.verbose)
                edis.append(ed)
                if self.verbose >= 1:
                    tip_str = "with tipper" if ed.has_tipper else "no tipper"
                    _log.info(
                        "Converted %-30s  Z=%d freq  %s",
                        label,
                        ed.Z.n_freq if ed.Z is not None else 0,
                        tip_str,
                    )
            except Exception as exc:  # noqa: BLE001
                msg = str(exc)
                failures.append(_FailRecord(source=label, error=msg))
                if not self.skip_errors:
                    raise RuntimeError(f"Conversion failed for {label}: {msg}") from exc
                _log.warning("Skipping %s — %s", label, msg)

        return TransformResult(
            collection=EDICollection(items=edis, verbose=0),
            failures=failures,
        )

    # ------------------------------------------------------------------
    # Single-record core
    # ------------------------------------------------------------------

    def _transform_one(
        self,
        record: Path | TSData,
        *,
        station_name: str | None = None,
    ) -> EDIFile:
        """Convert a single record to an impedance EDIFile."""
        from ..ts.convert import ts_to_edifile

        if self.verbose >= 2:
            _log.debug(
                "Processing time series: %s",
                self._label(record),
            )

        ed = ts_to_edifile(
            record,
            station=station_name,
            estimator=self.estimator,
            ridge=self.ridge,
            estimate_error=self.estimate_error,
            include_spectra=self.include_spectra,
            include_tseries=self.include_tseries,
            reader_kws=self.reader_kws,
            verbose=self.verbose,
            remote=self.remote,
            nfft=self.nfft,
            overlap=self.overlap,
            per_decade=self.per_decade,
            fmin=self.fmin,
            fmax=self.fmax,
            robust=self.robust,
            **self.spectra_kws,
        )
        # apply the suffix to the embedded station name
        # (setter also mirrors it to the MTSECT sectid)
        if station_name is None and self.station_suffix:
            base = ed.station or ""
            ed.station = f"{base}{self.station_suffix}"
        return ed

    # ------------------------------------------------------------------
    # Source resolution helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _label(record: Path | TSData) -> str:
        """Human-readable identifier for logs/failures."""
        if isinstance(record, TSData):
            return str(
                record.station or record.name or record.meta.get("source") or "TSData"
            )
        return str(record)

    def _resolve_sources(self, source: Any) -> list[Path | TSData]:
        """Normalise *source* to records (paths or TSData)."""
        if isinstance(source, TSData):
            return [source]

        if isinstance(source, (str, Path)):
            p = Path(source)
            if p.is_file():
                return [p.resolve()]
            if p.is_dir():
                found: list[Path] = []
                for pat in _TS_PATTERNS:
                    found.extend(f for f in p.rglob(pat) if f.is_file())
                return sorted(set(found))
            raise ValueError(f"{source!r} is neither a file nor a directory.")

        if isinstance(source, (list, tuple)):
            records: list[Path | TSData] = []
            for item in source:
                records.extend(self._resolve_sources(item))
            # keep TSData objects; de-duplicate paths only
            paths = sorted({r for r in records if isinstance(r, Path)})
            objs = [r for r in records if isinstance(r, TSData)]
            return objs + paths

        raise TypeError(
            "Cannot interpret source of type "
            f"{type(source).__name__!r}.  Pass a field "
            "file path, a directory, a TSData, or a list "
            "of such."
        )

    # ------------------------------------------------------------------
    # Write helper
    # ------------------------------------------------------------------

    @staticmethod
    def _write_edi(ed: EDIFile, out_dir: Path, verbose: int = 0) -> Path:
        """Write *ed* to *out_dir* and return the path."""
        try:
            name = f"{ed.station or 'site'}_ts.edi"
            out = ed.write(new_edifn=name, savepath=str(out_dir))
            if verbose >= 1:
                _log.info(
                    "  -> wrote %s",
                    Path(out).name if out else out_dir,
                )
            return Path(out) if out else out_dir
        except Exception as exc:  # noqa: BLE001
            _log.warning(
                "write failed for %s: %s",
                getattr(ed, "station", "?"),
                exc,
            )
            raise

    # ------------------------------------------------------------------
    # Convenience class-method constructors
    # ------------------------------------------------------------------

    @classmethod
    def with_errors(cls, **kw: Any) -> TStoEDI:
        """Return a transformer with ``estimate_error=True``."""
        return cls(estimate_error=True, **kw)

    @classmethod
    def with_remote(cls, remote: TSData, **kw: Any) -> TStoEDI:
        """Return a remote-reference transformer."""
        return cls(estimator="rr", remote=remote, **kw)
