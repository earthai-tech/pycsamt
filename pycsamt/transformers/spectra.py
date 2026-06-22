# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.transformers.spectra
============================

Transform SEG Spectra-format EDI files to MT-impedance EDI files.

The core operation is:

    Spectra EDI  →  Z (impedance) + Tipper  →  Impedance EDI

This module exposes :class:`SpectraToEDI`, the canonical transformer
for this conversion.  It accepts flexible input (single file, list of
files, directory, :class:`~pycsamt.seg.edi.EDIFile`, or
:class:`~pycsamt.seg.collection.EDICollection`) and always returns an
:class:`~pycsamt.seg.collection.EDICollection` — even for single-site
inputs.

Quick start
-----------
::

    from pycsamt.transformers import SpectraToEDI

    # --- single file ---
    col = SpectraToEDI().transform("site_spectra.edi")
    col[0].write(savepath="output/")

    # --- whole directory ---
    col = SpectraToEDI().transform("spectra_dir/", output_dir="imp_edis/")

    # --- with error propagation and tipper ---
    col = SpectraToEDI(estimate_error=True, verbose=1).transform(
        "spectra_dir/",
        output_dir="imp_edis/",
    )

    # --- inspect failures ---
    result = SpectraToEDI().transform_batch("spectra_dir/")
    for r in result.failures:
        print(r.source, r.error)

Spectra → Z math
-----------------
Per frequency *f*, the impedance tensor is estimated as:

    Z(f) = S_EH(f) · S_HH⁻¹(f)

where ``S_EH`` is the electric–magnetic cross-spectral matrix and
``S_HH`` is the magnetic auto-spectral block.  The tipper is:

    T(f) = S_ZH(f) · S_HH⁻¹(f)

Both are computed by :meth:`~pycsamt.seg.spectra.Spectra.to_Z`.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, List, Optional, Sequence, Tuple, Union

import numpy as np

from ..seg.collection import EDICollection
from ..seg.edi import EDIFile
from ._base import TransformerMixin

__all__ = ["SpectraToEDI", "TransformResult"]

_log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class _FailRecord:
    """A single failed conversion."""
    source: str
    error: str


@dataclass
class TransformResult:
    """Outcome of a batch spectra→impedance conversion.

    Attributes
    ----------
    collection : EDICollection
        Successfully converted EDI files.
    failures : list of _FailRecord
        Per-file failures (source path + error message).
    n_ok : int
        Number of successful conversions.
    n_fail : int
        Number of failed conversions.

    Examples
    --------
    ::

        result = SpectraToEDI().transform_batch("spectra_dir/")
        print(result.n_ok, "ok, ", result.n_fail, "failed")
        for f in result.failures:
            print(f.source, "->", f.error)
    """
    collection: EDICollection
    failures: List[_FailRecord] = field(default_factory=list)

    @property
    def n_ok(self) -> int:
        return len(self.collection)

    @property
    def n_fail(self) -> int:
        return len(self.failures)

    def __repr__(self) -> str:
        return (
            f"TransformResult(ok={self.n_ok}, fail={self.n_fail})"
        )


# ---------------------------------------------------------------------------
# SpectraToEDI
# ---------------------------------------------------------------------------

class SpectraToEDI(TransformerMixin):
    r"""Transform SEG spectra-EDI files to MT-impedance EDI files.

    The transformer calls :meth:`~pycsamt.seg.spectra.Spectra.from_file`
    on each input file, then :meth:`~pycsamt.seg.spectra.Spectra.to_edi`
    to recover ``Z`` and, when a vertical-magnetic channel is present,
    the tipper.  All successfully converted files are collected into a
    single :class:`~pycsamt.seg.collection.EDICollection`.

    Parameters
    ----------
    e_labels : tuple of str, default ("EX", "EY")
        Electric-channel type labels used to identify E-field channels
        in the spectra block.
    h_labels : tuple of str, default ("HX", "HY")
        Horizontal magnetic-channel type labels.
    ridge : float, optional
        Tikhonov regularisation added to the magnetic block before
        inversion.  Useful for numerically ill-conditioned spectra.
        ``None`` (default) applies no regularisation.
    estimate_error : bool, default False
        Propagate 1-σ impedance uncertainties into ``>ZXX.VAR`` /
        ``>ZXY.VAR`` / … blocks using first-order complex-Wishart
        error propagation.
    dof : float or ndarray, optional
        Effective degrees of freedom per frequency, forwarded to
        :meth:`~pycsamt.seg.spectra.Spectra.to_Z`.  When ``None``,
        the transformer tries to infer DoF from spectra metadata
        (``segnum``, ``avgt``, ``bw``).
    use_remote : bool, default False
        When both local and remote electric channels are present,
        set ``True`` to select the remote reference.
    station_suffix : str, default ""
        Appended to the station name of every converted EDI.  Use
        ``"_IMP"`` to reproduce the ``HBH03 → HBH03_IMP`` convention.
    skip_errors : bool, default True
        If ``True``, per-file conversion failures are caught, logged,
        and included in :attr:`TransformResult.failures` without
        aborting the batch.  Set to ``False`` to raise on first error.
    verbose : int, default 0
        Verbosity level.  ``1`` prints per-file progress; ``2`` adds
        debug detail.

    Examples
    --------
    Single file::

        col = SpectraToEDI().transform("HBH03.edi")
        ed  = col[0]
        ed.write(savepath="output/")

    Directory, write to disk::

        col = SpectraToEDI(estimate_error=True, station_suffix="_IMP").transform(
            "spectra_dir/",
            output_dir="imp_edis/",
        )
        print(len(col), "station(s) converted")

    Batch with failure report::

        result = SpectraToEDI(skip_errors=True).transform_batch(
            "spectra_dir/",
            output_dir="imp_edis/",
        )
        print(f"{result.n_ok} ok,  {result.n_fail} failed")

    Force a station-name override for a single file::

        col = SpectraToEDI().transform("site.edi", station_name="CUSTOM_01")

    See Also
    --------
    pycsamt.seg.spectra.Spectra.to_Z
        Underlying Z / tipper estimation.
    pycsamt.seg.spectra.Spectra.to_edi
        Single-file spectra → EDIFile conversion.
    pycsamt.transformers.AVGtoEDI
        Analogous transformer for Zonge AVG files.
    pycsamt.transformers.JtoEDI
        Analogous transformer for Jones J files.
    """

    def __init__(
        self,
        *,
        e_labels: Tuple[str, str] = ("EX", "EY"),
        h_labels: Tuple[str, str] = ("HX", "HY"),
        ridge: Optional[float] = None,
        estimate_error: bool = False,
        dof: Optional[Any] = None,
        use_remote: bool = False,
        station_suffix: str = "",
        skip_errors: bool = True,
        verbose: int = 0,
    ) -> None:
        super().__init__()
        self.e_labels      = tuple(e_labels)
        self.h_labels      = tuple(h_labels)
        self.ridge         = ridge
        self.estimate_error= estimate_error
        self.dof           = dof
        self.use_remote    = use_remote
        self.station_suffix= station_suffix
        self.skip_errors   = skip_errors
        self.verbose       = verbose

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def transform(
        self,
        source: Any,
        *,
        output_dir: Optional[Union[str, Path]] = None,
        station_name: Optional[str] = None,
    ) -> EDICollection:
        """Transform one or multiple spectra EDIs to impedance EDIs.

        Parameters
        ----------
        source : path-like, EDIFile, EDICollection, Spectra, or list
            Accepts:

            * A single ``.edi`` file path (str / Path).
            * A directory — all ``*.edi`` files inside are processed.
            * A :class:`~pycsamt.seg.edi.EDIFile` that holds a spectra
              section (its ``.path`` is used for the source metadata).
            * A :class:`~pycsamt.seg.collection.EDICollection`.
            * A pre-loaded :class:`~pycsamt.seg.spectra.Spectra` object.
            * A list of any of the above.
        output_dir : path-like, optional
            When given, each converted EDI is written to this directory.
            The directory is created if it does not exist.
        station_name : str, optional
            Override the station name for **single-file** input only.
            Ignored when *source* resolves to multiple files.

        Returns
        -------
        EDICollection
            Successfully converted impedance EDI files.

        Raises
        ------
        RuntimeError
            When ``skip_errors=False`` and any file fails to convert.
        ValueError
            When *source* resolves to zero convertible files.
        """
        result = self.transform_batch(
            source,
            output_dir=output_dir,
            station_name=station_name,
        )
        if result.n_fail and not self.skip_errors:
            msgs = "\n  ".join(f"{r.source}: {r.error}" for r in result.failures)
            raise RuntimeError(
                f"{result.n_fail} file(s) failed to convert:\n  {msgs}"
            )
        return result.collection

    def transform_batch(
        self,
        source: Any,
        *,
        output_dir: Optional[Union[str, Path]] = None,
        station_name: Optional[str] = None,
    ) -> TransformResult:
        """Like :meth:`transform` but always returns a :class:`TransformResult`.

        Useful when you want to inspect failures without raising.

        Parameters
        ----------
        source, output_dir, station_name
            Same as :meth:`transform`.

        Returns
        -------
        TransformResult
        """
        paths = self._resolve_sources(source)
        if not paths:
            raise ValueError(
                f"No spectra EDI files found in {source!r}.  "
                "Provide a .edi file or a directory containing .edi files."
            )

        out_dir = Path(output_dir) if output_dir is not None else None
        if out_dir is not None:
            out_dir.mkdir(parents=True, exist_ok=True)

        edis: List[EDIFile] = []
        failures: List[_FailRecord] = []

        single = len(paths) == 1

        for path in paths:
            sname = station_name if single else None
            try:
                ed = self._transform_one(path, station_name=sname)
                if out_dir is not None:
                    self._write_edi(ed, out_dir, verbose=self.verbose)
                edis.append(ed)
                if self.verbose >= 1:
                    tip_str = "with tipper" if ed.has_tipper else "no tipper"
                    _log.info(
                        "Converted %-30s  Z=%d freq  %s",
                        path.name,
                        ed.Z.n_freq if ed.Z is not None else 0,
                        tip_str,
                    )
            except Exception as exc:  # noqa: BLE001
                msg = str(exc)
                failures.append(_FailRecord(source=str(path), error=msg))
                if not self.skip_errors:
                    raise RuntimeError(
                        f"Conversion failed for {path}: {msg}"
                    ) from exc
                _log.warning("Skipping %s — %s", path.name, msg)

        return TransformResult(
            collection=EDICollection(items=edis, verbose=0),
            failures=failures,
        )

    # ------------------------------------------------------------------
    # Single-file core
    # ------------------------------------------------------------------

    def _transform_one(
        self,
        path: Path,
        *,
        station_name: Optional[str] = None,
    ) -> EDIFile:
        """Convert a single spectra EDI file to an impedance EDIFile.

        Parameters
        ----------
        path : Path
            Resolved path to the spectra EDI file.
        station_name : str, optional
            Station-name override.  When ``None``, the name embedded in
            the EDI ``>HEAD`` section is used, with
            :attr:`station_suffix` appended.

        Returns
        -------
        EDIFile
        """
        from ..seg.spectra import Spectra  # noqa: PLC0415

        if self.verbose >= 2:
            _log.debug("Reading spectra: %s", path)

        sp = Spectra.from_file(str(path))

        if sp.n_freq == 0:
            raise ValueError(f"No frequency blocks found in {path.name!r}")
        if sp.n_chan < 4:
            raise ValueError(
                f"{path.name!r} has only {sp.n_chan} channels; "
                "need at least HX, HY, EX, EY"
            )

        # Resolve station name
        if station_name is None:
            raw_name = getattr(sp, "name", None) or path.stem
            station_name = raw_name + self.station_suffix

        if self.verbose >= 2:
            _log.debug(
                "  to_edi(station=%s, estimate_error=%s)",
                station_name, self.estimate_error,
            )

        ed = sp.to_edi(
            source_edi=str(path),
            station_name=station_name,
            e_labels=self.e_labels,
            h_labels=self.h_labels,
            ridge=self.ridge,
            estimate_error=self.estimate_error,
            dof=self.dof,
        )
        return ed

    # ------------------------------------------------------------------
    # Source resolution helpers
    # ------------------------------------------------------------------

    def _resolve_sources(self, source: Any) -> List[Path]:
        """Normalise *source* to a sorted list of resolved ``.edi`` paths."""
        from ..seg.spectra import Spectra  # noqa: PLC0415

        if isinstance(source, Spectra):
            # Already a parsed Spectra — we need its backing file
            p = getattr(source, "_path", None) or getattr(source, "name", None)
            if p is not None:
                return [Path(p).resolve()]
            raise TypeError(
                "Spectra object has no path attribute; "
                "pass the file path directly instead."
            )

        if isinstance(source, EDIFile):
            if source.path:
                return [Path(source.path).resolve()]
            raise TypeError("EDIFile has no path; pass the file path directly.")

        if isinstance(source, EDICollection):
            paths = []
            for ed in source:
                p = getattr(ed, "path", None) or getattr(ed, "path_str", None)
                if p:
                    paths.append(Path(p).resolve())
            return sorted(set(paths))

        if isinstance(source, (str, Path)):
            p = Path(source)
            if p.is_dir():
                return sorted(p.glob("*.edi"))
            if p.suffix.lower() == ".edi":
                return [p.resolve()]
            raise ValueError(
                f"{source!r} is not a .edi file or a directory."
            )

        if isinstance(source, (list, tuple)):
            paths: List[Path] = []
            for item in source:
                paths.extend(self._resolve_sources(item))
            return sorted(set(paths))

        raise TypeError(
            f"Cannot interpret source of type {type(source).__name__!r}.  "
            "Pass a .edi path, a directory, EDIFile, EDICollection, "
            "Spectra, or a list of such."
        )

    # ------------------------------------------------------------------
    # Write helper
    # ------------------------------------------------------------------

    @staticmethod
    def _write_edi(ed: EDIFile, out_dir: Path, verbose: int = 0) -> Path:
        """Write *ed* to *out_dir* and return the output path."""
        try:
            out = ed.write(savepath=str(out_dir))
            if verbose >= 1:
                _log.info("  → wrote %s", Path(out).name if out else out_dir)
            return Path(out) if out else out_dir
        except Exception as exc:  # noqa: BLE001
            _log.warning("write failed for %s: %s", getattr(ed, "station", "?"), exc)
            raise

    # ------------------------------------------------------------------
    # Convenience class-method constructors
    # ------------------------------------------------------------------

    @classmethod
    def with_errors(cls, **kw: Any) -> "SpectraToEDI":
        """Return a transformer with ``estimate_error=True``."""
        return cls(estimate_error=True, **kw)

    @classmethod
    def with_tipper_suffix(cls, **kw: Any) -> "SpectraToEDI":
        """Return a transformer that appends ``'_IMP'`` to station names."""
        return cls(station_suffix="_IMP", **kw)
