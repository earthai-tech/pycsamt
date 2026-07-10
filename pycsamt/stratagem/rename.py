# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
stratagem.rename
================

EDI file renaming and rewriting utilities for Stratagem post-processing.

:class:`EDIRenamer`
    Rename a set of EDI files on disk according to a configurable basename,
    zero-pad width, and trailer, mirroring the ``watex.utils.rename_files``
    workflow used in the legacy pipeline.  Simultaneously updates the
    ``>HEAD`` ``DATAID`` and linked section ``SECTID`` fields so the EDI
    content stays consistent with the new filename.

:class:`EDIWriter`
    Write a list of in-memory :class:`~pycsamt.seg.edi.EDIFile` objects to
    a directory, with optional renaming and ``>HEAD`` metadata overrides
    (station name / DATAID, acquisition info, etc.).

Both classes follow the ``fit() → out()`` pattern used by all stratagem
classes.
"""

from __future__ import annotations

from pathlib import Path

from ..api.property import MetadataMixin, PyCSAMTObject
from ..exceptions import FileHandlingError, NotFittedError
from ..seg.edi import EDIFile

__all__ = ["EDIRenamer", "EDIWriter"]


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _make_new_name(basename: str, idx: int, zero_pad: int, trailer: str) -> str:
    """Build a new EDI filename from parts.

    Parameters
    ----------
    basename : str
        Prefix string, e.g. ``'T2.'``.
    idx : int
        Station index (0-based).
    zero_pad : int
        Width of the zero-padded index field.
    trailer : str
        Suffix appended after the index but before ``.edi``.

    Returns
    -------
    str
        E.g. ``'T2.000.edi'``, ``'T2.086.edi'``.
    """
    return f"{basename}{idx:0{zero_pad}d}{trailer}.edi"


def _new_dataid(basename: str, idx: int, zero_pad: int, trailer: str) -> str:
    """Build the DATAID / SECTID string matching the new filename."""
    return f"{basename}{idx:0{zero_pad}d}{trailer}"


# ---------------------------------------------------------------------------
# EDIRenamer
# ---------------------------------------------------------------------------

class EDIRenamer(PyCSAMTObject):
    """Rename Stratagem EDI files with a standardised naming convention.

    Reads each source EDI, updates ``>HEAD.DATAID`` and the linked
    ``SECTID`` fields to match the new name, then writes the result to
    *dst_path* (keeping the source files untouched).

    Parameters
    ----------
    basename : str, default ``'S'``
        Name prefix.  E.g. ``'T2.'`` produces ``T2.000.edi``,
        ``T2.001.edi``, …
    zero_pad : int, default 3
        Width of the zero-padded integer part (``'T2.000'`` has
        ``zero_pad=3``).
    trailer : str, default ``''``
        Optional string appended after the index (before ``.edi``).
    update_dataid : bool, default True
        When True, ``>HEAD.DATAID`` and all linked ``SECTID`` fields are
        updated to match the new filename stem.
    overwrite : bool, default False
        Overwrite existing files in *dst_path*.
    verbose : int, default 0

    Attributes
    ----------
    renamed_pairs_ : list of (Path, Path)
        ``(src, dst)`` path pairs for every file that was processed.
    skipped_ : list of Path
        Source files skipped because the destination already existed and
        ``overwrite=False``.

    Examples
    --------
    Rename processed EDIs to ``T2.000.edi`` … ``T2.082.edi``:

    >>> rn = EDIRenamer(basename="T2.", zero_pad=3)
    >>> rn.fit("2/2EDIP", "2/renamedEDIs")

    Or rename in-memory objects produced by the processing pipeline:

    >>> rn.fit(nr.edi_objects_, "2/renamedEDIs")
    """

    __repr_fields__ = ("basename", "zero_pad", "trailer", "n_renamed_")

    def __init__(
        self,
        *,
        basename: str = "S",
        zero_pad: int = 3,
        trailer: str = "",
        update_dataid: bool = True,
        overwrite: bool = False,
        verbose: int = 0,
    ) -> None:
        self.basename = basename
        self.zero_pad = zero_pad
        self.trailer = trailer
        self.update_dataid = update_dataid
        self.overwrite = overwrite
        self.verbose = verbose

    # ------------------------------------------------------------------
    def fit(
        self,
        source: str | Path | list[EDIFile] | list[Path],
        dst_path: str | Path,
    ) -> EDIRenamer:
        """Rename EDI files and write them to *dst_path*.

        Parameters
        ----------
        source : path-like, list of EDIFile, or list of Path
            Input EDI files.  Accepts:

            * A directory path — all ``.edi`` files in it (natural-sort
              order).
            * A list of :class:`~pycsamt.seg.edi.EDIFile` objects
              (e.g. from ``NoiseRemover.edi_objects_``).
            * A list of :class:`pathlib.Path` EDI paths.

        dst_path : path-like
            Output directory.  Created if absent.

        Returns
        -------
        self
        """
        out_dir = Path(dst_path).expanduser().resolve()
        out_dir.mkdir(parents=True, exist_ok=True)

        # resolve inputs to a list of (EDIFile | None, src_path | None)
        entries: list[tuple[EDIFile | None, Path | None]] = (
            self._resolve_source(source)
        )

        self.renamed_pairs_: list[tuple[Path, Path]] = []
        self.skipped_: list[Path] = []

        for i, (edi_obj, src_path) in enumerate(entries):
            new_fname = _make_new_name(
                self.basename, i, self.zero_pad, self.trailer
            )
            dst_file = out_dir / new_fname

            if dst_file.exists() and not self.overwrite:
                self.skipped_.append(dst_file)
                if self.verbose:
                    print(f"[EDIRenamer] skip existing: {new_fname}")
                continue

            # load from path if we only have a path, not an object
            if edi_obj is None and src_path is not None:
                try:
                    edi_obj = EDIFile(src_path, verbose=0)
                except Exception as exc:
                    if self.verbose:
                        print(f"[EDIRenamer] load failed {src_path.name}: {exc}")
                    continue

            if edi_obj is None:
                continue

            # update DATAID / SECTID
            if self.update_dataid:
                new_id = _new_dataid(self.basename, i, self.zero_pad, self.trailer)
                try:
                    edi_obj.station = new_id
                except Exception:
                    pass

            # write
            try:
                edi_obj.write(new_edifn=new_fname, savepath=str(out_dir))
                self.renamed_pairs_.append((
                    src_path or Path(new_fname),
                    dst_file,
                ))
            except Exception as exc:
                if self.verbose:
                    print(f"[EDIRenamer] write failed {new_fname}: {exc}")

        self.n_renamed_ = len(self.renamed_pairs_)

        if self.verbose:
            print(
                f"[EDIRenamer] renamed {self.n_renamed_} files → {out_dir} "
                f"({len(self.skipped_)} skipped)"
            )

        return self

    # ------------------------------------------------------------------
    def _resolve_source(
        self,
        source: str | Path | list[EDIFile] | list[Path],
    ) -> list[tuple[EDIFile | None, Path | None]]:
        """Normalise *source* to a list of (EDIFile | None, Path | None)."""
        import re as _re  # already imported at module level, but kept explicit

        if isinstance(source, (str, Path)):
            d = Path(source).expanduser().resolve()
            if not d.is_dir():
                raise FileHandlingError(f"Source directory not found: {d}")
            paths = sorted(
                d.glob("*.edi"),
                key=lambda p: tuple(
                    int(t) if t.isdigit() else t
                    for t in _re.split(r"(\d+)", p.stem.lower())
                ),
            )
            return [(None, p) for p in paths]

        if not source:
            return []

        first = source[0]
        if isinstance(first, EDIFile):
            return [(e, getattr(e, "path", None)) for e in source]
        if isinstance(first, Path):
            return [(None, p) for p in source]

        # fallback: try to treat as EDIFile-like
        return [(e, getattr(e, "path", None)) for e in source]

    # ------------------------------------------------------------------
    def dst_paths(self) -> list[Path]:
        """Return the list of written destination paths."""
        if not hasattr(self, "renamed_pairs_"):
            raise NotFittedError("Call fit() first.")
        return [dst for _, dst in self.renamed_pairs_]


# ---------------------------------------------------------------------------
# EDIWriter
# ---------------------------------------------------------------------------

class EDIWriter(PyCSAMTObject, MetadataMixin):
    """Write in-memory EDIFile objects to disk with optional HEAD overrides.

    Provides a thin, consistent wrapper around
    :meth:`~pycsamt.seg.edi.EDIFile.write` that also allows batch update
    of ``>HEAD`` fields (DATAID, ACQBY, DATAID prefix, etc.) before
    writing.

    Parameters
    ----------
    dataid_prefix : str, optional
        When set, each station's DATAID is overwritten with
        ``f"{dataid_prefix}{i:0{zero_pad}d}"``.  Useful for
        standardising station identifiers across a profile.
    zero_pad : int, default 3
        Zero-pad width used with *dataid_prefix*.
    overwrite : bool, default False
    verbose : int, default 0

    Attributes
    ----------
    written_ : list of Path
        Paths of successfully written files.
    failed_ : list of tuple(str, Exception)
        ``(filename, exc)`` for any file that could not be written.

    Examples
    --------
    Write the noise-corrected EDIs, keeping original file names:

    >>> wr = EDIWriter()
    >>> wr.fit(nr.edi_objects_, "2/final")
    >>> wr.written_

    Write with standardised DATAID ``S000`` … ``S082``:

    >>> wr = EDIWriter(dataid_prefix="S", zero_pad=3)
    >>> wr.fit(nr.edi_objects_, "2/final")
    """

    __repr_fields__ = ("dataid_prefix", "zero_pad", "n_written_")

    def __init__(
        self,
        *,
        dataid_prefix: str | None = None,
        zero_pad: int = 3,
        overwrite: bool = False,
        verbose: int = 0,
    ) -> None:
        self.dataid_prefix = dataid_prefix
        self.zero_pad = zero_pad
        self.overwrite = overwrite
        self.verbose = verbose

    # ------------------------------------------------------------------
    def fit(
        self,
        edi_objects: list[EDIFile],
        savepath: str | Path,
        *,
        head_overrides: dict[str, object] | None = None,
    ) -> EDIWriter:
        """Write *edi_objects* to *savepath*.

        Parameters
        ----------
        edi_objects : list of EDIFile
        savepath : path-like
            Output directory.
        head_overrides : dict, optional
            Key-value pairs applied to every EDI's ``>HEAD`` object before
            writing.  Keys must be valid HEAD attribute names (e.g.
            ``'acqby'``, ``'stdvers'``).

        Returns
        -------
        self
        """
        out_dir = Path(savepath).expanduser().resolve()
        out_dir.mkdir(parents=True, exist_ok=True)

        self.written_: list[Path] = []
        self.failed_: list[tuple[str, Exception]] = []

        for i, edi in enumerate(edi_objects):
            # determine output filename
            fname = (
                edi.path.name
                if getattr(edi, "path", None) is not None
                else f"station_{i:0{self.zero_pad}d}.edi"
            )

            out_path = out_dir / fname
            if out_path.exists() and not self.overwrite:
                self.written_.append(out_path)
                continue

            # apply DATAID prefix
            if self.dataid_prefix is not None:
                new_id = f"{self.dataid_prefix}{i:0{self.zero_pad}d}"
                try:
                    edi.station = new_id
                    fname = f"{new_id}.edi"
                    out_path = out_dir / fname
                except Exception:
                    pass

            # apply HEAD overrides
            if head_overrides:
                head = edi.get_section("head")
                if head is not None:
                    for k, v in head_overrides.items():
                        try:
                            setattr(head, k, v)
                        except Exception:
                            pass

            # write
            try:
                edi.write(new_edifn=fname, savepath=str(out_dir))
                self.written_.append(out_path)
            except Exception as exc:
                self.failed_.append((fname, exc))
                if self.verbose:
                    print(f"[EDIWriter] write failed {fname}: {exc}")

        self.n_written_ = len(self.written_)

        if self.verbose:
            print(
                f"[EDIWriter] wrote {self.n_written_} files → {out_dir}"
                + (f" ({len(self.failed_)} failed)" if self.failed_ else "")
            )

        return self
