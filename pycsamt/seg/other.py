# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
from __future__ import annotations

from pathlib import Path 
from typing import List, Optional, Dict, Any, Tuple

from ..log.logger import get_logger
from ..exceptions import EdIDataError
from .base import EDIComponentBase
from .validation import (
    _strip_norm,
    _to_int_or_none,
    _to_float_or_none,
    _split_comment,
    _is_tag,
)
from .validation import IsEdi 

logger = get_logger(__name__)


__all__ = ["OtherSECT", "OtherIO", "OtherMixin"]


class OtherSECT(EDIComponentBase):
    r"""
    Header container for ``>=OTHERSECT`` blocks.

    This lightweight container parses the header that opens an
    *Other Data Section* as defined by the SEG-EDI spec.  It
    collects canonical options, any extra key/value pairs,
    and the ordered list of measurement IDs that may follow
    the options list.

    Parameters
    ----------
    *args :
        Unused positional arguments. Present for MRO safety.
    verbose : int or bool, optional
        Verbosity level. Propagated by :class:`EDIComponentBase`.
    logger : object, optional
        Logger instance. If ``None``, a null logger is used.
    **kws :
        Field overrides to pre-populate attributes.

    Attributes
    ----------
    sectid : str or None
        Section identifier. Often mirrors ``DATAID``.
    nchan : int or None
        Number of channels in this section.
    nfreq : int or None
        Number of frequencies, if provided.
    maxblks : int or None
        Upper bound on the number of data blocks.
    ndipole : int or None
        EMAP-style dipole count, if present.
    type : str or None
        Free-form type tag found in the wild (e.g. ``FREE``).
    extra : dict
        Any unrecognized header keys are stored here.
    meas_ids : list of str
        Measurement IDs listed under the header.
    start_data_lines_num : int or None
        Absolute line index where the first ``>BLOCK`` begins.

    Methods
    -------
    from_file(path)
        Parse the first ``>=OTHERSECT`` in an EDI file.  The
        file is validated by :class:`~pycsamt.seg.validation.IsEdi`.
    write()
        Serialize the header back to EDI text lines.

    Notes
    -----
    *Unknown* header keys are preserved in :attr:`extra` so
    round-tripping does not lose information.  Measurement IDs
    are appended in the order they appear in the file.

    Examples
    --------
    >>> hdr = OtherSECT.from_file("site.edi")
    >>> hdr.sectid
    'B1'
    >>> hdr.meas_ids[:2]
    ['HX', 'HY']
    >>> lines = hdr.write()
    >>> print("".join(lines).splitlines()[0])
    >=OTHERSECT

    See Also
    --------
    OtherIO : Read and write generic ``>BLOCK`` data.
    OtherMixin : Convenience helpers for host classes.

    References
    ----------
    .. [1] SEG EDI Standard (MT/EMAP), 1987.  MTNet archive.
    """

    KEY_ORDER: List[str] = [
        "sectid",
        "nchan",
        "nfreq",
        "maxblks",
        "ndipole",
        "type",
    ]

    def __init__(
        self,
        *args: Any,
        verbose: int | bool = 0,
        logger=None,
        **kws: Any,
    ):
        super().__init__(verbose=verbose, logger=logger)
        self.sectid: Optional[str] = None
        self.nchan: Optional[int] = None
        self.nfreq: Optional[int] = None
        self.maxblks: Optional[int] = None
        self.ndipole: Optional[int] = None
        self.type: Optional[str] = None
        self.extra: Dict[str, Any] = {}
        self.meas_ids: List[str] = []
        self.start_data_lines_num: Optional[int] = None

        for k, v in kws.items():
            setattr(self, k, v)

    @classmethod
    def from_file(cls, edi_path: str) -> "OtherSECT":

        p = Path(edi_path)
        IsEdi._assert_edi(p, deep=True)
    
        lines = p.read_text(
            encoding="utf-8-sig", errors="replace"
        ).splitlines()

        start = None
        for i, ln in enumerate(lines):
            if ln.upper().lstrip().startswith(">=OTHERSECT"):
                start = i
                break
        if start is None:
            raise EdIDataError("No >=OTHERSECT found.")

        # stop at first next tag (> or >=)
        stop = len(lines)
        for j in range(start + 1, len(lines)):
            s = lines[j].lstrip()
            if s.startswith(">"):
                stop = j
                break

        inst = cls()
        for raw in lines[start + 1 : stop]:
            s = raw.strip()
            if not s:
                continue
            if s.startswith("//"):
                continue
            if "=" in s:
                k, v = s.split("=", 1)
                key = _strip_norm(k).lower()
                val = _strip_norm(v)
                if key == "sectid":
                    inst.sectid = val
                elif key == "nchan":
                    inst.nchan = _to_int_or_none(val)
                elif key == "nfreq":
                    inst.nfreq = _to_int_or_none(val)
                elif key == "maxblks":
                    inst.maxblks = _to_int_or_none(val)
                elif key == "ndipole":
                    inst.ndipole = _to_int_or_none(val)
                elif key == "type":
                    inst.type = val
                else:
                    inst.extra[key] = val
            else:
                before, _ = _split_comment(s)
                toks = [t for t in before.split() if t]
                if toks:
                    inst.meas_ids.extend(toks)

        inst.start_data_lines_num = stop
        return inst

    def write(self) -> List[str]:
        out: List[str] = [">=OTHERSECT\n"]

        vals: Dict[str, Any] = {
            "sectid": self.sectid,
            "nchan": self.nchan,
            "nfreq": self.nfreq,
            "maxblks": self.maxblks,
            "ndipole": self.ndipole,
            "type": self.type,
        }

        for key in self.KEY_ORDER:
            val = vals.get(key, None)
            if val in (None, "", "None"):
                continue
            out.append(
                f"  {key.upper()}={str(val).upper()}\n"
            )

        for k, v in sorted(self.extra.items()):
            if v in (None, "", "None"):
                continue
            out.append(f"  {k.upper()}={str(v).upper()}\n")

        if self.meas_ids:
            out.append(f"    // {len(self.meas_ids)}\n")
            for mid in self.meas_ids:
                out.append(f"     {mid}\n")

        return out


class _OtherBlock(EDIComponentBase):
    """Generic data block: >KEY <opts> <data_set>."""

    def __init__(
        self,
        *args: Any,
        verbose: int | bool = 0,
        logger=None,
        **kws: Any,
    ):
        super().__init__(verbose=verbose, logger=logger)
        self.keyword: Optional[str] = None
        self.options: Dict[str, Any] = {}
        self.nitems_hint: Optional[int] = None
        self.values: List[float] = []
        self.raw_lines: List[str] = []
        for k, v in kws.items():
            setattr(self, k, v)

    def _format_values(self) -> List[str]:
        kpl = int(self.PER_LINE)
        ffmt = str(self.FLOAT_FMT)
        out: List[str] = []
        buf: List[str] = []
        for v in self.values:
            buf.append(ffmt.format(v))
            if len(buf) == kpl:
                out.append("  " + " ".join(buf) + "\n")
                buf = []
        if buf:
            out.append("  " + " ".join(buf) + "\n")
        return out


class OtherIO(EDIComponentBase):
    r"""
    Reader/writer for generic ``>BLOCK`` data under OTHERSECT.

    The class iterates over consecutive ``>KEY`` data blocks
    that follow an ``>=OTHERSECT`` header, and stores each as
    a small record.  Numeric rows are collected into
    :attr:`_OtherBlock.values`.  Non-numeric rows are kept in
    :attr:`_OtherBlock.raw_lines` to preserve content that
    cannot be parsed as floats.

    Parameters
    ----------
    *args :
        Unused positional arguments. Present for MRO safety.
    verbose : int or bool, optional
        Verbosity level. Propagated by :class:`EDIComponentBase`.
    logger : object, optional
        Logger instance. If ``None``, a null logger is used.
    **kws :
        Field overrides to pre-populate attributes.

    Attributes
    ----------
    blocks : list of _OtherBlock
        Parsed data blocks in file order.

    Methods
    -------
    from_file(path, start_line=None, verbose=0, logger=None)
        Parse all ``>BLOCK`` entries.  If ``start_line`` is
        ``None``, the reader seeks the first ``>OTHER`` or the
        line after ``>=OTHERSECT``.  Raises
        :class:`EdIDataError` when no blocks are found.
    write()
        Serialize the parsed blocks to EDI text lines.

    Notes
    -----
    *Typing strategy*.  Option values are parsed with a
    best-effort rule: integers first, then floats, else kept
    as strings.  Numeric table rows are formatted using
    :attr:`~pycsamt.seg.base.Base.FLOAT_FMT` and wrapped
    using :attr:`~pycsamt.seg.base.Base.PER_LINE`.

    Examples
    --------
    >>> hdr = OtherSECT.from_file("site.edi")
    >>> io = OtherIO.from_file("site.edi",
    ...                        start_line=hdr.start_data_lines_num)
    >>> [b.keyword for b in io.blocks]
    ['>COH', '>ANNO']
    >>> out = io.write()
    >>> print(out[0].strip().split()[0])
    >COH

    See Also
    --------
    OtherSECT : Header that precedes the data blocks.
    OtherMixin : Helper methods for host classes.

    References
    ----------
    .. [1] SEG EDI Standard (MT/EMAP), 1987.  MTNet archive.
    """

    def __init__(
        self,
        *args: Any,
        verbose: int | bool = 0,
        logger=None,
        **kws: Any,
    ):
        super().__init__(verbose=verbose, logger=logger)
        self.blocks: List[_OtherBlock] = []
        for k, v in kws.items():
            setattr(self, k, v)


    @classmethod
    def from_file(
        cls,
        edi_path: str,
        start_line: Optional[int] = None,
        *,
        verbose: int | bool = 0,
        logger=None,
    ) -> "OtherIO":
        p = Path(edi_path)
        IsEdi._assert_edi(p, deep=True)
    
        lines = p.read_text(
            encoding="utf-8-sig", errors="replace"
        ).splitlines()

        if start_line is None:
            for i, ln in enumerate(lines):
                if ln.upper().lstrip().startswith(">OTHER"):
                    start_line = i
                    break
            if start_line is None:
                for i, ln in enumerate(lines):
                    if ln.upper().lstrip().startswith(
                        ">=OTHERSECT"
                    ):
                        start_line = i + 1
                        break
        if start_line is None:
            raise EdIDataError(
                "No OTHER data blocks found."
            )

        inst = cls(verbose=verbose, logger=logger)
        i = start_line
        n = len(lines)
        while i < n:
            u = lines[i].upper().lstrip()
            if u.startswith(">="):
                break
            if u.startswith(">END"):
                break
            if not u.startswith(">"):
                i += 1
                continue

            blk, nxt = cls._parse_block(
                lines, i, verbose=verbose, logger=logger
            )
            inst.blocks.append(blk)
            i = nxt

        # NEW: raise when no blocks were collected
        if not inst.blocks:
            raise EdIDataError(
                "No OTHER data blocks found."
            )

        return inst

    @staticmethod
    def _parse_block(
        lines: List[str],
        i: int,
        *,
        verbose: int | bool = 0,
        logger=None,
    ) -> Tuple[_OtherBlock, int]:
        head = lines[i].rstrip("\n")
        before, cmt = _split_comment(head)
        toks = before.split()
        if not toks:
            raise EdIDataError("Malformed OTHER block.")
        # e.g. >COH
        keyword = toks[0]
        opts = toks[1:]

        blk = _OtherBlock(verbose=verbose, logger=logger)
        blk.keyword = keyword

        if cmt is not None:
            try:
                blk.nitems_hint = int(float(cmt))
            except Exception:
                blk.nitems_hint = None

        for t in opts:
            if "=" not in t:
                continue
            k, v = t.split("=", 1)
            key = _strip_norm(k).lower()
            val = _strip_norm(v)
            # best-effort typing
            ival = _to_int_or_none(val)
            if ival is not None:
                blk.options[key] = ival
                continue
            fval = _to_float_or_none(val)
            blk.options[key] = fval if fval is not None else val

        j = i + 1
        while j < len(lines):
            s = lines[j]
            if _is_tag(s, ">") or _is_tag(s, ">="):
                break
            st = s.strip()
            if not st or st.startswith("//"):
                j += 1
                continue

            body, _ = _split_comment(st)
            toks = body.split()
            # try parse floats; if any token non-float, keep raw
            row_vals: List[float] = []
            numeric_row = True
            for tok in toks:
                try:
                    row_vals.append(float(tok))
                except Exception:
                    numeric_row = False
                    break
            if numeric_row and row_vals:
                blk.values.extend(row_vals)
            else:
                blk.raw_lines.append(s.rstrip("\n"))
            j += 1

        return blk, j

    def write(self) -> List[str]:
        out: List[str] = []
        for blk in self.blocks:
            if not blk.keyword:
                continue
            head = [blk.keyword]
            # stable order for common keys, then extras
            for k in ("id", "rot", "comp"):
                if k in blk.options:
                    head.append(f"{k.upper()}={blk.options[k]}")
            for k in sorted(blk.options.keys()):
                if k in {"id", "rot", "comp"}:
                    continue
                head.append(f"{k.upper()}={blk.options[k]}")

            n_hint = (
                blk.nitems_hint
                if blk.nitems_hint is not None
                else (
                    len(blk.values)
                    if blk.values
                    else len(blk.raw_lines)
                )
            )
            out.append(" ".join(head) + f" // {n_hint}\n")

            if blk.raw_lines and not blk.values:
                for ln in blk.raw_lines:
                    out.append(ln.rstrip("\n") + "\n")
            else:
                out.extend(blk._format_values())
        return out


class OtherMixin:
    r"""
    Convenience helpers for OTHER sections.

    This mixin supplies two small facade methods that delegate
    to :class:`OtherSECT` and :class:`OtherIO`.  It lets a host
    class add simple ``read_*`` helpers without duplicating
    plumbing.

    Methods
    -------
    read_other_header(edi_fn, verbose=0, logger=None)
        Return a parsed :class:`OtherSECT` from ``edi_fn``.
    read_other_blocks(edi_fn, verbose=0, logger=None)
        Return an :class:`OtherIO` built from the same file.
        The header is read first to locate the data blocks.

    Examples
    --------
    >>> class Host(OtherMixin): ...
    >>> hdr = Host.read_other_header("site.edi")
    >>> io = Host.read_other_blocks("site.edi")
    >>> len(io.blocks) >= 0
    True

    See Also
    --------
    OtherSECT, OtherIO

    References
    ----------
    .. [1] SEG EDI Standard (MT/EMAP), 1987.  MTNet archive.
    """
    @classmethod
    def read_other_header(
        cls,
        edi_fn: str,
        *,
        verbose: int | bool = 0,
        logger=None,
    ) -> OtherSECT:
        return OtherSECT.from_file(edi_fn)

    @classmethod
    def read_other_blocks(
        cls,
        edi_fn: str,
        *,
        verbose: int | bool = 0,
        logger=None,
    ) -> OtherIO:
        sect = OtherSECT.from_file(edi_fn)
        return OtherIO.from_file(
            edi_fn,
            start_line=sect.start_data_lines_num,
            verbose=verbose,
            logger=logger,
        )
