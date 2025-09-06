# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
from __future__ import annotations
from pathlib import Path
from typing import List, Optional, Dict, Any, Tuple
import os

from ..log.logger import get_logger
from ..exceptions import EdIDataError
from .base import EDIComponentBase
from .validation import ( 
    _strip_norm, 
    _to_int_or_none, 
    _to_float_or_none, 
    _split_comment, 
    IsEdi
 )

logger = get_logger(__name__)

__all__ = ["TSect", "TSIO", "TimeSeriesMixin"]


class TSect(EDIComponentBase):
    r"""
    Minimal container for the ``>=TSERIESSECT`` header block.
    It parses the section header and the ordered list of
    measurement IDs that follow the header.  The class keeps a
    pointer to where the first ``>TSERIES`` data block starts so
    downstream readers can jump straight to the data.
    
    Parameters
    ----------
    verbose : int or bool, optional
        Verbosity flag inherited from :class:`EDIComponentBase`.
    logger : object, optional
        Logger instance inherited from :class:`EDIComponentBase`.
    **kws
        Keyword overrides for any public attribute.  Unknown
        keys are ignored.
    
    Attributes
    ----------
    sectid : str or None
        Section identifier.  If absent in file it remains
        ``None``.
    nchan : int or None
        Number of channels declared in the header.
    nmeas : int or None
        Number of measurements declared in the header.
    npts : int or None
        Number of samples per trace if provided.
    maxblks : int or None
        Hint for the maximum number of data blocks.
    dt : float or None
        Sampling interval in seconds when present.
    meas_ids : list of str
        Ordered list of measurement IDs collected from the
        header tail.  One ID per line.
    extra : dict
        Any non standard key–value options preserved as strings.
    start_data_lines_num : int or None
        Absolute line index where the first ``>TSERIES`` block
        begins.  Useful for fast data scans.
    
    Methods
    -------
    from_file(edi_path)
        Parse a single ``>=TSERIESSECT`` from an EDI file.  The
        method validates the file structure with
        :meth:`validation.IsEdi._assert_edi` before parsing.
    write()
        Serialize the section back to EDI lines including the
        measurement ID list.
    
    Notes
    -----
    Parsing is tolerant.  Unknown keys are stored in ``extra``.
    Blank lines and comment lines beginning with ``//`` are
    ignored.  If multiple time-series sections exist, call
    :meth:`from_file` on the desired file view or use a higher
    level iterator to locate the right header first.
    
    Examples
    --------
    >>> sect = TSect.from_file("sound.edi")
    >>> sect.nchan, sect.dt
    (3, 0.01)
    >>> print("IDs:", sect.meas_ids[:2])
    IDs: ['HX', 'HY']
    
    See Also
    --------
    TSIO
        Reader and writer for ``>TSERIES`` data blocks.
    validation.IsEdi
        Lightweight EDI file validator used during reading.
    
    References
    ----------
    .. [1] SEG EDI MT/EMAP standard (1987).  MTNet.
           https://www.mtnet.info/docs/seg_mt_emap_1987.pdf
    """


    KEY_ORDER: List[str] = [
        "sectid",
        "nchan",
        "nmeas",
        "npts",
        "maxblks",
        "dt",
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
        self.nmeas: Optional[int] = None
        self.npts: Optional[int] = None
        self.maxblks: Optional[int] = None
        self.dt: Optional[float] = None
        self.meas_ids: List[str] = []
        self.extra: Dict[str, Any] = {}
        self.start_data_lines_num: Optional[int] = None

        for k, v in kws.items():
            setattr(self, k, v)

    @classmethod
    def from_file(cls, edi_path: str) -> "TSect":
        p = Path(edi_path)
        IsEdi._assert_edi(p, deep=True)
    
        lines = p.read_text(
            encoding="utf-8-sig", errors="replace"
        ).splitlines()
    
        start = None
        for i, ln in enumerate(lines):
            if ln.lstrip().upper().startswith(">=TSERIESSECT"):
                start = i
                break
        if start is None:
            raise EdIDataError("No >=TSERIESSECT found.")
    
        # stop at first >TSERIES, next >=..., or EOF
        stop = len(lines)
        for j in range(start + 1, len(lines)):
            u = lines[j].lstrip().upper()
            if u.startswith(">TSERIES") or u.startswith(">="):
                stop = j
                break
    
        inst = cls()
        for raw in lines[start + 1 : stop]:
            s = raw.strip()
            if not s or s.startswith("//"):
                continue
            if "=" in s:
                k, v = s.split("=", 1)
                key = _strip_norm(k).lower()
                val = _strip_norm(v)
                if key == "sectid":
                    inst.sectid = val
                elif key == "nchan":
                    inst.nchan = _to_int_or_none(val)
                elif key == "nmeas":
                    inst.nmeas = _to_int_or_none(val)
                elif key == "npts":
                    inst.npts = _to_int_or_none(val)
                elif key == "maxblks":
                    inst.maxblks = _to_int_or_none(val)
                elif key == "dt":
                    inst.dt = _to_float_or_none(val)
                else:
                    inst.extra[key] = val
            else:
                if s:
                    inst.meas_ids.append(_strip_norm(s))
    
        inst.start_data_lines_num = stop
        return inst

    def write(self) -> List[str]:
        out: List[str] = [">=TSERIESSECT\n"]
        vals: Dict[str, Any] = {
            "sectid": self.sectid,
            "nchan": self.nchan,
            "nmeas": self.nmeas,
            "npts": self.npts,
            "maxblks": self.maxblks,
            "dt": self.dt,
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
                out.append(f"     {str(mid)}\n")
        return out


class _TSBlock(EDIComponentBase):
    """
    Single >TSERIES block: flexible header + values.
    """

    def __init__(
        self,
        *args: Any,
        verbose: int | bool = 0,
        logger=None,
        **kws: Any,
    ):
        super().__init__(verbose=verbose, logger=logger)
        self.options: Dict[str, Any] = {}
        self.nvals_hint: Optional[int] = None
        self.values: List[float] = []

        # common aliases if present in header
        self.npts: Optional[int] = None
        self.dt: Optional[float] = None
        self.id: Optional[str] = None

        for k, v in kws.items():
            setattr(self, k, v)

    def apply_aliases(self) -> None:
        npts = self.options.get("npts", None)
        dt = self.options.get("dt", None)
        cid = self.options.get("id", None)
        self.npts = (
            _to_int_or_none(npts) if npts is not None else None
        )
        self.dt = (
            _to_float_or_none(dt) if dt is not None else None
        )
        self.id = str(cid) if cid not in (None, "") else None


class TSIO(EDIComponentBase):
    r"""
    Reader and writer for ``>TSERIES`` data blocks.  Each data
    block line starts with a flexible option list (e.g.
    ``ID=HX NPTS=4 DT=0.25``) followed by a ``// N`` hint and
    then one or more lines of numeric samples.
    
    Parameters
    ----------
    verbose : int or bool, optional
        Verbosity flag inherited from :class:`EDIComponentBase`.
    logger : object, optional
        Logger instance inherited from :class:`EDIComponentBase`.
    **kws
        Keyword overrides for public attributes.
    
    Attributes
    ----------
    blocks : list of _TSBlock
        Parsed time-series blocks in file order.  Every block
        exposes:
            
            - ``options`` : dict of parsed header options.
            - ``nvals_hint`` : int or None from the ``//`` count.
            - ``values`` : list[float] of samples.
            - ``id`` : str or None (alias of ``options['id']``).
            - ``npts`` : int or None (alias of ``options['npts']``).
            - ``dt`` : float or None (alias of ``options['dt']``).
    
    Methods
    -------
    from_file(edi_path, start_line=None, *, verbose=0, logger=None)
        Parse all ``>TSERIES`` blocks starting at ``start_line``.
        If ``start_line`` is ``None`` the first block is located
        automatically.  The method assumes the file already
        passed :meth:`validation.IsEdi._assert_edi` upstream.
    write(per_line=None, float_fmt=None)
        Serialize every block.  ``per_line`` controls how many
        samples are printed per line.  ``float_fmt`` controls the
        numeric format (e.g. ``"{: .6E}"``).
    
    Notes
    -----
    Header options are typed heuristically.  Integer-like tokens
    become integers.  Otherwise they are parsed as floats when
    possible, and finally left as strings.  The common aliases
    ``id``, ``npts`` and ``dt`` are mirrored onto block fields
    for convenience.
    
    Examples
    --------
    >>> sect = TSect.from_file("sound.edi")
    >>> io = TSIO.from_file("sound.edi",
    ...                     start_line=sect.start_data_lines_num)
    >>> len(io.blocks)
    2
    >>> io.blocks[0].id, io.blocks[0].dt
    ('HX', 0.25)
    >>> lines = io.write(per_line=5, float_fmt="{: .3E}")
    >>> print("".join(lines).splitlines()[0])
    >TSERIES ID=HX NPTS=4 DT=0.25 // 4
    
    See Also
    --------
    TSect
        Header reader for ``>=TSERIESSECT``.
    SpectraIO
        Similar reader for ``>SPECTRA`` blocks.
    
    References
    ----------
    .. [1] SEG EDI MT/EMAP standard (1987).  MTNet.
           https://www.mtnet.info/docs/seg_mt_emap_1987.pdf
    """

    def __init__(
        self,
        *args: Any,
        verbose: int | bool = 0,
        logger=None,
        **kws: Any,
    ):
        super().__init__(verbose=verbose, logger=logger)
        self.blocks: List[_TSBlock] = []
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
    ) -> "TSIO":
        if not os.path.isfile(edi_path):
            raise FileNotFoundError(
                f"{edi_path!r} is not a file."
            )
        with open(edi_path, "r", encoding="utf-8") as f:
            lines = f.readlines()

        if start_line is None:
            for i, ln in enumerate(lines):
                if ln.upper().lstrip().startswith(">TSERIES"):
                    start_line = i
                    break
        if start_line is None:
            raise EdIDataError("No >TSERIES blocks found.")

        inst = cls(verbose=verbose, logger=logger)
        i = start_line
        n = len(lines)
        while i < n:
            ln = lines[i].rstrip("\n")
            u = ln.upper().lstrip()
            if u.startswith(">="):
                break
            if not u.startswith(">TSERIES"):
                i += 1
                continue

            blk, nxt = cls._parse_block(
                lines, i, verbose=verbose, logger=logger
            )
            inst.blocks.append(blk)
            i = nxt

        return inst

    @staticmethod
    def _parse_block(
        lines: List[str],
        i: int,
        *,
        verbose: int | bool = 0,
        logger=None,
    ) -> Tuple[_TSBlock, int]:
        head = lines[i].rstrip("\n")
        body, cmt = _split_comment(head)
        toks = body.split()
        # toks[0] is ">TSERIES"
        opts = toks[1:]

        blk = _TSBlock(verbose=verbose, logger=logger)
        if cmt is not None:
            try:
                blk.nvals_hint = int(float(cmt))
            except Exception:
                blk.nvals_hint = None

        for t in opts:
            if "=" not in t:
                continue
            k, v = t.split("=", 1)
            key = _strip_norm(k).lower()
            val = _strip_norm(v)

            # best-effort typing: only ints for integer-like tokens
            vlow = val.lower()
            is_int_like = (
                vlow.isdigit()
                or (vlow.startswith(("+", "-")) and vlow[1:].isdigit())
            )
            
            if is_int_like:
                blk.options[key] = _to_int_or_none(val)
            else:
                fval = _to_float_or_none(val)
                blk.options[key] = fval if fval is not None else val


        blk.apply_aliases()

        j = i + 1
        while j < len(lines):
            s = lines[j].strip()
            if not s:
                j += 1
                continue
            if s.startswith(">"):
                break
            if s.startswith("//"):
                j += 1
                continue
            before, _ = _split_comment(s)
            for tok in before.split():
                try:
                    blk.values.append(float(tok))
                except Exception:
                    # tolerate non-numeric tokens
                    pass
            j += 1

        return blk, j

    def write(
        self,
        per_line: Optional[int] = None,
        float_fmt: Optional[str] = None,
    ) -> List[str]:
        kpl = self.PER_LINE if per_line is None else per_line
        ffmt = self.FLOAT_FMT if float_fmt is None else float_fmt

        out: List[str] = []
        for blk in self.blocks:
            head = [">TSERIES"]
            # keep deterministic order for common keys
            for k in ("id", "npts", "dt"):
                v = blk.options.get(k, None)
                if v is not None:
                    head.append(f"{k.upper()}={v}")
            # then any extras
            for k in sorted(blk.options.keys()):
                if k in {"id", "npts", "dt"}:
                    continue
                v = blk.options[k]
                head.append(f"{k.upper()}={v}")

            n_hint = (
                blk.nvals_hint
                if blk.nvals_hint is not None
                else len(blk.values)
            )
            out.append(" ".join(head) + f" // {n_hint}\n")

            vals: List[str] = []
            cnt = 0
            for v in blk.values:
                vals.append(ffmt.format(v))
                cnt += 1
                if cnt == kpl:
                    out.append("  " + " ".join(vals) + "\n")
                    vals = []
                    cnt = 0
            if vals:
                out.append("  " + " ".join(vals) + "\n")

        return out


class TimeSeriesMixin:
    r"""
    Convenience mixin that exposes two helpers so host classes
    can read time-series content without depending on concrete
    implementations.
    
    Methods
    -------
    read_tseries_header(edi_fn, *, verbose=0, logger=None)
        Return a :class:`TSect` parsed from ``edi_fn``.  The
        result holds the header fields and the position of the
        first data block.
    read_tseries_blocks(edi_fn, *, verbose=0, logger=None)
        Return a :class:`TSIO` built from the same file.  The
        method internally calls :class:`TSect` to find the first
        ``>TSERIES`` and then streams all blocks.
    
    Notes
    -----
    Use this mixin in higher level readers or project classes to
    offer a thin, stable API.  The methods only read data and do
    not modify files on disk.
    
    Examples
    --------
    >>> class Reader(TimeSeriesMixin):
    ...     pass
    >>> hdr = Reader.read_tseries_header("sound.edi")
    >>> ts = Reader.read_tseries_blocks("sound.edi")
    >>> hdr.nchan, len(ts.blocks)
    (2, 3)
    
    See Also
    --------
    TSect
        Header parser for time-series sections.
    TSIO
        Data block reader and writer.
    
    References
    ----------
    .. [1] SEG EDI MT/EMAP standard (1987).  MTNet.
           https://www.mtnet.info/docs/seg_mt_emap_1987.pdf
    """

    @classmethod
    def read_tseries_header(
        cls,
        edi_fn: str,
        *,
        verbose: int | bool = 0,
        logger=None,
    ) -> TSect:
        return TSect.from_file(edi_fn)

    @classmethod
    def read_tseries_blocks(
        cls,
        edi_fn: str,
        *,
        verbose: int | bool = 0,
        logger=None,
    ) -> TSIO:
        sect = TSect.from_file(edi_fn)
        return TSIO.from_file(
            edi_fn,
            start_line=sect.start_data_lines_num,
            verbose=verbose,
            logger=logger,
        )

