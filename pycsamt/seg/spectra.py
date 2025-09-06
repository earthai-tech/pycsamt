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
    )
from .validation import IsEdi 

logger = get_logger(__name__)

__all__ = ["SpectraSECT", "SpectraIO", "SpectraMixin"]


class SpectraSECT(EDIComponentBase):
    r"""
    Minimal container for the ``>=SPECTRASECT`` header.

    The class parses and serializes the spectra section
    header that precedes one or more ``>SPECTRA`` data
    blocks. It collects the option key/values and the
    ordered set of measurement IDs that the spectra apply
    to, as described by the SEG EDI convention [1]_.

    Parameters
    ----------
    verbose : int or bool, optional
        Verbosity flag propagated from :class:`Base`.
    logger : object, optional
        Logger instance to use. If ``None``, a default
        null-safe logger is attached.
    **kws :
        Additional field overrides. Keys may include
        any attribute listed below.

    Attributes
    ----------
    sectid : str or None
        Section identifier, often a site name. Some
        files omit this or use a numeric ID.
    nchan : int or None
        Number of channels in the spectra set.
    nfreq : int or None
        Number of frequencies expected in the section.
    maxblks : int or None
        Maximum number of blocks. Rarely used.
    meas_ids : list of str
        Ordered measurement ID list that follows the
        option lines in ``>=SPECTRASECT``.
    start_data_lines_num : int or None
        Line index in the EDI where the first
        ``>SPECTRA`` block begins. Set by
        :meth:`from_file`.

    Notes
    -----
    * Parsing is tolerant to case and extra whitespace.
    * Unknown header keys are ignored instead of raising.
    * The measurement ID list is collected from the
      header body once option lines end.
    * The start of the spectra data is detected by the
      first ``>SPECTRA`` tag, by the next ``>=...`` tag,
      or by end of file, whichever comes first.
    * For consistent processing, maintain the same
      frequency set across related data sections,
      as recommended in the EDI spec [1]_.

    See Also
    --------
    SpectraIO
        Reader/writer for the ``>SPECTRA`` data blocks.
    MTEMAP
        Header for ``>=MTSECT`` or ``>=EMAPSECT``. The
        spectra frequency set should match the MT set.
    TSect
        Header for ``>=TSERIESSECT`` (time series).

    Examples
    --------
    Read only the header and measurement IDs:

    >>> sect = SpectraSECT.from_file("site.edi")
    >>> sect.nfreq, sect.nchan
    (128, 5)
    >>> sect.meas_ids[:2]
    ['HX1', 'HY1']

    Serialize a header:

    >>> sect.nfreq = 3
    >>> sect.meas_ids = ["HX", "HY", "EX", "EY"]
    >>> lines = sect.write()
    >>> print("".join(lines).strip())  # doctest: +ELLIPSIS
    >=SPECTRASECT
      SECTID=...
      NCHAN=...
      NFREQ=3
      MAXBLKS=...
        // 4
         HX
         HY
         EX
         EY

    References
    ----------
    .. [1] SEG EDI standard, "Spectra Data Sections".
    """

    KEY_ORDER: List[str] = [
        "sectid",
        "nchan",
        "nfreq",
        "maxblks",
    ]

    def __init__(self, *args: Any, verbose:int=0, logger=None, **kws: Any):
        super().__init__(verbose=verbose, logger=logger)
        self.sectid: Optional[str] = None
        self.nchan: Optional[int] = None
        self.nfreq: Optional[int] = None
        self.maxblks: Optional[int] = None
        self.meas_ids: List[str] = []
        self.start_data_lines_num: Optional[int] = None

        for k, v in kws.items():
            setattr(self, k, v)

    @classmethod
    def from_file(cls, edi_path: str) -> "SpectraSECT":
        # Pseudo test with >Head missing
        
        #     if not os.path.isfile(edi_path):
        #         raise FileNotFoundError(
        #             f"{edi_path!r} is not a file."
        #         )
        #     with open(edi_path, "r", encoding="utf-8") as f:
        #         lines = f.readlines()
        
        p = Path(edi_path)
        IsEdi._assert_edi(p, deep=True)
    
        lines = p.read_text(
            encoding="utf-8-sig", errors="replace"
        ).splitlines()
    
        # find >=SPECTRASECT
        start = None
        for i, ln in enumerate(lines):
            if ln.lstrip().upper().startswith(">=SPECTRASECT"):
                start = i
                break
        if start is None:
            raise EdIDataError("No >=SPECTRASECT found.")
    
        # stop at first >SPECTRA, next >=..., or EOF
        stop = len(lines)
        for j in range(start + 1, len(lines)):
            u = lines[j].lstrip().upper()
            if u.startswith(">SPECTRA") or u.startswith(">="):
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
                elif key == "nfreq":
                    inst.nfreq = _to_int_or_none(val)
                elif key == "maxblks":
                    inst.maxblks = _to_int_or_none(val)
            else:
                if s:
                    inst.meas_ids.append(_strip_norm(s))
    
        inst.start_data_lines_num = stop
        return inst

    def write(self) -> List[str]:
        out: List[str] = [">=SPECTRASECT\n"]
        values: Dict[str, Any] = {
            "sectid": self.sectid,
            "nchan": self.nchan,
            "nfreq": self.nfreq,
            "maxblks": self.maxblks,
        }
        for key in self.KEY_ORDER:
            val = values.get(key, None)
            if val in (None, "", "None"):
                continue
            out.append(
                f"  {key.upper()}={str(val).upper()}\n"
            )
        if self.meas_ids:
            out.append(
                f"    // {len(self.meas_ids)}\n"
            )
            for mid in self.meas_ids:
                out.append(f"     {str(mid)}\n")
        return out


class _SpectraBlock(EDIComponentBase):
    """
    Single >SPECTRA block container.
    """

    def __init__(
        self, *args: Any, 
        verbose: int | bool =0 , 
        logger=None, **kws: Any
    ):
        super().__init__(verbose=verbose, logger=logger)
        self.freq: Optional[float] = None
        self.rotspec: Optional[int] = None
        self.bw: Optional[float] = None
        self.avgt: Optional[float] = None
        self.nvals_hint: Optional[int] = None
        self.options: Dict[str, Any] = {}
        self.values: List[float] = []

        for k, v in kws.items():
            setattr(self, k, v)

    def header_dict(self) -> Dict[str, Any]:
        d = dict(self.options)
        if self.freq is not None:
            d["freq"] = self.freq
        if self.rotspec is not None:
            d["rotspec"] = self.rotspec
        if self.bw is not None:
            d["bw"] = self.bw
        if self.avgt is not None:
            d["avgt"] = self.avgt
        return d


class SpectraIO(EDIComponentBase):
    r"""
    Read and write ``>SPECTRA`` data blocks.

    A spectra section contains one block per frequency.
    Each block begins with a ``>SPECTRA`` line that holds
    options such as frequency and bandwidth, optionally
    followed by a comment with the number of values,
    then one or more lines of numeric values.

    Known options are normalized:

    * ``FREQ`` : float
    * ``ROTSPEC`` : int
    * ``BW`` : float
    * ``AVGT`` : float

    Unrecognized options are preserved in a free-form
    mapping so that vendor-specific metadata is not lost.

    Parameters
    ----------
    verbose : int or bool, optional
        Verbosity flag propagated from :class:`Base`.
    logger : object, optional
        Logger instance to use. If ``None``, a default
        null-safe logger is attached.
    **kws :
        Additional field overrides.

    Attributes
    ----------
    blocks : list of _SpectraBlock
        Parsed spectra blocks, one per frequency. Each
        block stores header options, the optional value
        count hint, and the numeric values.

    Notes
    -----
    * :meth:`from_file` reads successive ``>SPECTRA``
      blocks starting from a given line or from the
      first match in the file.
    * Values are parsed as floats; non-numeric tokens
      in data lines are ignored rather than raising.
    * The writer orders known options first in the
      header line, then appends extra options sorted
      by key. Both option keys and values are written
      in upper case.
    * Line formatting uses the per-line and float
      format defaults from :class:`Base` unless you
      provide explicit overrides.

    See Also
    --------
    SpectraSECT
        Header container for spectra sections.
    TSIO
        Time-series counterpart for ``>TSERIES``.

    Examples
    --------
    Read all spectra blocks:

    >>> io = SpectraIO.from_file("site.edi")
    >>> len(io.blocks)
    128
    >>> b0 = io.blocks[0]
    >>> b0.freq, b0.bw  # doctest: +ELLIPSIS
    (..., ...)

    Build and serialize blocks:

    >>> from pycsamt.seg.spectra import _SpectraBlock
    >>> io = SpectraIO()
    >>> blk = _SpectraBlock()
    >>> blk.freq = 10.0
    >>> blk.rotspec = 1
    >>> blk.values = [0.1, 0.2, 0.3]
    >>> io.blocks.append(blk)
    >>> lines = io.write(per_line=2, float_fmt="{: .3E}")
    >>> print("".join(lines).strip())
    >SPECTRA FREQ=10.0 ROTSPEC=1 // 3
      1.000E-01  2.000E-01
      3.000E-01

    References
    ----------
    .. [1] SEG EDI standard, "Spectra Data Sections".
    """

    def __init__(
        self, *args: Any, 
        verbose: int | bool =0 , 
        logger=None, **kws: Any
    ):
        super().__init__(verbose=verbose, logger=logger)
        self.blocks: List[_SpectraBlock] = []
        for k, v in kws.items():
            setattr(self, k, v)

    # --------------------------
    # Load from file
    # --------------------------
    @classmethod
    def from_file(
        cls,
        edi_path: str,
        start_line: Optional[int] = None,
    ) -> "SpectraIO":
        if not os.path.isfile(edi_path):
            raise FileNotFoundError(
                f"{edi_path!r} is not a file."
            )
        with open(edi_path, "r", encoding="utf-8") as f:
            lines = f.readlines()

        if start_line is None:
            # find first >SPECTRA
            for i, ln in enumerate(lines):
                if ln.upper().lstrip().startswith(">SPECTRA"):
                    start_line = i
                    break
        if start_line is None:
            raise EdIDataError("No >SPECTRA blocks found.")

        inst = cls()
        i = start_line
        n = len(lines)
        while i < n:
            ln = lines[i].rstrip("\n")
            u = ln.upper().lstrip()
            if u.startswith(">="):
                break
            if not u.startswith(">SPECTRA"):
                i += 1
                continue

            blk, next_i = cls._parse_block(lines, i)
            inst.blocks.append(blk)
            i = next_i

        return inst

    # --------------------------
    # Parse one >SPECTRA block
    # --------------------------
    @staticmethod
    def _parse_block(
        lines: List[str],
        i: int,
    ) -> Tuple[_SpectraBlock, int]:
        head = lines[i].rstrip("\n")
        body, cmt = _split_comment(head)
        toks = body.split()
        # toks[0] is ">SPECTRA"
        opts = toks[1:]

        blk = _SpectraBlock()
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
            if key == "freq":
                blk.freq = _to_float_or_none(val)
            elif key == "rotspec":
                blk.rotspec = _to_int_or_none(val)
            elif key == "bw":
                blk.bw = _to_float_or_none(val)
            elif key == "avgt":
                blk.avgt = _to_float_or_none(val)
            else:
                # keep unknown options
                blk.options[key] = val

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
            # collect floats
            before, _ = _split_comment(s)
            for tok in before.split():
                try:
                    blk.values.append(float(tok))
                except Exception:
                    # tolerate bad tokens
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
            head = [">SPECTRA"]
            hd = blk.header_dict()
            # keep stable order in output
            for key in ("freq", "rotspec", "bw", "avgt"):
                if key in hd and hd[key] is not None:
                    head.append(
                        f"{key.upper()}={hd[key]}"
                    )
            # include any extra options
            for k, v in sorted(blk.options.items()):
                head.append(f"{k.upper()}={str(v).upper()}")

            n_hint = (
                blk.nvals_hint
                if blk.nvals_hint is not None
                else len(blk.values)
            )
            head_line = " ".join(head) + f" // {n_hint}\n"
            out.append(head_line)

            line_vals: List[str] = []
            cnt = 0
            for v in blk.values:
                line_vals.append(ffmt.format(v))
                cnt += 1
                if cnt == kpl:
                    out.append(
                        "  " + " ".join(line_vals) + "\n"
                    )
                    line_vals = []
                    cnt = 0
            if line_vals:
                out.append(
                    "  " + " ".join(line_vals) + "\n"
                )

        return out


class SpectraMixin:
    r"""
    Convenience facade for spectra access.

    This mixin exposes a compact API that host classes
    can reuse to discover and read spectra sections in
    an EDI file.

    Methods
    -------
    from_file(edi_fn)
        Return a :class:`SpectraSECT` parsed from the
        first ``>=SPECTRASECT`` header in ``edi_fn``.
    read_blocks(edi_fn)
        Return a :class:`SpectraIO` by scanning all
        subsequent ``>SPECTRA`` blocks that belong to the
        section discovered by :class:`SpectraSECT`.

    Notes
    -----
    Use this mixin in higher-level readers so spectra
    handling remains consistent and centralized. The
    method pair mirrors the design used for MT/EMAP
    headers and for time series sections.

    See Also
    --------
    SpectraSECT
        Header parsing and serialization.
    SpectraIO
        Data block reader/writer.
    MTEMAP
        MT/EMAP section header, often used alongside
        spectra for the same dataset.

    Examples
    --------
    >>> class Reader(SpectraMixin):
    ...     pass
    >>> sect = Reader.from_file("site.edi")
    >>> io = Reader.read_blocks("site.edi")
    >>> len(io.blocks) > 0
    True

    References
    ----------
    .. [1] SEG EDI standard, "Spectra Data Sections".
    """

    @classmethod
    def from_file(cls, edi_fn: str) -> SpectraSECT:
        return SpectraSECT.from_file(edi_fn)

    @classmethod
    def read_blocks(
        cls, edi_fn: str
    ) -> SpectraIO:
        sect = SpectraSECT.from_file(edi_fn)
        return SpectraIO.from_file(
            edi_fn, start_line=sect.start_data_lines_num
        )
