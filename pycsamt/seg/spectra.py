# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Dict, Any, Tuple
import os

import numpy as np 

from ..log.logger import get_logger
from ..exceptions import EdIDataError 
from ..z.base import BaseEM 

from .base import EDIComponentBase
from .validation import ( 
    _strip_norm, 
    _to_int_or_none, 
    _to_float_or_none, 
    _split_comment,
    )
from .validation import IsEdi 


logger = get_logger(__name__)

__all__ = ["SpectraSECT", "SpectraIO", "SpectraMixin", "Spectra"]

_EMPTY = 1.0e32

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

    def __iter__(self):
        return iter(self.blocks)

    def __len__(self):
        return len(self.blocks)

    def __getitem__(self, idx):
        return self.blocks[idx]

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

class Spectra(BaseEM):
    r"""
    Container for ``>SPECTRA`` blocks grouped per frequency.

    The class gathers one spectra record per frequency and
    exposes typed header fields (frequency, rotation flag,
    bandwidth, and averaging time) together with the numeric
    values stored in each block.  It is a compact, array-
    oriented view on top of :class:`~pycsamt.seg.spectra.
    SpectraSECT` and :class:`~pycsamt.seg.spectra.SpectraIO`.

    Parameters
    ----------
    name : str, optional
        Display name forwarded to :class:`BaseEM`.
    verbose : int, default 0
        Verbosity level, forwarded to :class:`BaseEM`.

    Attributes
    ----------
    freq : ndarray, shape ``(n_blk,)``
        Frequency (Hz) per block.  Missing values are set to
        ``np.nan``.
    rotspec : ndarray of int, shape ``(n_blk,)``
        Rotation specifier per block.  Missing values are set
        to ``-1``.
    bw : ndarray, shape ``(n_blk,)``
        Nominal bandwidth (Hz) per block or ``np.nan``.
    avgt : ndarray, shape ``(n_blk,)``
        Averaging time (s) per block or ``np.nan``.
    values : list of ndarray
        Numeric payload for each block.  Lengths may differ
        across blocks, as allowed by the SEG format.
    n_values : ndarray of int, shape ``(n_blk,)``
        Number of values in each block (as parsed or counted).

    Notes
    -----
    Blocks may contain vendor-specific options beyond the
    canonical ``FREQ``, ``ROTSPEC``, ``BW``, and ``AVGT``.
    Those options are preserved when round-tripping via
    :meth:`to_io`.  The class does **not** impose a common
    length across spectra vectors; if you require a 2-D array,
    pad the :attr:`values` list explicitly.

    The constructor itself does not read files.  Use
    :meth:`from_io` or :meth:`from_file` to populate an
    instance from sections and data blocks.

    Methods
    -------
    from_io(sect, io) : classmethod
        Build a :class:`Spectra` from :class:`SpectraSECT`
        and :class:`SpectraIO`.
    from_file(path) : classmethod
        Convenience that calls :class:`SpectraSECT.from_file`
        and :class:`SpectraIO.from_file`, then delegates to
        :meth:`from_io`.
    to_io()
        Serialize the current state to a fresh pair
        (:class:`SpectraSECT`, :class:`SpectraIO`) that can be
        written back to an EDI file.

    Examples
    --------
    Read, inspect, and serialize spectra::

        from pycsamt.seg.spectra import Spectra

        sp = Spectra.from_file("site.edi")
        f = sp.freq
        first = sp.values[0]

        sect2, io2 = sp.to_io()
        # writer can now combine sect2.write() and io2.write()

    See Also
    --------
    pycsamt.seg.spectra.SpectraSECT
        Header for ``>=SPECTRASECT`` sections.
    pycsamt.seg.spectra.SpectraIO
        Reader/writer for ``>SPECTRA`` blocks.
    pycsamt.seg.EDIFile
        High-level dispatcher that can attach spectra to an
        EDI session.

    References
    ----------
    .. [1] SEG EDI standard, *Spectra Data Sections*.
           Society of Exploration Geophysicists.
    .. [2] Chave, A. D., & Jones, A. G. (2012). *The
           Magnetotelluric Method: Theory and Practice*.
           Cambridge Univ. Press.
    """
    # Holds: freq(nf,), S(nf,nc,nc)
    # Hermitian, per-block meta.
    def __init__(
        self,
        name: Optional[str] = None,
        *,
        verbose: int = 0,
    ) -> None:
        super().__init__(name=name, verbose=verbose)
        self._freq = np.zeros(0, float)
        self._S = np.zeros((0, 0, 0), complex)

        self.bw = np.zeros(0, float)
        self.avgt = np.zeros(0, float)
        self.avgf = np.zeros(0, float)
        self.rotspec = np.zeros(0, float)
        self.segnum = np.zeros(0, int)
        self.band: List[str] = []

        self.chan_ids: List[str] = []

    # ---------------- basic props
    @property
    def freq(self) -> np.ndarray:
        return self._freq

    @property
    def S(self) -> np.ndarray:
        return self._S

    @property
    def n_freq(self) -> int:
        return int(self._freq.size)

    @property
    def n_chan(self) -> int:
        return int(self._S.shape[1]) if self._S.ndim else 0

    # -------------- pack/unpack helpers
    @staticmethod
    def _unpack(vals: np.ndarray, n: int, *, empty: float) -> np.ndarray:
        v = np.asarray(vals, float)
        need = n * n
        if v.size < need:
            raise EdIDataError("SPECTRA payload too short")
        M = v[:need].reshape(n, n)
        M = np.where(M == float(empty), 0.0, M)

        H = np.zeros((n, n), complex)
        for i in range(n):
            H[i, i] = M[i, i] + 0.0j
            for j in range(i + 1, n):
                re = M[j, i]
                im = M[i, j]
                H[i, j] = re + 1j * im
                H[j, i] = re - 1j * im
        return H

    @staticmethod
    def _pack(H: np.ndarray) -> np.ndarray:
        if H.ndim != 2 or H.shape[0] != H.shape[1]:
            raise ValueError("H must be square")
        n = H.shape[0]
        M = np.zeros((n, n), float)
        for i in range(n):
            M[i, i] = float(H[i, i].real)
            for j in range(i + 1, n):
                z = H[i, j]
                M[j, i] = float(z.real)
                M[i, j] = float(z.imag)
        return M.ravel()

    # -------------- build from IO
    @classmethod
    def from_io(
        cls,
        sect: "SpectraSECT",
        io: "SpectraIO",
        *,
        empty: float = 1.0e32,
        verbose: int = 0,
    ) -> "Spectra":
        # infer nchan robustly
        nc: Optional[int] = None
    
        # 1) prefer explicit header nchan
        if getattr(sect, "nchan", None) is not None:
            try:
                nc = int(sect.nchan)  # type: ignore[arg-type]
            except Exception:
                nc = None
    
        # 2) try channel ids length
        if not nc or nc <= 0:
            ids = getattr(sect, "meas_ids", None)
            if ids:
                try:
                    nc = int(len(ids))
                except Exception:
                    nc = None
    
        # 3) try block hint or values length of the first block
        if (not nc or nc <= 0) and getattr(io, "blocks", None):
            first = next(
                (b for b in io.blocks if getattr(
                    b, "values", None)), None)
            if first is not None:
                hint = getattr(first, "nvals_hint", None)
                if isinstance(hint, (int, float)) and hint > 0:
                    root = int(round(np.sqrt(float(hint))))
                    if root * root == int(hint):
                        nc = root
                if not nc or nc <= 0:
                    nv = len(np.asarray(
                        getattr(first, "values", []), float))
                    if nv > 0:
                        root = int(round(np.sqrt(nv)))
                        if root > 0 and root * root <= nv:
                            nc = root
    
        if not nc or nc <= 0:
            raise EdIDataError("bad SPECTRA header: cannot infer nchan")
            

        if not hasattr(io, "blocks") or not io.blocks:
            raise EdIDataError("no >SPECTRA blocks")
    
        self = cls(
            name=getattr(sect, "sectid", None),
            verbose=verbose,
        )
        self.chan_ids = list(sect.meas_ids or [])
    
        def _f(x, dv=0.0) -> float:
            try:
                v = float(x)
                return v
            except Exception:
                return float(dv)
    
        def _i(x, dv=0) -> int:
            try:
                return int(float(x))
            except Exception:
                return int(dv)
    
        def _opt(blk) -> dict:
            o = getattr(blk, "options", None)
            return o if isinstance(o, dict) else {}
    
        def _get_freq(blk) -> float | None:
            # attr first
            fv = getattr(blk, "freq", None)
            if fv is None:
                opts = _opt(blk)
                for k in ("freq", "cfreq", "f", "frequency"):
                    if k in opts:
                        fv = opts.get(k)
                        break
            try:
                v = float(fv)
                return v if np.isfinite(v) else None
            except Exception:
                return None
    
        def _get_float(
            blk, *names: str, default: float | None = None
        ) -> float | None:
            # attribute names then options
            for nm in names:
                v = getattr(blk, nm, None)
                if v is not None:
                    try:
                        vv = float(v)
                        return vv
                    except Exception:
                        pass
            o = _opt(blk)
            for nm in names:
                if nm in o:
                    try:
                        vv = float(o.get(nm))
                        return vv
                    except Exception:
                        pass
            return default
    
        def _get_int(
            blk, *names: str, default: int | None = None
        ) -> int | None:
            for nm in names:
                v = getattr(blk, nm, None)
                if v is not None:
                    try:
                        return int(float(v))
                    except Exception:
                        pass
            o = _opt(blk)
            for nm in names:
                if nm in o:
                    try:
                        return int(float(o.get(nm)))
                    except Exception:
                        pass
            return default
    
        freqs: list[float] = []
        mats: list[np.ndarray] = []
        bw: list[float] = []
        avgt: list[float] = []
        avgf: list[float] = []
        rots: list[float] = []
        segnum: list[int] = []
        band: list[str] = []
    
        for blk in io.blocks:
            f0 = _get_freq(blk)
            if f0 is None:
                continue
    
            vals = np.asarray(
                getattr(blk, "values", []), float
            )
            H = cls._unpack(vals, nc, empty=empty)
    
            freqs.append(float(f0))
            mats.append(H)
    
            bw.append(
                _get_float(blk, "bw", default=0.0) or 0.0
            )
            avgt.append(
                _get_float(blk, "avgt", default=1.0) or 1.0
            )
            avgf.append(
                _get_float(blk, "avgf", default=np.nan) or np.nan
            )
            rots.append(
                _get_float(blk, "rotspec", default=np.nan)
                or np.nan
            )
            sn = _get_int(blk, "segnum", default=0)
            segnum.append(0 if sn is None else int(sn))
            bo = _opt(blk).get("band", "")
            band.append(str(bo).upper() if bo else "")
    
        n = len(freqs)
        if n == 0:
            # tolerate empty after filtering
            self._freq = np.zeros(0, float)
            self._S = np.zeros((0, nc, nc), complex)
            self.bw = np.zeros(0, float)
            self.avgt = np.zeros(0, float)
            self.avgf = np.zeros(0, float)
            self.rotspec = np.zeros(0, float)
            self.segnum = np.zeros(0, int)
            self.band = []
            return self
    
        self._freq = np.asarray(freqs, float)
        self._S = np.stack(mats, axis=0)
        self.bw = np.asarray(bw, float)
        self.avgt = np.asarray(avgt, float)
        self.avgf = np.asarray(avgf, float)
        self.rotspec = np.asarray(rots, float)
        self.segnum = np.asarray(segnum, int)
        self.band = list(band)
    
        # ensure high→low order
        if n > 1 and self._freq[-1] > self._freq[0]:
            sl = slice(None, None, -1)
            self._freq = self._freq[sl]
            self._S = self._S[sl]
            self.bw = self.bw[sl]
            self.avgt = self.avgt[sl]
            self.avgf = self.avgf[sl]
            self.rotspec = self.rotspec[sl]
            self.segnum = self.segnum[sl]
            self.band = self.band[::-1]
    
        return self


    # -------------- round-trip to IO
    def to_io(self) -> Tuple["SpectraSECT", "SpectraIO"]:
        nc = self.n_chan
        nf = self.n_freq
        sect = SpectraSECT(
            sectid=self.name,
            nchan=nc,
            nfreq=nf,
        )
        sect.meas_ids = list(self.chan_ids)

        io = SpectraIO()
        for k in range(nf):
            blk = _SpectraBlock()
            blk.freq = float(self._freq[k])
            blk.rotspec = float(self.rotspec[k])
            blk.bw = float(self.bw[k])
            blk.avgt = float(self.avgt[k])
            blk.options["avgf"] = float(self.avgf[k])
            if int(self.segnum[k]) != 0:
                blk.options["segnum"] = int(self.segnum[k])
            if self.band[k]:
                blk.options["band"] = str(self.band[k])

            blk.nvals_hint = nc * nc
            blk.values = self._pack(self._S[k]).tolist()
            io.blocks.append(blk)

        return sect, io

    # -------------- conveniences
    def matrix(self, k: int) -> np.ndarray:
        return np.array(self._S[k], copy=True)

    def psd(self, idx: int) -> np.ndarray:
        return np.asarray(self._S[:, idx, idx].real)

    def cross(self, i: int, j: int) -> np.ndarray:
        return np.asarray(self._S[:, i, j])

    def rotate(
        self,
        theta_deg: float,
        *,
        pairs: Optional[List[Tuple[int, int]]] = None,
    ) -> None:
        if self.n_chan < 2:
            return
        th = float(theta_deg) * np.pi / 180.0
        R2 = np.array(
            [[np.cos(th), -np.sin(th)],
             [np.sin(th),  np.cos(th)]],
            float,
        )
        if pairs is None:
            pairs = []
            for b in range(0, self.n_chan - 1, 2):
                pairs.append((b, b + 1))

        for a, b in pairs:
            T = np.eye(self.n_chan, float)
            T[a : b + 1, a : b + 1] = R2
            # S' = T S Tᴴ  per frequency
            self._S = np.einsum(
                "ij,fjk,lk->fil", T, self._S, T, optimize=True
            )

