# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

from pathlib import Path
from typing import Any

# import os
import numpy as np

from ..exceptions import EdIDataError
from ..log.logger import get_logger
from ..utils.validation import has_read
from ..z.base import BaseEM
from ..z.tipper import Tipper
from ..z.z import Z
from .base import EDIComponentBase
from .ops import (
    compute_errors_from_S,
    effective_dof_from_meta,
    synthesize_spectra_from_z,
)
from .validation import (
    IsEdi,
    _split_comment,
    _strip_norm,
    _to_float_or_none,
    _to_int_or_none,
)

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

    KEY_ORDER: list[str] = [
        "sectid",
        "nchan",
        "nfreq",
        "maxblks",
    ]

    def __init__(self, *args: Any, verbose: int = 0, logger=None, **kws: Any):
        super().__init__(verbose=verbose, logger=logger)
        self.sectid: str | None = None
        self.nchan: int | None = None
        self.nfreq: int | None = None
        self.maxblks: int | None = None
        self.meas_ids: list[str] = []
        self.start_data_lines_num: int | None = None

        self.id_to_chtype: dict[str, str] = {}

        for k, v in kws.items():
            setattr(self, k, v)

    @classmethod
    def from_file(cls, edi_path: str) -> SpectraSECT:
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
        inst.id_to_chtype = cls._collect_id_to_chtype(lines)

        return inst

    def write(self) -> list[str]:
        out: list[str] = [">=SPECTRASECT\n"]
        values: dict[str, Any] = {
            "sectid": self.sectid,
            "nchan": self.nchan,
            "nfreq": self.nfreq,
            "maxblks": self.maxblks,
        }
        for key in self.KEY_ORDER:
            val = values.get(key, None)
            if val in (None, "", "None"):
                continue
            out.append(f"  {key.upper()}={str(val).upper()}\n")
        if self.meas_ids:
            out.append(f"    // {len(self.meas_ids)}\n")
            for mid in self.meas_ids:
                out.append(f"     {str(mid)}\n")
        return out

    @staticmethod
    def _collect_id_to_chtype(lines: list[str]) -> dict[str, str]:
        """
        Scan >HMEAS / >EMEAS lines and build {ID -> CHTYPE}.
        IDs are kept as strings; CHTYPE is upper-cased.
        """
        out: dict[str, str] = {}
        for raw in lines:
            s = raw.lstrip()
            u = s.upper()
            if not (u.startswith(">HMEAS") or u.startswith(">EMEAS")):
                continue
            toks = s[1:].split()
            kv: dict[str, str] = {}
            for t in toks[1:]:
                if "=" in t:
                    k, v = t.split("=", 1)
                    kv[_strip_norm(k).upper()] = _strip_norm(v)
            mid = kv.get("ID", None)
            cht = (
                kv.get("CHTYPE")
                or kv.get("CH")
                or kv.get("COMP")
                or kv.get("TYPE")
            )
            if mid and cht:
                out[str(mid)] = str(cht).upper()
        return out


class _SpectraBlock(EDIComponentBase):
    """
    Single >SPECTRA block container.
    """

    def __init__(
        self, *args: Any, verbose: int | bool = 0, logger=None, **kws: Any
    ):
        super().__init__(verbose=verbose, logger=logger)
        self.freq: float | None = None
        self.rotspec: int | None = None
        self.bw: float | None = None
        self.avgt: float | None = None
        self.nvals_hint: int | None = None
        self.options: dict[str, Any] = {}
        self.values: list[float] = []

        for k, v in kws.items():
            setattr(self, k, v)

    def header_dict(self) -> dict[str, Any]:
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
        self, *args: Any, verbose: int | bool = 0, logger=None, **kws: Any
    ):
        super().__init__(verbose=verbose, logger=logger)
        self.blocks: list[_SpectraBlock] = []
        for k, v in kws.items():
            setattr(self, k, v)

    # --------------------------
    # Load from file
    # --------------------------
    @classmethod
    def from_file(
        cls,
        edi_path: str,
        start_line: int | None = None,
    ) -> SpectraIO:
        # if not os.path.isfile(edi_path):
        #     raise FileNotFoundError(
        #         f"{edi_path!r} is not a file."
        #     )
        # with open(edi_path, "r", encoding="utf-8") as f:
        #     lines = f.readlines()

        p = Path(edi_path)
        IsEdi._assert_edi(p, deep=True)

        lines = p.read_text(
            encoding="utf-8-sig", errors="replace"
        ).splitlines()

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

    @staticmethod
    def _parse_block(
        lines: list[str],
        i: int,
    ) -> tuple[_SpectraBlock, int]:
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
        per_line: int | None = None,
        float_fmt: str | None = None,
    ) -> list[str]:
        kpl = self.PER_LINE if per_line is None else per_line
        ffmt = self.FLOAT_FMT if float_fmt is None else float_fmt

        out: list[str] = []
        for blk in self.blocks:
            head = [">SPECTRA"]
            hd = blk.header_dict()
            # keep stable order in output
            for key in ("freq", "rotspec", "bw", "avgt"):
                if key in hd and hd[key] is not None:
                    head.append(f"{key.upper()}={hd[key]}")
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

            line_vals: list[str] = []
            cnt = 0
            for v in blk.values:
                line_vals.append(ffmt.format(v))
                cnt += 1
                if cnt == kpl:
                    out.append("  " + " ".join(line_vals) + "\n")
                    line_vals = []
                    cnt = 0
            if line_vals:
                out.append("  " + " ".join(line_vals) + "\n")

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
    def read_blocks(cls, edi_fn: str) -> SpectraIO:
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
        name: str | None = None,
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
        self.band: list[str] = []

        self.chan_ids: list[str] = []
        self.id_to_chtype: dict[str, str] = {}

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

    @classmethod
    def from_io(
        cls,
        sect: SpectraSECT,
        io: SpectraIO,
        *,
        empty: float = 1.0e32,
        verbose: int = 0,
    ) -> Spectra:
        # infer nchan robustly
        nc: int | None = None

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
                (b for b in io.blocks if getattr(b, "values", None)), None
            )
            if first is not None:
                hint = getattr(first, "nvals_hint", None)
                if isinstance(hint, (int, float)) and hint > 0:
                    root = int(round(np.sqrt(float(hint))))
                    if root * root == int(hint):
                        nc = root
                if not nc or nc <= 0:
                    nv = len(np.asarray(getattr(first, "values", []), float))
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
        # carry DefineMeas mapping (if present on the header)
        idm = getattr(sect, "id_to_chtype", None)
        if isinstance(idm, dict) and idm:
            try:
                self.id_to_chtype = dict(idm)
            except Exception:
                pass

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

            vals = np.asarray(getattr(blk, "values", []), float)
            H = cls._unpack(vals, nc, empty=empty)

            freqs.append(float(f0))
            mats.append(H)

            # bw
            v = _get_float(blk, "bw", default=0.0)
            bw.append(0.0 if v is None else float(v))

            # avgt
            v = _get_float(blk, "avgt", default=1.0)
            avgt.append(1.0 if v is None else float(v))

            # avgf
            v = _get_float(blk, "avgf", default=np.nan)
            avgf.append(np.nan if v is None else float(v))

            # rotspec  (this was causing your failure)
            v = _get_float(blk, "rotspec", default=np.nan)
            rots.append(np.nan if v is None else float(v))

            # segnum
            sn = _get_int(blk, "segnum", default=0)
            segnum.append(0 if sn is None else int(sn))

            # band
            bo = _opt(blk).get("band", "")
            band.append(str(bo).upper() if bo else "")

        n = len(freqs)
        if n == 0:
            # tolerate empty after filtering
            self._freq = np.zeros(0, dtype=float)
            self._S = np.zeros((0, nc, nc), dtype=complex)
            self.bw = np.zeros(0, dtype=float)
            self.avgt = np.zeros(0, dtype=float)
            self.avgf = np.zeros(0, dtype=float)
            self.rotspec = np.zeros(0, dtype=float)
            self.segnum = np.zeros(0, dtype=int)
            self.band = []
            return self

        self._freq = np.asarray(freqs, dtype=float)
        self._S = np.stack(mats, axis=0)
        self.bw = np.asarray(bw, dtype=float)
        self.avgt = np.asarray(avgt, dtype=float)
        self.avgf = np.asarray(avgf, dtype=float)
        self.rotspec = np.asarray(rots, dtype=float)
        self.segnum = np.asarray(segnum, dtype=int)
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

    @classmethod
    def from_file(
        cls,
        path: str,
        *,
        empty: float = 1.0e32,
        verbose: int = 0,
    ) -> Spectra:
        """Read a :class:`Spectra` directly from an EDI file path.

        Convenience wrapper around :meth:`from_io` that calls
        :class:`SpectraSECT.from_file` and :class:`SpectraIO.from_file`
        internally.

        Parameters
        ----------
        path : str or Path
            Path to the EDI file containing ``>=SPECTRASECT`` and
            ``>SPECTRA`` blocks.
        empty : float
            Sentinel value for missing spectra entries.  Default 1e32.
        verbose : int
            Verbosity level forwarded to :meth:`from_io`.

        Returns
        -------
        Spectra
        """
        sect = SpectraSECT.from_file(str(path))
        sio = SpectraIO.from_file(
            str(path), start_line=sect.start_data_lines_num
        )
        return cls.from_io(sect, sio, empty=empty, verbose=verbose)

    def to_edi(
        self,
        source_edi=None,
        *,
        station_name: str | None = None,
        e_labels: tuple[str, str] = ("EX", "EY"),
        h_labels: tuple[str, str] = ("HX", "HY"),
        ridge: float | None = None,
        estimate_error: bool = False,
        dof: Any | None = None,
    ) -> Any:
        r"""Convert cross-spectra to an MT-impedance :class:`~pycsamt.seg.edi.EDIFile`.

        Calls :meth:`to_Z` and assembles a complete
        ``>=MTSECT`` / ``>FREQ`` / ``>ZXXR`` / ``>ZXYR`` / … EDI
        ready to be saved with :meth:`~pycsamt.seg.edi.EDIFile.write`.

        The structural sections (``>HEAD``, ``>INFO``, ``>=DEFINEMEAS``)
        are *re-used* from *source_edi* when provided, preserving all
        acquisition metadata; otherwise a minimal header is synthesised
        from the :class:`Spectra` metadata.

        According to the SEG EDI standard (§ 7.53, 12.1), an MT data
        section requires:

        .. code-block:: text

            >=MTSECT
              SECTID=...
              NFREQ=...
              HX=...  HY=...  HZ=...  EX=...  EY=...
            >FREQ   //N
              ...
            >ZXXR  ROT=ZROT  //N
              ...
            >ZXXI  ROT=ZROT  //N
              ...
            ...
            >END

        The measurement IDs for ``HX``, ``HY``, … in ``>=MTSECT`` are
        resolved from :attr:`id_to_chtype` (populated by
        :class:`SpectraSECT` from ``>HMEAS`` / ``>EMEAS`` lines).

        Parameters
        ----------
        source_edi : str, Path, or EDIFile, optional
            Spectra EDI file whose ``>HEAD``, ``>INFO``, and
            ``>=DEFINEMEAS`` sections are copied into the output.
            Pass the same path used with :meth:`from_file` to produce
            a fully metadata-rich result.  When ``None``, a minimal
            header is synthesised.
        station_name : str, optional
            Override for the ``DATAID`` in ``>HEAD`` and ``SECTID`` in
            ``>=MTSECT``.  Defaults to :attr:`name`.
        e_labels : tuple of str
            Electric channel type labels forwarded to :meth:`to_Z`.
        h_labels : tuple of str
            Horizontal magnetic channel type labels forwarded to :meth:`to_Z`.
        ridge : float, optional
            Tikhonov regularisation forwarded to :meth:`to_Z`.
        estimate_error : bool
            If ``True``, propagate 1-σ errors into ``>ZXX.VAR`` … blocks.
        dof : float or ndarray, optional
            Effective degrees of freedom forwarded to :meth:`to_Z`.

        Returns
        -------
        EDIFile
            Fully populated MT-impedance container.  Call
            :meth:`~pycsamt.seg.edi.EDIFile.write` to save.

        Raises
        ------
        EdIDataError
            If :meth:`to_Z` fails (channel types not resolved, singular
            magnetic block, etc.).

        Examples
        --------
        Convert and save::

            sp = Spectra.from_file("site.edi")
            ed = sp.to_edi("site.edi", estimate_error=False)
            out = ed.write(savepath="mt_output/")

        Convert with a custom station name and error propagation::

            ed = sp.to_edi(
                "site.edi",
                station_name="HBH03_imp",
                estimate_error=True,
                dof=24.0,
            )
            out = ed.write(savepath="mt_output/")

        Verify the round-trip::

            from pycsamt.seg.edi import EDIFile

            ed2 = EDIFile(out)
            assert ed2.Z.n_freq == sp.n_freq
        """
        # deferred import to avoid circular deps at module level
        from .edi import EDIFile as _EDIFile
        from .mtemap import MTEMAP as _MTEMAP

        # ── Step 1: recover Z and Tipper from cross-spectra ──────────────
        z_obj, tip = self.to_Z(
            e_labels=e_labels,
            h_labels=h_labels,
            ridge=ridge,
            estimate_error=estimate_error,
            dof=dof,
        )

        # ── Step 2: build or load the structural EDI template ─────────────
        if isinstance(source_edi, _EDIFile):
            ed = source_edi
        elif source_edi is not None:
            ed = _EDIFile(str(source_edi), verbose=0)
        else:
            # No source: create a blank container; caller must fill header
            ed = _EDIFile.__new__(_EDIFile)
            ed.path = None
            ed.verbose = 0
            ed._init_registry()
            ed._data_start = None
            from ..z.resphase import ResPhase as _RP
            from ..z.tipper import Tipper as _Tip
            from ..z.z import Z as _Z

            ed.Z = _Z(verbose=0)
            ed.Res = _RP(verbose=0)
            ed.Tip = _Tip()
            ed.block_size = 6
            ed.float_fmt = "{: .6E}"
            ed.header_tpl = ">!****{title}****!\n"

        # ── Step 3: build >=MTSECT section ───────────────────────────────
        # Invert id_to_chtype: CHTYPE → first matching measurement ID
        id_map: dict[str, str] = dict(getattr(self, "id_to_chtype", {}) or {})
        chtype_to_id: dict[str, str] = {}
        for mid, cht in id_map.items():
            key = "".join(c for c in str(cht).upper() if c.isalpha())
            if key and key not in chtype_to_id:
                chtype_to_id[key] = str(mid)

        # Also try chan_ids directly when they look like channel type labels
        if not chtype_to_id and self.chan_ids:
            for cid in self.chan_ids:
                key = "".join(c for c in str(cid).upper() if c.isalpha())
                if key in {"HX", "HY", "HZ", "EX", "EY", "RX", "RY"}:
                    chtype_to_id[key] = str(cid)

        mtsect_kw: dict[str, Any] = {}
        for attr in ("hx", "hy", "hz", "ex", "ey", "rx", "ry"):
            cht = attr.upper()
            if cht in chtype_to_id:
                mtsect_kw[attr] = chtype_to_id[cht]

        sid = station_name or self.name or "SITE"
        mtsect = _MTEMAP(
            sectid=str(sid),
            nfreq=int(z_obj.n_freq),
            **mtsect_kw,
        )
        ed.add_section("mtsect", mtsect)

        # ── Step 4: attach Z and Tipper ───────────────────────────────────
        ed.Z = z_obj
        if tip is not None:
            ed.Tip = tip

        # ── Step 5: propagate station name to >HEAD if possible ───────────
        if station_name:
            try:
                ed.station = str(station_name)
            except Exception:
                pass

        return ed

    # -------------- round-trip to IO
    def to_io(self) -> tuple[SpectraSECT, SpectraIO]:
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
        pairs: list[tuple[int, int]] | None = None,
    ) -> None:
        if self.n_chan < 2:
            return
        th = float(theta_deg) * np.pi / 180.0
        R2 = np.array(
            [[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]],
            float,
        )
        if pairs is None:
            pairs = []
            for b in range(0, self.n_chan - 1, 2):
                pairs.append((b, b + 1))

        for a, b in pairs:
            T = np.eye(self.n_chan, dtype=float)
            T[a : b + 1, a : b + 1] = R2
            # S' = T S Tᴴ  per frequency
            self._S = np.einsum("ij,fjk,lk->fil", T, self._S, T, optimize=True)

    def to_Z(
        self,
        *,
        id_to_chtype: dict[str, str] | None = None,
        e_labels: tuple[str, str] = ("EX", "EY"),
        h_labels: tuple[str, str] = ("HX", "HY"),
        use_remote: bool = False,
        ridge: float | None = None,
        estimate_error: bool = True,
        dof: float | np.ndarray | None = None,
    ):
        r"""
        Recover an impedance tensor ``Z`` and, if available, the
        tipper from cross-spectra stored in this :class:`Spectra`.

        The method resolves channel types, extracts the electric
        and magnetic sub-blocks, and computes per-frequency
        ``Z = S_EH @ inv(S_HH)``.  If a vertical magnetic channel
        is present it also computes the tipper
        ``T = S_ZH @ inv(S_HH)``.  Optional ridge regularization
        can be applied to stabilize the magnetic block.

        Parameters
        ----------
        id_to_chtype : dict of str to str, optional
            Mapping from measurement IDs (as recorded in
            ``>=SPECTRASECT`` or ``DefineMeas``) to channel
            types (``"HX"``, ``"HY"``, ``"HZ"``, ``"EX"``,
            ``"EY"``).  If omitted, the method uses
            ``self.id_to_chtype`` when available, otherwise it
            interprets ``self.chan_ids`` directly as labels.
        e_labels : tuple of str, default (``"EX"``, ``"EY"``)
            Labels that identify the two electric channels used
            for the ``E`` block.
        h_labels : tuple of str, default (``"HX"``, ``"HY"``)
            Labels that identify the two horizontal magnetic
            channels used for the ``H`` block.
        use_remote : bool, default ``False``
            When duplicate electric channels exist (e.g., local
            and remote), choose the second occurrence for the
            ``E`` block if ``True``; otherwise choose the first.
        ridge : float, optional
            Non-negative Tikhonov regularization added to
            ``S_HH`` prior to inversion, ``S_HH + ridge * I``.
        estimate_error : bool, default ``True``
            If ``True``, estimate per-component 1-sigma standard
            errors for ``Z`` (and tipper when available) using
            ``compute_errors_from_S`` and the degrees of freedom
            given by ``dof`` (or inferred; see Notes).
        dof : float or ndarray, optional
            Effective degrees of freedom per frequency.  If an
            array is provided it must broadcast to ``n_freq``.
            If ``None`` and ``estimate_error`` is ``True``, the
            method tries to infer DoF from metadata via
            ``effective_dof_from_meta`` using ``segnum``, or
            ``avgt * bw`` as a fallback.

        Returns
        -------
        z_obj : :class:`pycsamt.z.z.Z`
            Impedance object on the spectra frequency grid with
            ``z`` populated and, when estimated, ``z_err`` set.
        tip : :class:`pycsamt.z.tipper.Tipper` or ``None``
            Tipper on the same grid when ``HZ`` is available.
            When errors are estimated, tipper uncertainties are
            attached.

        Raises
        ------
        EdIDataError
            If spectra are empty, channel types cannot be
            resolved, or the stabilized magnetic block is
            singular.

        Notes
        -----
        Per frequency, ``Z`` is formed as ``Z = S_EH @ inv(S_HH)``,
        where ``S_EH`` is the cross-spectra between E and H, and
        ``S_HH`` is the magnetic auto/cross block. If a vertical
        magnetic channel is available, the tipper is computed as
        ``T = S_ZH @ inv(S_HH)``.

        Channel type resolution proceeds in this order:

        1. explicit ``id_to_chtype`` argument,
        2. ``self.id_to_chtype`` from the section header or
           ``DefineMeas``,
        3. direct interpretation of ``self.chan_ids``.

        When both local and remote electric channels are present,
        setting ``use_remote=True`` chooses the second occurrence
        as a simple heuristic.  Frequency ordering is preserved.

        Uncertainties are computed by first-order propagation
        under a complex-Wishart model and scale as
        ``1 / sqrt(DoF)``.  If DoF cannot be determined for every
        frequency, no per-frequency array is attached at all:
        ``z_err`` (and ``tip.tipper_err`` when applicable) come
        back as ``None`` rather than an array of NaN.

        Examples
        --------

        >>> Zhat, That = spectra.to_Z()
        >>> Zhat, _ = spectra.to_Z(use_remote=True, ridge=1e-6)
        >>> Zhat, _ = spectra.to_Z(dof=np.full(spectra.n_freq, 24.0))

        See Also
        --------
        Spectra.from_Z
            Inverse operation that synthesizes spectra.
        spectra_from_Z
            Functional wrapper for the inverse operation.
        effective_dof_from_meta
            Infer DoF from ``segnum``, ``avgt`` and ``bw``.
        compute_errors_from_S
            Per-frequency uncertainty estimator.

        References
        ----------
        .. [1] Chave, A. D., & Jones, A. G. (2012). *The
               Magnetotelluric Method: Theory and Practice*.
               Cambridge University Press.
        .. [2] Bendat, J. S., & Piersol, A. G. (2011). *Random
               Data: Analysis and Measurement Procedures*. Wiley.
        """
        has_read(
            self,
            msg="Spectra not populated; call from_file()/from_io() first.",
        )

        if self.n_freq == 0 or self.n_chan < 2:
            raise EdIDataError("Spectra is empty or too small.")

        # default to mapping collected on SpectraSECT (if any)
        if id_to_chtype is None:
            id_to_chtype = getattr(self, "id_to_chtype", None)

        def _norm(s: object) -> str:
            t = "".join(ch for ch in str(s).upper() if ch.isalpha())
            return t

        # Resolve per-channel types in the order used by S (self.chan_ids)
        # - If the file stored numeric IDs in SPECTRASECT, you *must* pass
        #   `id_to_chtype` (built from >HMEAS/>EMEAS). Otherwise, we accept
        #   chan_ids that are already like ["HX","HY","HZ","EX","EY",...].

        if id_to_chtype:
            kinds_raw = [
                id_to_chtype.get(str(mid), "") for mid in self.chan_ids
            ]
        else:
            kinds_raw = list(self.chan_ids)

        kinds = [_norm(x) for x in kinds_raw]

        # pick indices for H and E
        def _all_idx(lbl: str) -> list[int]:
            L = _norm(lbl)
            return [
                i for i, k in enumerate(kinds) if (k == L) or k.startswith(L)
            ]

        ex_all, ey_all = _all_idx(e_labels[0]), _all_idx(e_labels[1])
        hx_all, hy_all = _all_idx(h_labels[0]), _all_idx(h_labels[1])
        hz_all = _all_idx("HZ")

        def _choose(pair):
            # choose the first by default; second if remote requested
            if not pair:
                return None
            return pair[1] if (use_remote and len(pair) > 1) else pair[0]

        idx_ex = _choose(ex_all)
        idx_ey = _choose(ey_all)
        idx_hx = _choose(hx_all)
        idx_hy = _choose(hy_all)
        idx_hz = _choose(hz_all)  # optional

        if None in (idx_ex, idx_ey, idx_hx, idx_hy):
            raise EdIDataError(
                "Missing required channels to build Z "
                f"(EX={idx_ex}, EY={idx_ey}, HX={idx_hx}, HY={idx_hy})."
            )

        h_idx = [idx_hx, idx_hy]
        e_idx = [idx_ex, idx_ey]

        # Build Z(f) = S_EH @ inv(S_HH) per frequency; optionally
        # estimate errors using DoF.  Preallocate outputs.
        z_arr = np.zeros((self.n_freq, 2, 2), dtype=complex)
        z_err = (
            np.full((self.n_freq, 2, 2), np.nan, dtype=float)
            if estimate_error
            else None
        )

        tip_arr = (
            None
            if idx_hz is None
            else np.zeros((self.n_freq, 1, 2), dtype=complex)
        )
        tip_err = (
            None
            if (idx_hz is None or not estimate_error)
            else np.full((self.n_freq, 1, 2), np.nan, float)
        )

        I2 = np.eye(2, dtype=float)

        for k in range(self.n_freq):
            S = self._S[k]  # (nchan, nchan), Hermitian complex
            S_HH = S[np.ix_(h_idx, h_idx)].astype(complex)
            S_EH = S[np.ix_(e_idx, h_idx)].astype(complex)

            if ridge is not None and ridge > 0.0:
                S_HH = S_HH + float(ridge) * I2

            try:
                inv_SHH = np.linalg.inv(S_HH)
            except np.linalg.LinAlgError as exc:
                raise EdIDataError(f"S_HH singular at k={k}: {exc}") from exc

            z_arr[k] = S_EH @ inv_SHH

            if tip_arr is not None:
                S_ZH = S[np.ix_([idx_hz], h_idx)].astype(complex)
                tip_arr[k, 0, :] = S_ZH @ inv_SHH

            # error estimates (optional)
            if not estimate_error or z_err is None:
                continue

            # Resolve DoF M_k:
            # priority: explicit `dof` -> `segnum` -> avgt*bw.
            M_k = None
            if dof is not None:
                M_k = float(np.asarray(dof)[k])
            else:
                M_k = effective_dof_from_meta(
                    segnum=(
                        self.segnum[k] if hasattr(self, "segnum") else None
                    ),
                    avgt=(self.avgt[k] if hasattr(self, "avgt") else None),
                    bw=(self.bw[k] if hasattr(self, "bw") else None),
                )
                if M_k is not None:
                    M_k = float(M_k)

            if M_k and M_k > 0.0:
                z_e, t_e = compute_errors_from_S(
                    S=S,
                    e_idx=tuple(e_idx),
                    h_idx=tuple(h_idx),
                    hz_idx=idx_hz,
                    M=M_k,
                    ridge=ridge,
                )
                z_err[k] = z_e
                if tip_err is not None and t_e is not None:
                    tip_err[k, 0, :] = t_e

        # --- finalize/attach uncertainties with user-facing warnings
        valid_ze = (
            estimate_error
            and (z_err is not None)
            and np.all(np.isfinite(z_err))
            and np.all(z_err >= 0.0)
        )

        valid_te = (
            (tip_arr is not None)
            and estimate_error
            and (tip_err is not None)
            and np.all(np.isfinite(tip_err))
            and np.all(tip_err >= 0.0)
        )

        if estimate_error and not valid_ze:
            n_nan_ze = (
                0
                if z_err is None
                else int(np.count_nonzero(~np.isfinite(z_err)))
            )
            n_neg_ze = (
                0 if z_err is None else int(np.count_nonzero(z_err < 0.0))
            )
            logger.warning(
                "Z error estimation skipped: invalid entries "
                "(nan=%d, neg=%d) and/or missing DoF. "
                "Errors will not be attached. Consider passing "
                "`dof`, or disable with `estimate_error=False`.",
                n_nan_ze,
                n_neg_ze,
            )
            z_err = None

        if estimate_error and (tip_arr is not None) and not valid_te:
            n_nan_te = (
                0
                if tip_err is None
                else int(np.count_nonzero(~np.isfinite(tip_err)))
            )
            n_neg_te = (
                0 if tip_err is None else int(np.count_nonzero(tip_err < 0.0))
            )
            logger.warning(
                "Tipper error estimation skipped: invalid entries "
                "(nan=%d, neg=%d) and/or missing DoF. "
                "Errors will not be attached.",
                n_nan_te,
                n_neg_te,
            )
            tip_err = None

        # --- build objects (attach uncertainties only when valid)
        nm = getattr(self, "name", None)
        if z_err is not None:
            z_obj = Z(
                z_array=z_arr,
                freq=self.freq,
                name=nm,
                z_err_array=z_err,
                verbose=self.verbose,
            )
        else:
            z_obj = Z(
                z_array=z_arr,
                freq=self.freq,
                name=nm,
                verbose=self.verbose,
            )

        tip_obj = None
        if tip_arr is not None:
            tip_obj = Tipper()
            tip_obj._freq = np.array(self.freq, dtype=float)
            tip_obj._tipper = tip_arr
            if tip_err is not None:
                tip_obj._tipper_err = tip_err
            # compute derived quantities (no uncertainties needed)
            tip_obj.compute_amp_phase()
            tip_obj.compute_mag_direction()

        return z_obj, tip_obj

    def __has_read__(self) -> bool:
        """
        Lightweight check used by utils.validation.has_read().
        Returns True only if a non-empty spectral stack is present
        and shapes are coherent.
        """
        try:
            # freq present and 1-D
            if not isinstance(self._freq, np.ndarray) or self._freq.ndim != 1:
                return False
            if self._freq.size == 0:
                return False
            # spectral cube present and consistent
            S = getattr(self, "_S", None)
            if not isinstance(S, np.ndarray) or S.ndim != 3:
                return False
            nf, n1, n2 = S.shape
            if nf != self._freq.size or n1 != n2 or n1 < 2:
                return False
            return True
        except Exception:
            return False

    @classmethod
    def from_Z(
        cls,
        z_obj: Z,
        **kws: Any,
    ) -> Spectra:
        r"""
        Create a :class:`Spectra` from a transfer function
        :class:`~pycsamt.z.z.Z`.

        This class method is a thin, convenience wrapper
        around :func:`spectra_from_Z`.  It synthesizes a full
        Hermitian cross–spectral density tensor from the
        impedance tensor ``Z(f)`` and optional inputs that
        control magnetic power and tipper usage.

        Parameters
        ----------
        z_obj : :class:`~pycsamt.z.z.Z`
            Input impedance object.  The attributes
            ``z_obj.z`` (shape ``(n, 2, 2)``) and
            ``z_obj.freq`` (shape ``(n,)``) must be set.
        **kws : Any
            Forwarded to :func:`spectra_from_Z`.  See that
            function for the complete set of options such as
            ``S_HH``, ``H_psd``, ``tipper``, ``include_hz``,
            and ``chan_order``.

        Returns
        -------
        Spectra
            A spectra container on the same frequency grid as
            ``z_obj``.  Channel order follows the requested
            ``chan_order`` (default: ``HX, HY, EX, EY``).

        Raises
        ------
        EdIDataError
            If ``z_obj`` is incomplete (missing ``z`` or
            ``freq``).

        Notes
        -----
        Absolute spectral levels are not carried by the
        impedance tensor.  To obtain physically scaled
        spectra, provide magnetic spectra via ``S_HH`` or
        ``H_psd``.  If neither is given, a unit–power
        assumption is used (``S_HH = I``), which is suitable
        for tests but not for quantitative analysis.

        This method does not infer per–frequency metadata
        such as bandwidth or averaging time; those fields
        are initialized with zeros/NaNs.

        Examples
        --------
        >>> ed = EDIFile("site_imp.edi")
        >>> sp = Spectra.from_Z(
        ...     ed.Z,
        ...     H_psd=(np.ones(ed.Z.n_freq), np.ones(ed.Z.n_freq), None),
        ... )
        >>> sect, io = sp.to_io()
        >>> _ = ed.write_new_edi(
        ...     edi_fn="site_with_synth_spec.edi",
        ...     Spectra=sp,
        ... )

        See Also
        --------
        spectra_from_Z
            Functional API that performs the synthesis.
        pycsamt.seg.ops.synthesize_spectra_from_z
            Low–level array helper used under the hood.
        Spectra.to_Z
            Inverse operation (spectra → Z).

        References
        ----------
        .. [1] Chave, A. D., & Jones, A. G. (2012). *The
               Magnetotelluric Method: Theory and Practice*.
               Cambridge Univ. Press.
        .. [2] Bendat, J. S., & Piersol, A. G. (2011).
               *Random Data: Analysis and Measurement
               Procedures*. Wiley.
        .. [3] SEG EDI MT/EMAP standard (1987). MTNet.
        """

        return spectra_from_Z(z_obj=z_obj, **kws)


def spectra_from_Z(
    z_obj: Z,
    *,
    S_HH: np.ndarray | None = None,
    H_psd: tuple[
        np.ndarray,
        np.ndarray,
        np.ndarray | None,
    ]
    | None = None,
    tipper: Tipper | np.ndarray | None = None,
    include_hz: bool = False,
    chan_order: tuple[str, ...] = (
        "HX",
        "HY",
        "EX",
        "EY",
    ),
    e_noise: float | np.ndarray | None = None,
    h_noise: float | np.ndarray | None = None,
    name: str | None = None,
    verbose: int = 0,
) -> Spectra:
    r"""
    Synthesize a :class:`Spectra` from a transfer function
    :class:`~pycsamt.z.z.Z` and optional tipper.

    Given the relation ``E = Z H``, the cross–spectral
    blocks are

    * ``S_EH = Z S_HH``
    * ``S_EE = Z S_HH Z^H``

    and, if tipper ``T`` and HZ are supplied,

    * ``S_ZH = T S_HH``, ``S_ZZ = T S_HH T^H``,
    * ``S_EZ = Z S_HH T^H``.

    If magnetic spectra are not provided, a unit–power
    assumption (``S_HH = I``) is used.  The result is
    Hermitian by construction.

    Parameters
    ----------
    z_obj : :class:`~pycsamt.z.z.Z`
        Input impedance with ``z`` (``(n, 2, 2)``) and
        ``freq`` (``(n,)``) populated.
    S_HH : ndarray, optional
        Magnetic spectra per frequency, shape ``(n, 2, 2)``.
        Must be Hermitian.  Overrides ``H_psd`` if both are
        provided.
    H_psd : tuple of arrays, optional
        ``(Pxx, Pyy, Pxy)`` with ``Pxx, Pyy`` real
        ``(n,)`` and optional complex ``Pxy (n,)``.  Used
        to assemble ``S_HH`` when ``S_HH`` is not given.
    tipper : :class:`~pycsamt.z.tipper.Tipper` or ndarray, \
    optional
        Horizontal tipper as object or array.  Arrays must
        have shape ``(n, 1, 2)`` or ``(n, 2)``.  Required
        when ``include_hz`` is ``True`` or when ``HZ`` is
        in ``chan_order``.
    include_hz : bool, default False
        If ``True``, include HZ and all cross–terms.
        Requires ``tipper``.
    chan_order : tuple of str, default (\"HX\",\"HY\",\"EX\",\"EY\")
        Order of channels in the synthesized spectra.  If
        ``HZ`` is included, the tipper must be provided.
    e_noise : float or ndarray, optional
        Diagonal noise power added to the electric block
        ``S_EE``.  Scalar or ``(n,)``.
    h_noise : float or ndarray, optional
        Diagonal noise power added to the magnetic block
        ``S_HH``.  Scalar or ``(n,)``.
    name : str, optional
        Display name for the resulting :class:`Spectra`.
        Defaults to ``z_obj.name``.
    verbose : int, default 0
        Verbosity flag propagated to :class:`Spectra`.

    Returns
    -------
    Spectra
        A spectra container on ``z_obj.freq``.  The
        attribute ``chan_ids`` matches ``chan_order``.

    Raises
    ------
    EdIDataError
        If ``z_obj`` lacks ``z`` or ``freq``, or if HZ is
        requested without a tipper.
    ValueError
        If provided arrays have incompatible shapes.

    Notes
    -----
    Absolute scaling depends entirely on the magnetic
    spectra.  Without ``S_HH`` or ``H_psd``, the result
    is unit–power and useful mainly for testing or
    structural checks.

    Per–frequency metadata (e.g., bandwidth, averaging
    time) are not inferred and are initialized with
    zeros/NaNs.  Channel order is honored exactly; the
    resulting spectra are explicitly symmetrized to be
    Hermitian.

    Examples
    --------
    Minimal synthesis with unit power:

    >>> ed = EDIFile("site_imp.edi")
    >>> sp = spectra_from_Z(ed.Z)

    With magnetic PSD and tipper, including HZ:

    >>> Pxx = np.full(ed.Z.n_freq, 1e-2)
    >>> Pyy = np.full(ed.Z.n_freq, 1e-2)
    >>> sp = spectra_from_Z(
    ...     ed.Z,
    ...     H_psd=(Pxx, Pyy, None),
    ...     tipper=ed.Tip,
    ...     include_hz=True,
    ...     chan_order=("HX","HY","HZ","EX","EY"),
    ... )

    Writing synthesized spectra to EDI:

    >>> sect, io = sp.to_io()
    >>> _ = ed.write_new_edi(
    ...     edi_fn="site_with_synth_spec.edi",
    ...     Spectra=sp,
    ... )

    See Also
    --------
    Spectra.from_Z
        Classmethod wrapper returning a spectra from Z.
    pycsamt.seg.ops.synthesize_spectra_from_z
        Low–level array constructor used internally.
    Spectra.to_Z
        Inverse operation (spectra → Z and tipper).

    References
    ----------
    .. [1] Chave, A. D., & Jones, A. G. (2012). *The
           Magnetotelluric Method: Theory and Practice*.
           Cambridge Univ. Press.
    .. [2] Bendat, J. S., & Piersol, A. G. (2011).
           *Random Data: Analysis and Measurement
           Procedures*. Wiley.
    .. [3] SEG EDI MT/EMAP standard (1987). MTNet.
    """

    if (
        getattr(z_obj, "z", None) is None
        or getattr(z_obj, "freq", None) is None
    ):
        raise EdIDataError("Z object is incomplete (z or freq missing).")

    tip_arr = None
    if tipper is not None:
        tip_arr = getattr(tipper, "tipper", tipper)

    Sfull, order = synthesize_spectra_from_z(
        z_obj.z,
        S_HH=S_HH,
        H_psd=H_psd,
        tipper=tip_arr,
        include_hz=include_hz,
        chan_order=chan_order,
        e_noise=e_noise,
        h_noise=h_noise,
    )

    sp = Spectra(
        name=(name or getattr(z_obj, "name", None)),
        verbose=verbose,
    )
    sp._freq = np.asarray(z_obj.freq, float)
    sp._S = np.asarray(Sfull, complex)
    sp.chan_ids = list(order)

    nf = sp.n_freq
    sp.bw = np.zeros(nf, float)
    sp.avgt = np.zeros(nf, float)
    sp.avgf = np.full(nf, np.nan, float)
    sp.rotspec = np.full(nf, np.nan, float)
    sp.segnum = np.zeros(nf, int)
    sp.band = [""] * nf
    return sp
