# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
from __future__ import annotations

from typing import List, Optional, Dict, Union
import os
import re

import numpy as np

from ..log.logger import get_logger
from ..exceptions import EdIDataError
from .base import EdiComponentBase

logger = get_logger(__name__)


def _norm(s: str) -> str:
    """Trim and strip quotes."""
    return s.strip().strip('"').strip("'")


def _to_int(v: Optional[str]) -> Optional[int]:
    if v in (None, "", "None"):
        return None
    try:
        return int(float(v))
    except Exception:
        return None


def _to_float(v: Optional[str]) -> Optional[float]:
    if v in (None, "", "None"):
        return None
    try:
        return float(v)
    except Exception:
        return None


def _is_section_start(line_upper: str) -> bool:
    # Any new major block/section marker
    return (
        line_upper.startswith(">=")  # new section header
        or line_upper.startswith(">FREQ")
        or line_upper.startswith(">ZROT")
        or line_upper.startswith(">RHOROT")
        or line_upper.startswith(">!")  # comment/data intro
        or line_upper.startswith(">SPECTRA")  # next spectra block
    )


class SpectraBlock(EdiComponentBase):
    """
    One >SPECTRA block.

    Attributes
    ----------
    freq : float            # Hz (required)
    bw   : float            # Hz (required)
    nchan: int              # if absent in header, fallback to section.nchan
    rotspec : Optional[float]
    avgt    : Optional[float]  # can be non-integer in practice
    avgf    : Optional[float]
    band    : Optional[str]
    segnum  : Optional[int]    # default 0
    chksum  : Optional[float]
    matrix  : np.ndarray (nchan x nchan, complex)  # crosspower (Hermitian)
    """

    def __init__(
        self,
        *,
        freq: float,
        bw: float,
        nchan: int,
        matrix: np.ndarray,
        rotspec: Optional[float] = None,
        avgt: Optional[float] = 1.0,
        avgf: Optional[float] = 1.0,
        band: Optional[str] = None,
        segnum: Optional[int] = 0,
        chksum: Optional[float] = None,
        sectid: Optional[str] = None,
    ) -> None:
        super().__init__()
        self.sectid = sectid
        self.freq = float(freq)
        self.bw = float(bw)
        self.nchan = int(nchan)
        self.rotspec = None if rotspec is None else float(rotspec)
        self.avgt = None if avgt is None else float(avgt)
        self.avgf = None if avgf is None else float(avgf)
        self.band = band
        self.segnum = None if segnum is None else int(segnum)
        self.chksum = None if chksum is None else float(chksum)

        if not isinstance(matrix, np.ndarray):
            raise EdIDataError("SpectraBlock.matrix must be a numpy ndarray.")
        if matrix.shape != (self.nchan, self.nchan):
            raise EdIDataError(
                "SpectraBlock.matrix shape mismatch:"
                f" expected ({self.nchan},{self.nchan}), "
                f"got {matrix.shape}"
            )
        self.matrix = matrix.astype(np.complex128, copy=False)

    # Packing / Unpacking between Hermitian and 
    # SEG compressed real N×N 

    @staticmethod
    def pack_to_compressed_real(m: np.ndarray) -> np.ndarray:
        """
        Pack Hermitian NxN complex matrix into NxN real "compressed" array:
        - diagonal = real(diag)
        - lower-left (i>j) = real(m[i,j])
        - upper-right (i<j) = imag(m[i,j])  (sign preserved)

        Returns
        -------
        np.ndarray of shape (N, N), dtype=float64
        """
        if m.ndim != 2 or m.shape[0] != m.shape[1]:
            raise EdIDataError("Matrix must be square for spectra packing.")
        n = m.shape[0]
        out = np.zeros((n, n), dtype=np.float64)
        # Diagonal (real autos)
        out[np.diag_indices(n)] = m.diagonal().real
        # Lower-left: real parts
        il, jl = np.tril_indices(n, k=-1)
        out[il, jl] = m[il, jl].real
        # Upper-right: imaginary parts
        iu, ju = np.triu_indices(n, k=1)
        out[iu, ju] = m[iu, ju].imag
        return out

    @staticmethod
    def unpack_from_compressed_real(a: np.ndarray) -> np.ndarray:
        """
        Reconstruct Hermitian matrix from NxN real "compressed" array:
        - diag = real autos
        - real parts from lower-left
        - imaginary parts from upper-right (sign preserved)

        Returns complex Hermitian NxN matrix.
        """
        if a.ndim != 2 or a.shape[0] != a.shape[1]:
            raise EdIDataError("Compressed array must be square.")
        n = a.shape[0]
        m = np.zeros((n, n), dtype=np.complex128)
        # diagonal
        d = np.diag_indices(n)
        m[d] = a[d] + 0j
        # off-diags
        iu, ju = np.triu_indices(n, k=1)
        for i, j in zip(iu, ju):
            real = a[j, i]
            imag = a[i, j]
            m[i, j] = real + 1j * imag
            m[j, i] = real - 1j * imag
        return m


class Spectra(EdiComponentBase):
    """
    >=SPECTRASECT + multiple >SPECTRA blocks

    Section header (>=SPECTRASECT)
    ------------------------------
    sectid: Optional[str]
    nchan : int (>=1)
    nfreq : int (>=1)
    maxblks: Optional[int]
    chksum : Optional[float]
    chan_ids: List[str]  # ordered measurement IDs (length nchan)

    Blocks
    ------
    blocks: List[SpectraBlock]

    API
    ---
    - from_edi(edi_path: str, index: int = 0) -> Spectra
    - from_edi_all(edi_path: str) -> List[Spectra]
    - read(lines_for_one_section: List[str]) -> Spectra
    - write() -> List[str]  # returns lines with >=SPECTRASECT, IDs, and >SPECTRA blocks
    """

    # keys recognized in >=SPECTRASECT
    SECT_KEYS = ("sectid", "nchan", "nfreq", "maxblks", "chksum")

    def __init__(
        self,
        *,
        sectid: Optional[str] = None,
        nchan: Optional[int] = None,
        nfreq: Optional[int] = None,
        maxblks: Optional[int] = None,
        chksum: Optional[float] = None,
        chan_ids: Optional[List[str]] = None,
        blocks: Optional[List[SpectraBlock]] = None,
    ) -> None:
        super().__init__()
        self.sectid = sectid
        self.nchan = None if nchan is None else int(nchan)
        self.nfreq = None if nfreq is None else int(nfreq)
        self.maxblks = None if maxblks is None else int(maxblks)
        self.chksum = None if chksum is None else float(chksum)
        self.chan_ids = list(chan_ids) if chan_ids is not None else []
        self.blocks: List[SpectraBlock] = list(blocks) if blocks is not None else []

        # internal bookkeeping
        self._raw_section_lines: Optional[List[str]] = None


    @classmethod
    def from_edi_all(cls, edi_path: str) -> List["Spectra"]:
        """Parse all >=SPECTRASECT sections from an EDI file."""
        sections = cls._extract_all_section_blocks(edi_path)
        out: List[Spectra] = []
        for sec_lines in sections:
            inst = cls().read(sec_lines)
            out.append(inst)
        return out

    @classmethod
    def from_edi(cls, edi_path: str, index: int = 0) -> "Spectra":
        """Parse the `index`-th >=SPECTRASECT from an EDI file (0-based)."""
        sections = cls._extract_all_section_blocks(edi_path)
        if not sections:
            raise EdIDataError("No >=SPECTRASECT section found in EDI.")
        if index < 0 or index >= len(sections):
            raise IndexError(
                f"SPECTRASECT index out of range (0..{len(sections)-1}).")
        return cls().read(sections[index])

    @staticmethod
    def _extract_all_section_blocks(edi_path: str) -> List[List[str]]:
        """Return a list of 'one-section' line 
        lists (>=SPECTRASECT ... until next >= or EOF)."""
        if not os.path.isfile(edi_path):
            raise FileNotFoundError(f"{edi_path!r} is not a file.")

        with open(edi_path, "r", encoding="utf-8") as f:
            lines = [ln.rstrip("\n") for ln in f]

        start_idxs: List[int] = []
        for i, ln in enumerate(lines):
            if ln.upper().startswith(">=SPECTRASECT"):
                start_idxs.append(i)

        sections: List[List[str]] = []
        for k, si in enumerate(start_idxs):
            # section extends from si to next >=... or EOF
            sj = len(lines)
            for j in range(si + 1, len(lines)):
                if j in start_idxs:
                    sj = j
                    break
            sections.append(lines[si:sj])
        return sections

    # Parsing a single >=SPECTRASECT section (header + IDs + >SPECTRA blocks)
    def read(self, lines_for_section: List[str]) -> "Spectra":
        """
        Read one >=SPECTRASECT block and its >SPECTRA blocks.

        Parameters
        ----------
        lines_for_section : list of lines starting at '>=SPECTRASECT' up to
                            (but not including) the next '>=...' or EOF.
        """
        if not lines_for_section or not lines_for_section[0].upper(
                ).startswith(">=SPECTRASECT"):
            raise EdIDataError(
                "SPECTRASECT parsing requires a section"
                " beginning with '>=SPECTRASECT'.")

        self._raw_section_lines = list(lines_for_section)

        # 1) Parse header key-values (lines after the >=SPECTRASECT
        # line, until dataset '//N' or IDs start)
        i = 1
        while i < len(lines_for_section):
            lnu = lines_for_section[i].upper().lstrip()
            if lnu.startswith("//"):
                # channel count line (optional)
                break
            if lnu.startswith(">SPECTRA"):
                # means no channel-ID dataset; not compliant but tolerate
                break
            if lnu.startswith(">="):
                # next section already, stop
                break

            # key=val pairs possibly spaced
            if "=" in lines_for_section[i]:
                k, v = lines_for_section[i].split("=", 1)
                key, val = _norm(k).lower(), _norm(v)
                if key in self.SECT_KEYS:
                    if key in ("nchan", "nfreq", "maxblks"):
                        setattr(self, key, _to_int(val))
                    elif key == "chksum":
                        self.chksum = _to_float(val)
                    else:
                        setattr(self, key, val)
            i += 1

            # stop header if next line is IDs count marker or >SPECTRA
            if i < len(lines_for_section):
                nxt = lines_for_section[i].lstrip()
                if nxt.startswith("//") or nxt.upper().startswith(">SPECTRA"):
                    break

        if self.nchan is None or self.nfreq is None:
            raise EdIDataError("SPECTRASECT header must define NCHAN and NFREQ.")

        # 2) Parse optional //NCHAN line (count line)
        chan_count_from_slash: Optional[int] = None
        if i < len(lines_for_section) and lines_for_section[i].lstrip().startswith("//"):
            try:
                chan_count_from_slash = int(
                    lines_for_section[i].lstrip()[2:].strip())
            except Exception:
                chan_count_from_slash = None
            i += 1

        # 3) Read measurement IDs dataset: exactly NCHAN IDs (possibly split across lines)
        ids: List[str] = []
        while i < len(lines_for_section):
            ln = lines_for_section[i].strip()
            if not ln:
                i += 1
                continue
            u = ln.upper()
            if u.startswith(">SPECTRA") or u.startswith(">="):
                break
            # split by whitespace; these are IDs like 111.001
            ids.extend([_norm(tok) for tok in ln.split()])
            if len(ids) >= (self.nchan or 0):
                ids = ids[: self.nchan]
                i += 1
                break
            i += 1

        if len(ids) != self.nchan:
            raise EdIDataError(
                f"SPECTRASECT: expected {self.nchan}"
                f" channel IDs, got {len(ids)}.")
        if chan_count_from_slash is not None and chan_count_from_slash != self.nchan:
            logger.warning(
                "SPECTRASECT '//N' (%s) does not match NCHAN (%s). Proceeding.",
                chan_count_from_slash, self.nchan
            )
        self.chan_ids = ids

        # 4) Parse >SPECTRA blocks within this section
        self.blocks.clear()
        while i < len(lines_for_section):
            ln = lines_for_section[i]
            lnu = ln.upper().lstrip()

            if lnu.startswith(">="):  # next section
                break

            if not lnu.startswith(">SPECTRA"):
                i += 1
                continue

            # Parse header options in the >SPECTRA line
            # Example: >SPECTRA SECTID=DEMO-11 FREQ=144 BW=36 AVGF=12 AVGT=128 //25
            spec_header = ln
            header_pairs: Dict[str, str] = {}
            # Extract '//' trailing count, if present (not strictly needed)
            trailing_count: Optional[int] = None
            if "//" in spec_header:
                try:
                    trailing_count = int(spec_header.split("//", 1)[1].strip())
                except Exception:
                    trailing_count = None # noqa
                spec_header = spec_header.split("//", 1)[0]

            # Find key=val pairs
            for m in re.finditer(r"([A-Za-z]+)\s*=\s*([^ \t]+)", spec_header):
                header_pairs[m.group(1).lower()] = _norm(m.group(2))

            # Resolve header values with fallback
            sectid = header_pairs.get("sectid", self.sectid)
            freq = _to_float(header_pairs.get("freq"))
            bw = _to_float(header_pairs.get("bw"))
            if freq is None or bw is None:
                raise EdIDataError(">SPECTRA block requires FREQ and BW.")

            nchan_here = _to_int(header_pairs.get("nchan")) or self.nchan
            rotspec = _to_float(header_pairs.get("rotspec"))
            avgt = _to_float(header_pairs.get("avgt")) or 1.0
            avgf = _to_float(header_pairs.get("avgf")) or 1.0
            band = header_pairs.get("band")
            segnum = _to_int(header_pairs.get("segnum")) or 0
            chksum = _to_float(header_pairs.get("chksum"))

            # Read dataset: exactly nchan*nchan real values (packed NxN real array)
            i += 1
            needed = nchan_here * nchan_here
            vals: List[float] = []
            while i < len(lines_for_section) and len(vals) < needed:
                l = lines_for_section[i].strip()
                if not l:
                    i += 1
                    continue
                u2 = l.upper()
                if u2.startswith(">=") or u2.startswith(">SPECTRA"):
                    # Unexpected early end; stop to avoid infinite loop
                    break
                # collect all floats on the line
                for tok in l.split():
                    try:
                        vals.append(float(tok))
                    except Exception:
                        pass
                i += 1

            if len(vals) < needed:
                raise EdIDataError(
                    ">SPECTRA data incomplete: expected"
                    f" {needed} values, got {len(vals)}."
                )

            a = np.array(vals[:needed], dtype=np.float64).reshape(
                nchan_here, nchan_here)
            m = SpectraBlock.unpack_from_compressed_real(a)

            block = SpectraBlock(
                sectid=sectid,
                freq=freq,
                bw=bw,
                nchan=nchan_here,
                rotspec=rotspec,
                avgt=avgt,
                avgf=avgf,
                band=band,
                segnum=segnum,
                chksum=chksum,
                matrix=m,
            )
            self.blocks.append(block)

        return self

    def write(self, numbers_per_line: int = 6) -> List[str]:
        """
        Serialize one >=SPECTRASECT section 
        (header + IDs + all >SPECTRA blocks).

        Parameters
        ----------
        numbers_per_line : int
            How many numeric values per row when writing datasets.
        """
        if self.nchan is None or self.nfreq is None:
            raise EdIDataError("Cannot write SPECTRA without NCHAN and NFREQ.")
        if len(self.chan_ids) != self.nchan:
            raise EdIDataError(
                "Cannot write SPECTRA: chan_ids length"
                f" ({len(self.chan_ids)}) != NCHAN ({self.nchan})."
            )

        lines: List[str] = []
        lines.append(">=SPECTRASECT\n")
        # Header
        def _emit_kv(k: str, v: Union[str, int, float, None]) -> None:
            if v is None or v == "":
                return
            lines.append(f"  {k.upper()}={str(v).upper()}\n")

        _emit_kv("SECTID", self.sectid or "")
        _emit_kv("NCHAN", self.nchan)
        _emit_kv("NFREQ", self.nfreq)
        _emit_kv("MAXBLKS", self.maxblks)
        _emit_kv("CHKSUM", self.chksum)

        # IDs count marker
        lines.append(f"  //{self.nchan}\n")
        # Write IDs (as provided; 6 per line default)
        cur: List[str] = []
        for cid in self.chan_ids:
            cur.append(cid)
            if len(cur) == numbers_per_line:
                lines.append("  " + "\t".join(cur) + "\n")
                cur.clear()
        if cur:
            lines.append("  " + "\t".join(cur) + "\n")

        # Write blocks
        for blk in self.blocks:
            # header line
            header_items: List[str] = [">SPECTRA"]
            if blk.sectid:
                header_items.append(f"SECTID={blk.sectid}")
            header_items.append(f"FREQ={_fmt_num(blk.freq)}")
            header_items.append(f"BW={_fmt_num(blk.bw)}")
            if blk.nchan != self.nchan:
                header_items.append(f"NCHAN={blk.nchan}")
            if blk.rotspec is not None:
                header_items.append(f"ROTSPEC={_fmt_num(blk.rotspec)}")
            if blk.avgf is not None:
                header_items.append(f"AVGF={_fmt_num(blk.avgf)}")
            if blk.avgt is not None:
                header_items.append(f"AVGT={_fmt_num(blk.avgt)}")
            if blk.band:
                header_items.append(f"BAND={blk.band}")
            if blk.segnum is not None:
                header_items.append(f"SEGNUM={blk.segnum}")
            if blk.chksum is not None:
                header_items.append(f"CHKSUM={_fmt_num(blk.chksum)}")

            # trailing count is exactly nchan*nchan
            header = " ".join(header_items) + f"  //{blk.nchan * blk.nchan}\n"
            lines.append(header)

            # Pack matrix to NxN real array, emit in row-major
            a = SpectraBlock.pack_to_compressed_real(blk.matrix)
            flat = a.reshape(-1)
            chunk: List[str] = []
            for val in flat:
                chunk.append(_fmt_num(val))
                if len(chunk) == numbers_per_line:
                    lines.append("  " + "\t".join(chunk) + "\n")
                    chunk.clear()
            if chunk:
                lines.append("  " + "\t".join(chunk) + "\n")

        lines.append("\n")
        return lines


def _fmt_num(x: Union[int, float]) -> str:
    """
    Format numbers similar to common EDI examples (scientific for floats).
    Keep integers as int-like strings.
    """
    if isinstance(x, int) or (isinstance(x, float) and float(x).is_integer()):
        return f"{int(x)}"
    # 8 significant digits in E-format
    return f"{float(x):.8E}"
