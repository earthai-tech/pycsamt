# -*- coding: utf-8 -*-
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


# ------------------------- small helpers -------------------------

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


def _fmt_num(x: Union[int, float]) -> str:
    """Format like typical EDI numeric output."""
    if isinstance(x, int) or (isinstance(x, float) and float(x).is_integer()):
        return f"{int(x)}"
    return f"{float(x):.8E}"


# ------------------------- data blocks -------------------------

class TimeSeriesBlock(EdiComponentBase):
    """
    One >TSERIES block.

    Attributes
    ----------
    sectid : Optional[str]
    nchan  : int              # Number of channels; default from section if omitted
    npts   : int              # REQUIRED
    sr     : float            # Sample rate (Hz), REQUIRED
    mpx    : str              # "CHAN" (default) or "TIME"
    band   : Optional[str]
    chksum : Optional[float]
    data   : np.ndarray, shape (nchan, npts), dtype float64
    """

    def __init__(
        self,
        *,
        nchan: int,
        npts: int,
        sr: float,
        data: np.ndarray,
        sectid: Optional[str] = None,
        mpx: Optional[str] = "CHAN",
        band: Optional[str] = None,
        chksum: Optional[float] = None,
    ) -> None:
        super().__init__()
        self.sectid = sectid
        self.nchan = int(nchan)
        self.npts = int(npts)
        self.sr = float(sr)
        self.mpx = "CHAN" if (mpx is None) else str(mpx).upper()
        if self.mpx not in ("CHAN", "TIME"):
            raise EdIDataError("TSERIES MPX must be 'CHAN' or 'TIME'.")
        self.band = band
        self.chksum = None if chksum is None else float(chksum)

        if not isinstance(data, np.ndarray):
            raise EdIDataError("TimeSeriesBlock.data must be a numpy ndarray.")
        if data.shape != (self.nchan, self.npts):
            raise EdIDataError(
                f"TimeSeriesBlock.data shape mismatch:"
                f" expected ({self.nchan},{self.npts}), "
                f"got {data.shape}"
            )
        self.data = data.astype(np.float64, copy=False)

    # ---- flatten / unflatten helpers for MPX ----

    def flatten(self) -> np.ndarray:
        """
        Flatten data according to self.mpx into 1D vector of length nchan*npts.
        - MPX=CHAN: [ch0 series][ch1 series]...
        - MPX=TIME: [t0 of all chans][t1 of all chans]...
        """
        if self.mpx == "CHAN":
            return self.data.reshape(-1)  # (nchan * npts,)
        # TIME: interleave by time
        return self.data.T.reshape(-1)  # (npts * nchan,)

    @staticmethod
    def unflatten(vec: np.ndarray, nchan: int, npts: int, mpx: str
                  ) -> np.ndarray:
        """
        Inverse of flatten() for given MPX.
        Returns shape (nchan, npts).
        """
        if vec.size != nchan * npts:
            raise EdIDataError("TSERIES data length doesn't match NCHAN*NPTS.")
        mpx = str(mpx).upper()
        if mpx == "CHAN":
            return vec.reshape(nchan, npts)
        elif mpx == "TIME":
            return vec.reshape(npts, nchan).T
        else:
            raise EdIDataError("TSERIES MPX must be 'CHAN' or 'TIME'.")


class TimeSeries(EdiComponentBase):
    """
    >=TSERIESSECT + multiple >TSERIES blocks

    Section header (>=TSERIESSECT)
    ------------------------------
    sectid : Optional[str]
    nchan  : int (>=1)
    maxblks: Optional[int]
    chksum : Optional[float]
    chan_ids: List[str]  # ordered measurement IDs (length nchan)

    Blocks
    ------
    blocks: List[TimeSeriesBlock]

    API
    ---
    - from_edi_all(path) -> List[TimeSeries]
    - from_edi(path, index=0) -> TimeSeries
    - read(lines_for_section)
    - write(numbers_per_line=6) -> List[str]
    """

    SECT_KEYS = ("sectid", "nchan", "maxblks", "chksum")

    def __init__(
        self,
        *,
        sectid: Optional[str] = None,
        nchan: Optional[int] = None,
        maxblks: Optional[int] = None,
        chksum: Optional[float] = None,
        chan_ids: Optional[List[str]] = None,
        blocks: Optional[List[TimeSeriesBlock]] = None,
    ) -> None:
        super().__init__()
        self.sectid = sectid
        self.nchan = None if nchan is None else int(nchan)
        self.maxblks = None if maxblks is None else int(maxblks)
        self.chksum = None if chksum is None else float(chksum)
        self.chan_ids = list(chan_ids) if chan_ids is not None else []
        self.blocks: List[TimeSeriesBlock] = list(
            blocks) if blocks is not None else []

        self._raw_section_lines: Optional[List[str]] = None

    # ----------------- discovery & extraction -----------------

    @classmethod
    def from_edi_all(cls, edi_path: str) -> List["TimeSeries"]:
        sections = cls._extract_all_section_blocks(edi_path)
        out: List[TimeSeries] = []
        for sec_lines in sections:
            inst = cls().read(sec_lines)
            out.append(inst)
        return out

    @classmethod
    def from_edi(cls, edi_path: str, index: int = 0) -> "TimeSeries":
        sections = cls._extract_all_section_blocks(edi_path)
        if not sections:
            raise EdIDataError("No >=TSERIESSECT section found in EDI.")
        if index < 0 or index >= len(sections):
            raise IndexError(
                f"TSERIESSECT index out of range (0..{len(sections)-1}).")
        return cls().read(sections[index])

    @staticmethod
    def _extract_all_section_blocks(edi_path: str) -> List[List[str]]:
        """
        Return a list of 'one-section' line lists (>=TSERIESSECT ... until next >= or EOF).
        """
        if not os.path.isfile(edi_path):
            raise FileNotFoundError(f"{edi_path!r} is not a file.")

        with open(edi_path, "r", encoding="utf-8") as f:
            lines = [ln.rstrip("\n") for ln in f]

        start_idxs: List[int] = []
        for i, ln in enumerate(lines):
            if ln.upper().startswith(">=TSERIESSECT"):
                start_idxs.append(i)

        sections: List[List[str]] = []
        for si_idx, si in enumerate(start_idxs):
            sj = len(lines)
            for j in range(si + 1, len(lines)):
                if j in start_idxs:
                    sj = j
                    break
            sections.append(lines[si:sj])
        return sections

    # ----------------- parser -----------------

    def read(self, lines_for_section: List[str]) -> "TimeSeries":
        """
        Read one >=TSERIESSECT block and its >TSERIES blocks.
        """
        if not lines_for_section or not lines_for_section[0].upper(
                ).startswith(">=TSERIESSECT"):
            raise EdIDataError(
                "TSERIESSECT parsing requires a"
                " section beginning with '>=TSERIESSECT'.")

        self._raw_section_lines = list(lines_for_section)

        # Header KVs
        i = 1
        while i < len(lines_for_section):
            lnu = lines_for_section[i].upper().lstrip()
            if lnu.startswith("//") or lnu.startswith(
                    ">TSERIES") or lnu.startswith(">="):
                break

            if "=" in lines_for_section[i]:
                k, v = lines_for_section[i].split("=", 1)
                key, val = _norm(k).lower(), _norm(v)
                if key in self.SECT_KEYS:
                    if key == "nchan":
                        self.nchan = _to_int(val)
                    elif key == "maxblks":
                        self.maxblks = _to_int(val)
                    elif key == "chksum":
                        self.chksum = _to_float(val)
                    else:
                        setattr(self, key, val)
            i += 1

        if self.nchan is None:
            raise EdIDataError("TSERIESSECT header must define NCHAN.")

        # Optional //NCHAN
        chan_count_from_slash: Optional[int] = None
        if i < len(lines_for_section) and lines_for_section[i].lstrip().startswith("//"):
            try:
                chan_count_from_slash = int(lines_for_section[i].lstrip()[2:].strip())
            except Exception:
                chan_count_from_slash = None
            i += 1

        # Channel IDs dataset (exactly NCHAN)
        ids: List[str] = []
        while i < len(lines_for_section):
            ln = lines_for_section[i].strip()
            if not ln:
                i += 1
                continue
            u = ln.upper()
            if u.startswith(">TSERIES") or u.startswith(">="):
                break
            ids.extend([_norm(tok) for tok in ln.split()])
            if len(ids) >= (self.nchan or 0):
                ids = ids[: self.nchan]
                i += 1
                break
            i += 1

        if len(ids) != self.nchan:
            raise EdIDataError(
                f"TSERIESSECT: expected {self.nchan} channel IDs, got {len(ids)}.")
        if chan_count_from_slash is not None and chan_count_from_slash != self.nchan:
            logger.warning(
                "TSERIESSECT '//N' (%s) does not match NCHAN (%s). Proceeding.",
                chan_count_from_slash, self.nchan
            )
        self.chan_ids = ids

        # Parse >TSERIES blocks
        self.blocks.clear()
        while i < len(lines_for_section):
            ln = lines_for_section[i]
            lnu = ln.upper().lstrip()

            if lnu.startswith(">="):
                break

            if not lnu.startswith(">TSERIES"):
                i += 1
                continue

            # Header parse on >TSERIES line
            spec_header = ln
            trailing_count: Optional[int] = None
            if "//" in spec_header:
                try:
                    trailing_count = int(spec_header.split("//", 1)[1].strip())
                except Exception:
                    trailing_count = None # noqa
                spec_header = spec_header.split("//", 1)[0]

            header_pairs: Dict[str, str] = {}
            for m in re.finditer(r"([A-Za-z]+)\s*=\s*([^ \t]+)", spec_header):
                header_pairs[m.group(1).lower()] = _norm(m.group(2))

            sectid = header_pairs.get("sectid", self.sectid)
            nchan_here = _to_int(header_pairs.get("nchan")) or self.nchan
            npts = _to_int(header_pairs.get("npts"))
            sr = _to_float(header_pairs.get("sr"))
            mpx = header_pairs.get("mpx", "CHAN").upper()
            band = header_pairs.get("band")
            chksum = _to_float(header_pairs.get("chksum"))

            if npts is None or sr is None:
                raise EdIDataError(">TSERIES block requires NPTS and SR.")

            # Collect exactly nchan*npts values
            i += 1
            needed = nchan_here * npts
            vals: List[float] = []
            while i < len(lines_for_section) and len(vals) < needed:
                l = lines_for_section[i].strip()
                if not l:
                    i += 1
                    continue
                u2 = l.upper()
                if u2.startswith(">=") or u2.startswith(">TSERIES"):
                    # Unexpected early end; stop
                    break
                for tok in l.split():
                    try:
                        vals.append(float(tok))
                    except Exception:
                        pass
                i += 1

            if len(vals) < needed:
                raise EdIDataError(
                    f">TSERIES data incomplete: expected {needed} values, got {len(vals)}."
                )
            vec = np.array(vals[:needed], dtype=np.float64)
            data = TimeSeriesBlock.unflatten(vec, nchan_here, npts, mpx)

            block = TimeSeriesBlock(
                sectid=sectid,
                nchan=nchan_here,
                npts=npts,
                sr=sr,
                mpx=mpx,
                band=band,
                chksum=chksum,
                data=data,
            )
            self.blocks.append(block)

        return self

    def write(self, numbers_per_line: int = 6) -> List[str]:
        """
        Serialize one >=TSERIESSECT section (header + IDs + all >TSERIES blocks).
        """
        if self.nchan is None:
            raise EdIDataError("Cannot write TSERIES without NCHAN.")
        if len(self.chan_ids) != self.nchan:
            raise EdIDataError(
                f"Cannot write TSERIES: chan_ids length"
                f" ({len(self.chan_ids)}) != NCHAN ({self.nchan})."
            )

        lines: List[str] = []
        lines.append(">=TSERIESSECT\n")

        def _emit_kv(k: str, v: Union[str, int, float, None]) -> None:
            if v is None or v == "":
                return
            lines.append(f"  {k.upper()}={str(v).upper()}\n")

        _emit_kv("SECTID", self.sectid or "")
        _emit_kv("NCHAN", self.nchan)
        _emit_kv("MAXBLKS", self.maxblks)
        _emit_kv("CHKSUM", self.chksum)

        # IDs count marker and IDs themselves
        lines.append(f"  //{self.nchan}\n")
        cur: List[str] = []
        for cid in self.chan_ids:
            cur.append(cid)
            if len(cur) == numbers_per_line:
                lines.append("  " + "\t".join(cur) + "\n")
                cur.clear()
        if cur:
            lines.append("  " + "\t".join(cur) + "\n")

        # Write each >TSERIES block
        for blk in self.blocks:
            hdr: List[str] = [">TSERIES"]
            if blk.sectid:
                hdr.append(f"SECTID={blk.sectid}")
            if blk.nchan != self.nchan:
                hdr.append(f"NCHAN={blk.nchan}")
            hdr.append(f"NPTS={blk.npts}")
            hdr.append(f"SR={_fmt_num(blk.sr)}")
            if blk.mpx:
                hdr.append(f"MPX={blk.mpx}")
            if blk.band:
                hdr.append(f"BAND={blk.band}")
            if blk.chksum is not None:
                hdr.append(f"CHKSUM={_fmt_num(blk.chksum)}")

            hdr_line = " ".join(hdr) + f"  // {blk.nchan * blk.npts}\n"
            lines.append(hdr_line)

            vec = blk.flatten()  # honors MPX for output
            chunk: List[str] = []
            for val in vec:
                chunk.append(_fmt_num(val))
                if len(chunk) == numbers_per_line:
                    lines.append("  " + "\t".join(chunk) + "\n")
                    chunk.clear()
            if chunk:
                lines.append("  " + "\t".join(chunk) + "\n")

        lines.append("\n")
        return lines
