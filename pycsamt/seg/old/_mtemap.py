# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
from __future__ import annotations

from typing import List, Optional, Dict, Any
import os
import re

from ..log.logger import get_logger
from ..exceptions import EdIDataError
from .base import EDIComponentBase
from .utils import gather_measurement_key_value_with_str_parser

logger = get_logger(__name__)


def _strip_norm(s: str) -> str:
    return s.strip().strip('"').strip("'")


def _to_int_or_none(v: Any) -> Optional[int]:
    if v in (None, "", "None"):
        return None
    try:
        return int(float(v))
    except Exception:
        return None


class MTEMAP(EDIComponentBase):
    """
    MT/EMAP section header

    >=MTSECT
      SECTID=
      NFREQ=
      HX=
      HY=
      HZ=
      EX=
      EY=
      RX=
      RY=
      CHKSUM=

    >=EMAPSECT
      SECTID=
      NCHAN=
      MAXBLKS=
      NDIPOLE=
      NFREQ=
      HX=
      HY=
      (HZ=, CHKSUM= may appear in some variants)
      TYPE=

    Notes
    -----
    - We keep *measurement IDs* (hx, hy, hz, ex, ey, rx, ry) as strings to
      preserve formats like "10011.001".
    - Integers for counts: nfreq, maxblks, ndipole, chksum.
    - Presence of `ndipole` or `type` marks this as EMAP for writing.
    """

    KEY_ORDER: List[str] = [
        "sectid",
        "nfreq",
        "maxblks",
        "hx",
        "hy",
        "hz",
        "ex",
        "ey",
        "rx",
        "ry",
        "ndipole",
        "type",
        "chksum",
    ]

    # Parsed attributes
    sectid: Optional[str]
    nfreq: Optional[int]
    maxblks: Optional[int]
    hx: Optional[str]
    hy: Optional[str]
    hz: Optional[str]
    ex: Optional[str]
    ey: Optional[str]
    rx: Optional[str]
    ry: Optional[str]
    ndipole: Optional[int]
    type: Optional[str]
    chksum: Optional[int]

    # Carry-over info from legacy API
    start_data_lines_num: Optional[int] = None
    temp_sectid: Optional[str] = None

    def __init__(
            self, mt_or_emap_section_list: Optional[List[str]] = None, 
            verbose: int=0, logger =None, 
            **kwargs):
        super().__init__(verbose=verbose , logger=None)
        # defaults
        self.sectid = None
        self.nfreq = None
        self.maxblks = None
        self.hx = None
        self.hy = None
        self.hz = None
        self.ex = None
        self.ey = None
        self.rx = None
        self.ry = None
        self.ndipole = None
        self.type = None
        self.chksum = None

        # hold raw lines if provided
        self._section_lines: Optional[List[str]] = mt_or_emap_section_list

        # kwargs override
        for k, v in kwargs.items():
            setattr(self, k, v)

        # parse immediately if lines were given
        if self._section_lines is not None:
            self.read(self._section_lines)

    # Classmethod: extract section lines from EDI
    @classmethod
    def from_edi(cls, edi_path: str) -> "MTEMAP":
        """
        Extract >=MTSECT or >=EMAPSECT block from an EDI file
        and return an MTEMAP instance populated from it.
        """
        if not os.path.isfile(edi_path):
            raise FileNotFoundError(f"{edi_path!r} is not a file.")

        with open(edi_path, "r", encoding="utf-8") as f:
            lines = f.readlines()

        # capture DATAID early (used as fallback sectid)
        dataid = None
        for ln in lines:
            if re.match(r"^\s*DATAID\s*=", ln.upper()):
                try:
                    dataid = _strip_norm(ln.split("=", 1)[1])
                except Exception:
                    pass
                break

        # locate >=MTSECT or >=EMAPSECT
        start_idx = None
        for i, ln in enumerate(lines):
            u = ln.upper()
            if u.startswith(">=MTSECT") or u.startswith(">=EMAPSECT"):
                start_idx = i
                break

        if start_idx is None:
            raise EdIDataError("No >=MTSECT or >=EMAPSECT section found in EDI.")

        # find the end of this header block: stop at first data/comment block
        # like >! , >FREQ, >ZROT, >RHOROT, next ">=..." or end-of-file
        stop_idx = len(lines)
        for j in range(start_idx + 1, len(lines)):
            u = lines[j].upper().lstrip()
            if u.startswith(">!") or u.startswith(">FREQ") or u.startswith(">ZROT") \
               or u.startswith(">RHOROT") or u.startswith(">="):
                stop_idx = j
                break

        section_lines = [ln.rstrip("\n") for ln in lines[start_idx + 1: stop_idx]]

        # Normalize/flatten tokens using the shared parser util
        tokens = gather_measurement_key_value_with_str_parser(section_lines)

        inst = cls(mt_or_emap_section_list=tokens)
        inst.start_data_lines_num = stop_idx
        inst.temp_sectid = dataid
        return inst

    # ---------------------------------------------------------------------
    # Instance parsing / writing
    # ---------------------------------------------------------------------
    def read(self, mt_or_emap_section_list: Optional[List[str]] = None
             ) -> "MTEMAP":
        """
        Read tokens like ["SECTID=...", "NFREQ=...", "HX=..."] and populate attributes.
        """
        if mt_or_emap_section_list is not None:
            self._section_lines = mt_or_emap_section_list

        if not self._section_lines:
            raise EdIDataError("No MT/EMAP section lines to read.")

        # parse
        for raw in self._section_lines:
            if "=" not in raw:
                continue
            k, v = raw.split("=", 1)
            key = _strip_norm(k).lower()
            val = _strip_norm(v)

            if key == "sectid":
                self.sectid = val
            elif key == "nfreq":
                self.nfreq = _to_int_or_none(val)
            elif key == "maxblks":
                self.maxblks = _to_int_or_none(val)
            elif key in ("hx", "hy", "hz", "ex", "ey", "rx", "ry"):
                setattr(self, key, val if val not in ("", "NONE") else None)
            elif key == "ndipole":
                self.ndipole = _to_int_or_none(val)
            elif key == "type":
                self.type = val
            elif key == "chksum":
                self.chksum = _to_int_or_none(val)

        # Fallback: if sectid looks numeric or blank, use DATAID if captured
        if (self.sectid is None or self._looks_numeric(self.sectid)) and self.temp_sectid:
            self.sectid = self.temp_sectid

        return self

    def write(self, nfreq: Optional[int] = None) -> List[str]:
        """
        Serialize the section as a list of lines, including the opening tag
        (>=MTSECT or >=EMAPSECT) based on presence of EMAP fields (ndipole/type).

        Parameters
        ----------
        nfreq : Optional[int]
            If provided, overrides NFREQ in output.
        """
        is_emap = (self.ndipole is not None) or (self.type not in (None, ""))

        header = ">=EMAPSECT\n" if is_emap else ">=MTSECT\n"
        lines: List[str] = [header]

        # prepare values
        values: Dict[str, Any] = {
            "sectid": self.sectid,
            "nfreq": nfreq if nfreq is not None else self.nfreq,
            "maxblks": self.maxblks,
            "hx": self.hx,
            "hy": self.hy,
            "hz": self.hz,
            "ex": self.ex,
            "ey": self.ey,
            "rx": self.rx,
            "ry": self.ry,
            "ndipole": self.ndipole,
            "type": self.type,
            "chksum": self.chksum,
        }

        # EMAP: prior behavior removed EX/EY/HZ lines entirely
        omit_for_emap = {"ex", "ey", "hz"} if is_emap else set()

        for key in self.KEY_ORDER:
            if key in omit_for_emap:
                continue
            val = values.get(key, None)
            if val in (None, "", "None"):
                continue

            # Width / alignment similar to v1 for readability
            if key in ("hx", "hy", "hz", "ex", "ey", "rx", "ry"):
                fmt = "{:>12}"
            elif key in ("sectid", "nfreq"):
                fmt = "{:>7}"
            else:
                fmt = "{:>7}"

            lines.append(f"  {key.upper()}={fmt.format(str(val).upper())}\n")

        return lines

    # Helpers
    @staticmethod
    def _looks_numeric(s: str) -> bool:
        try:
            float(s)
            return True
        except Exception:
            return False

