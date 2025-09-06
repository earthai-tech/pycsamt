# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
from __future__ import annotations

from pathlib import Path
from typing import List, Dict, Optional

from .base import EdiComponentBase
from .he_meas import Hmeasurement, Emeasurement
from .utils import (
    strip_items,
    gather_measurement_key_value_with_str_parser,
    minimum_parser_to_write_edi,
)
from ..gis.utils import ( 
    assert_lat_value, 
    assert_lon_value, 
    assert_elevation_value, 
    convert_position_float2str, 
 )
from ..exceptions import (
    FileHandlingError,
    EdIDataError,
    HeaderError,
)
from .properties import IsEdi

from .heads import _KV_RE, _unquote, _norm_key, _is_tag, logger


class DefineMeasurement(EdiComponentBase):# In new structure, it is called EDIComponentBase
    r"""
    ``>=DEFINEMEAS`` section: sensor locations and per-channel
    parameters.

    Parameters
    ----------
    defineMeas_list : sequence of str, optional
        Raw lines of the ``>=DEFINEMEAS`` block (without the tag).
        If provided, they are parsed immediately.
    **kwargs
        Attribute overrides.

    Attributes
    ----------
    maxchan : int or None
        Maximum number of channels defined.
    maxrun : int
        Maximum number of runs (default 999).
    maxmeas : int
        Maximum number of measurements (default 9999).
    units : str
        Length unit (default ``"M"``).
    reflat : float or None
        Reference latitude (decimal degrees).
    reflong : float or None
        Reference longitude (decimal degrees).
    refelev : float or None
        Reference elevation.
    reftype : str
        Reference coordinate system (default ``"CART"``).
    define_measurement : list
        Normalized content of the block, a mix of ``KEY=VAL`` lines
        and dicts for channel lines (HMEAS/EMEAS).
    nchan : int
        Number of channels parsed.
    meas_ex, meas_ey, meas_hx, meas_hy, meas_hz : optional
        Populated `Emeasurement` / `Hmeasurement` instances per
        channel type when available.

    Notes
    -----
    - The parser accepts both DMS strings (e.g. ``26:35:11.0``) and
      decimal degrees for ``REFLAT`` / ``REFLONG``. Values are stored
      internally in decimal degrees and written back in DMS format.

    Examples
    --------
    >>> dm = DefineMeasurement.from_edi("E01.edi")
    >>> dm.reflat, dm.reflong, dm.refelev
    (26.586388..., 112.61115..., 78.0)
    >>> dm.meas_ex.id  # doctest: +SKIP
    10014.001
    """

    definemeasurementkeys = [
        "maxchan",
        "maxrun",
        "maxmeas",
        "units",
        "reftype",
        "reflat",
        "reflong",
        "refelev",
    ]

    definemeasurement_comment = (
        ">!***CHANNELS USING ORIGINAL SITE LAYOUT. "
        "FOR ROTATIONS SEE ZROT***!"
    )

    def __init__(self, defineMeas_list: Optional[List[str]] = None, **kwargs):
        super().__init__()
        self.define_measurement: Optional[List] = (
            list(defineMeas_list) if defineMeas_list else None
        )

        self.maxchan: Optional[int] = None
        self.maxrun: int = 999
        self.maxmeas: int = 9999
        self.units: str = "M"

        self._reflat: Optional[float] = None
        self._reflong: Optional[float] = None
        self._refelev: Optional[float] = None
        self.reftype: str = "CART"

        self.nchan: int = 0

        # Optional per-channel objects (created during read)
        self.meas_ex = None
        self.meas_ey = None
        self.meas_hx = None
        self.meas_hy = None
        self.meas_hz = None

        for k, v in kwargs.items():
            setattr(self, k, v)

        if self.define_measurement is not None:
            self.read(self.define_measurement)

    # ----------------------------
    # Reference coordinates props
    # ----------------------------
    @property
    def reflat(self) -> Optional[float]:
        return self._reflat

    @reflat.setter
    def reflat(self, val):
        self._reflat = None if val in (
            "", None) else assert_lat_value(val)

    @property
    def reflong(self) -> Optional[float]:
        return self._reflong

    @reflong.setter
    def reflong(self, val):
        self._reflong = None if val in (
            "", None) else assert_lon_value(val)

    @property
    def refelev(self) -> Optional[float]:
        return self._refelev

    @refelev.setter
    def refelev(self, val):
        self._refelev = None if val in (
            "", None) else assert_elevation_value(val)

    # ---------------
    # Class factories
    # ---------------
    @classmethod
    def from_edi(cls, edi_fn) -> "DefineMeasurement":
        """
        Extract and parse the ``>=DEFINEMEAS`` block from an EDI file.

        Parameters
        ----------
        edi_fn : str or Path-like
            Path to the EDI file.

        Returns
        -------
        DefineMeasurement
            Populated instance.
        """
        if edi_fn is None:
            raise FileHandlingError("No EDI path provided.")



        p = Path(edi_fn)
        IsEdi._assert_edi(p, deep=True)

        lines = p.read_text(encoding="utf-8-sig", errors="replace").splitlines()

        start, stop = None, None
        for i, ln in enumerate(lines):
            if _is_tag(ln, ">=DEFINEMEAS"):
                start = i
                break
        if start is None:
            raise EdIDataError("DEFINEMEAS block not found in EDI file.")

        # Stop at next section (>=MTSECT / >=EMAPSECT or any top-level tag)
        for j in range(start + 1, len(lines)):
            if _is_tag(lines[j], ">=MTSECT") or _is_tag(lines[j], ">=EMAPSECT"):
                stop = j
                break
            if lines[j].lstrip().startswith(">="):
                stop = j
                break

        if stop is None:
            stop = start + 1
            while stop < len(lines) and not lines[stop].lstrip().startswith(">"):
                stop += 1

        slice_ = [ln.rstrip("\n") for ln in lines[start + 1 : stop]]
        return cls(defineMeas_list=slice_)

    # ------------
    # Read content
    # ------------
    def read(self, define_measurement_list: Optional[List[str]] = None
             ) -> "DefineMeasurement":
        """
        Parse a ``>=DEFINEMEAS`` lines list and populate attributes.

        Parameters
        ----------
        define_measurement_list : sequence of str, optional
            Lines of the block (without the tag). If omitted, uses
            ``self.define_measurement``.

        Returns
        -------
        DefineMeasurement
            The instance (for chaining).
        """
        if define_measurement_list is not None:
            self.define_measurement = list(define_measurement_list)

        if not self.define_measurement:
            raise HeaderError("No DEFINEMEAS items to read.")

        items_out: List = []
        self.nchan = 0

        for raw in self.define_measurement:
            line = raw.strip()
            if not line:
                continue

            # Header KV (MAXCHAN=..., UNITS=..., REFLAT=..., ...)
            m = _KV_RE.match(line)
            if m and not line.upper().startswith((">HMEAS", ">EMEAS")):
                key = _norm_key(m.group("key"))
                val = _unquote(m.group("val"))
                key, = strip_items(key)
                val, = strip_items(val)

                # normalize synonyms
                if key == "reflon":
                    key = "reflong"

                if key not in self.definemeasurementkeys:
                    # keep unknowns as-is
                    items_out.append(f"  {key.upper()}={val}\n")
                    continue

                if key == "reflat":
                    self.reflat = val
                    items_out.append(
                        f"  REFLAT={convert_position_float2str(self.reflat)}\n"
                    )
                elif key == "reflong":
                    self.reflong = val
                    items_out.append(
                        f"  REFLONG={convert_position_float2str(self.reflong)}\n"
                    )
                elif key == "refelev":
                    self.refelev = val
                    items_out.append(f"  REFELEV={self.refelev}\n")
                elif key in {"maxchan", "maxrun", "maxmeas"}:
                    try:
                        setattr(self, key, int(val))
                    except Exception:
                        setattr(self, key, None)
                    finally:
                        outv = "" if getattr(self, key) is None else getattr(self, key)
                        items_out.append(f"  {key.upper()}={outv}\n")
                else:
                    setattr(self, key, val.upper() if isinstance(val, str) else val)
                    items_out.append(f"  {key.upper()}={val}\n")
                continue

            # Channel lines (HMEAS / EMEAS)
            up = line.upper()
            if up.startswith(">HMEAS") or up.startswith(">EMEAS"):
                self.nchan += 1
                tokens = line.split()
                tokens = gather_measurement_key_value_with_str_parser(tokens)
                mdict: Dict[str, str] = {}
                for tok in tokens:
                    if "=" in tok:
                        k, v = tok.split("=", 1)
                        k, = strip_items(k.lower())
                        v, = strip_items(v)
                        mdict[k] = v

                # keep the raw dict in the output list
                items_out.append(mdict)

                # instantiate measurement objects + attach
                chtype = mdict.get("chtype", "").lower()
                if chtype.startswith("h"):
                    obj = Hmeasurement(**mdict)
                    setattr(self, f"meas_{chtype}", obj)
                elif chtype.startswith("e"):
                    obj = Emeasurement(**mdict)
                    setattr(self, f"meas_{chtype}", obj)

                continue

            # comments or unknown lines: keep as-is
            items_out.append(line + ("\n" if not line.endswith("\n") else ""))

        logger.info('Found %d channel line(s) in DEFINEMEAS.', self.nchan)
        self.define_measurement = items_out
        
        return self

    # -------------
    # Write content
    # -------------
    def write(self, define_measurement_list: Optional[List] = None) -> List[str]:
        """
        Serialize the component to formatted lines.

        Parameters
        ----------
        define_measurement_list : list, optional
            Existing mixed list to normalize/rewrite. When omitted,
            the current object state is serialized.

        Returns
        -------
        list of str
            Lines including the leading ``>=DEFINEMEAS`` tag and a
            trailing blank line.
        """
        out: List[str] = [">=DEFINEMEAS\n"]

        # If an external list is provided, normalize and return
        if define_measurement_list is not None:
            if minimum_parser_to_write_edi(define_measurement_list) is None:
                raise EdIDataError(
                    "Malformed DEFINEMEAS content; expected a mix of "
                    "header 'KEY=VAL' lines and HMEAS/EMEAS dicts."
                )

            header_buf: List[str] = []
            chan_dicts: List[Dict[str, str]] = []

            for item in define_measurement_list:
                if isinstance(item, dict):
                    chan_dicts.append(item)
                else:
                    s = item.strip()
                    if not s or s.startswith(">"):
                        continue
                    m = _KV_RE.match(s)
                    if not m:
                        continue
                    k = _norm_key(m.group("key"))
                    v = _unquote(m.group("val"))
                    k, = strip_items(k)
                    v, = strip_items(v)
                    header_buf.append(f"  {k.upper()}={v}\n")

            out.extend(header_buf)
            out.append(self.definemeasurement_comment.upper() + "\n")

            for d in chan_dicts:
                tag = ">HMEAS" if d.get("chtype", "").lower().startswith("h") else ">EMEAS"
                out.append(tag + " ")
                for k, v in d.items():
                    vv = "" if v in ("", None, "None") else str(v).upper()
                    out.append(f"{k.upper()}={vv} ")
                out.append("\n")

            out.append("\n")
            return out

        # Otherwise, serialize from object state
        for k in self.definemeasurementkeys:
            val = getattr(self, k, None)
            if k == "reflat" and val is not None:
                val = convert_position_float2str(val)
            elif k == "reflong" and val is not None:
                val = convert_position_float2str(val)

            sval = "" if val in (None, "") else str(val).upper()
            out.append(f"  {k.upper()}={sval}\n")

        out.append(self.definemeasurement_comment.upper() + "\n")

        # Emit channel objects if present
        for chtype in ("ex", "ey", "hx", "hy", "hz"):
            obj = getattr(self, f"meas_{chtype}", None)
            if obj is None:
                continue
            is_e = chtype.startswith("e")
            out.append(">EMEAS " if is_e else ">HMEAS ")

            # order by class public keys, fallback to dict order
            keys = (
                getattr(Emeasurement, "emeasurementkey", None)
                if is_e
                else getattr(Hmeasurement, "hmeasurementkey", None)
            )
            if keys is None:
                # dump object's dict deterministically
                keys = sorted([k for k in obj.__dict__.keys() if not k.startswith("_")])

            for key in keys:
                val = getattr(obj, key, "")
                val = "" if val in (None, "") else str(val).upper()
                out.append(f"{key.upper()}={val} ")
            out.append("\n")

        out.append("\n")
        return out
