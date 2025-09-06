# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

from pathlib import Path
import datetime as _dt
import re
from typing import List, Optional, Sequence, Union

from .. import __version__ as _PKG_VERSION
from ..log.logger import get_logger
from ..z.z import Z as EMZ  # noqa: F401  
from ..site.location import Location
from ..gis.utils import dms_to_decimal, decimal_to_dms
from ..exceptions import (
    HeaderError,
    EdIDataError,
    FileHandlingError,
)
from .base import EdiComponentBase
from .properties import IsEdi, Source, Processing, Copyright


logger = get_logger(__name__)

# KEY=VALUE pairs, tolerant to spaces and optional quotes
_KV_RE = re.compile(
    r"""
    ^\s*                                   # leading spaces
    (?P<key>[A-Za-z][A-Za-z0-9_]*)         # key
    \s*=\s*                                # equals
    (?P<val>.+?)                           # value (lazy)
    \s*$                                   # trailing spaces
    """,
    re.VERBOSE,
)


def _unquote(s: str) -> str:
    """Remove one layer of surrounding quotes."""
    s = s.strip()
    if len(s) >= 2 and ((s[0] == s[-1]) and s[0] in ("'", '"')):
        return s[1:-1]
    return s


def _norm_key(k: str) -> str:
    """
    Normalize header keys to canonical names.

    SEG-EDI often mixes LON/LONG and case. Keep canonical names in lower.
    """
    k0 = k.strip().lower()
    if k0 == "lon":
        return "long"
    return k0


def _is_tag(line: str, tag: str) -> bool:
    """Case-insensitive, whitespace-tolerant check for EDI tag."""
    line = line.strip()
    if not line.startswith(">"):
        return False
    return line.upper().startswith(tag.upper())


class Head(EdiComponentBase):
    r"""
    EDI **HEAD** block.

    The HEAD section identifies the dataset, describes where/when/by whom
    data were acquired, and records how/when/by whom the EDI file was
    written (per SEG MT/EMAP standard).

    Parameters
    ----------
    edi_header_list : list of str, optional
        Raw lines belonging to the ``>HEAD`` block (up to the start of
        ``>INFO`` or the next section). If provided, they are parsed
        immediately.
    **kwargs
        Any HEAD attribute to override defaults (e.g., ``dataid``,
        ``acqby``, ...). Unknown keys are set as dynamic attributes.

    Attributes
    ----------
    Location : Location
        Location container (latitude, longitude, elevation).
    dataid, acqby, fileby : str or None
        Dataset id, acquisition contractor, file writer.
    acqdate, enddate, filedate : str or None
        Dates (as strings). ``filedate`` defaults to current UTC
        timestamp ``YYYY/MM/DD HH:MM:SS UTC``.
    country, state, county, prospect, loc : str or None
        Descriptive metadata.
    lat, long, elev : float or str
        Geographic coordinates; properties delegate into ``Location``.
        Readers accept either decimal degrees (float) or DMS strings;
        writers output DMS strings.
    units : {'m', 'ft'}
        Elevation units (defaults to 'm').
    stdvers : str
        EDI standard version; defaults to ``'SEG 1.0'``.
    progvers : str
        Program version string; defaults to ``'pyCSAMT <version>'``.
    progdate : str
        Program revision date ``YYYY/MM/DD``.
    coordsys : str
        Coordinate system description, default 'Geomagnetic North'.
    declination : float or None
        Geomagnetic declination, degrees.
    datum : str
        Geodetic datum (defaults to 'WGS84').
    maxsect : int or None
        Maximum data sections in file.
    bindata : str or None
        Optional tag for binary data file.
    project, survey : str or None
        Project/survey names.
    empty : float
        Missing-value sentinel (defaults to ``1.0e32``).
    edi_header : list of str or None
        Retained, normalized list of ``KEY=VALUE`` lines parsed from HEAD.

    Notes
    -----
    - Parsing is tolerant to case and single/double quotes.
    - For ``lat``/``long`` setters, DMS strings are converted to decimal
      degrees using the GIS helpers and stored in ``Location`` in
      decimal form. When writing, they are formatted back to DMS.

    Examples
    --------
    Read only the HEAD block from a file:

    >>> head = Head.get_header_list_from_edi("E01.edi")
    >>> head.dataid
    'E01'
    >>> head.lat, head.long, head.elev
    (.., .., ..)

    Create and write HEAD lines:

    >>> h = Head(dataid="E1_2", acqby="ULTREM", units="m")
    >>> lines = h.write_head_info()
    >>> print("".join(lines))
    """

    # Canonical keys that we know how to write in order
    head_keys: List[str] = [
        "dataid",
        "acqby",
        "fileby",
        "acqdate",
        "enddate",
        "filedate",
        "country",
        "state",
        "county",
        "prospect",
        "loc",
        "lat",
        "long",
        "elev",
        "declination",
        "datum",
        "units",
        "stdvers",
        "coordsys",
        "progvers",
        "progdate",
        "maxsect",
        "bindata",
        "project",
        "survey",
        "empty",
    ]

    # Keys that should be quoted in the EDI output
    _quoted_keys = {"dataid", "stdvers", "progvers"}

    def __init__(
        self, edi_header_list: Optional[Sequence[str]] = None, 
        **kwargs
        ):
        # Core containers
        self.Location = Location()

        # Defaults
        self.dataid: Optional[str] = None
        self.acqby: Optional[str] = None
        self.fileby: Optional[str] = None
        self.acqdate: Optional[str] = None
        self.enddate: Optional[str] = None
        self.filedate: str = _dt.datetime.utcnow(
            ).strftime("%Y/%m/%d %H:%M:%S UTC")

        self.country: Optional[str] = None
        self.state: Optional[str] = None
        self.county: Optional[str] = None
        self.prospect: Optional[str] = None
        self.loc: Optional[str] = None

        self.units: str = "m"
        self.stdvers: str = "SEG 1.0"
        self.progvers: str = f"pyCSAMT {_PKG_VERSION}"
        self.progdate: str = _dt.datetime.utcnow().strftime("%Y/%m/%d")
        self.coordsys: str = "Geomagnetic North"
        self.declination: Optional[float] = None
        self.datum: str = "WGS84"
        self.maxsect: Optional[int] = None
        self.bindata: Optional[str] = None
        self.project: Optional[str] = None
        self.survey: Optional[str] = None
        self.empty: float = 1.0e32

        self.edi_header: Optional[List[str]] = list(
            edi_header_list) if edi_header_list else None

        # Apply user overrides (including unknown keys)
        for k, v in kwargs.items():
            setattr(self, k, v)

        # Parse if provided
        if self.edi_header is not None:
            self.read_head(self.edi_header)

    # ----------------
    # Lat/Long/Elev IO
    # ----------------
    @property
    def lat(self) -> Optional[float]:
        return self.Location.latitude

    @lat.setter
    def lat(self, value: Union[str, float, int, None]) -> None:
        if value is None:
            self.Location.latitude = None
            return
        try:
            self.Location.latitude = float(value)
        except (TypeError, ValueError):
            self.Location.latitude = float(dms_to_decimal(str(value)))
            logger.info('Converted DMS string latitude to decimal degrees.')

    @property
    def long(self) -> Optional[float]:
        return self.Location.longitude

    @long.setter
    def long(self, value: Union[str, float, int, None]) -> None:
        if value is None:
            self.Location.longitude = None
            return
        try:
            self.Location.longitude = float(value)
        except (TypeError, ValueError):
            self.Location.longitude = float(dms_to_decimal(str(value)))
            logger.info('Converted DMS string longitude to decimal degrees.')

    @property
    def elev(self) -> Optional[float]:
        return self.Location.elevation

    @elev.setter
    def elev(self, value: Union[str, float, int, None]) -> None:
        if value is None or (isinstance(value, str) and not value.strip()):
            self.Location.elevation = None
        else:
            self.Location.elevation = float(value)

    # ------------------------------
    # Read HEAD block from EDI lines
    # ------------------------------
    @classmethod
    def from_edi(cls, edi_fn: Union[str, Path]) -> "Head":
        """
        Extract and parse the ``>HEAD`` block from an EDI file.

        Parameters
        ----------
        edi_fn : str or Path
            Path to the EDI file.

        Returns
        -------
        Head
            A ``Head`` instance populated from the file.

        Raises
        ------
        FileHandlingError
            If ``edi_fn`` is None or invalid.
        FileNotFoundError
            If the file does not exist.
        EdIDataError
            If file fails EDI validation or HEAD section is missing.
        """
        if edi_fn is None:
            raise FileHandlingError("No EDI path provided.")
        edi_path = Path(edi_fn)
        # Validate (deep)
        IsEdi._assert_edi(edi_path, deep=True)

        lines = edi_path.read_text(
            encoding="utf-8-sig", errors="replace").splitlines()

        # find >HEAD start and >INFO (or >=DEFINEMEAS) stop
        start, stop = None, None
        for i, ln in enumerate(lines):
            if _is_tag(ln, ">HEAD"):
                start = i
                break
        if start is None:
            raise EdIDataError("HEAD block not found in EDI file.")

        for j in range(start + 1, len(lines)):
            if _is_tag(lines[j], ">INFO") or _is_tag(lines[j], ">=DEFINEMEAS"):
                stop = j
                break
        if stop is None:
            # If INFO is genuinely missing, take contiguous lines until next tag
            stop = start + 1
            while stop < len(lines) and not lines[stop].lstrip().startswith(">"):
                stop += 1

        # Slice includes header content (excluding the >HEAD line itself)
        header_slice = [ln.strip().replace('"', "") for ln in lines[start + 1 : stop]]
        return cls(edi_header_list=header_slice)

    def read(self, edi_header_list: Optional[Sequence[str]] = None) -> "Head":
        """
        Parse a ``>HEAD`` key-value list and set attributes.

        Parameters
        ----------
        edi_header_list : sequence of str, optional
            Lines of the HEAD block (``KEY=VALUE``). If omitted, uses
            ``self.edi_header``.

        Returns
        -------
        Head
            The instance (for chaining).

        Raises
        ------
        HeaderError
            If nothing is available to read.
        """
        if edi_header_list is not None:
            self.edi_header = list(edi_header_list)

        if not self.edi_header:
            raise HeaderError("No HEAD items to read.")

        normalized: List[str] = []
        for raw in self.edi_header:
            line = raw.strip()
            if not line or line.startswith(">"):
                continue
            m = _KV_RE.match(line)
            if not m:
                # skip non KV lines quietly
                continue
            key = _norm_key(m.group("key"))
            val = _unquote(m.group("val"))

            # Set attribute, with special handling for coordsys alias
            if key == "coordsys":
                setattr(self, "coordsys", val)
            elif key == "lon":
                setattr(self, "long", val)
            else:
                setattr(self, key, val)
            normalized.append(f"{key.upper()}={val}")

        self.edi_header = normalized
        return self

    # --------------------
    # Write HEAD to lines
    # --------------------
    def write(
            self, head_list_infos: Optional[Sequence[str]] = None) -> List[str]:
        """
        Build formatted ``>HEAD`` lines ready to be written to an EDI file.

        When ``head_list_infos`` is provided, it is used as the source of
        ``KEY=VALUE`` pairs (after normalization). Otherwise, current
        attributes are serialized in the canonical order.

        Parameters
        ----------
        head_list_infos : sequence of str, optional
            Existing HEAD lines. If provided, they are normalized and
            returned; if omitted, the method serializes the current
            object state.

        Returns
        -------
        list of str
            Lines including the leading ``>HEAD`` and a trailing blank
            line for readability.
        """
        lines: List[str] = [">HEAD\n"]

        if head_list_infos is not None:
            # Normalize and rewrite
            cleaned = []
            for raw in head_list_infos:
                s = raw.strip()
                if not s or s.startswith(">"):
                    continue
                m = _KV_RE.match(s)
                if not m:
                    continue
                key = _norm_key(m.group("key"))
                val = _unquote(m.group("val"))
                # Formatting: quote some keys, uppercase right-hand when desired
                if key in self._quoted_keys:
                    out_val = f'"{val}"'
                else:
                    out_val = val.upper() if key not in {"lat", "long"} else val
                cleaned.append(f"  {key.upper()}={out_val}\n")
            lines.extend(cleaned)
            lines.append("\n")
            return lines

        # Serialize current attributes in canonical order
        for key in self.head_keys:
            val = getattr(self, key, None)
            if val in (None, "", "None"):
                continue

            # Lat/long as DMS strings if numeric
            if key in {"lat", "long"}:
                # Write as DMS text expected by many EDI tools
                vnum = getattr(self, key)
                try:
                    out_val = decimal_to_dms(float(vnum))
                except Exception:
                    out_val = str(vnum)
            elif key in {"elev"}:
                out_val = str(val)
            elif key in self._quoted_keys:
                out_val = f'"{val}"'
            else:
                out_val = str(val).upper()

            lines.append(f"  {key.upper()}={out_val}\n")

        lines.append("\n")
        return lines


class Info(EdiComponentBase):
    r"""
    EDI **INFO** block.

    Collects survey-level information and processing provenance.

    Parameters
    ----------
    edi_info_list : sequence of str, optional
        Raw lines of the ``>INFO`` block (without the ``>INFO`` tag).
        If provided, they are parsed immediately.
    **kwargs
        Any attribute overrides. Unknown keys are set as dynamic
        attributes.

    Attributes
    ----------
    maxinfo : int
        Maximum number of information lines (default 999).
    Source : Source
        Provenance of the dataset (project, survey, sitename, etc.).
    Processing : Processing
        Processing meta-information (software, processedby, etc.).
    Copyright : Copyright
        Optional copyright metadata container.
    filter : str or None
        Optional filter tag (e.g., for EMAP sections).
    ediinfo : list of str or None
        Normalized ``KEY=VALUE`` lines of the INFO block.

    Notes
    -----
    The following keys are routed to nested objects:

    - ``project``, ``survey``, ``creationdate``, ``sitename`` → ``Source``
    - ``processedby``, ``processingtag``, ``runlist``, ``remoteref``,
      ``remotesite``, ``signconvention`` → ``Processing``
    - ``processingsoftware`` → ``Processing.ProcessingSoftware.name``

    Examples
    --------
    Read the INFO block from a file:

    >>> info = Info.get_info_list_from_edi("E01.edi")
    >>> info.Source.project
    'SAMTEX'

    Create and serialize:

    >>> i = Info()
    >>> i.Source.project = "DemoProj"
    >>> i.Processing.processedby = "pyCSAMT"
    >>> print("".join(i.write()))
    """

    # canonical key serialization order
    infokeys = [
        "maxinfo",
        "project",
        "survey",
        "creationdate",
        "processedby",
        "processingsoftware",
        "processingtag",
        "sitename",
        "runlist",
        "remoteref",
        "remotesite",
        "signconvention",
    ]

    # keys to quote in output
    _quoted_keys = {"processedby", "processingsoftware"}

    def __init__(self, edi_info_list=None, **kwargs):
        super().__init__()
        self.ediinfo = list(edi_info_list) if edi_info_list else None
        self.filter = None  # used in some EMAP/processing contexts
        self.maxinfo = 999

        self.Source = Source()
        self.Processing = Processing()
        self.Copyright = Copyright()

        for k, v in kwargs.items():
            setattr(self, k, v)

        if self.ediinfo is not None:
            self.read(self.ediinfo)

    # ------------------------------
    # Read INFO block from EDI lines
    # ------------------------------
    @classmethod
    def from_edi(cls, edi_fn=None) -> "Info":
        """
        Extract and parse the ``>INFO`` block from an EDI file.

        Parameters
        ----------
        edi_fn : str or Path-like
            Path to the EDI file.

        Returns
        -------
        Info
            Populated ``Info`` instance.

        Raises
        ------
        FileHandlingError
            If path is invalid.
        FileNotFoundError
            If file does not exist.
        EdIDataError
            If EDI validation fails or INFO block is missing.
        """
        if edi_fn is None:
            raise FileHandlingError("No EDI path provided.")

        

        p = Path(edi_fn)
        IsEdi._assert_edi(p, deep=True)

        lines = p.read_text(encoding="utf-8-sig", errors="replace").splitlines()

        start, stop = None, None
        for i, ln in enumerate(lines):
            if _is_tag(ln, ">INFO"):
                start = i
                break
        if start is None:
            raise EdIDataError("INFO block not found in EDI file.")

        # Stop at next section tag (>=DEFINEMEAS, >=MTSECT, or any '>' line)
        for j in range(start + 1, len(lines)):
            if _is_tag(lines[j], ">=DEFINEMEAS") or _is_tag(lines[j], ">=MTSECT"):
                stop = j
                break
            # if any next top-level block appears
            if lines[j].lstrip().startswith(">"):
                stop = j
                break

        if stop is None:
            stop = start + 1
            while stop < len(lines) and not lines[stop].lstrip().startswith(">"):
                stop += 1

        info_slice = [ln.strip().replace('"', "") for ln in lines[start + 1 : stop]]
        return cls(edi_info_list=info_slice)

    def read(self, edi_info_list=None) -> "Info":
        """
        Parse an INFO key-value list and set attributes.

        Parameters
        ----------
        edi_info_list : sequence of str, optional
            Lines of the INFO block (``KEY=VALUE``). If omitted, uses
            ``self.ediinfo``.

        Returns
        -------
        Info
            The instance (for chaining).

        Raises
        ------
        HeaderError
            If nothing is available to read.
        """
        if edi_info_list is not None:
            self.ediinfo = list(edi_info_list)

        if not self.ediinfo:
            raise HeaderError("No INFO items to read.")

        normalized = []
        for raw in self.ediinfo:
            line = raw.strip()
            if not line or line.startswith(">"):
                continue
            # tolerate blank/flag lines (e.g., comments)
            m = _KV_RE.match(line)
            if not m:
                continue

            key = _norm_key(m.group("key"))
            val = _unquote(m.group("val"))

            # route keys to nested containers
            if key in {"project", "survey", "creationdate", "sitename"}:
                setattr(self.Source, key, val)
            elif key == "processingsoftware":
                # nested software name
                self.Processing.ProcessingSoftware.name = val
            elif key in {
                "processedby",
                "processingtag",
                "runlist",
                "remoteref",
                "remotesite",
                "signconvention",
            }:
                setattr(self.Processing, key, val)
            else:
                # accept other fields as direct attributes
                setattr(self, key, val)

            normalized.append(f"{key.upper()}={val}")

        self.ediinfo = normalized
        logger.debug("Parsed INFO: %s", self.ediinfo)
        return self

    # --------------------
    # Write INFO to lines
    # --------------------
    def write(self, edi_info_list=None):
        """
        Build formatted ``>INFO`` lines ready to write to an EDI file.

        When ``edi_info_list`` is provided, it is normalized and returned.
        Otherwise, the current object state is serialized in canonical
        order, and additional known fields (CREATINGSOFTWARE, FILTER) are
        appended when available.

        Parameters
        ----------
        edi_info_list : sequence of str, optional
            Existing INFO lines to normalize/rewrite.

        Returns
        -------
        list of str
            Lines including the leading ``>INFO`` and a trailing blank
            line.
        """
        out = [">INFO\n"]

        if edi_info_list is not None:
            cleaned = []
            for raw in edi_info_list:
                s = raw.strip()
                if not s or s.startswith(">"):
                    continue
                m = _KV_RE.match(s)
                if not m:
                    continue
                key = _norm_key(m.group("key"))
                val = _unquote(m.group("val"))

                # quote some keys, uppercase others for style
                if key in self._quoted_keys:
                    out_val = f'"{val}"'
                else:
                    out_val = val.upper()
                cleaned.append(f"  {key.upper()}={out_val}\n")

            out.extend(cleaned)
            out.append("\n")
            return out

        # serialize from current attributes
        for key in self.infokeys:
            if key in {"project", "survey", "creationdate", "sitename"}:
                val = getattr(self.Source, key, None)
            elif key == "processingsoftware":
                val = getattr(self.Processing.ProcessingSoftware, "name", None)
            elif key in {
                "processedby",
                "processingtag",
                "runlist",
                "remoteref",
                "remotesite",
                "signconvention",
            }:
                val = getattr(self.Processing, key, None)
            else:
                val = getattr(self, key, None)

            if val in (None, "", "None"):
                continue

            if key in self._quoted_keys:
                out_val = f'"{val}"'
            else:
                out_val = str(val).upper()

            out.append(f"  {key.upper()}={out_val}\n")

        # Additional metadata if available
        if getattr(self.Source, "creatingsoftware", None):
            out.append(
                f"  CREATINGSOFTWARE={getattr(self.Source, 'creatingsoftware')}\n"
            )
        if self.filter:
            out.append(f"  FILTER={str(self.filter).upper()}\n")

        out.append("\n")
        return out
