# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

import datetime as _dt
import math
import re
from collections.abc import Iterable, Sequence
from pathlib import Path

from .. import __version__ as _PKG_VERSION
from ..log.logger import get_logger

try:
    from ..gis.utils import decimal_to_dms, dms_to_decimal
except (ImportError, RuntimeError):
    def dms_to_decimal(s, hem=None):  # type: ignore[misc]
        try:
            return float(str(s).strip())
        except Exception:
            return float("nan")

    def decimal_to_dms(v, hem="N"):  # type: ignore[misc]
        return f"{float(v):.6f}"
from ..exceptions import (
    EdIDataError,
    FileHandlingError,
    HeaderError,
)

try:
    from ..loc import Location
except (ImportError, RuntimeError):
    class Location:  # type: ignore[no-redef]
        def __init__(self):
            self.latitude = None
            self.longitude = None
            self.elevation = None

from .base import EDIComponentBase
from .property import Copyright, Processing, Software, Source
from .validation import IsEdi

logger = get_logger(__name__)


__all__ = [
    "Head",
    "Info",
    "Heads",
    "HeadMixin",
    "InfoMixin",
]

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
    """Normalize header keys (e.g. LON => LONG);
    keep canonical lower-case."""
    k0 = k.strip().lower()
    if k0 == "lon":
        return "long"
    return k0


def _is_tag(line: str, tag: str) -> bool:
    """Case-insensitive, whitespace-tolerant
    check for EDI tag line."""
    line = line.strip()
    if not line.startswith(">"):
        return False
    return line.upper().startswith(tag.upper())


def _slice_section(
    lines: Sequence[str], start_tag: str,
    after_tags: Iterable[str]
    ) -> tuple[list[str], int, int]:
    """
    Slice section content between start_tag and
    the next tag in after_tags.

    Returns (payload_lines, i_start, i_stop).
    Payload excludes the start_tag line.
    Raises EdIDataError if start_tag not present.
    """
    i_start = None
    i_stop = None
    for i, ln in enumerate(lines):
        if _is_tag(ln, start_tag):
            i_start = i
            break
    if i_start is None:
        raise EdIDataError(f"{start_tag} block not found in EDI file.")

    for j in range(i_start + 1, len(lines)):
        if any(_is_tag(lines[j], tag) for tag in after_tags
               ) or lines[j].lstrip().startswith(">"):
            i_stop = j
            break

    if i_stop is None:
        i_stop = i_start + 1
        while i_stop < len(lines) and not lines[i_stop].lstrip(
                ).startswith(">"):
            i_stop += 1

    payload = [ln.strip() for ln in lines[i_start + 1 : i_stop]]
    return payload, i_start, i_stop


class Head(EDIComponentBase):
    r"""
    EDI ``>HEAD`` block container.

    The :class:`Head` class parses and serializes the header
    metadata of SEG-EDI files. It normalizes common variations
    in field names (for example *LON* → *LONG*) and supports
    latitude/longitude given either in decimal degrees or DMS
    strings (with or without cardinal letters).

    Parameters
    ----------
    edi_header_list : sequence of str, optional
        Raw lines that belong to the ``>HEAD`` section.  If
        provided, they are parsed immediately via
        :meth:`read`.
    **kwargs
        Attribute overrides (for example ``dataid``,
        ``acqby``, ``coordsys``). Unknown keys become
        dynamic attributes.

    Attributes
    ----------
    Location : :class:`~pycsamt.loc.Location`
        Container for geographic coordinates.  The
        :pyattr:`lat`, :pyattr:`long` and :pyattr:`elev`
        properties delegate to this object.
    dataid, acqby, fileby : str or None
        Dataset identifier, acquisition contractor, and file
        author.
    acqdate, enddate, filedate : str or None
        Acquisition start/end dates and the file timestamp.
        ``filedate`` defaults to current UTC.
    country, state, county, prospect, loc : str or None
        Descriptive location metadata.
    lat, long, elev : float or None
        Coordinates as decimal degrees and elevation in the
        current :pyattr:`units`.  DMS inputs are converted to
        decimal when reading; writers emit DMS strings.
    units : {'m', 'ft'}
        Elevation units (default ``'m'``).
    stdvers, progvers, progdate : str
        EDI standard version, program version, and revision
        date. ``progvers`` defaults to ``pyCSAMT <version>``.
    coordsys : str
        Coordinate system description (default
        ``'Geomagnetic North'``).
    declination : float or None
        Geomagnetic declination in degrees.
    datum : str
        Geodetic datum (default ``'WGS84'``).
    maxsect : int or None
        Maximum section count in file when present.
    bindata : str or None
        Optional tag for external binary payload.
    project, survey : str or None
        Project and survey names if supplied.
    empty : float
        Missing-value sentinel (default ``1.0e32``).
    edi_header : list of str or None
        Normalized ``KEY=VALUE`` lines retained after parse.

    Notes
    -----
    * Key names are normalized to lower case and stored in a
      canonical set (for example *lon* is exposed as
      :pyattr:`long`).
    * Latitude and longitude accept DMS strings of the form
      ``'DD:MM:SS'`` with optional decimals and optional
      cardinal letters.  See :mod:`pycsamt.gis.utils`.
    * Writers quote some text fields to match common EDI
      formatting tools.

    Examples
    --------
    Read and access coordinates::

        h = Head.from_file("E01.edi")
        (h.lat, h.long, h.elev)

    Create and serialize a header block::

        h = Head(dataid="E1_2", acqby="ULTREM", units="m")
        lines = h.write()
        print("".join(lines))

    See Also
    --------
    Info
        Companion container for the ``>INFO`` block.
    Heads
        Aggregator that bundles :class:`Head` and
        :class:`Info`.
    pycsamt.gis.utils.dms_to_decimal
        Robust DMS → decimal converter.
    pycsamt.gis.utils.decimal_to_dms
        Decimal → DMS formatter.
    pycsamt.seg.properties.IsEdi
        File validator used by :meth:`from_file`.

    References
    ----------
    .. [1] SEG MT/EMAP EDI specification (1987/2006).  MTNet:
           https://www.mtnet.info/
    """
    # prog details
    _PROGVERS = f"pyCSAMT {_PKG_VERSION}"
    _PROGDATE = _dt.datetime.utcnow().strftime("%Y/%m/%d")
    _FILEDATE = _dt.datetime.utcnow().strftime("%Y/%m/%d %H:%M:%S UTC")

    # Canonical keys order for serialization
    head_keys: list[str] = [
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
        "chainage",
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
        self, edi_header_list: Sequence[str] | None = None,
        verbose =0, logger =None, **kwargs
        ):
        super().__init__(verbose = verbose, logger=logger)
        # Location container
        self.Location = Location()

        # Defaults
        self.dataid: str | None = None
        self.acqby: str | None = None
        self.fileby: str | None = None
        self.acqdate: str | None = None
        self.enddate: str | None = None
        self.filedate: str = self._FILEDATE

        self.country: str | None = None
        self.state: str | None = None
        self.county: str | None = None
        self.prospect: str | None = None
        self.loc: str | None = None

        self.units: str = "M".lower() if False else "m"
        self.stdvers: str = "SEG 1.0"
        self.progvers: str = self._PROGVERS
        self.progdate: str = self._PROGDATE
        self.coordsys: str = "Geomagnetic North"
        self.declination: float | None = None
        self.datum: str = "WGS84"
        self.maxsect: int | None = None
        self.bindata: str | None = None
        self.project: str | None = None
        self.survey: str | None = None
        self.empty: float = 1.0e32

        self.chainage: float | None = None

        self.edi_header: list[str] | None = list(
            edi_header_list) if edi_header_list else None

        # User overrides (unknown keys allowed)
        for k, v in kwargs.items():
            setattr(self, k, v)

        # Parse if provided
        if self.edi_header is not None:
            self.read(self.edi_header)

    # --------- Geographic IO via Location ----------
    @property
    def lat(self) -> float | None:
        return self.Location.latitude

    @lat.setter
    def lat(self, value: str | float | int | None) -> None:
        if value is None:
            self.Location.latitude = None
            return
        try:
            self.Location.latitude = float(value)
        except (TypeError, ValueError):
            self.Location.latitude = float(dms_to_decimal(str(value)))
            if self.verbose:
                logger.info(
                    "Converted DMS latitude to decimal degrees."
                )

    @property
    def long(self) -> float | None:
        return self.Location.longitude

    @long.setter
    def long(self, value: str | float | int | None) -> None:
        if value is None:
            self.Location.longitude = None
            return
        try:
            self.Location.longitude = float(value)
        except (TypeError, ValueError):
            self.Location.longitude = float(dms_to_decimal(str(value)))
            if self.verbose:
                logger.info(
                    "Converted DMS longitude to decimal degrees."
                    )

    @property
    def lon(self) -> float | None:
        return self.long

    @lon.setter
    def lon(self, value: str | float | int | None) -> None:
        self.long = value

    @property
    def elev(self) -> float | None:
        return self.Location.elevation

    @elev.setter
    def elev(self, value: str | float | int | None) -> None:
        if value is None or (isinstance(value, str) and not value.strip()):
            self.Location.elevation = None
        else:
            self.Location.elevation = float(value)


    # --------- API ----------
    @classmethod
    def from_file(cls, edi_fn: str | Path) -> Head:
        """
        Extract and parse the >HEAD block from an EDI file path.
        """
        if edi_fn is None:
            raise FileHandlingError("No EDI path provided.")
        edi_path = Path(edi_fn)

        # Validate EDI (deep check)
        IsEdi._assert_edi(edi_path, deep=True)

        lines = edi_path.read_text(
            encoding="utf-8-sig", errors="replace").splitlines()
        payload, _, _ = _slice_section(
            lines, ">HEAD", after_tags=[">INFO", ">=DEFINEMEAS"])
        # drop quotes in payload for tolerant parsing
        payload = [ln.replace('"', "") for ln in payload]
        return cls(edi_header_list=payload)

    def read(self, edi_header_list: Sequence[str] | None = None) -> Head:
        """
        Parse HEAD KV lines and set attributes.
        """
        if edi_header_list is not None:
            self.edi_header = list(edi_header_list)
        if not self.edi_header:
            raise HeaderError("No HEAD items to read.")

        normalized: list[str] = []
        for raw in self.edi_header:
            line = raw.strip()
            if not line or line.startswith(">"):
                continue
            m = _KV_RE.match(line)
            if not m:
                continue
            key = _norm_key(m.group("key"))
            val = _unquote(m.group("val"))

            if key == "chainage":
                try:
                    self.chainage = float(val)
                except Exception:
                    self.chainage = None
                normalized.append(f"{key.upper()}={val}")
                continue

            if key == "coordsys":
                self.coordsys = val
            elif key == "lon":
                self.long = val
            else:
                setattr(self, key, val)

            normalized.append(f"{key.upper()}={val}")

        self.edi_header = normalized
        return self

    def write(
        self, head_list_infos: Sequence[str] | None = None,
        stamp: bool=True,
        ) -> list[str]:
        """
        Build formatted >HEAD lines (including trailing blank line).
        """
        lines: list[str] = [">HEAD\n"]

        # If given explicit lines, normalize & echo them
        if head_list_infos is not None:
            cleaned: list[str] = []
            for raw in head_list_infos:
                s = raw.strip()
                if not s or s.startswith(">"):
                    continue
                m = _KV_RE.match(s)
                if not m:
                    continue
                key = _norm_key(m.group("key"))
                val = _unquote(m.group("val"))
                if key in self._quoted_keys:
                    out_val = f'"{val}"'
                else:
                    out_val = val.upper() if key not in {"lat", "long"} else val
                cleaned.append(f"  {key.upper()}={out_val}\n")
            lines.extend(cleaned)
            lines.append("\n")
            return lines

        # Serialize internal state in canonical order
        for key in self.head_keys:
            if stamp and key == "progvers":
                val = self._PROGVERS  # This gets the current default value
            elif stamp and key == "progdate":
                val = self._PROGDATE  # This gets the current default value
            elif stamp and key=="filedate":
                val= self._FILEDATE
            else:
                val = getattr(self, key, None)

            if val in (None, "", "None"):
                continue

            if key in {"lat", "long"}:
                vnum = getattr(self, key)
                try:
                    out_val = decimal_to_dms(float(vnum))
                except Exception:
                    out_val = str(vnum)
            elif key == "elev":
                out_val = str(val)
            elif key in self._quoted_keys:
                out_val = f'"{val}"'

            elif key == "chainage":
                try:
                    out_val = str(float(val))
                except:
                    continue

            else:
                out_val = str(val).upper()

            lines.append(f"  {key.upper()}={out_val}\n")

        lines.append("\n")
        return lines

    # Convenience helpers
    def as_dict(self) -> dict:
        """Export known fields as a plain dict."""
        d = {k: getattr(self, k, None) for k in self.head_keys}
        # inject decimals for lat/long (if present)
        d["lat"] = self.lat
        d["long"] = self.long
        d["elev"] = self.elev
        return d

    def update(self, **kwargs) -> Head:
        """Update fields (accepts both known and unknown keys)."""
        for k, v in kwargs.items():
            setattr(self, k, v)
        return self

    def compute_chainage(
        self,
        origin: tuple[float, float],
        azimuth: float,
        *,
        set_attr: bool = True,
    ) -> float:
        """
        Compute chainage (meters) along a profile defined by an
        origin and azimuth, using a local flat metric.

        Chainage is positive in the forward profile direction and
        negative behind the origin.

        Parameters
        ----------
        origin : tuple of float
            (lat0, lon0) in decimal degrees (WGS84 assumed).
        azimuth : float
            Profile azimuth in degrees. North = 0, East = 90.
        set_attr : bool, optional
            If True, store the computed value into
            ``self.chainage``.

        Returns
        -------
        float
            Chainage in meters. Returns ``nan`` if coordinates
            are missing.
        """

        la0, lo0 = origin
        if self.lat is None or self.long is None:
            ch = float("nan")
            if set_attr:
                self.chainage = ch
            return ch

        # local flat approximation
        dy = (float(self.lat) - float(la0)) * 111_000.0
        dx = ((float(self.long) - float(lo0)) * 111_000.0
              * math.cos(math.radians(float(la0))))

        az = math.radians(float(azimuth))
        # north=0 => dx*sin + dy*cos
        ch = dx * math.sin(az) + dy * math.cos(az)

        if set_attr:
            try:
                self.chainage = float(ch)
            except:
                self.chainage = None
        return float(ch)

    def __repr__(self) -> str:  # compact identity for debugging
        return f"<Head dataid={self.dataid!r} lat={self.lat} long={self.long}>"


class Info(EDIComponentBase):
    r"""
    EDI ``>INFO`` block container.

    The :class:`Info` class collects survey-level metadata and
    processing provenance.  It supports **both** classic
    ``KEY=VALUE`` blocks and spectra-style **free-text** INFO
    blocks.  Non-KV lines are preserved in
    :pyattr:`info_text` and are written back unmodified.

    Parameters
    ----------
    edi_info_list : sequence of str, optional
        Raw ``KEY=VALUE`` lines for the INFO block.  If
        provided, they are parsed immediately via
        :meth:`read`.  Non-KV lines should be passed via
        :pyattr:`info_text` or supplied through
        :meth:`from_file`.
    verbose : int, optional
        Verbosity level forwarded to the base component.
    logger : logging.Logger, optional
        Logger instance to use for diagnostics.
    info_text : sequence of str, optional
        Free-text lines from spectra-style INFO blocks that
        should be preserved.
    **kwargs
        Attribute overrides; unknown keys are accepted.

    Attributes
    ----------
    maxinfo : int
        Maximum number of textual entries (default ``999``).
    Source : :class:`~pycsamt.seg.property.Source`
        Survey provenance (project, survey, site name, and
        creation date).
    Processing : :class:`~pycsamt.seg.property.Processing`
        Processing meta-information (software, run list,
        sign convention, etc.).
    Copyright : :class:`~pycsamt.seg.property.Copyright`
        Optional copyright container.
    filter : str or None
        Optional filter tag.
    ediinfo : list of str
        Normalized KV lines retained after parse.
    info_text : list of str
        Free-form text preserved and written back.

    Notes
    -----
    * Known keys are routed into nested containers (for
      example ``processedby`` → :pyattr:`Processing`).
    * If an INFO block contains **no** KV lines, parse is
      still successful; fields remain at defaults and
      :pyattr:`info_text` carries the original text.

    Examples
    --------
    Parse INFO from file with mixed content::

        i = Info.from_file("15125A_spe.edi")
        i.Source.project
        len(i.info_text)  # free-text lines preserved

    Build and write a KV-only INFO block::

        i = Info()
        i.Source.project = "Demo"
        i.Processing.processedby = "pyCSAMT"
        print("".join(i.write()))

    See Also
    --------
    Head
        Header metadata companion.
    Heads
        Aggregator that bundles :class:`Head` and
        :class:`Info`.
    pycsamt.seg.property.Source
    pycsamt.seg.property.Processing
    pycsamt.seg.property.Copyright

    References
    ----------
    .. [1] SEG MT/EMAP EDI specification (1987/2006).  MTNet:
           https://www.mtnet.info/
    """
    _PROCESSEDBY = f"pyCSAMT v{_PKG_VERSION}"
    _PROCESSINGSOFTWARE =Software(name= "PYCSAMT")

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

    def __init__(
        self,
        edi_info_list: Sequence[str] | None = None,
        *,
        verbose: int = 0,
        logger=None,
        info_text: Sequence[str] | None = None,
        **kwargs,
    ):
        super().__init__(verbose=verbose, logger=logger)

        # Raw KV lines we could parse (normalized later)
        self.ediinfo: list[str] | None = list(
            edi_info_list
        ) if edi_info_list else None
        # Free-text payload (lines that are *not* KEY=VALUE)
        self.info_text: list[str] = list(info_text) if info_text else []

        self.filter: str | None = None  # used in some EMAP/processing contexts
        self.maxinfo: int = 999

        # nested containers
        self.Source = Source()
        self.Processing = Processing(
            processedby=self._PROCESSEDBY,
            ProcessingSoftware = self._PROCESSINGSOFTWARE,
        )
        self.Copyright = Copyright()

        for k, v in kwargs.items():
            setattr(self, k, v)

        # If KV lines were passed, parse immediately
        if self.ediinfo is not None:
            self.read(self.ediinfo)


    @classmethod
    def from_file(cls, edi_fn: str | Path) -> Info:
        """
        Extract and parse the ``>INFO`` block from an EDI file.

        Works for both KV-style INFO and free-text INFO blocks.
        """
        if edi_fn is None:
            raise FileHandlingError("No EDI path provided.")

        p = Path(edi_fn)
        IsEdi._assert_edi(p, deep=True)

        lines = p.read_text(encoding="utf-8-sig", errors="replace").splitlines()

        # Slice INFO payload until next section (>=DEFINEMEAS, >=MTSECT, or any '>' tag)
        start, stop = None, None
        for i, ln in enumerate(lines):
            if _is_tag(ln, ">INFO"):
                start = i
                break
        if start is None:
            raise EdIDataError("INFO block not found in EDI file.")

        for j in range(start + 1, len(lines)):
            if _is_tag(lines[j], ">=DEFINEMEAS") or _is_tag(lines[j], ">=MTSECT"):
                stop = j
                break
            if lines[j].lstrip().startswith(">"):
                stop = j
                break
        if stop is None:
            stop = start + 1
            while stop < len(lines) and not lines[stop].lstrip().startswith(">"):
                stop += 1

        payload = lines[start + 1 : stop]

        # Split payload into KV lines vs free text; strip outer quotes on KV.
        kv_lines: list[str] = []
        text_lines: list[str] = []
        for raw in payload:
            s = raw.strip()
            if not s:
                # keep spacing inside text block; treat as text line
                text_lines.append(raw)
                continue
            m = _KV_RE.match(s)
            if m:
                key = _norm_key(m.group("key"))
                val = _unquote(m.group("val"))
                kv_lines.append(f"{key}={val}")
            else:
                text_lines.append(raw)

        obj = cls(edi_info_list=kv_lines, info_text=text_lines)
        # If we had KV lines, __init__ already parsed them. If not, that's fine.
        return obj

    def read(
            self, edi_info_list: Sequence[str] | None = None,
        ) -> Info:
        """
        Parse an INFO key-value list and set attributes.

        If there are no KEY=VALUE items (e.g., spectra-style free-text INFO),
        this method simply leaves metadata at defaults and returns without error.
        Any non-KV lines should be provided via ``info_text``
        (from_file handles this).
        """
        if edi_info_list is not None:
            self.ediinfo = list(edi_info_list)

        # Nothing to read is a valid state for spectra-style INFO; just return.
        if not self.ediinfo:
            self.ediinfo = []
            return self

        normalized: list[str] = []
        for raw in self.ediinfo:
            line = raw.strip()
            if not line or line.startswith(">"):
                continue
            m = _KV_RE.match(line)
            if not m:
                # KV list given but line isn't KV; treat as free text
                self.info_text.append(raw)
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
            elif key == "maxinfo":
                # normalize to int when possible
                try:
                    self.maxinfo = int(float(val))
                except Exception:
                    self.maxinfo = 999
            else:
                # accept other fields as direct attributes
                setattr(self, key, val)

            normalized.append(f"{key.upper()}={val}")

        self.ediinfo = normalized
        logger.debug("Parsed INFO: %s", self.ediinfo)
        return self


    def write(self, edi_info_list: Sequence[str] | None = None) -> list[str]:
        """
        Build formatted ``>INFO`` lines ready to write to an EDI file.

        - If ``edi_info_list`` is provided, it is normalized and written.
        - Otherwise, current attributes are serialized in canonical order,
          followed by any preserved free-text lines stored in ``info_text``.
        """
        out: list[str] = [">INFO\n"]

        # If caller supplied lines, normalize KV and pass through non-KV as text.
        if edi_info_list is not None:
            cleaned: list[str] = []
            text_passthru: list[str] = []
            for raw in edi_info_list:
                s = raw.strip()
                if not s or s.startswith(">"):
                    continue
                m = _KV_RE.match(s)
                if not m:
                    text_passthru.append(raw.rstrip("\n"))
                    continue
                key = _norm_key(m.group("key"))
                val = _unquote(m.group("val"))

                if key in self._quoted_keys:
                    out_val = f'"{val}"'
                else:
                    out_val = val.upper()
                cleaned.append(f"  {key.upper()}={out_val}\n")

            out.extend(cleaned)
            # Preserve any non-KV lines we were given
            for t in text_passthru:
                out.append(f"{t.rstrip()}\n")
            out.append("\n")
            return out

        # Serialize KV items from current attributes/containers
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

        # Append preserved free-text INFO lines (spectra-style blocks, etc.)
        if self.info_text:
            for ln in self.info_text:
                # write as-is (no extra '  ' if original spacing is desired)
                out.append(f"{ln.rstrip()}\n")

        out.append("\n")
        return out

    # -------- convenience --------
    def as_dict(self) -> dict:
        """Flattened INFO dictionary (KV items only)."""
        d = {k: getattr(self, k, None) for k in self.infokeys}
        d.update(
            {
                "project": getattr(self.Source, "project", None),
                "survey": getattr(self.Source, "survey", None),
                "creationdate": getattr(self.Source, "creationdate", None),
                "sitename": getattr(self.Source, "sitename", None),
                "processingsoftware": getattr(
                    self.Processing.ProcessingSoftware, "name", None
                ),
            }
        )
        return d

    def update(self, **kwargs) -> Info:
        """Update INFO fields (routes to nested containers where appropriate)."""
        for k, v in kwargs.items():
            key = _norm_key(k)
            if key in {"project", "survey", "creationdate", "sitename"}:
                setattr(self.Source, key, v)
            elif key == "processingsoftware":
                self.Processing.ProcessingSoftware.name = v
            elif key in {
                "processedby",
                "processingtag",
                "runlist",
                "remoteref",
                "remotesite",
                "signconvention",
            }:
                setattr(self.Processing, key, v)
            elif key == "info_text" and isinstance(v, (list, tuple)):
                self.info_text = list(v)
            else:
                setattr(self, key, v)
        return self

    def __repr__(self) -> str:
        return (
            f"<Info project={self.Source.project!r} "
            f"processedby={self.Processing.processedby!r} "
            f"text_lines={len(self.info_text)}>"
        )



class Heads(EDIComponentBase):
    r"""
    Convenience aggregator for ``>HEAD`` and ``>INFO``.

    The :class:`Heads` helper wraps a parsed
    :class:`Head` and :class:`Info` pair and provides simple
    I/O helpers to extract or write both blocks together.  It
    is useful when a workflow wants to treat the *top matter*
    of an EDI file as a single unit.

    Parameters
    ----------
    head : :class:`Head`, optional
        Parsed header.  One will be created by
        :meth:`from_file` when reading from disk.
    info : :class:`Info`, optional
        Parsed INFO.  One will be created by
        :meth:`from_file` when reading from disk.

    Attributes
    ----------
    head : :class:`Head`
        Header container.
    info : :class:`Info`
        INFO container.

    Notes
    -----
    * :meth:`from_file` validates the input via
      :class:`~pycsamt.seg.properties.IsEdi`, then slices the
      two blocks and delegates to :class:`Head` and
      :class:`Info`.
    * :meth:`write` serializes in canonical order
      (``>HEAD`` first, then ``>INFO``), preserving free-text
      INFO content when present.

    Examples
    --------
    Read both blocks at once and re-emit them::

        hs = Heads.from_file("000CSA_csamt.edi")
        lines = hs.write()
        open("out.edi", "w", encoding="utf-8").writelines(lines)

    Access nested objects directly::

        hs.head.dataid, hs.info.Source.project

    See Also
    --------
    Head, Info
        Underlying block containers.
    pycsamt.seg.properties.IsEdi
        Validator used at read time.

    References
    ----------
    .. [1] SEG MT/EMAP EDI specification (1987/2006).  MTNet:
           https://www.mtnet.info/
    """


    def __init__(
            self, head: Head | None = None, info: Info | None = None,
            verbose: int = 0, logger=None ):
        super().__init__(verbose=verbose, logger=logger )
        self.head: Head = head if head is not None else Head()
        self.info: Info = info if info is not None else Info()

    @classmethod
    def from_file(cls, edi_fn: str | Path) -> Heads:
        """
        Load >HEAD and >INFO from EDI path and return an aggregate container.
        """
        p = Path(edi_fn)
        if not p:
            raise FileHandlingError("No EDI path provided.")
        IsEdi._assert_edi(p, deep=True)

        lines = p.read_text(encoding="utf-8-sig", errors="replace").splitlines()

        # HEAD
        head_payload, _, stop_h = _slice_section(
            lines, ">HEAD", after_tags=[">INFO", ">=DEFINEMEAS"])
        head_payload = [ln.replace('"', "") for ln in head_payload]
        head = Head(edi_header_list=head_payload).read()

        # INFO (optional in theory, but standard requires it right after HEAD)
        try:
            info_payload, _, _ = _slice_section(
                lines, ">INFO", after_tags=[">=DEFINEMEAS", ">=MTSECT"])
            info_payload = [ln.replace('"', "") for ln in info_payload]
            info = Info(edi_info_list=info_payload).read()
        except EdIDataError:
            # Keep an empty info; callers may enrich later
            logger.warning(
                "INFO block not found immediately after HEAD; created empty Info()")
            info = Info()

        return cls(head=head, info=info)

    def read(self, text_or_lines: str | Sequence[str]) -> Heads:
        """
        Parse >HEAD and >INFO from provided content (string or list of lines).
        """
        if isinstance(text_or_lines, str):
            lines = text_or_lines.splitlines()
        else:
            lines = list(text_or_lines)

        # HEAD
        head_payload, _, _ = _slice_section(
            lines, ">HEAD", after_tags=[">INFO", ">=DEFINEMEAS"])
        self.head = Head(edi_header_list=[ln.replace(
            '"', "") for ln in head_payload]).read()

        # INFO
        try:
            info_payload, _, _ = _slice_section(
                lines, ">INFO", after_tags=[">=DEFINEMEAS", ">=MTSECT"])
            self.info = Info(edi_info_list=[
                ln.replace('"', "") for ln in info_payload]).read()
        except EdIDataError:
            self.info = Info()  # leave empty
        return self

    def write(self) -> list[str]:
        """
        Serialize >HEAD + >INFO blocks back-to-back.
        """
        out: list[str] = []
        out.extend(self.head.write())
        out.extend(self.info.write())
        return out

    # Simple helpers
    def to_text(self) -> str:
        """Return concatenated HEAD+INFO text."""
        return "".join(self.write())

    def __repr__(self) -> str:
        return (
            f"<Heads dataid={self.head.dataid!r}"
            f" project={self.info.Source.project!r}>"
            )

class HeadMixin:
    r"""
    Mixin that exposes :class:`Head` helpers on a host class.

    ``HeadMixin`` adds convenience constructors for reading
    the ``>HEAD`` block without forcing the host to inherit
    from :class:`Head` directly.  This is primarily used by
    higher-level SEG components that *contain* a header but
    are not themselves :class:`Head` objects.

    Methods
    -------
    from_file(edi_fn)
        Return a new :class:`Head` parsed from *edi_fn*.
    read(lines)
        Shortcut to :meth:`Head.read`.
    write(...)
        Shortcut to :meth:`Head.write`.

    Notes
    -----
    The mixin does not own any state; it simply forwards to
    :class:`Head`.

    Examples
    --------
    Use the mixin inside a host container::

        class Host(HeadMixin):
            pass

        h = Host.from_file("E01.edi")
        print(h.dataid)

    See Also
    --------
    Head, InfoMixin
    """


    head: Head

    # -- class-level convenience --
    @classmethod
    def from_file(cls, edi_fn: str | Path) -> Head:
        """Return a parsed Head instance from EDI path."""
        return Head.from_file(edi_fn)

    # -- instance-level --
    def read(self, edi_header_list: Sequence[str] | None = None) -> Head:
        """Parse HEAD KV lines into `self.head` and return it."""
        if not hasattr(self, "head") or self.head is None:
            self.head = Head()
        return self.head.read(edi_header_list)

    def write(self, head_list_infos: Sequence[str] | None = None) -> list[str]:
        """Serialize the host's `head` as >HEAD lines."""
        if not hasattr(self, "head") or self.head is None:
            self.head = Head()
        return self.head.write(head_list_infos)


class InfoMixin:
    r"""
    Mixin that exposes :class:`Info` helpers on a host class.

    ``InfoMixin`` adds convenience constructors for reading
    the ``>INFO`` block while keeping the host class focused
    on its domain logic.  It forwards calls to
    :class:`Info` and returns an :class:`Info` instance.

    Methods
    -------
    from_file(edi_fn)
        Return a new :class:`Info` parsed from *edi_fn*.
    read(lines)
        Shortcut to :meth:`Info.read`.
    write(...)
        Shortcut to :meth:`Info.write`.

    Notes
    -----
    Free-text INFO blocks (for example spectra-style RUN
    INFORMATION) are supported transparently by
    :class:`Info`.  The mixin needs no special handling.

    Examples
    --------
    Mix in the helpers::

        class Host(InfoMixin):
            pass

        info = Host.from_file("15125A_spe.edi")
        len(info.info_text)  # preserved text lines

    See Also
    --------
    Info, HeadMixin
    """


    info: Info

    # -- class-level convenience --
    @classmethod
    def from_file(cls, edi_fn: str | Path) -> Info:
        """Return a parsed Info instance from EDI path."""
        return Info.from_file(edi_fn)

    # -- instance-level --
    def read(self, edi_info_list: Sequence[str] | None = None) -> Info:
        """Parse INFO KV lines into `self.info` and return it."""
        if not hasattr(self, "info") or self.info is None:
            self.info = Info()
        return self.info.read(edi_info_list)

    def write(self, edi_info_list: Sequence[str] | None = None) -> list[str]:
        """Serialize the host's `info` as >INFO lines."""
        if not hasattr(self, "info") or self.info is None:
            self.info = Info()
        return self.info.write(edi_info_list)


