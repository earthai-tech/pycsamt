# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

from pathlib import Path
from typing import Any, List, Optional, Sequence, Union

from datetime import datetime 

from ..exceptions import JParseError
from .base import JComponentBase
from .config import (
    RE_DATATYPE_UNITS,
    ENCODING_DEFAULT,
    RE_BLANK,
    RE_COMMENT,
    RE_BANNER, 
    RE_INFO,
    RE_NPOINTS,
    RE_STATION,
    UNITS_CANONICAL,
)
from .utils import DataType, iter_lines, parse_datatype_units
from .property import JSiteProperty


__all__ = [
    "Head",
    "Info",
    "Heads",
    "HeadMixin",
    "InfoMixin",
    "Banner"
]


class Banner(JComponentBase):
    r"""
    Parse and serialize the top provenance comment line.

    This helper targets lines such as::

        #WRITTEN BY GEOTOOLS: kb0-s001 10/06/95 RAW RECS

    It extracts the producer software name, an optional station
    hint, a free-form date token, and an optional trailing note.
    The writer defaults the software field to ``PYSCAMT`` when
    none is provided.

    Parameters
    ----------
    top_lines : sequence of str, optional
        An initial slice of file lines. The constructor scans
        these lines and parses the first matching banner.
    software : str, optional
        Producer software name (e.g., ``'GEOTOOLS'``). When
        omitted and the banner is not present, the writer uses
        ``'PYSCAMT'`` as a default.
    station : str, optional
        Station identifier hint found in the banner. This does
        not replace the station found in the :class:`Head`.
    date : str, optional
        Free-form date token as found in the banner line.
    note : str, optional
        Extra text following the date token on the same line.
    verbose : int, default=0
        Verbosity for warnings during parsing.

    Attributes
    ----------
    software : str or None
        Producer software parsed from the banner or provided by
        the user.
    station_hint : str or None
        Station token parsed from the banner; may differ from the
        station found in the :class:`Head`.
    date : str or None
        Banner date token as text; not parsed into a datetime.
    note : str or None
        Trailing free-form text after the date token.
    path : pathlib.Path or None
        Source path when constructed via :meth:`from_file`.
    encoding : str
        Text encoding used by :meth:`from_file`.
    verbose : int
        Verbosity level inherited from :class:`JComponentBase`.

    Methods
    -------
    from_file(j_fn, *, verbose=0)
        Read a file, scan the top lines, and return a parsed
        :class:`Banner`.
    from_lines(lines, *, verbose=0)
        Build a :class:`Banner` from an in-memory sequence of
        lines.
    read(top_lines)
        Parse the first matching banner from the given lines and
        mark the instance as read.
    write()
        Render a banner line. Defaults the software field to
        ``PYSCAMT`` if missing.

    Notes
    -----
    The banner parser is intentionally liberal: it ignores
    leading whitespace, is case-insensitive on the ``WRITTEN BY``
    marker, and preserves the raw text of the trailing note. It
    does not validate dates or enforce a specific date format.

    Examples
    --------
    >>> b = Banner().read([  # doctest: +NORMALIZE_WHITESPACE
    ...   '#WRITTEN BY GEOTOOLS: kb0-s001 10/06/95 RAW RECS'
    ... ])
    >>> b.software, b.station_hint, b.date
    ('GEOTOOLS', 'kb0-s001', '10/06/95')
    >>> Banner().write()[0].startswith('#WRITTEN BY PYSCAMT:')
    True

    See Also
    --------
    Head : Station / data-type / count header triple.
    Info : Site information (``>KEY=VALUE`` records).
    Heads : Minimal container that includes a :class:`Banner`.

    References
    ----------
    .. [1] A. G. Jones (1994). Magnetotelluric data file J-format,
       version 2.0.
    .. [2] MTNet. "J format documentation".
    """
    def __init__(
        self,
        top_lines: Optional[Sequence[str]] = None,
        *,
        software: Optional[str] = None,
        station: Optional[str] = None,
        date: Optional[str] = None,
        note: Optional[str] = None,
        verbose: int = 0,
    ) -> None:
        super().__init__(verbose=verbose)
       
        self.station_hint = station
        self.date = date
        self.note = note
        self._software = software
        self._raw: Optional[str] = None
        if top_lines is not None:
            self.read(top_lines)

    @classmethod
    def from_file(
        cls, j_fn: Union[str, Path], *, verbose: int = 0
    ) -> "Banner":
        lines = list(iter_lines(j_fn, encoding=ENCODING_DEFAULT))
        return cls.from_lines(lines, verbose=verbose)

    @classmethod
    def from_lines(
        cls, lines: Sequence[str], *, verbose: int = 0
    ) -> "Banner":
        inst = cls(verbose=verbose)
        inst.read(lines)
        return inst

    def read(
        self, top_lines: Optional[Sequence[str]] = None
    ) -> "Banner":
        if top_lines is None:
            raise ValueError("top_lines is required")
        for ln in top_lines:
            m = RE_BANNER.match(ln.rstrip("\n"))
            if m:
                self._software = m.group("software").strip()
                self.station_hint = m.group("station").strip()
                self.date = m.group("date").strip()
                self.note = m.group("note").strip()
                self._raw = ln.rstrip("\n")
                break
        self._mark_read(True)
        return self

    def write(
        self,
        *,
        new: bool = True,
        include_origin: bool = False,
    ) -> list[str]:
        st = self.station_hint or ""
        lines: list[str] = []
    
        if new:
            sw = "PYCSAMT"
            dt = datetime.now().strftime("%d/%m/%y")
            nt = (self.note or "RAW RECS").strip().upper()
            lines.append(
                f"#WRITTEN BY {sw}: {st} {dt} {nt}".rstrip()
            )
        else:
            sw = self.software or "PYCSAMT"
            dt = (
                self.date
                if self.date is not None
                else datetime.now().strftime("%d/%m/%y")
            )
            nt = (self.note or "RAW RECS").strip().upper()
            lines.append(
                f"#WRITTEN BY {sw}: {st} {dt} {nt}".rstrip()
            )
    
        if include_origin and self._raw:
            # keep original line verbatim, but mark as provenance
            lines.append(
                "#FROM "
                + self._raw.lstrip("#").strip().lstrip("WRITTEN BY")
            )
    
        return lines

    @property
    def software(self) -> str | None:
        return ( 
            None if self._software is None 
            else self._software.upper()
        )
    
    @property
    def date_parsed(self) -> datetime | None:
        if not self.date:
            return None
        txt = self.date.strip()
        fmts = ("%d/%m/%y", "%Y-%m-%d", "%d-%m-%Y")
        for f in fmts:
            try:
                return datetime.strptime(txt, f)
            except Exception:
                pass
        return None


class Info(JComponentBase):
    r"""
    Parse and serialize the J-format information block.

    This component collects ``>KEY=VALUE`` records along with
    leading comment lines. It exposes a small set of convenience
    properties derived from :class:`JSiteProperty`.

    Parameters
    ----------
    j_info_list : sequence of str, optional
        A sequence containing the comment and information
        records. The parser stops at the first non-info,
        non-comment, non-blank line.
    verbose : int, default=0
        Verbosity for warnings during parsing.
    **kwargs
        Accepted for API forwards-compatibility; ignored here.

    Attributes
    ----------
    items : dict of (str -> str)
        Mapping of upper-cased keys to unmodified string values.
    comments : list of str
        Preserved leading comment lines (starting with ``#``).
    site : JSiteProperty
        Lazily parsed view providing normalized latitude,
        longitude, azimuth and elevation (see properties below).
    latitude : float or None
        Decimal degrees; hemisphere and DMS tolerated.
    longitude : float or None
        Decimal degrees in ``[-180, 180)``.
    azimuth : float or None
        Site X-axis azimuth (degrees, true north).
    elevation : float or None
        Elevation in metres.
    path : pathlib.Path or None
        Source path when constructed via :meth:`from_file`.
    encoding : str
        Text encoding used by :meth:`from_file`.
    verbose : int
        Verbosity level inherited from :class:`JComponentBase`.

    Methods
    -------
    from_file(j_fn, *, verbose=0)
        Read only the comment/info header from a file.
    from_lines(j_info_list, *, verbose=0)
        Build from an in-memory sequence.
    read(j_info_list)
        Parse the header and mark the instance as read.
    write(j_info_list=None)
        Render comments and ``>KEY = VALUE`` lines.

    Notes
    -----
    Unknown keys are preserved verbatim. Coordinate and azimuth
    values are normalized via :class:`JSiteProperty`, which
    handles ranges, hemispheres and DMS.

    Examples
    --------
    >>> info = Info.from_file('data/j/kb0-s001.txt')
    >>> info.latitude, info.longitude
    (41.9782, 140.8958)
    >>> lines = info.write()
    >>> Info.from_lines(lines).azimuth == info.azimuth
    True

    See Also
    --------
    Head : One data-block header (station / dtype / count).
    Heads : Combined view with convenience properties.
    JSiteProperty : Robust parsing used by this class.

    References
    ----------
    .. [1] A. G. Jones (1994). Magnetotelluric data file
       J-format, version 2.0.
    """
    def __init__(
        self,
        j_info_list: Optional[Sequence[str]] = None,
        *,
        verbose: int = 0,
        **kwargs: Any,
    ) -> None:
        super().__init__(verbose=verbose)
        self.comments: List[str] = []
        self.items: dict[str, str] = {}
        self._site_cache: JSiteProperty | None = None
        if j_info_list is not None:
            self.read(j_info_list)


    @classmethod
    def from_file(
        cls, j_fn: Union[str, Path], *, verbose: int = 0
    ) -> "Info":
        lines = list(iter_lines(j_fn, encoding=ENCODING_DEFAULT))
        info_list = _extract_info_header_list(lines)
        return cls(info_list, verbose=verbose)

    @classmethod
    def from_lines(
            cls,
            j_info_list: Sequence[str] | None = None,
            *,
            verbose: int = 0,
        ) -> "Info":
        if j_info_list is None:
            raise ValueError("j_info_list is required")
        seq = list(j_info_list)
        info_list = _extract_info_header_list(seq)
        
        return cls(info_list, verbose=verbose)

    def write(
        self, j_info_list: Optional[Sequence[str]] = None
    ) -> List[str]:
        if j_info_list is not None:
            self.read(j_info_list)
        # out: List[str] = []
        # for c in self.comments:
        #     out.append(c.rstrip("\n"))
        # for k, v in self.items.items():
        #     out.append(f">{k} = {v}")
        # return out
        out: list[str] = []
        for c in self.comments:
            # skip any original “#WRITTEN BY …” banner lines
            if RE_BANNER.match(c):
                continue
            out.append(c.rstrip("\n"))
        for k, v in self.items.items():
            out.append(f">{k} = {v}")
        return out
    
    def read(
        self, j_info_list: Optional[Sequence[str]] = None
    ) -> "Info":
        if j_info_list is None:
            raise ValueError("j_info_list is required")
        self.comments = []
        self.items = {}
        self._site_cache = None
        for ln in j_info_list:
            if _is_comment(ln):
                self.comments.append(ln)
                continue
            if _is_blank(ln):
                continue
            if _is_info(ln):
                k, v = ln.split("=", 1)
                key = k.strip().lstrip(">").strip().upper()
                self.items[key] = v.strip()
        self._mark_read(True)
        return self

    @property
    def site(self) -> JSiteProperty:
        if self._site_cache is None:
            self._site_cache = JSiteProperty.from_lines(
                self.write(), verbose=self.verbose
            )
        return self._site_cache

    @property
    def latitude(self) -> float | None:
        return self.site.latitude

    @property
    def longitude(self) -> float | None:
        return self.site.longitude

    @property
    def azimuth(self) -> float | None:
        return self.site.azimuth

    @property
    def elevation(self) -> float | None:
        return self.site.elevation


    def get(self, key: str, default: str | None = None
            ) -> str | None:
        return self.items.get(key.upper(), default)
    
    def keys(self) -> tuple[str, ...]:
        return tuple(self.items.keys())
    
    def values(self) -> tuple[str, ...]:
        return tuple(self.items.values())
    
    def items_map(self) -> dict[str, str]:
        # explicit copy to avoid accidental mutation
        return dict(self.items)
    
    @property
    def lat(self) -> float | None:
        return self.latitude
    
    @property
    def lon(self) -> float | None:
        return self.longitude

    
    def __contains__(self, key: str) -> bool:
        return key.upper() in self.items
    
    def __getitem__(self, key: str) -> str:
        return self.items[key.upper()]
    
    def __str__(self) -> str:
        n = len(self.items)
        return f"Info(items={n})"

    def __repr__(self) -> str:
        return ( 
            f"Info(items={len(self.items)},"
            " cmts={len(self.comments)})"
            )


class Head(JComponentBase):
    r"""
    Parse and serialize a single J-format *head* triple.

    This lightweight component represents one data-block header
    consisting of the station line, the data-type line, and the
    row count. It provides a consistent constructor and the
    ``read`` / ``write`` round-trip expected by the package.

    Parameters
    ----------
    j_header_list : sequence of str, optional
        Three logical header lines. The sequence may be a clean
        triple ``[station, dtype, n]`` or a longer slice from a
        file; the parser will extract the first valid triple.
    verbose : int, default=0
        Verbosity for warnings during parsing.
    **kwargs
        Accepted for API forwards-compatibility; ignored here.

    Attributes
    ----------
    station : str or None
        Upper-cased station identifier. Trailing azimuth tokens
        on the same line are supported and ignored here (but see
        ``az_hint``).
    dtype : DataType or None
        Parsed data-type token (e.g., ``ZXY``, ``RTE``) with
        optional normalized units and a TE/TM tensor hint.
    n : int or None
        Number of data rows that follow this header.
    az_hint : float or None
        Optional azimuth extracted from the station line when
        present (e.g., ``KB0001  -30``). This is *not* the same
        as site ``AZIMUTH`` from the info block.
    path : pathlib.Path or None
        Source path when constructed via :meth:`from_file`.
    encoding : str
        Text encoding used by :meth:`from_file`.
    verbose : int
        Verbosity level inherited from :class:`JComponentBase`.

    Methods
    -------
    from_file(j_fn, *, verbose=0)
        Read a file and return a :class:`Head`. Only the first
        valid triple is parsed.
    from_lines(j_header_list, *, verbose=0)
        Build a :class:`Head` from an in-memory sequence. If the
        sequence is longer than three lines, the first valid
        triple is extracted.
    read(j_header_list)
        Parse the header triple into attributes and mark the
        instance as read (see ``__has_read__``).
    write(head_list_infos=None)
        Return a list of three strings representing the header
        triple. If ``head_list_infos`` is given, it is parsed
        first and then serialized.

    Notes
    -----
    The station regex accepts common field variants, including
    hyphens and underscores, and an optional numeric azimuth
    token after the station id. The data-type parser is liberal
    with whitespace and recognizes ``SI`` (``S.I.``) and
    ``FIELD`` units.

    Examples
    --------
    >>> Head().read(['KB0001  -30', 'ZXY SI', '29'])
    Head(station='KB0001', n=29)
    >>> h = Head.from_file('data/j/kb0-s001.txt')
    >>> h.station
    'KB0001'

    See Also
    --------
    Info : Site information (``>KEY=VALUE`` lines).
    Heads : Container combining one :class:`Head` and one
        :class:`Info`.
    JSiteProperty : Robust latitude/longitude parsing.
    DataType : Structured token for kind / comp / units.

    References
    ----------
    .. [1] A. G. Jones (1994). Magnetotelluric data file
       J-format, version 2.0.
    .. [2] MTNet. "J format documentation".
    """
    _repr_keys = ["station", "n"]

    def __init__(
        self,
        j_header_list: Optional[Sequence[str]] = None,
        *,
        verbose: int = 0,
        **kwargs: Any,
    ) -> None:
        super().__init__(verbose=verbose)
        self.station: str | None = None
        self.dtype: DataType | None = None
        self.n: int | None = None
        if j_header_list is not None:
            self.read(j_header_list)

    @classmethod
    def from_file(
        cls, j_fn: Union[str, Path], *, verbose: int = 0
    ) -> "Head":
        lines = list(iter_lines(j_fn, encoding=ENCODING_DEFAULT))
        head_list = _extract_first_head_list(lines)
        return cls(head_list, verbose=verbose)
    
    @classmethod
    def from_lines(
            cls,
            j_header_list: Sequence[str] | None = None,
            *,
            verbose: int = 0,
        ) -> "Head":
        if j_header_list is None:
            raise ValueError("j_header_list is required")
        seq = list(j_header_list)
        # Be tolerant: if seq already looks like a triple, use it
        if not (
                len(seq) >= 3 and RE_STATION.match(seq[0]) and (
                RE_DATATYPE_UNITS.match(seq[1]) or True
                )
            ):
            seq = _extract_first_head_list(seq)
        return cls(seq, verbose=verbose)


    def read(
        self, j_header_list: Optional[Sequence[str]] = None
    ) -> "Head":
        if j_header_list is None:
            raise ValueError("j_header_list is required")
        i = 0
        nln = len(j_header_list)
        while i < nln and (
            _is_blank(j_header_list[i]) or _is_comment(j_header_list[i]) or
            _is_info(j_header_list[i])
        ):
            i += 1
        
        m = RE_STATION.match(j_header_list[i])
        if not m:
            raise ValueError("expected station line in header list")
            
        token = m.group("station").strip()
        self.station = token.upper()
        try:
            parts = j_header_list[i].split()
            if len(parts) > 1:
                self.az_hint = float(parts[1])
        except Exception:
            self.az_hint = None
            
        # if i >= nln or not RE_STATION.match(j_header_list[i]):
        #     raise ValueError("expected station line in header list")
        # station = j_header_list[i].strip().upper()
        i += 1
        while i < nln and _is_blank(j_header_list[i]):
            i += 1
        if i >= nln:
            raise ValueError("missing data-type line")
        try:
            dtype = parse_datatype_units(j_header_list[i])
        except JParseError as exc:
            raise ValueError(
                "expected data-type line in header list"
            ) from exc
        i += 1
        while i < nln and _is_blank(j_header_list[i]):
            i += 1
        if i >= nln or not RE_NPOINTS.match(j_header_list[i]):
            raise ValueError("expected count line in header list")
        n = int(RE_NPOINTS.match(j_header_list[i]).group("n"))  # type: ignore
        self.dtype, self.n = dtype, n
        self._mark_read(True)
        return self

    def write(
        self, head_list_infos: Optional[Sequence[str]] = None
    ) -> List[str]:
        if head_list_infos is not None:
            self.read(head_list_infos)
        if self.station is None or self.dtype is None or self.n is None:
            raise ValueError("head not fully defined")
        return [f"{self.station}", f"{_fmt_dtype(self.dtype)}", f"{self.n}"]

    @property
    def kind(self) -> str | None:
        return ( None if self.dtype is None 
                else self.dtype.kind
            )
    
    @property
    def comp(self) -> str | None:
        return ( 
            None if self.dtype is None 
            else self.dtype.comp
        )
    
    @property
    def units(self) -> str | None:
        return ( 
            None if self.dtype is None 
            else self.dtype.units
        )
    
    @property
    def tensor_hint(self) -> str | None:
        return (
            None if self.dtype is None 
            else self.dtype.tensor_hint
        )
    
    @property
    def header(self) -> tuple[str, str, int] | None:
        if self.station is None or self.dtype is None:
            return None
        if self.n is None:
            return None
        return (self.station, f"{self.kind}{self.comp}", self.n)

    def __str__(self) -> str:
        return f"Head(station={self.station!r}, n={self.n})"

    def __repr__(self) -> str:
        return self.__str__()


class Heads(JComponentBase):
    r"""
    Minimal container for one :class:`Head` and one
    :class:`Info`.

    The class provides convenient accessors for station and
    site-level properties, and a small banner recorder for the
    top provenance line (``#WRITTEN BY ...``). It is intended as
    the lightest useful representation of a single J header
    section.

    Parameters
    ----------
    head : Head, optional
        Existing head to attach. If omitted, an empty
        :class:`Head` is created.
    info : Info, optional
        Existing info to attach. If omitted, an empty
        :class:`Info` is created.
    verbose : int, default=0
        Verbosity for warnings during parsing.

    Attributes
    ----------
    head : Head
        The parsed header triple.
    info : Info
        The parsed site information block.
    banner : Banner
        Parsed top comment provenance. The writer defaults the
        ``software`` field to ``PYSCAMT`` if missing.
    n : int
        ``0`` or ``1`` depending on whether a head has been
        parsed.
    station : str or None
        Shortcut to ``head.station``.
    latitude, longitude, elevation, azimuth : float or None
        Shortcuts to values provided by ``info``. ``azimuth``
        falls back to ``head.az_hint`` when the info block does
        not contain ``AZIMUTH``.

    Methods
    -------
    from_file(j_fn, *, verbose=0)
        Read the file then delegate to :meth:`read`.
    from_lines(lines, *, verbose=0)
        Build from an in-memory line sequence.
    read(text_or_lines)
        Extract info + head lists and parse both.
    write()
        Serialize banner, head and info back-to-back.

    Notes
    -----
    This class does *not* parse the subsequent data rows. It is
    focused on fast, dependency-free header discovery to support
    scanning tasks and metadata extraction.

    Examples
    --------
    >>> h = Heads.from_file('data/j/kb0-s001.txt')
    >>> h.station, h.latitude, h.software
    ('KB0001', 41.9782, 'GEOTOOLS')

    See Also
    --------
    Head : Single data-block header.
    Info : Site-level header.
    Banner : Top provenance line handler.

    References
    ----------
    .. [1] A. G. Jones (1994). Magnetotelluric data file
       J-format, version 2.0.
    """
    _repr_keys = ["n"]

    def __init__(
        self,
        head: Optional[Head] = None,
        info: Optional[Info] = None,
        *,
        verbose: int = 0,
    ) -> None:
        super().__init__(verbose=verbose)
        self.banner = Banner(verbose=verbose) 
        self.head: Head = head if head is not None else Head(
            verbose=verbose
        )
        self.info: Info = info if info is not None else Info(
            verbose=verbose
        )


    def read(
        self, text_or_lines: Union[str, Sequence[str]]
    ) -> "Heads":
        if isinstance(text_or_lines, str):
            lines = text_or_lines.splitlines()
        else:
            lines = list(text_or_lines)
        self.banner.read(lines)
        info_list = _extract_info_header_list(lines)
        head_list = _extract_first_head_list(lines)
        self.info.read(info_list)
        self.head.read(head_list)
        self._mark_read(True)
        return self

    def write(self, include_origin=False) -> List[str]:
        out: list[str] = []
        out.extend(self.banner.write(
            new=True, include_origin=include_origin))
        out.extend(self.head.write())
        out.extend(self.info.write())
        return out

    @classmethod
    def from_lines(
        cls, lines: Sequence[str], *, verbose: int = 0
        ) -> "Heads":
        inst = cls(verbose=verbose)
        inst.read(lines)
        return inst
    
    
    @classmethod
    def from_file(
        cls, j_fn: Union[str, Path], *, verbose: int = 0
        ) -> "Heads":
        lines = list(iter_lines(j_fn, encoding=ENCODING_DEFAULT))
        return cls.from_lines(lines, verbose=verbose)


    @property
    def n(self) -> int:
        return 0 if self.head.n is None else 1

    @property
    def station(self) -> str | None:
        return self.head.station

    @property
    def latitude(self) -> float | None:
        return self.info.latitude

    @property
    def longitude(self) -> float | None:
        return self.info.longitude

    @property
    def elevation(self) -> float | None:
        return self.info.elevation

    @property
    def azimuth(self) -> float | None:
        return (
            self.info.azimuth if self.info.azimuth is not None
            else self.head.az_hint
        )

    @property
    def software(self) -> str | None:
        return self.banner.software

    def __str__(self) -> str:
        return f"Heads(n={self.n})"

    def __repr__(self) -> str:
        return self.__str__()
    

class HeadMixin:
    r"""
    Mixin that provides a :class:`Head` on host classes.

    The mixin offers class-level and instance-level helpers that
    delegate to :class:`Head` while keeping a single copy of the
    header on the host.

    Attributes
    ----------
    head : Head
        Lazily created and cached on first use.

    Methods
    -------
    from_file(edi_fn)
        Class method that returns ``Head.from_file(edi_fn)``.
    read(j_header_list=None)
        Ensure a ``Head`` exists on the host, parse into it, and
        return it.
    write(head_list_infos=None)
        Serialize the host's ``Head`` to three header lines.

    Examples
    --------
    >>> class Host(HeadMixin):  # doctest: +SKIP
    ...     pass
    >>> Host.from_file('file.j')  # doctest: +SKIP
    Head(...)
    """
    
    head: Head

    @classmethod
    def from_file(cls, edi_fn: Union[str, Path]) -> Head:
        return Head.from_file(edi_fn)

    def read(
        self, j_header_list: Optional[Sequence[str]] = None
    ) -> Head:
        if not hasattr(self, "head") or self.head is None:
            self.head = Head()
        return self.head.read(j_header_list)

    def write(
        self, head_list_infos: Optional[Sequence[str]] = None
    ) -> List[str]:
        if not hasattr(self, "head") or self.head is None:
            self.head = Head()
        return self.head.write(head_list_infos)


class InfoMixin:
    r"""
    Mixin that provides an :class:`Info` on host classes.

    The mixin mirrors :class:`HeadMixin` but for the information
    block. It centralizes parsing and writing of ``>KEY=VALUE``
    lines.

    Attributes
    ----------
    info : Info
        Lazily created and cached on first use.

    Methods
    -------
    from_file(edi_fn)
        Class method that returns ``Info.from_file(edi_fn)``.
    read(j_info_list=None)
        Ensure an ``Info`` exists on the host, parse into it, and
        return it.
    write(j_info_list=None)
        Serialize the host's ``Info`` to comment and info lines.

    Examples
    --------
    >>> class Host(InfoMixin):  # doctest: +SKIP
    ...     pass
    >>> Host().read(['>LATITUDE=10'])  # doctest: +SKIP
    Info(items=1)
    """
    info: Info

    @classmethod
    def from_file(cls, edi_fn: Union[str, Path]) -> Info:
        return Info.from_file(edi_fn)

    def read(
        self, j_info_list: Optional[Sequence[str]] = None
    ) -> Info:
        if not hasattr(self, "info") or self.info is None:
            self.info = Info()
        return self.info.read(j_info_list)

    def write(
        self, j_info_list: Optional[Sequence[str]] = None
    ) -> List[str]:
        if not hasattr(self, "info") or self.info is None:
            self.info = Info()
        return self.info.write(j_info_list)


def _extract_info_header_list(lines: Sequence[str]) -> List[str]:
    header: List[str] = []
    for ln in lines:
        if _is_comment(ln) or _is_info(ln) or _is_blank(ln):
            header.append(ln)
            continue
        break
    return header



def _extract_first_head_list(lines: Sequence[str]) -> List[str]:
    i, n = 0, len(lines)
    while i < n and (
        _is_comment(lines[i]) or _is_info(lines[i]) or _is_blank(lines[i])
    ):
        i += 1
    while i < n and not RE_STATION.match(lines[i]):
        i += 1
    if i >= n:
        raise ValueError("no head triple found")
    s = i

    # first non-blank after station
    j = i + 1
    while j < n and _is_blank(lines[j]):
        j += 1
    if j >= n:
        raise ValueError("missing header after station")

    # case A: standard order (dtype then count)
    try:
        parse_datatype_units(lines[j])
        k = j + 1
        while k < n and _is_blank(lines[k]):
            k += 1
        if k >= n or not RE_NPOINTS.match(lines[k]):
            raise ValueError("missing count after data-type")
        return [lines[s], lines[j], lines[k]]
    except JParseError:
        pass  # not a dtype here

    # case B: count first, then dtype later
    if RE_NPOINTS.match(lines[j]):
        k = j
        # scan forward for first dtype before next station/EOF
        j = k + 1
        while j < n:
            if RE_STATION.match(lines[j]):
                break
            if _is_blank(lines[j]):
                j += 1
                continue
            try:
                parse_datatype_units(lines[j])
                return [lines[s], lines[j], lines[k]]
            except JParseError:
                # not a dtype; likely data row — keep scanning a bit
                j += 1
                continue
        raise ValueError("no data-type found after count")
    else:
        raise ValueError("bad header after station")


def _is_comment(s: str) -> bool:
    return bool(RE_COMMENT.match(s))


def _is_blank(s: str) -> bool:
    return bool(RE_BLANK.match(s))


def _is_info(s: str) -> bool:
    return bool(RE_INFO.match(s))


def _fmt_dtype(dt: DataType) -> str:
    base = f"{dt.kind}{dt.comp}"
    if dt.units:
        for k, v in UNITS_CANONICAL.items():
            if v == dt.units:
                return f"{base} {k}"
    return base
