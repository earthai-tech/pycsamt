# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

import copy
import math
import re
from collections import namedtuple
from collections.abc import Iterator, Sequence
from os import PathLike
from pathlib import Path
from typing import (
    Any,
    Callable,
)

import numpy as np

from ..core.base import ensure_station
from ..seg.collection import EDICollection
from ..seg.edi import EDIFile
from ..seg.heads import Head

__all__ = [
    "is_pathlike",
    "is_edi_file",
    "is_edi_collection",
    "iter_edifiles",
    "as_edicollection",
    "station_name",
    "set_station_name",
    "get_coords",
    "set_coords",
    "maybe_copy",
    "apply_inplace",
    "get_freq",
    "freq_match",
    "freq_select",
    "select_by_name",
    "match_name",
    "wrap_azimuth",
    "deg_to_mrad",
    "mrad_to_deg",
]

_Coord = namedtuple("_Coord", ["lat", "lon", "elev"])


def is_pathlike(x: Any) -> bool:
    r"""
    Return True if *x* looks like a filesystem path.

    Parameters
    ----------
    x : Any
        Object to test.

    Returns
    -------
    bool
        True when *x* is a `str` or `pathlib.Path`, otherwise
        False.

    Notes
    -----
    This helper is intentionally conservative. Containers that
    happen to be iterable but are not paths will return False.

    Examples
    --------
    >>> from pathlib import Path
    >>> from pycsamt.site.utils import is_pathlike
    >>> is_pathlike("data/line01.edi")
    True
    >>> is_pathlike(Path("data/line01.edi"))
    True
    >>> is_pathlike(["a", "b"])
    False

    See Also
    --------
    pathlib.Path
    """

    return isinstance(x, (str, Path))


def is_edi_file(x: Any) -> bool:
    r"""
    Heuristically detect an EDI file object.

    Parameters
    ----------
    x : Any
        Object to test.

    Returns
    -------
    bool
        True if *x* exposes both ``get_section`` and an attribute
        named ``Z`` (the impedance section). Otherwise False.

    Notes
    -----
    This check is duck-typed and does not import the backend
    class. It is compatible with multiple EDI implementations.

    Examples
    --------
    >>> class Dummy:
    ...     def get_section(self, *_): ...
    ...     Z = object()
    ...
    >>> from pycsamt.site.utils import is_edi_file
    >>> is_edi_file(Dummy())
    True
    >>> is_edi_file(object())
    False

    See Also
    --------
    is_edi_collection, iter_edifiles
    """

    if x is None:
        return False
    return hasattr(x, "get_section") and hasattr(x, "Z")


def is_edi_collection(x: Any) -> bool:
    r"""
    Detect whether *x* is an iterable collection of EDI objects.

    Parameters
    ----------
    x : Any
        Candidate collection. May be an `EDICollection` instance,
        or any iterable containing EDI-like objects.

    Returns
    -------
    bool
        True if *x* looks like a collection of EDI files.

    Notes
    -----
    The check is non-destructive. It peeks at the first element
    of the iterable to decide. Pathlike objects are not treated
    as iterables here, even though `str` is iterable.

    Examples
    --------
    >>> from pycsamt.site.utils import is_edi_collection, is_edi_file
    >>> class E:
    ...     def get_section(self, *_): ...
    ...     Z = object()
    ...
    >>> is_edi_collection([E(), E()])
    True
    >>> is_edi_collection("folder/*.edi")
    False
    >>> is_edi_collection(None)
    False
    """

    if x is None:
        return False
    if isinstance(x, EDICollection):
        return True
    if hasattr(x, "__iter__") and not is_pathlike(x):
        try:
            it = iter(x)  # type: ignore
            first = next(it, None)
        except Exception:
            return False
        return is_edi_file(first)
    return False


def iter_edifiles(edic: Any) -> Iterator[EDIFile]:
    r"""
    Iterate over EDI files contained in *edic*.

    Parameters
    ----------
    edic : Any
        Single EDI object or an iterable of EDI objects. Pathlike
        inputs are ignored (yield nothing).

    Yields
    ------
    pycsamt.seg.edi.EDIFile
        Each detected EDI-like object from *edic*.

    Notes
    -----
    This is a safe iterator. Non EDI elements are skipped. A
    single EDI object will be yielded once.

    Examples
    --------
    >>> from pycsamt.site.utils import iter_edifiles
    >>> class E:
    ...     def get_section(self, *_): ...
    ...     Z = object()
    ...
    >>> list(iter_edifiles(E()))
    [<E object at ...>]
    >>> lst = [E(), object(), E()]
    >>> len(list(iter_edifiles(lst)))
    2

    See Also
    --------
    is_edi_file, is_edi_collection, as_edicollection
    """

    if is_edi_file(edic):
        yield edic  # type: ignore
        return
    if hasattr(edic, "__iter__") and not is_pathlike(edic):
        for it in edic:  # type: ignore
            if is_edi_file(it):
                yield it  # type: ignore


def _is_pathlike(obj: Any) -> bool:
    return isinstance(obj, (str, bytes, Path, PathLike))


def _is_seq_of_pathlike(x: Any) -> bool:
    if isinstance(x, (str, bytes, Path, PathLike)):
        return False
    try:
        xs = list(x)
    except Exception:
        return False
    return len(xs) > 0 and all(_is_pathlike(t) for t in xs)


def as_edicollection(
    edic: Any,
    *,
    # Discovery / parsing knobs (used when edic is path-like)
    recursive: bool = True,
    strict: bool = False,
    on_dup: str = "replace",
    verbose: int = 0,
) -> EDICollection | None:
    r"""
    Coerce *edic* into an :class:`EDICollection` when possible.

    Parameters
    ----------
    edic : Any
        Single/sequence of EDI objects, an existing EDICollection,
        or path-like (file/dir/glob) input.
    recursive : bool, default ``True``
        When ``edic`` is path-like (or a sequence of), recurse into
        directories while discovering files.
    strict : bool, default ``False``
        If ``True`` and ``edic`` is path-like, propagate read errors;
        otherwise errors are captured in the parser.
    on_dup : {'replace', 'keep'}, default ``'replace'``
        Duplicate policy during path-like discovery. See
        :meth:`EDICollection.from_sources`.
    verbose : int, default ``0``
        Verbosity forwarded to underlying readers.

    Returns
    -------
    EDICollection or None
        A new collection containing the input items, or ``None`` if
        nothing EDI-like was found.

    Notes
    -----
    * If *edic* is path-like → use
      :meth:`EDICollection.from_sources(...)`.
    * If *edic* is already an EDICollection → return it.
    * Else → collect EDI-like items via ``iter_edifiles(edic)`` and
      build a collection.

    Examples
    --------
    >>> from pycsamt.site.utils import as_edicollection
    >>> coll = as_edicollection("data/*.edi", recursive=True)
    >>> coll is not None
    >>>
    >>> class E:
    ...     def get_section(self, *_): ...
    ...     Z = object()
    ...
    >>> coll = as_edicollection([E(), E()])
    >>> coll is not None
    True
    >>> as_edicollection([]) is None
    True

    See Also
    --------
    iter_edifiles, is_edi_collection

    References
    ----------
    .. [1] Robust factory patterns for heterogeneous inputs.

    True
    """

    # 1) Path-like discovery path
    if _is_pathlike(edic) or _is_seq_of_pathlike(edic):
        return EDICollection.from_sources(
            edic,
            recursive=recursive,
            strict=strict,
            on_dup=on_dup,
            verbose=verbose,
        )

    # 2) Already a collection
    if isinstance(edic, EDICollection):
        return edic

    # 3) Try to iterate EDI-like objects
    items = list(iter_edifiles(edic))
    if not items:
        return None

    try:
        return EDICollection(items=items, verbose=verbose)  # type: ignore
    except TypeError:
        return EDICollection(items, verbose=verbose)  # type: ignore


def station_name(ed: Any) -> str:
    r"""
    Return a best-effort station name for an EDI-like object.

    The lookup order is:

    1. Attribute ``station`` on the object itself (if truthy).
    2. Header field ``dataid`` (via an internal head accessor).
    3. Fallback attributes on the object: ``name``, ``site``,
       or ``dataid`` (first non-empty one).
    4. Empty string when nothing suitable is found.

    Parameters
    ----------
    ed : Any
        EDI-like object exposing either a header section or plain
        attributes for station naming.

    Returns
    -------
    str
        The resolved station name, possibly an empty string.

    Notes
    -----
    This helper is intentionally tolerant so it can work across
    different EDI backends. It prefers explicit ``station`` on
    the object, then the header ``dataid`` which is commonly used
    as the unique site identifier.

    Examples
    --------
    >>> from types import SimpleNamespace
    >>> from pycsamt.site.utils import station_name
    >>> ed = SimpleNamespace(station="E01")
    >>> station_name(ed)
    'E01'
    >>> ed = SimpleNamespace(name="Alt01")
    >>> station_name(ed)
    'Alt01'

    See Also
    --------
    set_station_name : Synchronously update object and header names.
    get_coords : Read latitude, longitude and elevation from header.
    """

    if hasattr(ed, "station") and ed.station:
        return str(ed.station)
    h = _get_head(ed)
    if h is not None and getattr(h, "dataid", None):
        return str(h.dataid)
    for k in ("name", "site", "dataid"):
        v = getattr(ed, k, None)
        if v:
            return str(v)
    return ""


def set_station_name(
    ed: Any,
    name: str | None = None,
    *,
    station_id: str | int | None = None,
    policy: Any = None,
    inplace: bool = True,
) -> Any:
    r"""
    Set the station name on an EDI-like object and its header.

    This function updates multiple common identifiers so that the
    object-level name and the header ``dataid`` stay in sync.

    Parameters
    ----------
    ed : Any
        EDI-like object.
    name : str, optional
        New station name. If not provided, it is derived from
        ``station_id`` or a ``policy`` callable.
    station_id : str or int, optional
        Alternate identifier used to build the new name when
        ``name`` is not explicitly given.
    policy : Any, optional
        A callable or policy object that can convert an input
        identifier to a station string, e.g. ``lambda s: s.upper()``.
    inplace : bool, default True
        If True, mutate *ed* in place; otherwise return a modified
        copy when possible.

    Returns
    -------
    Any
        The same object (when ``inplace=True``) or a best-effort
        copy with the new identifiers applied.

    Notes
    -----
    Internally, the header object is created if missing. Both
    ``ed.station`` and ``head.dataid`` are written when available.
    Errors from non-writable fields are suppressed.

    Examples
    --------
    >>> from types import SimpleNamespace
    >>> from pycsamt.site.utils import set_station_name, station_name
    >>> ed = SimpleNamespace()
    >>> set_station_name(ed, name="X01")
    ... # doctest: +ELLIPSIS
    <...>
    >>> station_name(ed)
    'X01'

    See Also
    --------
    station_name : Read the current station name.
    maybe_copy : Lightweight deep-copy helper used when not in place.
    """

    tgt = ed if inplace else maybe_copy(ed)
    nm = ensure_station(name, station_id, policy=policy)
    try:
        tgt.station = nm
    except Exception:
        pass
    try:
        h = _ensure_head(tgt)
        if h is not None:
            h.dataid = nm
    except Exception:
        pass
    return tgt


def get_coords(ed: Any) -> _Coord:
    r"""
    Read latitude, longitude and elevation from an EDI header.

    Parameters
    ----------
    ed : Any
        EDI-like object holding a header section.

    Returns
    -------
    _Coord
        A small tuple-like object with fields ``lat``, ``lon`` and
        ``elev``. Missing values are returned as ``nan``.

    Notes
    -----
    Both ``lon`` and legacy ``long`` header attribute names are
    supported. If the header is missing or fields cannot be read,
    all values are ``nan``.

    Examples
    --------
    >>> from types import SimpleNamespace
    >>> from pycsamt.site.utils import get_coords
    >>> head = SimpleNamespace(lat=10.0, long=20.0, elev=300.0)
    >>> ed = SimpleNamespace(get_section=lambda *_: head)
    >>> c = get_coords(ed)
    >>> (c.lat, c.lon, c.elev)
    (10.0, 20.0, 300.0)

    See Also
    --------
    set_coords : Write latitude, longitude and elevation to header.
    """

    h = _get_head(ed)
    la = getattr(h, "lat", float("nan")) if h else float("nan")
    lo = getattr(h, "long", float("nan")) if h else float("nan")
    ev = getattr(h, "elev", float("nan")) if h else float("nan")
    try:
        return _Coord(float(la), float(lo), float(ev))
    except Exception:
        return _Coord(float("nan"), float("nan"), float("nan"))


def set_coords(
    ed: Any,
    *,
    lat: float | None = None,
    lon: float | None = None,
    elev: float | None = None,
    inplace: bool = True,
) -> Any:
    r"""
    Write latitude, longitude and/or elevation into the header.

    Parameters
    ----------
    ed : Any
        EDI-like object.
    lat : float, optional
        Latitude in degrees. When None, the value is left unchanged.
    lon : float, optional
        Longitude in degrees. When None, the value is left unchanged.
    elev : float, optional
        Elevation in meters. When None, the value is left unchanged.
    inplace : bool, default True
        If True, mutate *ed* in place; otherwise apply the update to
        a copy and return it.

    Returns
    -------
    Any
        The updated object (in place) or a best-effort copy.

    Notes
    -----
    Both ``lon`` and ``long`` attribute names are supported on the
    header. If neither exists, a small adapter header may be created
    and attached back to the object when possible. Failures to write
    individual fields are ignored to maximize portability.

    Examples
    --------
    >>> from types import SimpleNamespace
    >>> from pycsamt.site.utils import set_coords, get_coords
    >>> head = SimpleNamespace(lat=0.0, long=0.0, elev=0.0)
    >>> ed = SimpleNamespace(get_section=lambda *_: head)
    >>> set_coords(ed, lat=12.5, lon=1.2, elev=350.0)
    ... # doctest: +ELLIPSIS
    <...>
    >>> tuple(map(float, get_coords(ed)))
    (12.5, 1.2, 350.0)

    See Also
    --------
    get_coords : Read header coordinates.
    maybe_copy : Helper used when not writing in place.
    """

    tgt = ed if inplace else maybe_copy(ed)
    h = _ensure_head(tgt)
    if h is None:
        return tgt
    if lat is not None:
        try:
            h.lat = float(lat)
        except Exception:
            pass
    if lon is not None:
        ok = False
        try:
            h.lon = float(lon)
            ok = True
        except Exception:
            pass
        try:
            h.long = float(lon)
            ok = True
        except Exception:
            pass
        if not ok:
            # replace with a tiny adapter if needed
            try:
                nh = type("Head", (), {})()
                for k in ("lat", "long", "elev", "dataid"):
                    if hasattr(h, k):
                        setattr(nh, k, getattr(h, k))
                nh.lon = float(lon)
                tgt.set_section("head", nh)  # type: ignore
            except Exception:
                pass
    if elev is not None:
        try:
            h.elev = float(elev)
        except Exception:
            pass
    return tgt


def maybe_copy(x: Any) -> Any:
    r"""
    Attempt a deep copy of *x*, falling back to identity.

    Parameters
    ----------
    x : Any
        Object to copy.

    Returns
    -------
    Any
        A deep-copied object when possible; otherwise *x* itself.

    Notes
    -----
    Some EDI backends contain objects that are not deepcopyable.
    This helper catches such cases and simply returns the original
    object to avoid raising.

    Examples
    --------
    >>> from pycsamt.site.utils import maybe_copy
    >>> maybe_copy({"a": [1, 2]})
    {'a': [1, 2]}
    """

    try:
        return copy.deepcopy(x)
    except Exception:
        return x


def apply_inplace(
    x: Any,
    fn: Callable[[Any], Any],
    *,
    inplace: bool = False,
) -> Any:
    r"""
    Apply a function to an object with optional in-place semantics.

    Parameters
    ----------
    x : Any
        Input object.
    fn : callable
        Function ``fn(obj) -> obj`` that mutates and/or returns the
        object.
    inplace : bool, default False
        If True, call ``fn`` directly on *x*. Otherwise a deep copy
        is attempted first.

    Returns
    -------
    Any
        The result of applying ``fn`` to the chosen object.

    Notes
    -----
    This utility centralizes a common pattern used across the
    package: honor an ``inplace`` flag while being resilient to
    objects that cannot be deep-copied.

    Examples
    --------
    >>> from pycsamt.site.utils import apply_inplace
    >>> def inc(d):
    ...     d["n"] = d.get("n", 0) + 1
    ...     return d
    ...
    >>> data = {}
    >>> r1 = apply_inplace(data, inc, inplace=False)
    >>> r1 is data
    False
    >>> r2 = apply_inplace(data, inc, inplace=True)
    >>> r2 is data
    True
    """

    if inplace:
        return fn(x)
    y = maybe_copy(x)
    return fn(y)


def get_freq(ed: Any) -> np.ndarray:
    r"""
    Return the frequency vector from the impedance section.

    Parameters
    ----------
    ed : Any
        EDI-like object exposing a ``Z`` section with either a
        public ``freq`` array or a private ``_freq`` array.

    Returns
    -------
    numpy.ndarray
        A 1-D float array sorted in ascending order. Returns an
        empty array when no frequency information is available.

    Notes
    -----
    This is a read-only accessor. When both ``freq`` and ``_freq``
    are absent or unreadable, an empty array is returned. Non-finite
    values are preserved; only ordering is enforced.

    Examples
    --------
    >>> import numpy as np
    >>> from types import SimpleNamespace
    >>> from pycsamt.site.utils import get_freq
    >>> Z = SimpleNamespace(freq=np.array([2.0, 1.0]))
    >>> ed = SimpleNamespace(Z=Z)
    >>> get_freq(ed)
    array([1., 2..])

    See Also
    --------
    freq_match : Find indices of exact matches within a tolerance.
    freq_select : Build index selections from ranges or scalars.
    """

    z = getattr(ed, "Z", None)
    if z is None:
        return np.asarray([], float)
    for k in ("freq", "_freq"):
        try:
            f = getattr(z, k)
            if f is None:
                continue
            a = np.asarray(f, float).ravel()
            if a.size:
                a = a[np.argsort(a)]
            return a
        except Exception:
            pass
    return np.asarray([], float)


def freq_match(
    f: np.ndarray,
    target: float | Sequence[float],
    *,
    tol: float = 1e-6,
) -> np.ndarray:
    r"""
    Indices where frequency values match targets within tolerance.

    Parameters
    ----------
    f : numpy.ndarray
        Frequency array to search. It is cast to ``float``.
    target : float or Sequence[float]
        One or more target frequencies to match against ``f``.
    tol : float, default 1e-6
        Absolute tolerance for a match. A value ``v`` in ``f`` is
        considered a match to ``t`` when
        :math:`|v - t| \\leq \\mathrm{tol}`.

    Returns
    -------
    numpy.ndarray
        Sorted integer indices into ``f`` where matches occur. May be
        empty when no value matches within ``tol``.

    Notes
    -----
    Non-finite values in ``f`` never match. Duplicates in ``f`` are
    all returned if they fall within the tolerance window.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.site.utils import freq_match
    >>> f = np.array([10.0, 20.0, 20.0000004, 30.0])
    >>> freq_match(f, 20.0, tol=5e-7)
    array([1, 2])
    >>> freq_match(f, [5.0, 30.0])
    array([3])

    See Also
    --------
    freq_select : More general selection by ranges or slices.
    get_freq : Read the frequency vector from an EDI object.
    """

    a = np.asarray(f, float)
    t = np.atleast_1d(np.asarray(target, float))
    mask = np.zeros(a.shape, dtype=bool)
    for v in t:
        mask |= np.isfinite(a) & (np.abs(a - v) <= tol)
    idx = np.where(mask)[0]
    return idx.astype(int)


def freq_select(
    f: np.ndarray,
    sel: slice | float | Sequence[float] | tuple[float, float],
    *,
    tol: float = 1e-6,
) -> np.ndarray:
    r"""
    Build an index selection from ranges, slices or scalars.

    Parameters
    ----------
    f : numpy.ndarray
        Frequency vector to index. Cast to ``float``.
    sel : slice or float or Sequence[float] or tuple(float, float)
        Selection specifier:
          * ``slice(lo, hi)`` inclusive on both ends (with ``tol``).
          * ``(lo, hi)`` tuple uses the same semantics as above.
          * single ``float`` selects exact matches (with ``tol``).
          * sequence of floats selects those exact targets.
    tol : float, default 1e-6
        Absolute tolerance used for end inclusions and exact
        matching.

    Returns
    -------
    numpy.ndarray
        Sorted integer indices into ``f`` that satisfy the selection.
        Returns an empty array when nothing matches.

    Notes
    -----
    Bounds are treated as inclusive with a tolerance, i.e.
    :math:`\\mathrm{lo} - \\mathrm{tol} \\le f \\le
    \\mathrm{hi} + \\mathrm{tol}`. For exact targets, matching is
    performed by :func:`freq_match`.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.site.utils import freq_select
    >>> f = np.array([1.0, 10.0, 100.0])
    >>> freq_select(f, slice(5.0, 100.0))
    array([1, 2])
    >>> freq_select(f, (0.9, 1.1))
    array([0])
    >>> freq_select(f, [10.0, 999.0])
    array([1])

    See Also
    --------
    freq_match : Matching by exact targets with tolerance.
    get_freq : Frequency accessor for EDI objects.
    """

    a = np.asarray(f, float)
    if isinstance(sel, slice):
        lo = float(sel.start) if sel.start else -np.inf
        hi = float(sel.stop) if sel.stop else np.inf
        m = (a >= lo - tol) & (a <= hi + tol)
        return np.where(m)[0].astype(int)
    if isinstance(sel, (tuple, list)) and len(sel) == 2:
        lo, hi = float(sel[0]), float(sel[1])
        m = (a >= lo - tol) & (a <= hi + tol)
        return np.where(m)[0].astype(int)
    if isinstance(sel, (float, int)):
        return freq_match(a, float(sel), tol=tol)
    try:
        return freq_match(a, list(sel), tol=tol)  # type: ignore
    except Exception:
        return np.asarray([], int)


def match_name(
    pat: str | re.Pattern[str] | Callable[[str], bool],
    name: str,
) -> bool:
    r"""
    Case-insensitive station-name matcher with flexible patterns.

    Parameters
    ----------
    pat : str or Pattern or callable
        Matching strategy. If callable, it is invoked as
        ``bool(pat(name))``. If a compiled regex, ``pat.search`` is
        used. If a string, the function tries in order:
        glob-like wildcards (``*``, ``?``, ``[]``), regex-looking
        strings, and finally a case-insensitive literal compare.
    name : str
        The candidate station name.

    Returns
    -------
    bool
        ``True`` if the pattern matches ``name``, ``False`` on any
        error or no match.

    Notes
    -----
    Glob-like patterns are translated to regular expressions with
    ``re.IGNORECASE``. Strings that contain typical regex meta-
    characters (e.g. ``. ^ $ + | ( ) \\``) are treated as regexes.

    Examples
    --------
    >>> import re
    >>> from pycsamt.site.utils import match_name
    >>> match_name("E*", "E01")
    True
    >>> match_name(re.compile(r"^X\\d+$"), "X123")
    True
    >>> match_name(lambda s: s.endswith("99"), "A99")
    True

    See Also
    --------
    select_by_name : Collect items from an iterable using a pattern.
    """

    try:
        if callable(pat):
            return bool(pat(name))
        if isinstance(pat, re.Pattern):
            return bool(pat.search(name))
        # glob-like?
        if any(c in pat for c in ("*", "?", "[")):
            rx = re.compile(
                "^"
                + re.escape(pat).replace("\\*", ".*").replace("\\?", ".")
                + "$",
                flags=re.IGNORECASE,
            )
            return bool(rx.match(name))
        # regex-looking string?
        if any(c in pat for c in (".", "^", "$", "+", "|", "(", ")", "\\")):
            rx = re.compile(pat, flags=re.IGNORECASE)
            return bool(rx.search(name))
        return name.upper() == str(pat).upper()
    except Exception:
        return False


def select_by_name(
    edic: Any,
    pat: str | re.Pattern[str] | Callable[[str], bool],
) -> list[EDIFile]:
    r"""
    Collect EDI-like objects whose station names match a pattern.

    Parameters
    ----------
    edic : Any
        Iterable of EDI-like objects. It may also be a single EDI
        object, in which case a list with zero or one element is
        returned depending on the match.
    pat : str or Pattern or callable
        Pattern passed to :func:`match_name`.

    Returns
    -------
    list of pycsamt.seg.edi.EDIFile
        Items whose names match ``pat`` in iteration order.

    Examples
    --------
    >>> from types import SimpleNamespace
    >>> from pycsamt.site.utils import select_by_name
    >>> eds = [SimpleNamespace(station="E01"),
    ...        SimpleNamespace(station="N10")]
    >>> keep = select_by_name(eds, "E*")
    >>> [getattr(x, "station") for x in keep]
    ['E01']

    See Also
    --------
    match_name : Flexible single-name predicate.
    iter_edifiles : Safe iterator over EDI-like inputs.
    """

    out: list[EDIFile] = []
    for ed in iter_edifiles(edic):
        nm = station_name(ed)
        if match_name(pat, nm):
            out.append(ed)
    return out


def wrap_azimuth(az: float) -> float:
    r"""
    Wrap an azimuth angle in degrees to the half-open range [0, 360).

    Parameters
    ----------
    az : float
        Angle in degrees, possibly outside the canonical range.

    Returns
    -------
    float
        Angle in degrees in the interval ``[0, 360)``.

    Notes
    -----
    The operation is equivalent to ``az % 360`` with care for
    negative inputs.

    Examples
    --------
    >>> from pycsamt.site.utils import wrap_azimuth
    >>> wrap_azimuth(370.0)
    10.0
    >>> wrap_azimuth(-10.0)
    350.0
    """
    a = float(az) % 360.0
    return a if a >= 0 else a + 360.0


def deg_to_mrad(x: float | np.ndarray) -> np.ndarray:
    r"""
    Convert degrees to milliradians.

    Parameters
    ----------
    x : float or numpy.ndarray
        Angle(s) in degrees.

    Returns
    -------
    numpy.ndarray
        Converted values in milliradians.

    Notes
    -----
    The conversion uses:
    :math:`\\mathrm{mrad} = \\mathrm{deg} \\times
    \\pi / 180 \\times 1000`.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.site.utils import deg_to_mrad
    >>> float(deg_to_mrad(180.0))
    3141.592653589793
    >>> deg_to_mrad(np.array([0.0, 90.0]))
    array([   0.        , 1570.79632679])
    """

    a = np.asarray(x, float)
    return a * (math.pi / 180.0) * 1000.0


def mrad_to_deg(x: float | np.ndarray) -> np.ndarray:
    r"""
    Convert milliradians to degrees.

    Parameters
    ----------
    x : float or numpy.ndarray
        Angle(s) in milliradians.

    Returns
    -------
    numpy.ndarray
        Converted values in degrees.

    Notes
    -----
    The conversion uses:
    :math:`\\mathrm{deg} = \\mathrm{mrad} \\times
    180 / (\\pi \\times 1000)`.

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.site.utils import mrad_to_deg
    >>> float(mrad_to_deg(3141.592653589793))
    180.0
    >>> mrad_to_deg(np.array([0.0, 1570.79632679]))
    array([ 0., 90.])
    """

    a = np.asarray(x, float)
    return a * (180.0 / math.pi) / 1000.0


def _get_head(ed: EDIFile) -> Any:
    try:
        return ed.get_section("head")  # type: ignore
    except Exception:
        return None


def _ensure_head(ed: EDIFile) -> Any:
    h = _get_head(ed)
    if h is not None:
        return h
    try:
        h = Head()
    except Exception:

        class _H:
            pass

        h = _H()
    try:
        ed.set_section("head", h)  # type: ignore
    except Exception:
        try:
            ed.Head = h
        except Exception:
            pass
    return h
