# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

import math
import re
from collections.abc import Iterable
from typing import Any, Callable

import numpy as np

from .utils import (
    get_coords,
    get_freq,
    iter_edifiles,
    match_name,
    station_name,
)

__all__ = [
    "by_names",
    "by_index",
    "by_chainage",
    "by_freq",
    "by_bbox",
    "by_predicate",
    "keep_finite_z",
    "mask_large_phase_err",
    "drop_empty",
]


def by_names(
    sites: Any,
    patterns: Iterable[Any] | Any,
    *,
    case: bool = False,
):
    r"""
    Select sites by matching station names against one or more
    patterns.

    This is a flexible name-based selector that accepts several
    pattern types:

    * string with optional glob-like wildcards ``*`` and ``?``
    * compiled regular expression (``re.Pattern``)
    * callable ``fn(name)->bool``
    * iterable of any mix of the above

    A site is kept if **any** pattern matches its station name.
    Matching is stable: the relative order of kept sites is the
    same as in the input.

    Parameters
    ----------
    sites : Any
        A :class:`~pycsamt.site.base.Sites` instance, an
        ``EDICollection``, a sequence of ``EDIFile`` objects,
        or any iterable yielding EDI-like objects.
    patterns : Iterable[Any] or Any
        One pattern or an iterable of patterns. See the list
        above for supported pattern types.
    case : bool, optional
        If ``True``, perform case-sensitive matching. If
        ``False`` (default) names and string patterns are
        upper-cased before comparison.

    Returns
    -------
    pycsamt.site.base.Sites
        A new ``Sites`` wrapper containing only the matched
        EDI items. The original container is not modified.

    Notes
    -----
    String patterns support a minimal glob syntax. ``*`` matches
    any sequence (possibly empty) and ``?`` matches any single
    character. If you need full regular expressions, pass a
    compiled ``re.Pattern``.

    When multiple patterns are given, the match is an OR over all
    patterns. Matching uses the station name as returned by
    ``station_name(ed)`` which reflects header normalization.

    Examples
    --------
    >>> from pycsamt.site.base import Sites
    >>> from pycsamt.site.selection import by_names
    >>> sites = Sites([e1, e2, e3])  # EDIFile objects
    >>> out = by_names(sites, "K*")  # glob: all names starting K
    >>> [s.name for s in out]
    ['K01', 'K02']

    >>> import re
    >>> rx = re.compile(r"^S0[1-3]$")
    >>> out = by_names(sites, rx)
    >>> [s.name for s in out]
    ['S01', 'S02', 'S03']

    >>> out = by_names(sites, lambda n: n.endswith("A"))
    >>> [s.name for s in out]
    ['X1A', 'X2A']

    See Also
    --------
    pycsamt.site.selection.by_index :
        Select by numeric positions.
    pycsamt.site.selection.by_predicate :
        Keep sites for which a boolean predicate returns True.
    pycsamt.site.selection.by_freq :
        Keep sites that contain data within a frequency window.

    References
    ----------
    .. [1] Python re module documentation.
    .. [2] Unix shell-style wildcards (glob) convention.
    """

    s = _to_sites(sites)
    pats = patterns if isinstance(patterns, (list, tuple)) else [patterns]
    keep = []
    for ed in iter_edifiles(s.edic):
        nm = station_name(ed)
        ok = any(_name_matches(nm, p, case) for p in pats)
        if ok:
            keep.append(ed)
    return _new_sites(s, keep)


def by_index(sites: Any, indices: Iterable[int] | int):
    r"""
    Select sites by zero-based numeric indices, supporting negative
    indices.

    Indices are normalized exactly like Python sequence indexing:
    ``-1`` addresses the last item, ``-2`` the one before last,
    and so on. Out-of-range or non-integer indices are ignored.
    The resulting subset preserves the original ordering of the
    selected items, not the order in which indices are provided.

    Parameters
    ----------
    sites : Any
        A :class:`~pycsamt.site.base.Sites` instance, an
        ``EDICollection``, a sequence of ``EDIFile`` objects,
        or any iterable yielding EDI-like objects.
    indices : Iterable[int] or int
        One index or an iterable of indices. Negative indices
        are supported and mapped to their Python equivalents.

    Returns
    -------
    pycsamt.site.base.Sites
        A new ``Sites`` wrapper containing only items at the
        requested positions. If no valid indices remain after
        normalization, an empty ``Sites`` is returned.

    Notes
    -----
    Duplicate indices are de-duplicated in the output since the
    selection is implemented as a membership test over the set
    of normalized indices. The relative order of kept items is
    the same as in the original sequence.

    Examples
    --------
    >>> from pycsamt.site.base import Sites
    >>> from pycsamt.site.selection import by_index
    >>> sites = Sites([e1, e2, e3])  # names: A, B, C
    >>> out = by_index(sites, [0, -1])   # first and last
    >>> [s.name for s in out]
    ['A', 'C']

    >>> out = by_index(sites, 1)  # single integer
    >>> [s.name for s in out]
    ['B']

    >>> out = by_index(sites, [10, -10])  # both invalid -> empty
    >>> len(out)
    0

    See Also
    --------
    pycsamt.site.selection.by_names :
        Name-based matching using strings, regex, or callables.
    pycsamt.site.selection.by_chainage :
        Select by stored chainage range when available.
    pycsamt.site.base.Sites.by_index :
        Random access to a single site by position.
    """

    s = _to_sites(sites)
    idxs = indices if isinstance(indices, (list, tuple)) else [indices]
    # coerce to valid indices with Python-like negative support
    all_items = list(iter_edifiles(s.edic))
    n = len(all_items)
    norm: list[int] = []
    for i in idxs:
        if not isinstance(i, int):
            continue
        j = i if i >= 0 else (n + i)
        if 0 <= j < n:
            norm.append(j)
    if not norm:
        return _new_sites(s, [])
    keep = [e for k, e in enumerate(all_items) if k in set(norm)]
    return _new_sites(s, keep)


def by_chainage(sites: Any, smin: float, smax: float):
    r"""
    Select sites whose stored chainage falls within a closed
    interval.

    This helper reads the chainage value first from the EDI
    ``HEAD`` section (``head.chainage``) and, if missing, from a
    top-level attribute ``edi.chainage``. Sites for which a
    numeric chainage cannot be determined are silently skipped.

    Parameters
    ----------
    sites : Any
        A :class:`~pycsamt.site.base.Sites` instance, an
        ``EDICollection``, a sequence of ``EDIFile`` objects, or
        any iterable yielding EDI-like items.
    smin : float
        Minimum chainage (inclusive).
    smax : float
        Maximum chainage (inclusive).

    Returns
    -------
    pycsamt.site.base.Sites
        A new ``Sites`` wrapper containing only the EDI items
        whose chainage :math:`c` satisfies
        :math:`smin \\le c \\le smax`.

    Notes
    -----
    Chainage is a linear reference commonly used along profiles
    or lines, typically measured in meters from a chosen origin.
    If chainage is not present on a site, that site is excluded.
    The original order of sites is preserved among the kept
    items.

    Examples
    --------
    >>> from pycsamt.site.base import Sites
    >>> from pycsamt.site.selection import by_chainage
    >>> s = Sites([e1, e2, e3])  # EDIFile objects
    >>> out = by_chainage(s, smin=100.0, smax=300.0)
    >>> [t.name for t in out]
    ['L02', 'L03']

    See Also
    --------
    pycsamt.site.selection.by_index :
        Select by zero-based positions with negative support.
    pycsamt.site.selection.by_names :
        Select by station names using glob, regex, or callables.
    pycsamt.site.selection.by_freq :
        Keep sites that contain data in a frequency window.
    pycsamt.site.base.Sites.to_profile :
        Build a profile or ordered view along a line.

    References
    ----------
    .. [1] Linear referencing and chainage in civil engineering.
    """

    s = _to_sites(sites)
    keep = []
    for ed in iter_edifiles(s.edic):
        h = getattr(ed, "get_section", lambda *_: None)("head")
        ch = getattr(h, "chainage", None) if h else None
        if ch is None:
            ch = getattr(ed, "chainage", None)
        try:
            v = float(ch)
        except:
            continue
        if smin <= v <= smax:
            keep.append(ed)
    return _new_sites(s, keep)


def by_freq(sites: Any, fmin: float, fmax: float):
    r"""
    Select sites that contain at least one data row with
    frequency inside a closed interval.

    A site is kept if its frequency array ``f`` contains any
    finite value satisfying :math:`fmin \le f \le fmax`. Sites
    with empty or non-finite frequency arrays are skipped.

    Parameters
    ----------
    sites : Any
        A :class:`~pycsamt.site.base.Sites` instance, an
        ``EDICollection``, a sequence of ``EDIFile`` objects, or
        any iterable yielding EDI-like items.
    fmin : float
        Minimum frequency (inclusive).
    fmax : float
        Maximum frequency (inclusive).

    Returns
    -------
    pycsamt.site.base.Sites
        A new ``Sites`` wrapper containing only the EDI items with
        at least one finite frequency in the requested window.

    Notes
    -----
    Frequencies are obtained via
    ``pycsamt.site.utils.get_freq(ed)``. The check is membership
    based (any row in range), not a full slicing or resampling.
    Use :func:`pycsamt.site.edit.select_freq` to actually subset
    rows by frequency.

    Examples
    --------
    >>> from pycsamt.site.base import Sites
    >>> from pycsamt.site.selection import by_freq
    >>> s = Sites([e1, e2, e3])  # EDIFile objects
    >>> out = by_freq(s, fmin=0.5, fmax=2.0)
    >>> [t.name for t in out]
    ['A02', 'A03']

    See Also
    --------
    pycsamt.site.selection.by_names :
        Name-based selection using glob, regex, or callables.
    pycsamt.site.selection.by_chainage :
        Select by stored chainage range when available.
    pycsamt.site.edit.select_freq :
        Subset frequency rows within sites.
    """

    s = _to_sites(sites)
    keep = []
    for ed in iter_edifiles(s.edic):
        f = get_freq(ed)
        if f.size == 0:
            continue
        m = np.isfinite(f) & (f >= fmin) & (f <= fmax)
        if bool(m.any()):
            keep.append(ed)
    return _new_sites(s, keep)


def by_bbox(
    sites: Any,
    minlat: float,
    minlon: float,
    maxlat: float,
    maxlon: float,
):
    r"""
    Select sites that fall inside an axis-aligned geographic box.

    The selection is performed in latitude/longitude degrees and
    assumes a geographic CRS (WGS84-like). A site is kept if its
    stored coordinates satisfy

    .. math::

       minlat \le lat \le maxlat \;\;\text{and}\;\;
       minlon \le lon \le maxlon .

    Parameters
    ----------
    sites : Any
        A :class:`~pycsamt.site.base.Sites` instance, an
        ``EDICollection``, a sequence of ``EDIFile`` objects, or
        any iterable yielding EDI-like items.
    minlat, minlon, maxlat, maxlon : float
        Latitude and longitude bounds in degrees. Bounds are
        inclusive.

    Returns
    -------
    pycsamt.site.base.Sites
        A new ``Sites`` wrapper with only the items whose coords
        are inside the box.

    Notes
    -----
    This is a simple axis-aligned test in lat/lon and does not
    handle antimeridian wrapping. If longitudes cross the
    antimeridian (for example, from 170 to -170 deg), split the
    selection into two boxes and union the results. Coordinates
    are retrieved via
    :func:`pycsamt.site.utils.get_coords`.

    Examples
    --------
    >>> from pycsamt.site.base import Sites
    >>> from pycsamt.site.selection import by_bbox
    >>> s = Sites([e1, e2, e3])  # EDIFile objects
    >>> out = by_bbox(s, 24.0, 9.0, 27.0, 11.0)
    >>> [site.name for site in out]
    ['S01', 'S03']

    See Also
    --------
    pycsamt.site.selection.by_freq :
        Keep sites with at least one frequency inside a window.
    pycsamt.site.selection.by_chainage :
        Select by stored chainage range.
    pycsamt.site.selection.by_predicate :
        Arbitrary user-defined filtering.
    pycsamt.site.base.Sites.closest :
        Find the closest site to a target coordinate.

    References
    ----------
    .. [1] Snyder, J. P., "Map Projections: A Working Manual",
           USGS Professional Paper 1395.
    """

    s = _to_sites(sites)
    keep = []
    for ed in iter_edifiles(s.edic):
        c = get_coords(ed)
        if _in_box(c.lat, c.lon, minlat, minlon, maxlat, maxlon):
            keep.append(ed)
    return _new_sites(s, keep)


def by_predicate(sites: Any, pred: Callable[[Any], bool]):
    r"""
    Select sites using a user-supplied predicate function.

    The predicate is called for each EDI-like object and should
    return ``True`` to keep the site. Any exception raised by the
    predicate is caught and treated as a ``False`` (site is not
    kept). This makes bulk filtering robust against occasional
    data issues.

    Parameters
    ----------
    sites : Any
        A :class:`~pycsamt.site.base.Sites` instance, an
        ``EDICollection``, a sequence of ``EDIFile`` objects, or
        any iterable yielding EDI-like items.
    pred : Callable[[Any], bool]
        Function receiving a single EDI-like object and returning
        a boolean.

    Returns
    -------
    pycsamt.site.base.Sites
        A new ``Sites`` wrapper containing only the sites for
        which ``pred(site)`` returned ``True``.

    Notes
    -----
    The objects passed to ``pred`` are the raw EDI containers, not
    the :class:`~pycsamt.site.base.Site` wrapper. If you prefer
    the wrapper API, wrap the object inside the predicate:

    ``lambda ed: Site(ed).has_component("Zxy")``.

    Examples
    --------
    Keep sites that have at least 3 frequency rows:

    >>> from pycsamt.site.base import Sites, Site
    >>> from pycsamt.site.selection import by_predicate
    >>> s = Sites([e1, e2, e3])
    >>> out = by_predicate(
    ...     s, lambda ed: (Site(ed).freq is not None and
    ...                    len(Site(ed).freq) >= 3)
    ... )
    >>> [t.name for t in out]
    ['A01', 'A03']

    Keep sites whose name matches a rule:

    >>> import re
    >>> from pycsamt.site.utils import station_name
    >>> rule = re.compile(r"^X_")
    >>> out = by_predicate(s, lambda ed: bool(rule.search(
    ...     station_name(ed))))
    >>> [t.name for t in out]
    ['X_E01', 'X_E02']

    See Also
    --------
    pycsamt.site.selection.by_names :
        Name-based selection with glob or regex patterns.
    pycsamt.site.selection.drop_empty :
        Remove sites with no usable data arrays.
    pycsamt.site.base.Sites.select :
        Selection API on the wrapper.

    References
    ----------
    .. [1] Gamble, T. D. et al., "Magnetotellurics with a remote
           reference", Geophysics, 44(1), 53-68, 1979.
    """

    s = _to_sites(sites)
    keep = []
    for ed in iter_edifiles(s.edic):
        try:
            if pred(ed):
                keep.append(ed)
        except Exception:
            pass
    return _new_sites(s, keep)


def keep_finite_z(sites: Any):
    r"""
    Keep sites that contain at least one finite impedance value.

    A site is considered to have finite data if either of the
    following is true:

    1. The impedance tensor array (``Z.z`` or ``Z._z``) contains
       any finite real or imaginary entry.
    2. If the tensor is not present, a resistivity array
       (``Z._resistivity`` or ``Z.rho``) exists and has at least
       one finite value.

    Parameters
    ----------
    sites : Any
        A :class:`~pycsamt.site.base.Sites` instance, an
        ``EDICollection``, a sequence of ``EDIFile`` objects, or
        any iterable yielding EDI-like items.

    Returns
    -------
    pycsamt.site.base.Sites
        A new ``Sites`` wrapper with only the sites that contain
        finite impedance (or resistivity) values.

    Notes
    -----
    This function is intended as a coarse pre-filter to remove
    empty placeholders and fully invalid sites before more costly
    processing. It does not inspect errors or phases, and it does
    not modify the data. If a site has a ``Z`` container but all
    arrays are missing or fully non-finite, the site is dropped.

    Examples
    --------
    >>> from pycsamt.site.base import Sites
    >>> from pycsamt.site.selection import keep_finite_z
    >>> s = Sites([e1, e2, e3])
    >>> out = keep_finite_z(s)
    >>> [t.name for t in out]
    ['MT01', 'MT03']

    See Also
    --------
    pycsamt.site.selection.drop_empty :
        Remove sites with empty frequency or missing Z section.
    pycsamt.site.edit.fill_missing :
        Allocate arrays and replace invalid entries.
    pycsamt.site.compute.res_at_freq :
        Compute apparent resistivity at a specific frequency.
    """

    s = _to_sites(sites)
    keep = [ed for ed in iter_edifiles(s.edic) if _any_finite_z(ed)]
    return _new_sites(s, keep)


def mask_large_phase_err(sites: Any, thresh: float):
    r"""
    Filter out sites whose maximum phase-error exceeds a threshold.

    For each site, the function inspects the phase-error array
    when present (common attribute names are tried, e.g.
    ``_phase_err`` or ``phase_err``). If no phase-error array is
    found, the site is conservatively **kept**. Otherwise, the
    site is kept only when the maximum finite phase-error is less
    than or equal to ``thresh``.

    Parameters
    ----------
    sites : Any
        A :class:`~pycsamt.site.base.Sites` object, an
        ``EDICollection``, a sequence of ``EDIFile`` objects, or
        any iterable yielding EDI-like items.
    thresh : float
        Threshold on phase-error (same units as stored by the
        processing pipeline, usually degrees).

    Returns
    -------
    pycsamt.site.base.Sites
        New wrapper containing only sites that pass the phase
        error test.

    Notes
    -----
    The check uses a "best effort" attribute lookup and ignores
    non-finite values during the maximum computation. If the
    phase-error array is entirely missing, the site is kept.
    This behavior makes the filter robust when some sites did not
    store uncertainty products.

    Examples
    --------
    >>> from pycsamt.site.base import Sites
    >>> from pycsamt.site.selection import mask_large_phase_err
    >>> s = Sites([e1, e2, e3])
    >>> out = mask_large_phase_err(s, thresh=10.0)
    >>> [t.name for t in out]
    ['E01', 'E03']

    See Also
    --------
    pycsamt.site.selection.keep_finite_z :
        Keep sites that contain at least one finite impedance.
    pycsamt.site.selection.drop_empty :
        Remove sites with no usable arrays.
    pycsamt.site.edit.fill_missing :
        Allocate arrays and replace invalid entries.

    References
    ----------
    .. [1] Gamble, T. D., Goubau, W. M., Clarke, J., "Magneto-
           tellurics with a remote reference", Geophysics,
           44(1), 53-68, 1979.
    """

    s = _to_sites(sites)
    keep = []
    for ed in iter_edifiles(s.edic):
        m = _max_phase_err(ed)
        if not math.isfinite(m):
            keep.append(ed)
        elif m <= float(thresh):
            keep.append(ed)
    return _new_sites(s, keep)


def drop_empty(sites: Any):
    r"""
    Drop sites that are effectively empty.

    A site is considered empty when either:

    * The frequency vector is missing or has zero length.
    * The impedance container ``Z`` is missing.
    * The ``Z`` container is present but holds no usable arrays
      (for example, no ``z`` and no resistivity surrogate).

    Parameters
    ----------
    sites : Any
        A :class:`~pycsamt.site.base.Sites` object, an
        ``EDICollection``, a sequence of ``EDIFile`` objects, or
        any iterable yielding EDI-like items.

    Returns
    -------
    pycsamt.site.base.Sites
        New wrapper that excludes empty sites.

    Notes
    -----
    This is a coarse, fast filter that checks structural
    presence and basic array availability. It does **not** test
    for NaN-only content; for that, consider
    :func:`pycsamt.site.selection.keep_finite_z`.

    Examples
    --------
    >>> from pycsamt.site.base import Sites
    >>> from pycsamt.site.selection import drop_empty
    >>> s = Sites([e1, e2, e3])
    >>> out = drop_empty(s)
    >>> [t.name for t in out]
    ['MT01', 'MT02']

    See Also
    --------
    pycsamt.site.selection.keep_finite_z :
        Keep only sites with finite impedance or resistivity.
    pycsamt.site.selection.by_freq :
        Keep sites that touch a target frequency window.
    """

    s = _to_sites(sites)
    keep = [ed for ed in iter_edifiles(s.edic) if not _is_empty_site(ed)]
    return _new_sites(s, keep)


# ------------- Internal helpers --------------------------


def _to_sites(x: Any):
    """Coerce any edi-like into a Sites wrapper."""
    from .base import _to_sites as __to_sites

    return __to_sites(x)


def _new_sites(src: Any, items: list[Any]):
    """Build a new Sites preserving wrapper semantics."""
    from .base import Sites  # lazy import

    try:
        return Sites(items)
    except TypeError:
        return Sites(edic=items)


def _in_box(
    lat: float,
    lon: float,
    mnla: float,
    mnlo: float,
    mxla: float,
    mxlo: float,
) -> bool:
    if not (math.isfinite(lat) and math.isfinite(lon)):
        return False
    return (mnla <= lat <= mxla) and (mnlo <= lon <= mxlo)


def _name_matches(nm: str, pat: Any, case: bool) -> bool:
    # Prefer utils.match_name for DRY; keep a tiny fallback
    try:
        if not case and isinstance(pat, str):
            # match_name is already case-insensitive
            return match_name(pat, nm)
        if callable(pat) or isinstance(pat, re.Pattern):
            return match_name(pat, nm)
        return match_name(str(pat), nm)
    except Exception:
        # very small fallback
        s = str(pat)
        if not case:
            nm = nm.upper()
            s = s.upper()
        return nm == s


def _get_attr_any(obj: Any, *names: str):
    for nm in names:
        v = getattr(obj, nm, None)
        if v is not None:
            return v
    return None


def _any_finite_z(ed: Any) -> bool:
    z = getattr(ed, "Z", None)
    if z is None:
        return False
    arr = _get_attr_any(z, "_z", "z")
    if arr is None:
        r = _get_attr_any(z, "_resistivity", "resistivity", "rho")
        if r is None:
            return False
        a = np.asarray(r, float)
        return np.isfinite(a).any()
    a = np.asarray(arr, complex)
    return np.isfinite(a.real).any() or np.isfinite(a.imag).any()


def _max_phase_err(ed: Any) -> float:
    z = getattr(ed, "Z", None)
    if z is None:
        return float("nan")
    pe = _get_attr_any(z, "_phase_err", "phase_err")
    if pe is None:
        return float("nan")
    a = np.asarray(pe, float)
    if a.size == 0:
        return float("nan")
    with np.errstate(all="ignore"):
        return float(np.nanmax(a))


def _is_empty_site(ed: Any) -> bool:
    f = get_freq(ed)
    if f.size == 0:
        return True
    z = getattr(ed, "Z", None)
    if z is None:
        return True

    # If Z present, consider "empty" when there is no usable
    # impedance nor valid resistivity (sentinel-aware).
    arr = _get_attr_any(z, "_z", "z")
    if arr is not None:
        a = np.asarray(arr)
        with np.errstate(all="ignore"):
            if a.size == 0:
                return True
            return not np.isfinite(a).any()

    r = (
        getattr(z, "_resistivity", None)
        or getattr(z, "rho", None)
        or getattr(z, "res", None)
        or getattr(z, "resistivity", None)
    )
    if r is None:
        return True
    a = np.asarray(r, float)
    if a.size == 0:
        return True
    with np.errstate(all="ignore"):
        fin = np.isfinite(a)
        if not fin.any():
            return True
        # only huge sentinels -> empty
        return (np.abs(a[fin]) >= 1.0e30).all()
