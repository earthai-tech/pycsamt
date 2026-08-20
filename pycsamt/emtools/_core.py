from __future__ import annotations

import warnings
from collections.abc import Iterable
from pathlib import Path
from typing import (
    Any,
    Callable,
)

import numpy as np

from ..seg.collection import EDICollection
from ..site.base import Sites


def _axes_list(axes: Any, n: int, *, label: str = "axes") -> list[Any] | None:
    """Return *n* flattened matplotlib axes, or ``None`` if not supplied."""
    if axes is None:
        return None
    if isinstance(axes, np.ndarray):
        out = list(axes.ravel())
    elif isinstance(axes, (list, tuple)):
        out = list(np.asarray(axes, dtype=object).ravel())
    else:
        out = [axes]
    if len(out) < n:
        raise ValueError(
            f"{label} must provide at least {n} axes; got {len(out)}."
        )
    return out[:n]


def hide_polar_radius_labels(ax: Any) -> Any:
    """Hide radial tick labels on a polar axes while keeping grid rings."""
    if ax is None or getattr(ax, "name", "") != "polar":
        return ax
    try:
        ax.set_yticklabels([])
    except Exception:
        pass
    try:
        ax.tick_params(axis="y", labelleft=False, labelright=False)
    except Exception:
        pass
    return ax


def _wrap_one(ed):
    # For Sites-wrapped Site objects, unwrap to the underlying EDI item so
    # that ensure_sites can build a proper single-item Sites around it.
    edi = getattr(ed, "edi", None)
    item = (
        edi
        if (edi is not None and getattr(edi, "Z", None) is not None)
        else ed
    )
    return ensure_sites([item], recursive=False, strict=False)


def _sites_from_items(items, verbose=0):
    coll = EDICollection(items=items, verbose=verbose)
    try:
        return Sites(coll)
    except TypeError:
        return Sites(edic=coll)


def _apply_each(
    sites,
    fn: Callable[[Sites], Sites],
    *,
    inplace: bool,
    verbose: int = 0,
) -> Sites:
    # Collect EDI-level items so that _sites_from_items can rebuild a Sites.
    # For Sites-wrapped Site objects (site.edi = original EDI), _unwrap gets
    # the raw EDI so that EDICollection.add can find site.station / site.path.
    items_in = list(_iter_items(sites))

    if not inplace:
        # ``fn`` scales Z.z *in place* on the underlying EDI item. To honour
        # inplace=False we correct independent copies of the raw EDI objects
        # and leave the caller's sites pristine — otherwise the mutation
        # aliases back onto the input, so a "before" view re-read from it
        # shows the corrected data and before/after look identical. Copy at
        # the unwrapped-EDI level (deep-copying the Site wrapper breaks the
        # EDICollection rebuild, dropping every station).
        import copy

        raw = [copy.deepcopy(_unwrap(ed)) for ed in items_in]
        for r in raw:
            fn(_wrap_one(r))
        return _sites_from_items(raw, verbose=verbose)

    for ed in items_in:
        Si = _wrap_one(ed)
        fn(Si)  # fn modifies Z.z in-place on the underlying EDI item
    return sites


def _iter_items(sites: Any) -> Iterable[Any]:
    try:
        for it in sites:
            yield it
        return
    except Exception:
        pass
    items = getattr(sites, "items", None)
    if isinstance(items, dict):
        for _, it in items.items():
            yield it


def _unwrap(ed: Any) -> Any:
    """Return the raw EDI object behind a Site wrapper when possible."""
    try:
        from ..site.base import to_edis

        raw = to_edis(ed, strict=False)
        if isinstance(raw, list):
            return raw[0] if raw else ed
        return raw if raw is not None else ed
    except Exception:
        edi = getattr(ed, "edi", None)
        return edi if edi is not None else ed


def _name(ed: Any, i: int) -> str:
    # For Sites-wrapped Site objects the real station name lives in ed.edi.
    # Check that path first so the generic Site.name = "site" is not picked up.
    edi = getattr(ed, "edi", None)
    if edi is not None:
        for k in ("station", "id"):
            v = getattr(edi, k, None)
            if isinstance(v, str) and v:
                return v
    for k in ("station", "name", "site", "id"):
        v = getattr(ed, k, None)
        if isinstance(v, str) and v:
            return v
    return f"site_{i}"


def _first_attr(obj: Any, names: tuple[str, ...]) -> Any:
    for n in names:
        try:
            v = getattr(obj, n)
        except Exception:
            v = None
        if v is not None:
            return v
    return None


def _as_1d_float(x: Any) -> np.ndarray | None:
    try:
        a = np.asarray(x, dtype=float).ravel()
    except Exception:
        return None
    if a.ndim != 1 or a.size == 0:
        return None
    if not np.isfinite(a).any():
        return None
    return a


def _as_cmplx_nd(x: Any, shape0: int | None = None) -> np.ndarray | None:
    try:
        a = np.asarray(x, dtype=np.complex128)
    except Exception:
        return None
    if a.ndim < 1:
        return None
    if shape0 is not None and a.shape[0] != shape0:
        # allow later trimming; keep as-is
        return a
    return a


def _trim_to_min(*arrs: np.ndarray) -> tuple[np.ndarray, ...]:
    n = min(a.shape[0] for a in arrs if a is not None)
    out = []
    for a in arrs:
        out.append(a[:n] if a is not None else None)
    return tuple(out)


def _get_z_block(
    ed: Any,
    *,
    with_errors: bool = False,
) -> tuple:
    # Prefer uppercase-Z wrapper (EDIFile, _FakeZ).  If not found, try the
    # Sites-wrapped Site pattern where the real wrapper lives at ed.edi.Z.
    Z = _first_attr(ed, ("Z",))
    if Z is None:
        edi = getattr(ed, "edi", None)
        if edi is not None:
            Z = _first_attr(edi, ("Z",))
    if Z is None:
        Z = _first_attr(ed, ("z",))  # last resort: raw array
    if Z is None:
        return (
            (None, None, None) if not with_errors else (None, None, None, None)
        )
    z = _first_attr(Z, ("z", "Z"))
    fr = _first_attr(Z, ("freq",))
    if fr is None:
        fr = _first_attr(ed, ("freq",))
    ze = _first_attr(Z, ("z_err", "err_z", "zerr", "errors"))
    z = _as_cmplx_nd(z)
    fr = _as_1d_float(fr)
    if z is None or fr is None:
        return (
            (None, None, None) if not with_errors else (None, None, None, None)
        )
    # enforce (n,2,2) when possible
    if z.ndim == 3 and z.shape[1:3] == (2, 2):
        pass
    elif z.ndim == 2 and z.shape == (2, 2):
        z = z[None, ...]
    else:
        return (
            (None, None, None) if not with_errors else (None, None, None, None)
        )
    # trim all to min length
    if isinstance(ze, np.ndarray):
        ze = np.asarray(ze)
        if ze.ndim == 3 and ze.shape[1:3] != (2, 2):
            ze = None
    if with_errors and isinstance(ze, np.ndarray):
        z, fr, ze = _trim_to_min(z, fr, ze)
        return Z, z, fr, ze
    z, fr = _trim_to_min(z, fr)
    return Z, z, fr


def _get_t_block(
    ed: Any,
    *,
    with_errors: bool = False,
) -> tuple:
    # Only probe the exact-cased wrapper names here (mirrors
    # _get_z_block's "Z"-only first pass). Site-wrapped objects expose
    # a lowercase ``tipper`` property that already returns the final
    # raw array rather than a wrapper with its own ``.tipper``/``.freq``
    # sub-attributes; including it here would short-circuit the
    # ``ed.edi`` fallback below and get treated as an unwrapped array.
    T = _first_attr(ed, ("Tipper", "Tip"))
    if T is None:
        edi = getattr(ed, "edi", None)
        if edi is not None:
            T = _first_attr(edi, ("Tipper", "tipper", "Tip"))
    if T is None:
        # pycsamt.airborne.site.AirborneSite has no .Tipper/.Tip
        # wrapper and no .edi bridge at all (that bridge is exactly
        # what it exists to avoid) -- it exposes .tipper/.freq
        # directly off its EMTF document instead. Use the site
        # itself as the "T" handle so downstream mutation code (e.g.
        # ztem.mask_outside_ztem_band) can still write back through
        # AirborneSite.tipper's setter.
        if _first_attr(ed, ("tipper",)) is not None:
            T = ed
    if T is None:
        return (
            (None, None, None) if not with_errors else (None, None, None, None)
        )
    t = _first_attr(T, ("tipper", "T", "tx_ty"))
    fr = _first_attr(T, ("freq",))
    if fr is None:
        fr = _first_attr(ed, ("freq",))
    te = _first_attr(
        T,
        ("tipper_err", "err_tipper", "terr", "errors"),
    )
    t = _as_cmplx_nd(t)
    fr = _as_1d_float(fr)
    if t is None or fr is None:
        return (
            (None, None, None) if not with_errors else (None, None, None, None)
        )
    # enforce (n,2).  v2 EDI tipper objects commonly store the array as
    # (n_freq, 1, 2), mirroring impedance tensor dimensions.
    if t.ndim == 2 and t.shape[1] == 2:
        pass
    elif t.ndim == 3 and t.shape[1] == 1 and t.shape[2] == 2:
        t = t[:, 0, :]
    elif t.ndim == 1 and t.size == 2:
        t = t[None, ...]
    else:
        return (
            (None, None, None) if not with_errors else (None, None, None, None)
        )
    if isinstance(te, np.ndarray):
        te = np.asarray(te)
        if te.ndim == 3 and te.shape[1] == 1 and te.shape[2] == 2:
            te = te[:, 0, :]
        if te.ndim == 2 and te.shape[1] != 2:
            te = None
    if with_errors and isinstance(te, np.ndarray):
        t, fr, te = _trim_to_min(t, fr, te)
        return T, t, fr, te
    t, fr = _trim_to_min(t, fr)
    return T, t, fr


_GEO_OFFSET_CACHE_ATTR = "_pycsamt_utm_offset_m"
_WGS84_ELLIPSOID_ID = 23  # ELLIPSOIDS id used across pycsamt.gis


def _station_latlon_offsets(
    eds,
) -> list[tuple[float, float] | None]:
    """Return real (easting_m, northing_m) offsets from EDI HEAD lat/lon.

    Uses the dependency-free UTM projection in
    :func:`pycsamt.gis.utils.ll_to_utm` (no GDAL/pyproj required), so it
    gives true metric station coordinates rather than a hardcoded or
    reference-relative approximation. Each result is cached on the
    underlying EDI object (keyed by the exact lat/lon it was derived
    from) so repeated calls -- e.g. from successive
    ``plot_confidence_profile``/``station_confidence_table`` calls on the
    same survey, which each get a freshly re-wrapped ``Sites`` instance --
    do not re-parse the EDI header or redo the projection every time. The
    cache self-invalidates whenever a station's coordinates change.

    Returns one ``(easting, northing)`` (or ``None``) per input station,
    all in the *same* UTM zone. If the finite coordinates straddle more
    than one UTM zone, eastings would not be directly comparable, so
    every entry is returned as ``None`` and the caller should fall back
    to uniform spacing instead.
    """
    from ..gis.utils import ll_to_utm
    from ..site.utils import get_coords

    offsets: list[tuple[float, float] | None] = []
    zones: list[str | None] = []
    for ed in eds:
        # Unchanged primary path: ground Site (or any bare EDI-like
        # object exposing get_section()/HEAD lat-lon, as used by
        # e.g. this module's own test fixtures) resolves through
        # get_coords() exactly as before.
        edi = getattr(ed, "edi", ed)
        cache_key = edi
        try:
            c = get_coords(edi)
            lat, lon = float(c.lat), float(c.lon)
        except Exception:
            lat, lon = float("nan"), float("nan")

        if not (np.isfinite(lat) and np.isfinite(lon)):
            # AirborneSite (pycsamt.airborne.site) has no .edi at all
            # for get_coords() to resolve anything from -- it exposes
            # (lat, lon, elev) directly as .coords instead, the same
            # fallback tf.py's own _station_xy already uses for
            # exactly this reason. Only tried as a fallback here, so
            # every existing get_coords()-resolvable case (including
            # bare test fixtures with no .coords at all) is untouched.
            coords = getattr(ed, "coords", None)
            if coords is not None:
                try:
                    lat, lon = float(coords[0]), float(coords[1])
                    cache_key = ed
                except (TypeError, ValueError, IndexError):
                    lat, lon = float("nan"), float("nan")

        if not (np.isfinite(lat) and np.isfinite(lon)):
            offsets.append(None)
            zones.append(None)
            continue

        cached = getattr(cache_key, _GEO_OFFSET_CACHE_ATTR, None)
        if cached is not None and cached[0] == lat and cached[1] == lon:
            _, _, east, north, zone = cached
        else:
            try:
                zone, east, north = ll_to_utm(_WGS84_ELLIPSOID_ID, lat, lon)
            except Exception:
                offsets.append(None)
                zones.append(None)
                continue
            try:
                setattr(
                    cache_key,
                    _GEO_OFFSET_CACHE_ATTR,
                    (lat, lon, east, north, zone),
                )
            except Exception:
                pass
        offsets.append((east, north))
        zones.append(zone)

    if len({z for z in zones if z is not None}) > 1:
        return [None] * len(offsets)
    return offsets


def _station_positions(
    eds,
    spacing_m: float = 200.0,
    *,
    force_spacing: bool = False,
) -> np.ndarray:
    """
    Return station positions [m] projected along the survey line.

    By default, computes real inter-station distances: it tries
    ``east``/``north`` (or ``x``/``y``, ``easting``/``northing``)
    attributes on each EDI object first, then falls back to a local-metre
    projection of the EDI HEAD ``lat``/``lon`` coordinates (cached per EDI
    object, see :func:`_station_latlon_offsets`). Either way, positions are
    the projection of each station onto the bearing from the first to the
    last valid station -- i.e. real chainage, not a hardcoded step.

    ``spacing_m`` is only used as a uniform per-station fallback, for
    stations that have no usable coordinates (or, if none at all do, for
    the whole line). Set ``force_spacing=True`` to skip coordinate lookup
    entirely and always lay stations out at uniform ``spacing_m`` steps --
    e.g. when the available coordinates are known to be unreliable and a
    user-supplied spacing should be trusted instead.
    """
    eds = list(eds)
    if force_spacing:
        return np.arange(len(eds), dtype=float) * spacing_m

    coords = []
    for ed in eds:
        e = None
        for attr in ("east", "x", "easting"):
            v = getattr(ed, attr, None)
            if v is not None:
                try:
                    e = float(v)
                    break
                except (TypeError, ValueError):
                    pass
        n = None
        for attr in ("north", "y", "northing"):
            v = getattr(ed, attr, None)
            if v is not None:
                try:
                    n = float(v)
                    break
                except (TypeError, ValueError):
                    pass
        coords.append((e, n) if e is not None and n is not None else None)

    valid = [(i, c) for i, c in enumerate(coords) if c is not None]
    if len(valid) < 2:
        latlon_offsets = _station_latlon_offsets(eds)
        latlon_valid = [
            (i, c) for i, c in enumerate(latlon_offsets) if c is not None
        ]
        if len(latlon_valid) >= 2:
            coords = latlon_offsets
            valid = latlon_valid

    if len(valid) < 2:
        return np.arange(len(coords), dtype=float) * spacing_m

    pts = np.array([list(c) for _, c in valid], dtype=float)
    origin = pts[0]
    direction = pts[-1] - origin
    norm = float(np.linalg.norm(direction))
    if norm < 1.0:
        return np.arange(len(coords), dtype=float) * spacing_m

    direction /= norm
    positions = np.arange(len(coords), dtype=float) * spacing_m
    for idx, c in valid:
        positions[idx] = float(
            np.dot(np.array(list(c), dtype=float) - origin, direction)
        )
    return positions


def ensure_sites(
    sites: Any,
    *,
    recursive: bool = True,
    on_dup: str = "replace",
    order_by: str | None = None,
    strict: bool = False,
    verbose: int = 0,
):
    r"""
    Normalize arbitrary user input to a ``Sites`` object.

    This is the single entry-point validator for all ``emtools``
    public APIs. It guarantees that downstream code receives a
    :class:`~pycsamt.site.base.Sites` instance, no matter whether the
    caller passed a path, an ``EDIFile``/``EDICollection``,
    a single ``Site``, an existing ``Sites``, or an iterable of
    EDI-like items.

    Parameters
    ----------
    sites : Any
        User-provided input. Accepts filesystem paths, glob patterns,
        EDI-like objects (``EDIFile``, ``EDICollection``), ``Site``,
        ``Sites``, or iterables of such.
    recursive : bool, default=True
        When walking directories, recurse into subfolders if the
        lower-level coercion utility supports it.
    on_dup : {"replace", "keep_first", "keep_last", "raise"}, \
default="replace"
        Duplicate site-name policy. See
        :func:`pycsamt.seg.base.to_sites` for semantics.
    order_by : {"auto", "chainage", "input", "station", "latitude", \
"longitude"}, optional
        Site ordering policy. ``None`` uses the package-wide
        :data:`pycsamt.api.PYCSAMT_ORDERING` setting. Automatic mode uses
        coordinate-derived profile chainage only when the coordinates pass a
        conservative single-line geometry check; otherwise input order is
        preserved.
    strict : bool, default=False
        If ``True``, raise when no items can be resolved.
    verbose : int, default=0
        Verbosity level. ``>0`` emits warnings about duplicates and
        coercion steps.

    Returns
    -------
    pycsamt.site.base.Sites
        Canonical ``Sites`` wrapper suitable for all processing tools.

    Raises
    ------
    ValueError
        If ``sites`` is ``None`` or, in ``strict`` mode, nothing could
        be resolved into EDI-like items.
    TypeError
        If the result is not a ``Sites`` instance (indicates a broken
        installation or import cycle).

    Notes
    -----
    All ``emtools`` module functions should start with::

        S = ensure_sites(sites, ...)

    to guarantee API consistency across the package.
    """
    if sites is None:
        raise ValueError(
            "ensure_sites: got None. Provide a path, EDI-like object, "
            "Site/Sites, or an iterable of such."
        )

    from pycsamt.site.base import Sites
    from pycsamt.site.base import to_sites as _to_sites

    # Try rich signature first; fall back for older versions.
    try:
        S = _to_sites(
            sites,
            recursive=recursive,
            on_dup=on_dup,
            strict=strict,
            verbose=verbose,
        )
    except TypeError:
        if verbose > 0:
            warnings.warn(
                "ensure_sites: 'to_sites' does not accept the extended "
                "signature in this environment; falling back to "
                "legacy call.",
                RuntimeWarning,
                stacklevel=2,
            )
        S = _to_sites(sites)

    if not isinstance(S, Sites):  # very defensive
        raise TypeError(
            "ensure_sites: expected a Sites instance after coercion, "
            f"got {type(S)!r}. Check your installation/imports."
        )

    if strict and len(S) == 0:
        raise ValueError(
            "ensure_sites(strict=True): no sites were resolved from the "
            "given input."
        )

    # ``Sites.ordered`` resolves ``None`` through the process-wide ordering
    # configuration, including its spatial-validation thresholds.  Keeping
    # that resolution in one place also preserves compatibility with custom
    # Sites-like implementations that expose the simple ``ordered(by)`` API.
    return S.ordered(order_by)


def _looks_edi_incompatible(path: Any) -> bool | None:
    """Peek at one EMTF-XML file to see if it needs ``AirborneSites``.

    Returns ``True`` when the first matching file parses and has no
    impedance transfer function -- the one condition
    :func:`~pycsamt.emtf.converters.edi.emtf_to_edi` unconditionally
    refuses (see :mod:`pycsamt.airborne.site`'s module docstring),
    whatever the *other* transfer functions on the document are: a
    ZTEM-style tipper, an AFMAG interstation tensor, or a bare
    AFMAG scalar tilt with no transfer function at all besides
    itself. Returns ``False`` when it parses and carries an
    impedance TF (ordinary ground MT). Returns ``None`` when nothing
    could be determined (not a path, no matching file, or the file
    does not parse), in which case the caller should fall back to
    the existing ground-``Sites`` behavior rather than guess.
    """
    p = Path(path)
    if p.is_file():
        candidate = p
    elif p.is_dir():
        found = sorted(p.rglob("*.xml")) or sorted(p.glob("*.xml"))
        candidate = found[0] if found else None
    else:
        return None
    if candidate is None:
        return None
    try:
        from ..emtf.document import EMTF

        doc = EMTF.from_xml(candidate)
    except Exception:
        return None
    return doc.impedance is None


def ensure_any_sites(
    sites: Any,
    *,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    r"""
    Coerce to :class:`~pycsamt.site.base.Sites` or
    :class:`~pycsamt.airborne.site.AirborneSites`, whichever fits.

    Tipper-only technologies (ZTEM, AFMAG) can be delivered either
    as ground-style EDI/EMTF-XML carrying an impedance transfer
    function too (routes through the existing :func:`ensure_sites`),
    or as genuinely tipper-only EMTF-XML with no electric channel at
    all -- the case ``Sites``/``Site`` cannot represent (see
    :mod:`pycsamt.airborne.site`'s module docstring). This function
    is the single entry point tipper-based ``emtools`` functions
    (:mod:`pycsamt.emtools.ztem`, and any future ``emtools`` module
    facing the same choice) should coerce their ``sites`` parameter
    through instead of calling :func:`ensure_sites` directly, so that
    parameter stays a single, plainly-named argument accepting either
    container rather than gaining a second, redundant one.

    Parameters
    ----------
    sites : Any
        Anything :func:`ensure_sites` accepts, plus anything
        :func:`pycsamt.airborne.site.ensure_asites` accepts
        (:class:`~pycsamt.airborne.AirborneEMLine`/``AirborneEMDataset``,
        an existing
        :class:`~pycsamt.airborne.site.AirborneSites`/``AirborneSite``,
        or a path). A bare path/directory is routed to the airborne
        side only when peeking at its first EMTF-XML file shows no
        impedance transfer function at all -- whatever else the
        document carries instead (a ZTEM-style tipper, an AFMAG
        interstation tensor, a bare AFMAG scalar tilt); anything
        else (including an ambiguous or unreadable peek) goes
        through the unchanged, existing :func:`ensure_sites` path.
    recursive, on_dup, strict, verbose
        Forwarded to whichever of :func:`ensure_sites`/
        :func:`~pycsamt.airborne.site.ensure_asites` is selected.

    Returns
    -------
    pycsamt.site.base.Sites or pycsamt.airborne.site.AirborneSites

    See Also
    --------
    ensure_sites : The ground-only entry point this wraps.
    pycsamt.airborne.site.ensure_asites : The airborne-only entry
        point this wraps.
    """
    from ..airborne.base import AirborneEMDataset, AirborneEMLine
    from ..airborne.site import AirborneSite, AirborneSites, ensure_asites

    _airborne_types = (
        AirborneSites, AirborneEMLine, AirborneEMDataset, AirborneSite,
    )

    def _to_airborne():
        return ensure_asites(
            sites,
            recursive=recursive,
            on_dup=on_dup,
            strict=strict,
            verbose=verbose,
        )

    if isinstance(sites, _airborne_types):
        return _to_airborne()

    if isinstance(sites, (str, Path)):
        if _looks_edi_incompatible(sites):
            return _to_airborne()
        return ensure_sites(
            sites,
            recursive=recursive,
            on_dup=on_dup,
            strict=strict,
            verbose=verbose,
        )

    if isinstance(sites, (list, tuple)) and any(
        isinstance(it, _airborne_types) for it in sites
    ):
        return _to_airborne()

    return ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
