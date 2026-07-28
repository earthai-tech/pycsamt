from __future__ import annotations

import warnings
from collections.abc import Iterable
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
        raise ValueError(f"{label} must provide at least {n} axes; got {len(out)}.")
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
    item = edi if (edi is not None and getattr(edi, "Z", None) is not None) else ed
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
        return (None, None, None) if not with_errors else (None, None, None, None)
    z = _first_attr(Z, ("z", "Z"))
    fr = _first_attr(Z, ("freq",))
    if fr is None:
        fr = _first_attr(ed, ("freq",))
    ze = _first_attr(Z, ("z_err", "err_z", "zerr", "errors"))
    z = _as_cmplx_nd(z)
    fr = _as_1d_float(fr)
    if z is None or fr is None:
        return (None, None, None) if not with_errors else (None, None, None, None)
    # enforce (n,2,2) when possible
    if z.ndim == 3 and z.shape[1:3] == (2, 2):
        pass
    elif z.ndim == 2 and z.shape == (2, 2):
        z = z[None, ...]
    else:
        return (None, None, None) if not with_errors else (None, None, None, None)
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
        return (None, None, None) if not with_errors else (None, None, None, None)
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
        return (None, None, None) if not with_errors else (None, None, None, None)
    # enforce (n,2).  v2 EDI tipper objects commonly store the array as
    # (n_freq, 1, 2), mirroring impedance tensor dimensions.
    if t.ndim == 2 and t.shape[1] == 2:
        pass
    elif t.ndim == 3 and t.shape[1] == 1 and t.shape[2] == 2:
        t = t[:, 0, :]
    elif t.ndim == 1 and t.size == 2:
        t = t[None, ...]
    else:
        return (None, None, None) if not with_errors else (None, None, None, None)
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


def _station_positions(eds, spacing_m: float = 200.0) -> np.ndarray:
    """
    Return station positions [m] projected along the survey line.

    Tries ``east``/``north`` (or ``x``/``y``, ``easting``/``northing``)
    attributes on each EDI object and projects onto the bearing from the
    first to the last valid station.  Falls back to integer spacing when
    no coordinate metadata is found.
    """
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
            "ensure_sites(strict=True): no sites were resolved from the " "given input."
        )

    # ``Sites.ordered`` resolves ``None`` through the process-wide ordering
    # configuration, including its spatial-validation thresholds.  Keeping
    # that resolution in one place also preserves compatibility with custom
    # Sites-like implementations that expose the simple ``ordered(by)`` API.
    return S.ordered(order_by)
