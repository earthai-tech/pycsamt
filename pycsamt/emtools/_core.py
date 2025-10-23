
from __future__ import annotations

from typing import Any, Iterable, Optional, Tuple 
import warnings

import numpy as np 

from typing import Callable
from pycsamt.seg.collection import EDICollection
from pycsamt.site.base import Sites


def _wrap_one(ed):
    from ._input import ensure_sites
    return ensure_sites([ed], recursive=False, strict=False)


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
    if inplace:
        for ed in _iter_items(sites):
            Si = _wrap_one(ed)
            fn(Si)  # fn must honor inplace=True internally
        return sites
    out_items = []
    for ed in _iter_items(sites):
        Si = _wrap_one(ed)
        So = fn(Si)
        it = next(_iter_items(So))
        out_items.append(it)
    return _sites_from_items(out_items, verbose=verbose)



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


def _name(ed: Any, i: int) -> str:
    for k in ("station", "name", "site", "id"):
        v = getattr(ed, k, None)
        if isinstance(v, str) and v:
            return v
    return f"site_{i}"


def _first_attr(obj: Any, names: Tuple[str, ...]) -> Any:
    for n in names:
        try:
            v = getattr(obj, n)
        except Exception:
            v = None
        if v is not None:
            return v
    return None


def _as_1d_float(x: Any) -> Optional[np.ndarray]:
    try:
        a = np.asarray(x, dtype=float).ravel()
    except Exception:
        return None
    if a.ndim != 1 or a.size == 0:
        return None
    if not np.isfinite(a).any():
        return None
    return a


def _as_cmplx_nd(x: Any, shape0: Optional[int] = None) -> Optional[np.ndarray]:
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


def _trim_to_min(*arrs: np.ndarray) -> Tuple[np.ndarray, ...]:
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
    Z = _first_attr(ed, ("Z", "z"))
    if Z is None:
        return (None, None, None) if not with_errors else (
            None, None, None, None
        )
    z = _first_attr(Z, ("z", "Z"))
    fr = _first_attr(Z, ("freq",))
    if fr is None:
        fr = _first_attr(ed, ("freq",))
    ze = _first_attr(Z, ("z_err", "err_z", "zerr", "errors"))
    z = _as_cmplx_nd(z)
    fr = _as_1d_float(fr)
    if z is None or fr is None:
        return (None, None, None) if not with_errors else (
            None, None, None, None
        )
    # enforce (n,2,2) when possible
    if z.ndim == 3 and z.shape[1:3] == (2, 2):
        pass
    elif z.ndim == 2 and z.shape == (2, 2):
        z = z[None, ...]
    else:
        return (None, None, None) if not with_errors else (
            None, None, None, None
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
    T = _first_attr(ed, ("Tipper", "tipper"))
    if T is None:
        return (None, None, None) if not with_errors else (
            None, None, None, None
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
        return (None, None, None) if not with_errors else (
            None, None, None, None
        )
    # enforce (n,2)
    if t.ndim == 2 and t.shape[1] == 2:
        pass
    elif t.ndim == 1 and t.size == 2:
        t = t[None, ...]
    else:
        return (None, None, None) if not with_errors else (
            None, None, None, None
        )
    if isinstance(te, np.ndarray):
        te = np.asarray(te)
        if te.ndim == 2 and te.shape[1] != 2:
            te = None
    if with_errors and isinstance(te, np.ndarray):
        t, fr, te = _trim_to_min(t, fr, te)
        return T, t, fr, te
    t, fr = _trim_to_min(t, fr)
    return T, t, fr


def ensure_sites(
    sites: Any,
    *,
    recursive: bool = True,
    on_dup: str = "replace",
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

    from pycsamt.site.base import Sites,  to_sites as _to_sites 

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

    return S
