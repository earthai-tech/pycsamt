# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

import math
from collections.abc import Iterable
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd

from .utils import (
    _ensure_head,
    get_freq,
    iter_edifiles,
    maybe_copy,
    set_station_name,
    station_name,
)
from .utils import (
    set_coords as _set_coords,
)

__all__ = [
    "rotate",
    "select_freq",
    "rename",
    "set_coords",
    "fill_missing",
    "recompute_res_phase",
    "rotate_all",
    "select_freq_all",
    "rename_all",
    "set_coords_all",
    "set_coords_from_table",
    "set_coords_from_en"
]


def rotate(site: Any, angle_deg: float, *,
           inplace: bool = False) -> Any:
    r"""
    Rotate impedance tensor Z (and tipper T, if present) by an
    azimuthal angle in degrees.

    The rotation is applied in the horizontal plane using the
    similarity transform :math:`Z' = R Z R^{-1}`, where :math:`R`
    is the 2x2 rotation matrix built from ``angle_deg``. When a
    tipper is available (either on the EDI object as ``T`` /
    ``TIP`` / ``Tip`` or attached to ``Z``), the 2-component
    vector is rotated consistently.

    Parameters
    ----------
    site : Any
        An EDI-like object (e.g., ``pycsamt.seg.edi.EDIFile``) or
        wrapper exposing a ``Z`` section compatible with a
        complex 2x2 impedance array and, optionally, a tipper.
    angle_deg : float
        Rotation angle in degrees. Positive values rotate the
        measurement axes according to the internal convention of
        this package. If your acquisition system defines the sign
        oppositely, use the negative of your desired angle.
    inplace : bool, optional
        If ``True``, mutate ``site`` in place. Otherwise, work on
        a shallow copy and return that copy. Default is
        ``False``.

    Returns
    -------
    Any
        The rotated object. If ``inplace`` is ``True``, this is
        the same object as ``site``; otherwise a new object.

    Notes
    -----
    - Error arrays (``z_error`` or aliases) are rotated with a
      magnitude-only scheme (using absolute values of the
      rotation matrices) as a pragmatic best-effort. This is a
      common, but approximate, practice.
    - Only arrays with shapes consistent with MT tensors
      (``(n, 2, 2)`` for Z, ``(n, 2)`` for T) are rotated. Other
      shapes are ignored silently.
    - The function is no-throw by design. If a section is not
      present or incompatible, it is skipped.

    Examples
    --------
    >>> from pycsamt.seg.edi import EDIFile
    >>> from pycsamt.site.edit import rotate
    >>> ed = EDIFile("path/to/station.edi")
    >>> ed_rot = rotate(ed, 30.0)       # copy by default
    >>> ed_rot2 = rotate(ed, -45.0, inplace=True)

    See Also
    --------
    select_freq : Subset rows by frequency range or indices.
    rename : Rename the station using a policy or explicit name.
    pycsamt.site.edit.rotate_all : Broadcast rotation over a
        collection.

    References
    ----------
    .. [1] Simpson, F. and Bahr, K. (2005). Practical
       Magnetotellurics. Cambridge Univ. Press.
    .. [2] SEG EDI format usage notes for MT tensors.
    """

    ed = _to_mutable(site, inplace=inplace)

    z = getattr(ed, "Z", None)
    if z is not None:
        zz = _get_attr_any(z, "z", "impedance", "_z")
        if zz is not None:
            r, rinv = _rotm(angle_deg)
            a = np.asarray(zz, complex)
            if a.ndim == 3 and a.shape[-2:] == (2, 2):
                newz = r[None, :, :] @ a @ rinv[None, :, :]
                _set_attr_first(z, ("z", "impedance", "_z"),
                                newz)
            ze = _get_attr_any(
                z, "z_error", "z_err", "impedance_err", "_z_err"
            )
            if ze is not None:
                ar = np.abs(r)
                ari = np.abs(rinv)
                e = np.asarray(ze, float)
                if e.ndim == 3 and e.shape[-2:] == (2, 2):
                    newe = ar[None, :, :] * e * ari[None, :, :]
                    _set_attr_first(
                        z, ("z_error", "z_err", "impedance_err",
                            "_z_err"),
                        newe,
                    )

    # tipper may live in ed.T / ed.TIP or Z.tipper
    tip_obj = None
    for cand in ("T", "TIP", "Tip"):
        tip_obj = getattr(ed, cand, None)
        if tip_obj is not None:
            break
    if tip_obj is not None:
        t = _get_attr_any(tip_obj, "tipper", "_tipper")
        if t is not None:
            _, rinv = _rotm(angle_deg)
            a = np.asarray(t, complex)
            if a.ndim == 2 and a.shape[1] == 2:
                _set_attr_first(
                    tip_obj, ("tipper", "_tipper"), a @ rinv.T
                )
    else:
        if z is not None:
            t = _get_attr_any(z, "tipper", "tip", "_tipper")
            if t is not None:
                _, rinv = _rotm(angle_deg)
                a = np.asarray(t, complex)
                if a.ndim == 2 and a.shape[1] == 2:
                    _set_attr_first(
                        z, ("tipper", "tip", "_tipper"),
                        a @ rinv.T,
                    )

    return ed


def select_freq(
    site: Any,
    *,
    fmin: float | None = None,
    fmax: float | None = None,
    keep: Iterable[int] | np.ndarray | None = None,
    inplace: bool = False,
) -> Any:
    r"""
    Subset the dataset along frequency by range or explicit
    indices, keeping all affected arrays aligned.

    This applies the selection to every frequency-indexed array
    found on the object, including Z, Z errors, derived
    resistivity/phase, and any tipper arrays resident on either
    the root EDI object or under its ``Z`` section.

    Parameters
    ----------
    site : Any
        An EDI-like object (e.g., ``pycsamt.seg.edi.EDIFile``)
        with a discoverable frequency vector.
    fmin : float, optional
        Keep rows with ``freq >= fmin``. Ignored if ``None``.
    fmax : float, optional
        Keep rows with ``freq <= fmax``. Ignored if ``None``.
    keep : Iterable[int] or numpy.ndarray, optional
        Explicit indices or a boolean mask to keep. If provided,
        ``fmin`` and ``fmax`` are ignored. Use integer indices
        for exact row picks, or a boolean mask of the same length
        as the frequency vector.
    inplace : bool, optional
        If ``True``, mutate ``site`` in place. Otherwise, operate
        on a shallow copy. Default is ``False``.

    Returns
    -------
    Any
        The object after selection. If ``inplace`` is ``True``,
        this is the same object as ``site``; otherwise a new
        object.

    Notes
    -----
    - All known aliases are sliced consistently (e.g., frequency,
      Z, Z error, resistivity, phase, tipper, and related
      per-row arrays) to preserve alignment.
    - If the frequency vector is empty or missing, the call is a
      no-op.
    - The function is no-throw by design; incompatible shapes are
      skipped silently.

    Examples
    --------
    Keep only rows with frequency >= 10 Hz:

    >>> from pycsamt.seg.edi import EDIFile
    >>> from pycsamt.site.edit import select_freq
    >>> ed = EDIFile("path/to/station.edi")
    >>> ed_hi = select_freq(ed, fmin=10.0)

    Keep the first and last rows explicitly:

    >>> sel = [0, -1]
    >>> ed_edge = select_freq(ed, keep=sel)

    Apply in place:

    >>> _ = select_freq(ed, fmin=1.0, fmax=100.0, inplace=True)

    See Also
    --------
    rotate : Rotate Z and tipper by an azimuth angle.
    rename : Rename the station identifiers.
    pycsamt.site.edit.select_freq_all : Broadcast selection over
        a collection.

    References
    ----------
    .. [1] SEG EDI format usage notes for frequency-indexed MT
       arrays.
    """

    ed = _to_mutable(site, inplace=inplace)
    f = get_freq(ed)

    if keep is not None:
        sel = np.asarray(list(keep))
        if getattr(ed, "Z", None) is not None:
            _slice_fields(ed.Z, sel)
        for cand in ("T", "TIP", "Tip"):
            obj = getattr(ed, cand, None)
            if obj is not None:
                _slice_fields(obj, sel)
        return ed

    if f.size == 0:
        return ed
    m = np.isfinite(f)
    if fmin is not None:
        m &= (f >= float(fmin))
    if fmax is not None:
        m &= (f <= float(fmax))
    if getattr(ed, "Z", None) is not None:
        _slice_fields(ed.Z, m)
    for cand in ("T", "TIP", "Tip"):
        obj = getattr(ed, cand, None)
        if obj is not None:
            _slice_fields(obj, m)
    return ed


def rename(
    site: Any,
    name: str | None = None,
    policy: Callable[[str], str] | None = None,
    *,
    inplace: bool = False,
) -> Any:
    r"""
    Rename a station using an explicit name or a policy function.

    This updates common station identifiers across the EDI header
    and attempts to keep them in sync so that downstream code
    resolves the new name consistently.

    Parameters
    ----------
    site : Any
        An EDI-like object (e.g., ``pycsamt.seg.edi.EDIFile``) or
        wrapper with a modifiable HEAD section.
    name : str, optional
        Explicit new station name to set. If provided, this takes
        precedence over ``policy``.
    policy : Callable[[str], str], optional
        A function mapping the current station name to a new one.
        Ignored when ``name`` is provided.
    inplace : bool, optional
        If ``True``, mutate ``site`` in place. Otherwise, operate
        on a shallow copy. Default is ``False``.

    Returns
    -------
    Any
        The object with updated identifiers. If ``inplace`` is
        ``True``, this is the same object as ``site``; otherwise a
        new object.

    Notes
    -----
    - The function writes multiple header fields when present
      (e.g., ``dataid``, ``station``, and other common aliases)
      and also mirrors the name to ``edi.name``. This increases
      the chance that name resolution remains stable across
      different readers.
    - The rename operation does not touch on-disk filenames. File
      paths remain unchanged unless you later write out using a
      template that depends on the station name.
    - If both ``name`` and ``policy`` are given, ``name`` wins.

    Examples
    --------
    Policy-based rename:

    >>> from pycsamt.seg.edi import EDIFile
    >>> from pycsamt.site.edit import rename
    >>> ed = EDIFile("path/to/station.edi")
    >>> ed2 = rename(ed, policy=lambda n: f"X_{n}")

    Explicit name, in place:

    >>> _ = rename(ed, name="ST123A", inplace=True)

    See Also
    --------
    rotate : Rotate Z and tipper by an azimuth angle.
    select_freq : Subset rows by frequency range or indices.
    pycsamt.site.edit.rename_all : Broadcast rename over a
        collection.
    pycsamt.site.base.Site : Wrapper that resolves a stable site
        name for indexing.

    References
    ----------
    .. [1] SEG EDI format field naming and common aliases for
       station identifiers.
    """

    ed = _to_mutable(site, inplace=inplace)
    cur = station_name(ed)

    if name is None and policy is None:
        return ed

    new = name if name is not None else cur
    if policy is not None:
        try:
            new = str(policy(cur))
        except Exception:
            new = str(cur)

    # Best-effort: update all usual identifiers and edi.name.
    set_station_name(ed, new)
    try:
        h = _ensure_head(ed)
        h.dataid = str(new or "")
        if not getattr(h, "station", None):
            h.station = str(new or "")
        # also mirror common alternates when present
        for k in ("sitename", "name", "STATION"):
            try:
                setattr(h, k, str(new))
            except Exception:
                pass
    except Exception:
        pass
    try:
        ed.name = str(new)
    except Exception:
        pass
    return ed

def set_coords(
    site: Any,
    *,
    lat: float | None = None,
    lon: float | None = None,
    elev: float | None = None,
    inplace: bool = False,
) -> Any:
    r"""
    Set geographic coordinates on the EDI header.

    Only the values explicitly provided are updated. The call
    delegates to the same coordinate writer used by the Site
    API, so downstream tools see consistent lat, lon, and elev
    fields.

    Parameters
    ----------
    site : Any
        An EDI-like object (e.g., ``pycsamt.seg.edi.EDIFile``)
        or wrapper exposing a mutable HEAD section.
    lat : float, optional
        Latitude in degrees. If ``None``, the field is left
        unchanged.
    lon : float, optional
        Longitude in degrees. If ``None``, the field is left
        unchanged.
    elev : float, optional
        Elevation in meters. If ``None``, the field is left
        unchanged.
    inplace : bool, optional
        If ``True``, mutate ``site`` in place. Otherwise, work
        on a shallow copy and return that copy. Default is
        ``False``.

    Returns
    -------
    Any
        The updated object. If ``inplace`` is ``True``, this is
        the same object as ``site``; otherwise a new object.

    Notes
    -----
    - The function validates numeric types and writes to common
      HEAD attribute names (``lat``, ``lon``, ``elev``).
    - If a field is not present in the header, a best-effort
      attribute is created.

    Examples
    --------
    >>> from pycsamt.seg.edi import EDIFile
    >>> from pycsamt.site.edit import set_coords
    >>> ed = EDIFile("path/to/station.edi")
    >>> ed2 = set_coords(ed, lat=35.1, lon=12.8, elev=1234.0)
    >>> _ = set_coords(ed, lat=36.0, inplace=True)

    See Also
    --------
    pycsamt.site.base.Site.set_coords : Object-oriented wrapper.
    pycsamt.site.edit.set_coords_all : Broadcast over a
        collection.

    References
    ----------
    .. [1] SEG EDI format, HEAD section fields for station
       coordinates.
    """

    ed = _to_mutable(site, inplace=inplace)
    _set_coords(ed, lat=lat, lon=lon, elev=elev, inplace=True)
    return ed


def fill_missing(
    site: Any,
    *,
    how: str = "zero",
    components: Iterable[str] = ("Z", "Tip"),
    inplace: bool = False,
) -> Any:
    r"""
    Replace missing or non-finite values in Z and/or tipper
    arrays with zeros or NaNs.

    The operation preserves shapes and alignment across arrays.
    If an array is absent, a new one is allocated with the
    correct shape inferred from the frequency vector.

    Parameters
    ----------
    site : Any
        An EDI-like object with a discoverable frequency vector
        and optionally Z and tipper sections.
    how : {"zero", "nan"}, optional
        Replacement policy. Use ``"zero"`` to fill with numeric
        zeros. Use ``"nan"`` to fill with NaN for all non-finite
        entries. Default is ``"zero"``.
    components : Iterable[str], optional
        Which components to process. Accepts items like
        ``"Z"`` or ``"Tip"`` (case-insensitive). Default is
        ``("Z", "Tip")``.
    inplace : bool, optional
        If ``True``, mutate in place. Otherwise, operate on a
        shallow copy and return that copy. Default is ``False``.

    Returns
    -------
    Any
        The object after filling. If ``inplace`` is ``True``,
        this is the same object as ``site``; otherwise a new
        object.

    Notes
    -----
    - Z arrays are expected as shape ``(n, 2, 2)`` and tipper as
      ``(n, 2)``. Only arrays with the expected shapes are
      modified or allocated.
    - When Z exists, the function also fills common aliases for
      errors and derived quantities, such as ``z_error``,
      ``rho`` (resistivity), and ``phase`` arrays, if present.
    - The frequency vector length ``n`` defines the number of
      rows used for any new arrays.
    - This function is no-throw by design. Incompatible or absent
      pieces are skipped silently.

    Examples
    --------
    Fill Z and tipper with zeros where values are non-finite:

    >>> from pycsamt.seg.edi import EDIFile
    >>> from pycsamt.site.edit import fill_missing
    >>> ed = EDIFile("path/to/station.edi")
    >>> ed2 = fill_missing(ed, how="zero")

    Fill only Z with NaNs, in place:

    >>> _ = fill_missing(ed, how="nan", components=("Z",),
    ...                  inplace=True)

    See Also
    --------
    select_freq : Subset rows by frequency while keeping arrays
        aligned.
    rotate : Rotate Z and tipper by an azimuth angle.

    References
    ----------
    .. [1] Simpson, F. and Bahr, K. (2005). Practical
       Magnetotellurics. Cambridge Univ. Press.
    .. [2] SEG EDI format usage notes for MT arrays and tipper.
    """

    how = str(how).lower()
    if how not in {"zero", "nan"}:
        raise ValueError("how must be 'zero' or 'nan'")
    ed = _to_mutable(site, inplace=inplace)
    f = get_freq(ed)
    n = int(f.size)

    if "Z" in set(map(str.upper, components)):
        z = getattr(ed, "Z", None)
        if z is not None:
            # Z 2x2 blocks
            zcur = _get_attr_any(z, "z", "impedance", "_z")
            znew = (
                np.zeros((n, 2, 2), complex)
                if zcur is None else
                _fill_array(zcur, (n, 2, 2), how)
            )
            _set_attr_first(z, ("z", "impedance", "_z"), znew)

            ze = _get_attr_any(
                z, "z_error", "z_err", "impedance_err", "_z_err"
            )
            if ze is not None:
                _set_attr_first(
                    z,
                    ("z_error", "z_err", "impedance_err",
                     "_z_err"),
                    _fill_array(ze, (n, 2, 2), how),
                )

            rr = _get_attr_any(
                z, "rho", "res", "resistivity", "_resistivity"
            )
            if rr is not None:
                _set_attr_first(
                    z,
                    ("rho", "res", "resistivity",
                     "_resistivity"),
                    _fill_array(rr, (n, 2, 2), how),
                )

            ph = _get_attr_any(z, "phase", "phi", "_phase")
            if ph is not None:
                _set_attr_first(
                    z,
                    ("phase", "phi", "_phase"),
                    _fill_array(ph, (n, 2, 2), how),
                )

            phe = _get_attr_any(z, "phase_err", "_phase_err")
            if phe is not None:
                _set_attr_first(
                    z,
                    ("phase_err", "_phase_err"),
                    _fill_array(phe, (n, 2, 2), how),
                )

    if "TIP" in set(map(str.upper, components)):
        tp = None
        for cand in ("T", "TIP", "Tip"):
            tp = getattr(ed, cand, None)
            if tp is not None:
                break
        if tp is not None:
            tcur = _get_attr_any(tp, "tipper", "_tipper")
            tnew = (
                np.zeros((n, 2), complex)
                if tcur is None else
                _fill_array(tcur, (n, 2), how)
            )
            _set_attr_first(tp, ("tipper", "_tipper"), tnew)

            te = _get_attr_any(tp, "tipper_err", "_tipper_err")
            if te is not None:
                _set_attr_first(
                    tp,
                    ("tipper_err", "_tipper_err"),
                    _fill_array(te, (n, 2), how),
                )
    return ed



def rotate_all(sites: Any, angle_deg: float, *,
               inplace: bool = False):
    r"""
    Rotate every site in a collection by an azimuthal angle in
    degrees.

    This is the broadcast variant of :func:`rotate`. It accepts
    either a ``Sites`` wrapper or any iterable of EDI-like
    objects and returns a new ``Sites`` unless ``inplace`` is
    requested.

    Parameters
    ----------
    sites : Any
        A ``pycsamt.site.base.Sites`` instance or any iterable
        of EDI-like objects.
    angle_deg : float
        Rotation angle in degrees, passed through to
        :func:`rotate`.
    inplace : bool, optional
        If ``True``, attempt to apply the rotation in place on
        the given container. Otherwise, return a new ``Sites``.
        Default is ``False``.

    Returns
    -------
    pycsamt.site.base.Sites
        A sites collection holding the rotated items, or the
        original container when mutated in place.

    Notes
    -----
    - The function preserves the input order of items.
    - If some items lack compatible arrays, they are skipped
      silently.

    Examples
    --------
    Using a list of EDI files:

    >>> from pycsamt.seg.edi import EDIFile
    >>> from pycsamt.site.edit import rotate_all
    >>> eds = [EDIFile("A.edi"), EDIFile("B.edi")]
    >>> ro = rotate_all(eds, 15.0)

    Using a Sites wrapper:

    >>> from pycsamt.site.base import Sites
    >>> sites = Sites(eds)
    >>> ro2 = rotate_all(sites, -30.0)

    See Also
    --------
    rotate : Single-site rotation.
    select_freq_all : Broadcast selection by frequency.
    """

    items = []
    for ed, _ in _each_site(sites):
        items.append(rotate(ed, angle_deg, inplace=False))
    return _wrap_output(sites, items, inplace=inplace)


def select_freq_all(
    sites: Any,
    *,
    fmin: float | None = None,
    fmax: float | None = None,
    keep: Iterable[int] | np.ndarray | None = None,
    inplace: bool = False,
):
    r"""
    Subset all sites in a collection along frequency, keeping
    arrays aligned.

    This is the broadcast variant of :func:`select_freq`. It
    accepts a ``Sites`` wrapper or any iterable of EDI-like
    objects and returns a new ``Sites`` unless ``inplace`` is
    requested.

    Parameters
    ----------
    sites : Any
        A ``pycsamt.site.base.Sites`` instance or any iterable
        of EDI-like objects.
    fmin : float, optional
        Keep rows with ``freq >= fmin``. Ignored if ``None``.
    fmax : float, optional
        Keep rows with ``freq <= fmax``. Ignored if ``None``.
    keep : Iterable[int] or numpy.ndarray, optional
        Explicit indices or a boolean mask to keep. If provided,
        ``fmin`` and ``fmax`` are ignored.
    inplace : bool, optional
        If ``True``, attempt to modify the given container in
        place. Otherwise, return a new ``Sites``. Default is
        ``False``.

    Returns
    -------
    pycsamt.site.base.Sites
        A sites collection with selection applied, or the
        original container when mutated in place.

    Notes
    -----
    - All known frequency-indexed arrays are sliced in sync for
      each site (Z, errors, derived quantities, and tipper).
    - Items missing a frequency vector are left unchanged.

    Examples
    --------
    >>> from pycsamt.seg.edi import EDIFile
    >>> from pycsamt.site.base import Sites
    >>> from pycsamt.site.edit import select_freq_all
    >>> eds = [EDIFile("A.edi"), EDIFile("B.edi")]
    >>> sites = Sites(eds)
    >>> out = select_freq_all(sites, fmin=1.0, fmax=100.0)

    Use explicit row indices for all sites:

    >>> out2 = select_freq_all(eds, keep=[0, -1])

    See Also
    --------
    select_freq : Single-site selection by frequency.
    rotate_all : Broadcast rotation over a collection.

    References
    ----------
    .. [1] SEG EDI format usage notes for frequency-indexed MT
       arrays.
    """

    items = []
    for ed, _ in _each_site(sites):
        items.append(
            select_freq(
                ed,
                fmin=fmin,
                fmax=fmax,
                keep=keep,
                inplace=False,
            )
        )
    return _wrap_output(sites, items, inplace=inplace)


def rename_all(
    sites: Any,
    *,
    policy: Callable[[str], str] | None = None,
    name_fn: Callable[[Any], str] | None = None,
    inplace: bool = False,
):
    r"""
    Batch rename a collection of sites.

    This is the broadcast variant of :func:`rename`. It accepts a
    ``Sites`` wrapper or any iterable of EDI-like objects and
    produces a new ``Sites`` (unless ``inplace`` is ``True``).

    Parameters
    ----------
    sites : Any
        A :class:`pycsamt.site.base.Sites` instance or any
        iterable of EDI-like objects.
    policy : Callable[[str], str], optional
        Function mapping the *current* station name to a new one,
        e.g. ``lambda n: f"X_{n}"``. Ignored if ``name_fn`` is
        provided.
    name_fn : Callable[[Any], str], optional
        Function mapping each EDI object to a new name, e.g. from
        the file stem. This takes precedence over ``policy`` and
        is useful to guarantee uniqueness across sites.
    inplace : bool, optional
        If ``True``, attempt to modify the given container in
        place. Otherwise, return a new ``Sites``. Default is
        ``False``.

    Returns
    -------
    pycsamt.site.base.Sites
        A collection holding the renamed items, or the original
        container when mutated in place.

    Notes
    -----
    - The rename updates common header identifiers and mirrors
      the name to ``edi.name`` for robust downstream resolution.
    - If multiple sites share the same original name, a simple
      ``policy`` like ``lambda n: "X_" + n`` may produce duplicate
      outputs. Prefer ``name_fn`` that uses a unique attribute
      (e.g., the file stem) to avoid collisions.
    - The function preserves input order and is no-throw; items
      that cannot be renamed are skipped.

    Examples
    --------
    Policy-based rename for all items:

    >>> from pycsamt.site.edit import rename_all
    >>> out = rename_all(eds, policy=lambda n: f"X_{n}")

    Unique names from file stems:

    >>> from pathlib import Path
    >>> out = rename_all(
    ...     eds,
    ...     name_fn=lambda e: f"X_{Path(getattr(e,'path','')).stem}",
    ... )

    See Also
    --------
    rename : Single-site rename helper.
    set_coords_all : Batch coordinate assignment.
    pycsamt.site.base.Sites : Collection wrapper used here.

    References
    ----------
    .. [1] SEG EDI format field naming conventions for station
       identifiers.
    """

    items = []
    for ed, _ in _each_site(sites):
        cur = station_name(ed)
        new = name_fn(ed) if name_fn else (policy(cur)
                                           if policy else cur)
        items.append(rename(ed, name=str(new), inplace=False))
    return _wrap_output(sites, items, inplace=inplace)


def set_coords_all(
    sites: Any,
    src: Any,
    *,
    inplace: bool = False,
):
    r"""
    Batch set coordinates for a collection of sites.

    Coordinates can be provided by a callable, a mapping keyed by
    station name, or an object exposing a ``.frame`` attribute
    with tabular data.

    Parameters
    ----------
    sites : Any
        A :class:`pycsamt.site.base.Sites` instance or any
        iterable of EDI-like objects.
    src : Any
        One of:
          * ``callable(edi) -> (lat, lon, elev)``
          * ``mapping[name] -> (lat, lon, elev)``
          * object with ``.frame`` DataFrame-like with columns:
            ``station`` plus either ``lat``/``lon`` or
            ``latitude``/``longitude``, and optionally
            ``elev``/``elevation``.
    inplace : bool, optional
        If ``True``, attempt to modify the given container in
        place. Otherwise, return a new ``Sites``. Default is
        ``False``.

    Returns
    -------
    pycsamt.site.base.Sites
        A collection with updated coordinates, or the original
        container when mutated in place.

    Notes
    -----
    - The lookup order is: ``callable`` then ``mapping`` by
      station name, then the optional ``.frame`` table. The first
      source that returns a non-``None`` triple is used.
    - Latitude and longitude are expected in degrees; elevation
      in meters.
    - The function preserves input order and is no-throw; items
      without a matching entry are left unchanged.

    Examples
    --------
    From a mapping keyed by names:

    >>> from pycsamt.site.edit import set_coords_all
    >>> coords = {"S01": (35.1, 12.8, 1234.0)}
    >>> out = set_coords_all(eds, coords)

    From a callable using file stems:

    >>> from pathlib import Path
    >>> def pick(edi):
    ...     stem = Path(getattr(edi, "path", "")).stem
    ...     return (10.0, 20.0, 100.0) if stem == "S01" else \
    ...            (11.0, 21.0, 110.0)
    >>> out = set_coords_all(eds, pick)

    From a pandas DataFrame holder:

    >>> class Holder:
    ...     def __init__(self, frame):
    ...         self.frame = frame
    >>> out = set_coords_all(eds, Holder(df))

    See Also
    --------
    set_coords : Single-site coordinate update.
    rename_all : Batch rename helper.
    pycsamt.site.base.Site.set_coords : OOP variant per site.

    References
    ----------
    .. [1] SEG EDI HEAD section fields for station coordinates.
    """

    def _lookup_by_frame(name: str):
        fr = getattr(src, "frame", None)
        if fr is None:
            return None
        try:
            col = ("station" if "station" in fr.columns else None)
            if not col:
                return None
            row = fr[fr[col] == name]
            if row is None or len(row) == 0:
                return None
            la = float(row["lat"].iloc[0] if "lat" in row else
                       row["latitude"].iloc[0])
            if "lon" in row:
                lo = float(row["lon"].iloc[0])
            else:
                lo = float(row["longitude"].iloc[0])
            ev = float(row.get("elev", row.get("elevation",
                                               0.0)).iloc[0])
            return la, lo, ev
        except Exception:
            return None

    items = []
    for ed, _ in _each_site(sites):
        nm = station_name(ed)
        lat = lon = elv = None
        if callable(src):
            try:
                v = src(ed)
                if v is not None:
                    lat, lon, elv = v
            except Exception:
                pass
        elif hasattr(src, "get"):
            try:
                v = src.get(nm)
                if v is not None:
                    lat, lon, elv = v
            except Exception:
                pass
        else:
            v = _lookup_by_frame(nm)
            if v is not None:
                lat, lon, elv = v
        items.append(
            set_coords(
                ed,
                lat=lat,
                lon=lon,
                elev=elv,
                inplace=False,
            )
        )
    return _wrap_output(sites, items, inplace=inplace)

def recompute_res_phase(
    site: Any,
    *,
    inplace: bool = False,
) -> Any:
    r"""
    Recompute apparent resistivity and phase from the impedance
    tensor Z for a single site.

    The function looks for a ``Z`` section and, if present, calls
    its ``compute_resistivity_phase()`` method. The operation is
    best-effort and suppresses exceptions.

    Parameters
    ----------
    site : Any
        An EDI-like object (``EDIFile``) or a wrapper exposing a
        ``Z`` section with the expected API.
    inplace : bool, optional
        If ``True``, mutate the given object in place. Otherwise
        work on a shallow copy and return it. Default is
        ``False``.

    Returns
    -------
    Any
        The mutated object (in place) or a new object (copy).

    Notes
    -----
    - Derived quantities are typically written under common
      aliases (e.g., ``rho`` or ``resistivity`` for apparent
      resistivity in ohm-m, and ``phase`` in degrees).
    - The method assumes Z has shape ``(n, 2, 2)`` and that a
      frequency vector is present. If these are missing or
      incompatible, nothing is changed.

    Examples
    --------
    Single-site:

    >>> from pycsamt.seg.edi import EDIFile
    >>> from pycsamt.site.edit import recompute_res_phase
    >>> ed = EDIFile("path/to/station.edi")
    >>> ed2 = recompute_res_phase(ed)

    In place:

    >>> _ = recompute_res_phase(ed, inplace=True)

    See Also
    --------
    select_freq : Keep a subset of rows by frequency.
    fill_missing : Ensure arrays are allocated and finite before
        recomputation.
    pycsamt.site.base.Site.to_dataframe : Inspect derived arrays.

    References
    ----------
    .. [1] Simpson, F. and Bahr, K. (2005). Practical
       Magnetotellurics. Cambridge Univ. Press.
    .. [2] SEG EDI format usage notes for derived MT quantities.
    """

    ed = _to_mutable(site, inplace=inplace)
    try:
        z = getattr(ed, "Z", None)
        if z is None:
            return ed
        fn = getattr(z, "compute_resistivity_phase", None)
        if callable(fn):
            try:
                # All args are optional in ResPhase API.
                fn()
            except TypeError:
                # Some forks might require explicit Nones.
                try:
                    fn(None, None, None)
                except Exception:
                    pass
    except:
        # Best-effort: never raise here
        pass
    return ed

def set_coords_from_table(
    sites: Any,
    table: Any,
    *,
    columns: dict | None = None,
    crs_from: str | None = None,
    to_crs: str = "EPSG:4326",
    inplace: bool = False,
):
    r"""
    Set site coordinates for many EDI files from a table.

    This high-level helper accepts a wide range of table-like
    objects (CSV path, pandas DataFrame, numpy structured array,
    or list of dicts / tuples). It normalizes column names,
    optionally projects easting/northing to lon/lat, builds a
    mapping ``{station: (lat, lon, elev)}``, and delegates to
    :func:`set_coords_all`.

    Parameters
    ----------
    sites : Any
        A :class:`~pycsamt.site.base.Sites` instance, an iterable
        of ``EDIFile`` objects, or anything accepted by
        :func:`set_coords_all`.
    table : Any
        One of:
          - Path to a CSV or whitespace-delimited text file.
          - A ``pandas.DataFrame``.
          - A numpy structured array.
          - A list of dicts or a list of tuples.
    columns : dict, optional
        Explicit column mapping. Keys are canonical names and
        values are the actual column names present in ``table``.
        Supported canonical keys are:
        ``'station'``, ``'lat'``, ``'lon'``, ``'elev'``,
        ``'easting'``, ``'northing'``. Matching is case-insensitive.
    crs_from : str, optional
        Source CRS used when the table provides ``easting`` and
        ``northing`` instead of ``lat`` and ``lon``. Required in
        that case. Example: ``'EPSG:32631'`` for UTM 31N.
    to_crs : str, default "EPSG:4326"
        Target CRS for output coordinates. The default is WGS84
        lon/lat.
    inplace : bool, default False
        If ``True``, mutate the given collection and return it.
        Otherwise return a new :class:`~pycsamt.site.base.Sites`
        with updated EDI objects.

    Returns
    -------
    Any
        The same semantics as :func:`set_coords_all`: either the
        mutated input (``inplace=True``) or a new
        :class:`~pycsamt.site.base.Sites` instance.

    Raises
    ------
    ValueError
        If the ``station`` column cannot be resolved, or if neither
        (``lat``, ``lon``) nor (``easting``, ``northing``) can be
        resolved, or if ``crs_from`` is required but missing.
    ImportError
        If projection is needed (``easting``/``northing`` present)
        but ``pyproj`` is not installed.

    Notes
    -----
    Column detection is case-insensitive and understands common
    aliases:

    - ``station``: ``station``, ``name``, ``site``, ``id``.
    - ``lat``: ``lat``, ``latitude``.
    - ``lon``: ``lon``, ``long``, ``longitude``.
    - ``elev``: ``elev``, ``elevation``, ``z``.
    - ``easting``: ``easting``, ``x``.
    - ``northing``: ``northing``, ``y``.

    When both geographic and projected fields are present, the
    geographic pair (``lat``, ``lon``) is preferred. If only
    projected fields are present, a valid ``crs_from`` must be
    provided and :mod:`pyproj` will be used for projection.

    Examples
    --------
    Load from a CSV path with standard columns::

        >>> from pycsamt.site.edit import set_coords_from_table
        >>> from pycsamt.site.base import Sites
        >>> edis = Sites([...])  # your EDI files
        >>> out = set_coords_from_table(
        ...     edis, "coords.csv", inplace=False
        ... )
        >>> isinstance(out, Sites)
        True

    Pass a DataFrame with aliases and an explicit mapping::

        >>> import pandas as pd
        >>> df = pd.DataFrame({
        ...     "name": ["S01", "S02"],
        ...     "latitude": [35.1, 35.2],
        ...     "long": [12.8, 12.9],
        ...     "elevation": [120.0, 130.0],
        ... })
        >>> out = set_coords_from_table(
        ...     edis,
        ...     df,
        ...     columns={"station": "name",
        ...              "lat": "latitude",
        ...              "lon": "long",
        ...              "elev": "elevation"},
        ...     inplace=False,
        ... )

    Use easting/northing with a source CRS::

        >>> df = pd.DataFrame({
        ...     "station": ["S10"],
        ...     "easting": [400000.0],
        ...     "northing": [5750000.0],
        ... })
        >>> out = set_coords_from_table(
        ...     edis,
        ...     df,
        ...     crs_from="EPSG:32631",  # UTM 31N
        ...     inplace=False,
        ... )

    See Also
    --------
    set_coords_all : Broadcast setting of coordinates using a
        mapping.
    set_coords : Single-site coordinate update helper.

    References
    ----------
    .. [1] EPSG Geodetic Parameter Registry,
           https://epsg.org/
    .. [2] pyproj documentation, https://pyproj4.github.io/pyproj/
    """

    df, cols = _maybe_df(table, columns=columns)
    mapping = _frame_to_mapping(
        df, cols, crs_from=crs_from, to_crs=to_crs
    )
    return set_coords_all(sites, mapping, inplace=inplace)


def set_coords_from_en(
    site: Any,
    easting: float,
    northing: float,
    *,
    crs_from: str,
    elev: float | None = None,
    to_crs: str = "EPSG:4326",
    inplace: bool = False,
) -> Any:
    r"""
    Project (easting, northing) to lon/lat and set a site's coords.

    This convenience wraps projection and assignment for a single
    EDI site. It projects the provided easting/northing from
    ``crs_from`` to ``to_crs`` (default WGS84 lon/lat), then calls
    :func:`set_coords`.

    Parameters
    ----------
    site : Any
        A single EDI-like object compatible with
        :func:`set_coords`.
    easting : float
        Easting value in meters for the source CRS.
    northing : float
        Northing value in meters for the source CRS.
    crs_from : str
        The EPSG or PROJ string that identifies the source CRS,
        for example ``'EPSG:32631'`` for UTM 31N.
    elev : float, optional
        Elevation in meters to store in the EDI header. If not
        provided, the previous elevation is kept (if any).
    to_crs : str, default "EPSG:4326"
        Target CRS for output coordinates. The default is WGS84
        lon/lat.
    inplace : bool, default False
        If ``True``, mutate the given object and return it.
        Otherwise work on a copy and return the copy.

    Returns
    -------
    Any
        The mutated site (``inplace=True``) or a new site object
        with updated coordinates.

    Raises
    ------
    ImportError
        If :mod:`pyproj` is not installed.
    Exception
        Any errors raised by the underlying projection engine or
        by :func:`set_coords` may propagate.

    Notes
    -----
    Projection uses :mod:`pyproj` with ``always_xy=True`` so that
    the axis order is interpreted as (lon, lat). Units are assumed
    to be meters for the easting and northing values.

    The returned object follows the same semantics as
    :func:`set_coords`. If you need to update many sites, prefer
    :func:`set_coords_from_table` or :func:`set_coords_all`.

    Examples
    --------
    Update a single site from UTM 31N (EPSG:32631) coordinates::

        >>> from pycsamt.site.edit import set_coords_from_en
        >>> site = ...  # an EDIFile
        >>> site2 = set_coords_from_en(
        ...     site,
        ...     easting=400000.0,
        ...     northing=5750000.0,
        ...     crs_from="EPSG:32631",
        ...     elev=250.0,
        ...     inplace=False,
        ... )

    Do the update in place and keep the existing elevation::

        >>> _ = set_coords_from_en(
        ...     site,
        ...     easting=500000.0,
        ...     northing=4600000.0,
        ...     crs_from="EPSG:32630",
        ...     inplace=True,
        ... )

    See Also
    --------
    set_coords : Assign lat, lon, elev on a single site.
    set_coords_from_table : Batch update from tabular input.
    set_coords_all : Broadcast update using a mapping.

    References
    ----------
    .. [1] EPSG Geodetic Parameter Registry,
           https://epsg.org/
    .. [2] pyproj documentation, https://pyproj4.github.io/pyproj/
    """

    lo, la = _project_en_to_lonlat(
        np.asarray([easting], float),
        np.asarray([northing], float),
        crs_from,
        to_crs=to_crs,
    )
    return set_coords(
        site,
        lat=float(la[0]),
        lon=float(lo[0]),
        elev=(None if elev is None else float(elev)),
        inplace=inplace,
    )


def _to_mutable(ed: Any, *, inplace: bool) -> Any:
    return ed if inplace else maybe_copy(ed)


def _rotm(angle_deg: float) -> tuple[np.ndarray, np.ndarray]:
    th = float(angle_deg) * (math.pi / 180.0)
    c = math.cos(th)
    s = math.sin(th)
    r = np.array([[c, s], [-s, c]], dtype=float)
    rinv = np.array([[c, -s], [s, c]], dtype=float)
    return r, rinv


def _slice_fields(obj: Any, sl: Any) -> None:
    groups = [
        ("freq", ("freq", "frequency", "_freq")),
        ("z", ("z", "impedance", "_z")),
        ("z_error", ("z_error", "z_err", "impedance_err",
                     "_z_err")),
        ("rho", ("rho", "res", "resistivity", "_resistivity")),
        ("phase", ("phase", "phi", "_phase")),
        ("tipper", ("tipper", "tip", "_tipper")),
        ("tipper_err", ("tipper_err", "_tipper_err")),
        ("amp", ("amp", "_amp")),
        ("azimuth", ("azimuth", "_azimuth")),
    ]
    for _, names in groups:
        for nm in names:
            a = getattr(obj, nm, None)
            if a is None:
                continue
            try:
                arr = np.asarray(a)
                setattr(obj, nm, arr[sl])
            except Exception:
                # keep going to other fields
                pass

def _maybe_df(table: Any, columns: dict | None = None):
    """
    Coerce many table-like inputs to a small DataFrame with
    columns for station and coordinates. Returns (df, cols)
    where cols is a dict with resolved canonical names:
      {'station','lat','lon','elev'} or {'station','easting',
       'northing','elev'} when lat/lon are absent.
    """
    cols = {"station": None, "lat": None, "lon": None,
            "elev": None, "easting": None, "northing": None}

    def _pick(name: str, choices: list[str]) -> str | None:
        for c in choices:
            if c in df.columns:
                return c
        return None

    # 1) import pandas lazily

    # 2) Build a DataFrame
    if isinstance(table, (str, Path)):
        p = Path(table)
        try:
            df = pd.read_csv(p)
        except Exception:
            # whitespace-delimited fallback
            df = pd.read_csv(p, delim_whitespace=True)
    elif pd is not None and isinstance(table, pd.DataFrame):
        df = table.copy()
    else:
        # try numpy structured or list[dict]/list[tuple]
        try:
            df = pd.DataFrame(table)
        except Exception as exc:
            raise TypeError(
                "Unsupported table-like input for coordinates."
            ) from exc

    # normalize column names to lower for matching
    df.columns = [str(c).strip() for c in df.columns]
    lowmap = {c.lower(): c for c in df.columns}
    df.columns = list(lowmap.keys())

    # 3) If explicit mapping provided, honor it
    if columns:
        for k, v in columns.items():
            k0 = str(k).lower()
            v0 = str(v).lower()
            if v0 in df.columns and k0 in cols:
                cols[k0] = v0

    # 4) Auto-detect columns when not provided
    if cols["station"] is None:
        cols["station"] = _pick("station", ["station", "name",
                                            "site", "id"])
    if cols["lat"] is None and cols["easting"] is None:
        cols["lat"] = _pick("lat", ["lat", "latitude"])
    if cols["lon"] is None and cols["easting"] is None:
        cols["lon"] = _pick("lon", ["lon", "long", "longitude"])
    if cols["easting"] is None and cols["lat"] is None:
        cols["easting"] = _pick("easting", ["easting", "x"])
    if cols["northing"] is None and cols["lat"] is None:
        cols["northing"] = _pick("northing", ["northing", "y"])
    if cols["elev"] is None:
        cols["elev"] = _pick("elev", ["elev", "elevation", "z"])

    if cols["station"] is None:
        raise ValueError(
            "Could not resolve 'station' column. Provide "
            "`columns={'station':'...'} or standard name."
        )

    # keep only the relevant columns
    keep = [v for v in cols.values() if v is not None]
    df = df.loc[:, sorted(set(keep))].copy()

    return df, cols


def _project_en_to_lonlat(
    e: np.ndarray,
    n: np.ndarray,
    crs_from: str,
    to_crs: str = "EPSG:4326",
) -> tuple[np.ndarray, np.ndarray]:
    """
    Project easting/northing to lon/lat using pyproj if present.
    """
    try:
        from pyproj import Transformer  # type: ignore
    except Exception as exc:
        raise ImportError(
            "pyproj is required for easting/northing projection. "
            "Install pyproj or pre-convert to lon/lat."
        ) from exc

    tr = Transformer.from_crs(crs_from, to_crs, always_xy=True)
    lon, lat = tr.transform(e, n)
    return np.asarray(lon, float), np.asarray(lat, float)

def _frame_to_mapping(
    df: Any,
    cols: dict,
    *,
    crs_from: str | None = None,
    to_crs: str = "EPSG:4326",
) -> dict[str, tuple[float, float, float]]:
    """
    Return mapping: name -> (lat, lon, elev) from a frame.
    """

    namec = cols["station"]
    elevc = cols["elev"]

    if cols["lat"] and cols["lon"]:
        latv = np.asarray(df[cols["lat"]], float)
        lonv = np.asarray(df[cols["lon"]], float)
    elif cols["easting"] and cols["northing"]:
        if not crs_from:
            raise ValueError(
                "Found easting/northing but `crs_from` not "
                "provided. Set crs_from='EPSG:32632' (example) or "
                "pre-convert to lon/lat."
            )
        ev = np.asarray(df[cols["easting"]], float)
        nv = np.asarray(df[cols["northing"]], float)
        lonv, latv = _project_en_to_lonlat(ev, nv, crs_from,
                                           to_crs=to_crs)
    else:
        raise ValueError(
            "Could not resolve lat/lon nor easting/northing."
        )

    if elevc and elevc in df.columns:
        elv = np.asarray(df[elevc], float)
    else:
        elv = np.full(latv.size, np.nan, float)

    names = [str(x) for x in df[namec].astype(str).tolist()]
    out: dict[str, tuple[float, float, float]] = {}
    for nm, la, lo, ev in zip(names, latv, lonv, elv):
        out[nm] = (float(la), float(lo), float(ev))
    return out


def _get_attr_any(obj: Any, *names: str) -> Any:
    for n in names:
        v = getattr(obj, n, None)
        if v is not None:
            return v
    return None


def _set_attr_first(obj: Any, names: Iterable[str], val: Any):
    for n in names:
        try:
            setattr(obj, n, val)
            return
        except Exception:
            pass
    # fallback: attach anyway on the first name
    names = list(names)
    if names:
        try:
            setattr(obj, names[0], val)
        except Exception:
            pass


def _each_site(container: Any):
    try:
        from .base import Sites  # lazy import
        if isinstance(container, Sites):
            for i, s in enumerate(getattr(container, "_items",
                                          [])):
                yield getattr(s, "edi", s), i
            return
    except Exception:
        pass
    for i, ed in enumerate(iter_edifiles(container)):
        yield ed, i


def _wrap_output(
    src: Any, items: list[Any], *, inplace: bool
) -> Any:
    if inplace:
        try:
            from .base import Sites
            if isinstance(src, Sites):
                for (s, _), ne in zip(_each_site(src), items):
                    # preserve Site wrappers, swap underlying EDI
                    s.edi = ne
                return src
        except Exception:
            pass
    try:
        from .base import Sites
        return Sites(items)
    except Exception:
        return items


def _fill_array(a: Any, shape: tuple[int, ...], how: str) -> np.ndarray:
    if a is None:
        if how == "zero":
            return np.zeros(shape, float)
        return np.full(shape, np.nan, float)
    out = np.array(a, copy=True)
    if how == "zero":
        out = np.nan_to_num(
            out, nan=0.0, posinf=0.0, neginf=0.0
        )
    else:
        out[~np.isfinite(out)] = np.nan
    return out
