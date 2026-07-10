# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

from collections.abc import Iterator
from typing import Any

from .base import CoreObject, TFBundle, ensure_station, to_edi

__all__ = [
    "bundle_from_edi",
    "BundleMixin",
    "BundleContainerMixin",
]


def _get(obj: Any, *names: str) -> Any:
    for n in names:
        if hasattr(obj, n):
            return getattr(obj, n)
    return None


def bundle_from_edi(obj: Any) -> TFBundle:
    r"""
    Build a :class:`TFBundle` from an EDI-like object.

    The function probes common attribute names found in EDI objects
    and populates a neutral bundle. It tolerates missing fields and
    does not enforce array types.

    Parameters
    ----------
    obj : Any
        EDI-like object. Attributes such as ``freq``, ``z``,
        ``z_err``, ``tipper``, ``rho``, ``phase``, ``station``,
        ``lat``, ``lon``, ``elev``, ``azimuth`` are searched with
        several conventional aliases.

    Returns
    -------
    TFBundle
        A bundle with fields filled where data could be found.

    Notes
    -----
    If the tipper has shape ``(n, 1, 2)`` it is normalized to
    ``(n, 2)``. Unknown shapes are left unchanged.

    Examples
    --------
    >>> from pycsamt.core.mixins import bundle_from_edi
    >>> class E:
    ...     def __init__(self):
    ...         self.freq = [1.0, 2.0]
    ...         self.phase = [45.0, 50.0]
    ...         self.rho = [100.0, 120.0]
    ...
    >>> b = bundle_from_edi(E())
    >>> b.freq, b.phase[0]
    ([1.0, 2.0], 45.0)
    """
    freq = _get(obj, "freq", "frequency")
    z = _get(obj, "z", "Z")
    z_err = _get(obj, "z_err", "Z_err", "z_error")
    tip = _get(obj, "tipper", "Tipper", "tip")
    tip_err = _get(obj, "tipper_err", "tip_err")
    rho = _get(obj, "resistivity", "rho")
    pha = _get(obj, "phase", "phi", "phase_deg")
    name = _get(obj, "station", "site", "name")
    lat = _get(obj, "lat", "latitude")
    lon = _get(obj, "lon", "longitude")
    elev = _get(obj, "elev", "elevation")
    az = _get(obj, "azimuth", "az")

    # normalize tipper shape to (n, 2)
    try:
        if tip is not None and getattr(tip, "ndim", 1) == 3:
            if tip.shape[1:] == (1, 2):
                tip = tip[:, 0, :]
    except Exception:
        pass

    return TFBundle(
        freq=freq,
        z=z,
        z_err=z_err,
        tipper=tip,
        tipper_err=tip_err,
        rho=rho,
        phase=pha,
        station=name,
        lat=lat,
        lon=lon,
        elev=elev,
        azimuth=az,
    )


def _looks_collection(obj: Any) -> bool:
    r"""
    Heuristically decide if ``obj`` behaves like a collection.

    Criteria are conservative: basic sequences return True, while
    strings and bytes return False. Any iterable that does not
    expose a ``z`` attribute is considered a collection.

    Parameters
    ----------
    obj : Any
        Candidate object.

    Returns
    -------
    bool
        True when it looks like a collection of items.

    Notes
    -----
    This is a soft check used by :meth:`BundleMixin.from_edi`
    to support collections without importing container types.
    """

    if isinstance(obj, (list, tuple, set)):
        return True
    if isinstance(obj, (str, bytes)):
        return False
    if hasattr(obj, "__iter__") and not hasattr(obj, "z"):
        return True
    return False


class BundleMixin(CoreObject):
    r"""
    Mixin for single items convertible to and from bundles.

    Subclasses implement ``to_bundle`` and ``from_bundle`` to map
    between their internal representation and the neutral
    :class:`TFBundle`.

    Notes
    -----
    ``to_edi`` delegates to
    :func:`pycsamt.core.base.to_edi`, which uses the adapter
    registry to pick the right converter.

    Examples
    --------
    >>> from pycsamt.core.mixins import BundleMixin
    >>> from pycsamt.core.base import TFBundle
    >>>
    >>> class Item(BundleMixin):
    ...     def __init__(self, b): self._b = b
    ...     def to_bundle(self): return self._b
    ...     @classmethod
    ...     def from_bundle(cls, b): return cls(b)
    ...
    >>> b = TFBundle(freq=[1.0], rho=[100.], phase=[45.])
    >>> it = Item.from_bundle(b)
    >>> isinstance(it.to_bundle(), TFBundle)
    True
    """

    def to_bundle(self) -> TFBundle:  # noqa: D401
        r"""
        Return a :class:`TFBundle` representation.
        Must be implemented by subclasses.
        """
        raise NotImplementedError

    @classmethod
    def from_bundle(cls, bundle: TFBundle):  # noqa: D401
        r"""
        Construct an instance from a :class:`TFBundle`.

        Returns a new instance of the implementing class.

        Must be implemented by subclasses.
        """
        raise NotImplementedError

    @staticmethod
    def ensure_station_name(
        name: str | None,
        station_id: str | int | None,
    ) -> str:
        r"""
        Validate or synthesize a station name via policy.

        Parameters
        ----------
        name : str or None
            Preferred station/site name if available.
        station_id : str, int or None
            Identifier used to synthesize a name when needed.

        Returns
        -------
        str
            Validated or synthetic name.

        See Also
        --------
        pycsamt.core.base.ensure_station
            Underlying helper with policy lookup.
        """

        return ensure_station(name, station_id)

    def to_edi(self, *, key: str | None = None, **k: Any) -> Any:
        r"""
        Convert ``self`` to EDI or EDICollection via dispatch.

        Parameters
        ----------
        key : str, optional
            Adapter key override (e.g. ``'avg'``, ``'j'``).
        **k : Any
            Extra keyword arguments forwarded to the adapter.

        Returns
        -------
        Any
            EDI object or EDICollection, depending on adapter.

        Notes
        -----
        Dispatch is handled by
        :func:`pycsamt.core.base.to_edi`. A suitable adapter must
        have been registered, otherwise a runtime error is raised.
        """

        return to_edi(self, key=key, **k)

    @classmethod
    def from_edi(cls, edi_obj: Any):
        r"""
        Construct instance(s) from an EDI object or a collection.

        If ``edi_obj`` looks like a collection, every element is
        converted to a :class:`TFBundle` and fed to
        ``from_bundle``. Otherwise the single object is converted
        and fed to ``from_bundle``.

        Parameters
        ----------
        edi_obj : Any
            EDI object or an iterable of EDI objects.

        Returns
        -------
        Any
            Instance or list of instances of the implementing class.

        Examples
        --------
        >>> from pycsamt.core.mixins import BundleMixin, bundle_from_edi
        >>> class I(BundleMixin):
        ...     def __init__(self, b): self._b = b
        ...     def to_bundle(self): return self._b
        ...     @classmethod
        ...     def from_bundle(cls, b): return cls(b)
        ...
        >>> e = type('E', (), {'freq':[1.0], 'rho':[100.], 'phase':[45.]})()
        >>> obj = I.from_edi(e)  # doctest: +SKIP
        """

        if _looks_collection(edi_obj):
            out = []
            for it in edi_obj:  # type: ignore
                b = bundle_from_edi(it)
                out.append(cls.from_bundle(b))
            return out
        b = bundle_from_edi(edi_obj)
        return cls.from_bundle(b)


class BundleContainerMixin(BundleMixin):
    r"""
    Mixin for containers of items that implement ``to_bundle``.

    This mixin provides iteration over bundles and a convenience
    method to convert all items to an EDI collection using the
    adapter dispatch.

    Notes
    -----
    Detection of children is permissive: if the container has an
    ``items()`` method, values are inspected; otherwise the
    container is iterated directly.

    Examples
    --------
    >>> from pycsamt.core.mixins import BundleContainerMixin
    >>> from pycsamt.core.base import TFBundle
    >>>
    >>> class Bag(dict, BundleContainerMixin):
    ...     pass
    ...
    >>> _ = Bag({0: type('X',(),{'to_bundle':lambda s: TFBundle()})()})
    >>> list(_.iter_bundles())[0].__class__.__name__
    'TFBundle'
    """

    def iter_bundles(self) -> Iterator[TFBundle]:
        """ Yield contained items as :class:`TFBundle` objects."""
        if hasattr(self, "items"):
            for _, v in self.items():  # type: ignore
                if hasattr(v, "to_bundle"):
                    yield v.to_bundle()
        elif hasattr(self, "__iter__"):
            for v in self:  # type: ignore
                if hasattr(v, "to_bundle"):
                    yield v.to_bundle()

    def to_edi_collection(
        self,
        *,
        key: str | None = None,
        **k: Any,
    ) -> list[Any]:
        r"""
        Convert all child bundles through adapter dispatch.

        Parameters
        ----------
        key : str, optional
            Adapter key override for all children.
        **k : Any
            Extra keyword arguments forwarded to the adapter.

        Returns
        -------
        list of Any
            A list of EDI objects. Callers may wrap this into an
            actual :class:`~pycsamt.seg.collection.EDICollection`.

        Examples
        --------
        >>> from pycsamt.core.mixins import BundleContainerMixin
        >>> # Requires registered adapters to actually run:
        >>> # objs = MyContainer(...).to_edi_collection()  # doctest:+SKIP
        """

        edis: list[Any] = []
        for b in self.iter_bundles():
            proxy = _ProxyFromBundle(b)
            edis.append(to_edi(proxy, key=key, **k))
        return edis


class _ProxyFromBundle:
    r"""
    Tiny proxy that exposes ``to_bundle`` for a stored bundle.

    Used internally to feed :func:`pycsamt.core.base.to_edi`
    with a minimal object that satisfies the dispatch contract.
    """
    def __init__(self, b: TFBundle) -> None:
        self._b = b

    def to_bundle(self) -> TFBundle:
        return self._b
