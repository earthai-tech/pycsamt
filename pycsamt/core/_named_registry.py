# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Generic name/alias registry bookkeeping shared across pyCSAMT.

This module carries no domain knowledge about transfer functions,
EDI/EMTF formats, or airborne technologies. It exists only to hold the
alias-collision bookkeeping and detector-then-extension matching
algorithm that :mod:`pycsamt.io.formats` (serialization formats) and
:mod:`pycsamt.airborne.registry` (native delivery formats and
technologies) each independently implemented, so that logic now has one
implementation instead of several near-identical copies. Domain-specific
item types, error classes, and lookup/detection policy stay in each
caller.
"""

from __future__ import annotations

from typing import Any, Generic, TypeVar

__all__ = [
    "NamedRegistryError",
    "NamedRegistry",
    "normalize_key",
    "normalize_extension",
]

T = TypeVar("T")


class NamedRegistryError(ValueError):
    """Default error raised by :class:`NamedRegistry` bookkeeping.

    Callers that need a domain-specific exception type should pass
    ``error_cls`` to :class:`NamedRegistry` instead of catching this.
    """


def normalize_key(value: Any) -> str:
    """Return a canonical, case-insensitive registry key.

    Returns an empty string for empty/whitespace-only input rather than
    raising, since some callers use this purely for tag/name comparison
    on values that may legitimately be blank (e.g. an unset ``subtype``).
    :meth:`NamedRegistry.register` rejects an empty canonical name on its
    own, where "empty" genuinely is invalid.
    """
    return str(value).strip().lower().replace("-", "_").replace(" ", "_")


def normalize_extension(value: Any) -> str:
    """Return a canonical, dot-prefixed, lowercase filename extension."""
    text = str(value).strip().lower()
    if not text:
        return ""
    return text if text.startswith(".") else f".{text}"


class NamedRegistry(Generic[T]):
    """Register items under a canonical name plus optional aliases.

    Parameters
    ----------
    kind : str
        Human-readable noun used in error messages, e.g. ``"format"``.
    error_cls : type of Exception, default=NamedRegistryError
        Exception type raised on an unresolved name/alias collision.
        Callers pass their own domain-specific error type so this
        primitive never dictates the caller's public exception surface.

    Notes
    -----
    Detection helpers (:meth:`match_by_detector`,
    :meth:`match_by_extension`) expect each registered item to expose
    ``.detector`` (an optional callable) and ``.extensions`` (a tuple of
    strings) attributes; they place no other requirement on the item
    type. :meth:`register`/:meth:`get`/:meth:`all` place no requirement
    on the item type at all.
    """

    def __init__(
        self,
        *,
        kind: str,
        error_cls: type[Exception] = NamedRegistryError,
    ) -> None:
        self._kind = str(kind)
        self._error_cls = error_cls
        self._items: dict[str, T] = {}
        self._aliases: dict[str, str] = {}

    def register(
        self,
        name: str,
        item: T,
        *,
        aliases: tuple[str, ...] = (),
        replace: bool = False,
    ) -> str:
        """Register *item* under *name* and *aliases*.

        Returns
        -------
        str
            The canonical (normalized) name *item* was registered under.

        Raises
        ------
        Exception
            An instance of ``error_cls`` if *name* is already registered
            and ``replace`` is ``False``, or if any alias already
            belongs to a different canonical name.
        """
        canonical = normalize_key(name)
        if not canonical:
            raise ValueError("name must be non-empty")
        alias_keys = tuple(normalize_key(alias) for alias in aliases)
        all_keys = (canonical, *alias_keys)

        if canonical in self._items and not replace:
            raise self._error_cls(
                f"{self._kind} already registered: {canonical!r}"
            )
        if not replace:
            for key in all_keys:
                owner = self._aliases.get(key)
                if owner is not None and owner != canonical:
                    raise self._error_cls(
                        f"{self._kind} alias {key!r} already belongs "
                        f"to {owner!r}"
                    )
        elif canonical in self._items:
            for key, owner in tuple(self._aliases.items()):
                if owner == canonical:
                    del self._aliases[key]

        self._items[canonical] = item
        for key in all_keys:
            self._aliases[key] = canonical
        return canonical

    def get(self, name: str) -> T | None:
        """Return the item registered under *name* or one of its
        aliases, or ``None`` if unregistered."""
        canonical = self._aliases.get(normalize_key(name))
        return None if canonical is None else self._items.get(canonical)

    def all(self) -> tuple[T, ...]:
        """Return all registered items in registration order."""
        return tuple(self._items.values())

    def names(self) -> tuple[str, ...]:
        """Return all canonical names, sorted."""
        return tuple(sorted(self._items))

    def match_by_detector(self, source: Any) -> list[str]:
        """Return canonical names whose ``.detector(source)`` is truthy.

        Detector exceptions (``OSError``, ``TypeError``, ``ValueError``)
        are swallowed per-item so one broken detector cannot block
        detection for the rest of the registry.
        """
        matches: list[str] = []
        for name, item in self._items.items():
            detector = getattr(item, "detector", None)
            if detector is None:
                continue
            try:
                if detector(source):
                    matches.append(name)
            except (OSError, TypeError, ValueError):
                continue
        return matches

    def match_by_extension(self, suffix: str) -> list[str]:
        """Return canonical names whose ``.extensions`` contains
        *suffix*."""
        return [
            name
            for name, item in self._items.items()
            if suffix in getattr(item, "extensions", ())
        ]
