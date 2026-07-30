from __future__ import annotations

import dataclasses as _dc
import sys as _sys
from collections.abc import Iterable, Iterator
from itertools import zip_longest as _zip_longest
from typing import Any, Callable, TypeVar

__all__ = ["PY310_PLUS", "dc", "DATACLASS_SLOTS", "zip_strict"]

PY310_PLUS: bool = _sys.version_info >= (3, 10)

_T = TypeVar("_T")
_U = TypeVar("_U")

if PY310_PLUS:

    def zip_strict(
        a: Iterable[_T], b: Iterable[_U]
    ) -> Iterator[tuple[_T, _U]]:
        """``zip(a, b, strict=True)``, back-ported for Python 3.9
        (which raises ``TypeError: zip() takes no keyword arguments``
        for the ``strict`` kwarg)."""
        return zip(a, b, strict=True)

else:

    def zip_strict(
        a: Iterable[_T], b: Iterable[_U]
    ) -> Iterator[tuple[_T, _U]]:
        """``zip(a, b, strict=True)``, back-ported for Python 3.9
        (which raises ``TypeError: zip() takes no keyword arguments``
        for the ``strict`` kwarg)."""
        sentinel = object()
        for x, y in _zip_longest(a, b, fillvalue=sentinel):
            if x is sentinel or y is sentinel:
                raise ValueError(
                    "zip_strict() argument lengths did not match."
                )
            yield x, y


# Convenience dict for the SO-style decorator usage:
#   @dataclasses.dataclass(**DATACLASS_SLOTS)
DATACLASS_SLOTS: dict[str, Any] = {"slots": True} if PY310_PLUS else {}


def _filter_kwargs(kwargs: dict[str, Any]) -> dict[str, Any]:
    """Remove dataclass kwargs not supported on this Python."""
    if not PY310_PLUS:
        # 'slots', 'kw_only', and 'match_args' landed in 3.10.
        kwargs.pop("slots", None)
        kwargs.pop("kw_only", None)
        kwargs.pop("match_args", None)
    return kwargs


def dc(**kwargs: Any) -> Callable[[type], type]:
    """Version-safe dataclass decorator factory.

    Usage:
        from pycsamt.compat.dataclasses import dc

        @dc(slots=True, frozen=True)
        class MyData:
            ...
    """
    return _dc.dataclass(**_filter_kwargs(dict(kwargs)))
