# -*- coding: utf-8 -*-
from __future__ import annotations

import dataclasses as _dc
import sys as _sys
from typing import Any, Callable, Dict

__all__ = ["PY310_PLUS", "dc", "DATACLASS_SLOTS"]

PY310_PLUS: bool = _sys.version_info >= (3, 10)

# Convenience dict for the SO-style decorator usage:
#   @dataclasses.dataclass(**DATACLASS_SLOTS)
DATACLASS_SLOTS: Dict[str, Any] = {"slots": True} if \
    PY310_PLUS else {}


def _filter_kwargs(kwargs: Dict[str, Any]) -> Dict[str, Any]:
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
