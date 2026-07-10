# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Minimal bases for the Jones (J-format) subsystem.
"""
from __future__ import annotations

import io
from collections.abc import Iterable
from pathlib import Path
from typing import Any, ClassVar

from ..utils.validation import has_read
from .config import ENCODING_DEFAULT
from .utils import iter_lines

__all__ = ["BaseJones", "JComponentBase"]

class BaseJones:
    """Root base for all Jones objects.

    Provides a common interface and light I/O hooks. Subclasses
    must implement ``from_lines``, ``read`` and ``write``.
    """

    _repr_keys: ClassVar[list[str]] = ["path", "encoding"]

    def __init__(self, *, verbose: int | bool = 0) -> None:
        self.verbose: int | bool = verbose
        self.path: Path | None = None
        self.encoding: str = ENCODING_DEFAULT
        self._has_read: bool = False


    @classmethod
    def from_file(
        cls,
        obj: io.TextIOBase | str | Path,
        *,
        encoding: str = ENCODING_DEFAULT,
        verbose: int | bool = 0,
        **kws: Any,
    ) -> BaseJones:
        lines = iter_lines(obj, encoding=encoding)
        inst = cls.from_lines(lines, verbose=verbose, **kws)
        if isinstance(obj, (str, Path)):
            inst.path = Path(obj)
        inst.encoding = encoding
        return inst

    @classmethod
    def from_lines(
        cls,
        lines: Iterable[str],
        *,
        verbose: int | bool = 0,
        **kws: Any,
    ) -> BaseJones:
        raise NotImplementedError(
            f"{cls.__name__}.from_lines() not implemented"
        )


    def read(self, *args: Any, **kws: Any) -> BaseJones:
        raise NotImplementedError(
            f"{self.__class__.__name__}.read() not implemented"
        )

    def write(
        self,
        obj: io.TextIOBase | str | Path,
        *,
        encoding: str | None = None,
        **kws: Any,
    ) -> None:
        raise NotImplementedError(
            f"{self.__class__.__name__}.write() not implemented"
        )


    def __has_read__(self) -> bool:
        return bool(self._has_read)

    def _mark_read(self, state: bool = True) -> None:
        self._has_read = bool(state)

    def require_read(
        self,
        *,
        attributes: str | list[str] | None = None,
        msg: str | None = None,
    ) -> None:
        has_read(self, attributes=attributes, msg=msg)


    def to_dict(self) -> dict[str, Any]:
        out: dict[str, Any] = {}
        for k, v in vars(self).items():
            if k.startswith("_"):
                continue
            if callable(v):
                continue
            out[k] = v
        return out

    def summary(self) -> str:
        parts = []
        for k in self._repr_keys:
            if hasattr(self, k):
                v = getattr(self, k)
                if v is not None:
                    parts.append(f"{k}={v!r}")
        return f"{self.__class__.__name__}({', '.join(parts)})"

    def __str__(self) -> str:
        return self.summary()

    def __repr__(self) -> str:
        s = self.summary()
        if s.endswith("())"):
            return f"{self.__class__.__name__}({self.to_dict()})"
        return s


class JComponentBase(BaseJones):
    """Common base for concrete components.

    All components must accept ``verbose`` in ``__init__`` and
    must implement ``from_file``, ``from_lines``, ``read`` and
    ``write``.
    """

    _repr_keys: ClassVar[list[str]] = [
        "name",
        "station",
        "kind",
        "comp",
        "units",
        "n",
        "shape",
    ]

