# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Base class shared by all ModEM model objects."""

from __future__ import annotations

from pathlib import Path
from typing import Union

from ...log.logger import get_logger as _get_logger

PathLike = Union[str, Path]

__all__ = ["ModEmBase"]


class ModEmBase:
    """Lightweight base for all ModEM objects.

    Provides
    --------
    verbose : int
        Verbosity level (0 = silent, 1 = info, 2 = debug).
    logger  : logging.Logger
        Per-instance logger named after the concrete class.
    """

    def __init__(self, *args, verbose: int | bool = 0, logger=None, **kwargs):
        try:
            super().__init__(*args, **kwargs)
        except Exception:
            pass

        self.verbose: int = int(verbose)
        name = f"{self.__class__.__module__}.{self.__class__.__qualname__}"
        self.logger = logger if logger is not None else _get_logger(name)

    def __repr__(self) -> str:
        cls = self.__class__.__name__
        attrs = {
            k: v
            for k, v in self.__dict__.items()
            if not k.startswith("_") and k not in ("logger", "verbose")
        }
        pairs = ", ".join(f"{k}={v!r}" for k, v in list(attrs.items())[:4])
        return f"{cls}({pairs})"

    def __str__(self) -> str:
        return self.__repr__()
