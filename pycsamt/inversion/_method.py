# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Shared method-specific workflow helpers."""

from __future__ import annotations

from typing import Any

from .config import InversionConfig
from .workflow import InversionWorkflow

__all__ = ["ConfiguredInversionWorkflow"]


class ConfiguredInversionWorkflow(InversionWorkflow):
    """Workflow with method, dimension, and backend defaults.

    Method namespace classes should inherit this helper and set
    ``method``, ``dimension``, and ``backend``. This keeps convenience
    wrappers thin and avoids repeating constructor logic across MT, AMT,
    EMAP, and TDEM modules.
    """

    method: str = "mt"
    dimension: str = "1d"
    backend: str = "builtin"

    def __init__(self, data: Any = None, **kwargs: Any) -> None:
        kwargs.setdefault("method", self.method)
        kwargs.setdefault("dimension", self.dimension)
        kwargs.setdefault("backend", self.backend)
        if data is not None:
            kwargs.setdefault("data", data)
        super().__init__(InversionConfig(**kwargs))
