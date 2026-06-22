# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Base classes for inversion backends."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any

from ..api.property import PyCSAMTObject

from .config import InversionConfig
from .data import EMData
from .results import InversionResult

__all__ = ["BaseInversionBackend"]


class BaseInversionBackend(PyCSAMTObject, ABC):
    """Abstract interface implemented by inversion backends."""

    name = "base"
    supports: tuple[tuple[str, str], ...] = ()

    def __init__(self, config: InversionConfig) -> None:
        self.config = config

    def prepare_data(self, data: Any | None = None) -> EMData:
        """Coerce input data into :class:`EMData` when possible."""
        raw = self.config.data if data is None else data
        return EMData.coerce(raw, method=self.config.method)

    def check_supported(self) -> None:
        """Raise if this backend does not support the configured method/dim."""
        if not self.supports:
            return
        key = (self.config.method, self.config.dimension)
        wild = (self.config.method, "*")
        if key not in self.supports and wild not in self.supports:
            raise NotImplementedError(
                f"{self.name!r} backend does not support "
                f"{self.config.method}/{self.config.dimension}."
            )

    @abstractmethod
    def run(self, data: Any | None = None) -> InversionResult:
        """Run or prepare the inversion and return an :class:`InversionResult`."""
