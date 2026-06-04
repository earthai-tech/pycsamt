# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""High-level inversion workflow entry point."""

from __future__ import annotations

from typing import Any

from ..api.property import PyCSAMTObject

from .config import InversionConfig
from .results import InversionResult

__all__ = ["InversionWorkflow", "run_inversion"]


class InversionWorkflow(PyCSAMTObject):
    """Execute a configured EM inversion backend."""

    def __init__(self, config: InversionConfig | dict[str, Any] | None = None, **kw: Any):
        if config is None:
            config = InversionConfig(**kw)
        elif isinstance(config, dict):
            values = dict(config)
            values.update(kw)
            config = InversionConfig(**values)
        elif kw:
            config = config.clone(**kw)
        self.config = config
        self.backend = config.to_backend()

    def run(self, data: Any | None = None) -> InversionResult:
        """Run the selected backend."""
        return self.backend.run(data=data)


def run_inversion(config: InversionConfig | dict[str, Any] | None = None, **kw: Any) -> InversionResult:
    """Convenience function equivalent to ``InversionWorkflow(...).run()``."""
    return InversionWorkflow(config, **kw).run()
