# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli.commands.config
============================

Click command group for reading and writing pyCSAMT configuration.

Configuration is persisted in ``~/.pycsamt.toml`` and applied to all
10 API singletons at every CLI start-up via :func:`load_all_config`.

Sub-command modules
-------------------
_base      TOML I/O, key parsing, per-section appliers, singleton getters.
config     ``pycsamt config`` root group and all sub-commands:
             list / get / set / unset / reset / show / env /
             style / interp / agent.
"""

from ._base import (
    load_all_config,  # noqa: F401  (re-exported for cli/_base.py)
)
from .config import config  # noqa: F401

__all__ = ["config", "load_all_config"]
