"""Backward-compatible import path for PyCSAMT typing helpers.

The v2 typing definitions live in :mod:`pycsamt.api.typing`.
This module remains as a thin shim for older modules that still import
``pycsamt._typing`` during the transition.
"""

from __future__ import annotations

from .api.typing import *  # noqa: F403

