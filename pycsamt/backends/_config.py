# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Backend configuration — global state for the active AI backend.

The active backend is resolved with this priority chain:

1. Explicit call to :func:`set_backend` in the current session.
2. ``PYCSAMT_AI_BACKEND`` environment variable (``torch`` / ``tensorflow``
   / ``auto``).
3. ``~/.pycsamt/config.json`` key ``"ai_backend"``.
4. Auto-detection: first available framework in ``['torch', 'tensorflow']``.

No ML framework is imported here.  Actual framework objects are loaded
lazily by the backend implementations.
"""
from __future__ import annotations

import json
import os
import threading
from pathlib import Path
from typing import Optional

__all__ = ["BackendConfig"]

_VALID_NAMES = frozenset({"torch", "tensorflow", "auto", "none"})
_LOCK = threading.Lock()


class BackendConfig:
    """Thread-safe global backend configuration (singleton)."""

    _instance: Optional["BackendConfig"] = None

    def __new__(cls) -> "BackendConfig":
        with _LOCK:
            if cls._instance is None:
                obj = object.__new__(cls)
                obj._backend_name: Optional[str] = None  # None = not yet resolved
                cls._instance = obj
        return cls._instance

    # ─── public API ───────────────────────────────────────────────────────

    @property
    def backend_name(self) -> str:
        """
        Return the resolved backend name (``'torch'`` or ``'tensorflow'``).

        Resolves lazily on first access.
        """
        with _LOCK:
            if self._backend_name is None or self._backend_name == "auto":
                self._backend_name = self._resolve()
            return self._backend_name

    def set(self, name: str) -> None:
        """
        Set the active backend.

        Parameters
        ----------
        name : {'torch', 'tensorflow', 'auto'}

        Raises
        ------
        ValueError
            If *name* is not recognised.
        """
        name = name.lower().strip()
        if name not in _VALID_NAMES:
            raise ValueError(
                f"Unknown backend {name!r}.  "
                f"Choose from: {sorted(_VALID_NAMES)}"
            )
        with _LOCK:
            if name == "auto":
                self._backend_name = self._resolve()
            else:
                self._backend_name = name

    def reset(self) -> None:
        """Clear the resolved backend (forces re-detection on next access)."""
        with _LOCK:
            self._backend_name = None

    # ─── resolution chain ─────────────────────────────────────────────────

    def _resolve(self) -> str:
        """Walk the priority chain and return a concrete backend name."""
        # 1. Environment variable
        env = os.environ.get("PYCSAMT_AI_BACKEND", "").strip().lower()
        if env and env in _VALID_NAMES and env != "auto":
            return env

        # 2. Config file
        cfg_name = self._read_config_file()
        if cfg_name and cfg_name in _VALID_NAMES and cfg_name != "auto":
            return cfg_name

        # 3. Auto-detect
        return self._auto_detect()

    @staticmethod
    def _read_config_file() -> Optional[str]:
        """Read ai_backend from ~/.pycsamt/config.json if it exists."""
        fpath = Path.home() / ".pycsamt" / "config.json"
        try:
            with open(fpath) as f:
                d = json.load(f)
            return str(d.get("ai_backend", "")).lower() or None
        except (FileNotFoundError, json.JSONDecodeError, PermissionError):
            return None

    @staticmethod
    def _auto_detect() -> str:
        """Return the first available framework, or ``'none'`` if none found."""
        from ._detect import detect_available_backends
        available = detect_available_backends()
        if not available:
            return "none"
        return available[0]

    def write_config_file(self, name: str) -> None:
        """
        Persist the backend choice to ``~/.pycsamt/config.json``.

        Parameters
        ----------
        name : str
            Backend name to save.
        """
        name = name.lower()
        if name not in _VALID_NAMES:
            raise ValueError(f"Cannot persist unknown backend {name!r}.")
        cfg_dir = Path.home() / ".pycsamt"
        cfg_dir.mkdir(parents=True, exist_ok=True)
        cfg_file = cfg_dir / "config.json"
        data: dict = {}
        if cfg_file.exists():
            try:
                with open(cfg_file) as f:
                    data = json.load(f)
            except (json.JSONDecodeError, PermissionError):
                data = {}
        data["ai_backend"] = name
        with open(cfg_file, "w") as f:
            json.dump(data, f, indent=2)
