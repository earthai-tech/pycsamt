"""Global configuration for pyCSAMT API view wrappers."""

from __future__ import annotations

import copy
import os
from collections.abc import Callable, Generator
from contextlib import contextmanager
from typing import Any

__all__ = [
    "APIViewConfig",
    "PYCSAMT_API_VIEW",
    "configure_api_view",
    "reset_api_view",
]

ViewWrapper = Callable[..., Any]

_PYCSAMT_BACKENDS = {"api", "apiframe", "pycsamt", "view"}
_PANDAS_BACKENDS = {"false", "none", "off", "pandas", "raw"}


def _normalise_backend(value: Any) -> str | ViewWrapper:
    if callable(value):
        return value
    if isinstance(value, bool):
        return "pycsamt" if value else "pandas"
    if value is None:
        return "pandas"
    name = str(value).strip().lower()
    if name in _PYCSAMT_BACKENDS:
        return "pycsamt"
    if name in _PANDAS_BACKENDS:
        return "pandas"
    msg = (
        "api view backend must be 'pycsamt', 'pandas'/False, "
        "or a callable wrapper."
    )
    raise ValueError(msg)


class APIViewConfig:
    """Package-wide policy for dataframe-like public API views.

    The backend controls what :func:`pycsamt.api.view.wrap_frame` returns:

    * ``"pycsamt"`` / ``True``: return the built-in ``APIFrame``.
    * ``"pandas"`` / ``False``: return the original dataframe-like object.
    * callable: call ``wrapper(data, **metadata)`` and return its result.
    """

    def __init__(self) -> None:
        self._init_defaults()

    def _init_defaults(self) -> None:
        self.backend: str | ViewWrapper = _normalise_backend(
            os.environ.get("PYCSAMT_API_VIEW", "pycsamt")
        )

    def configure(self, **kw: Any) -> None:
        """Configure the view policy."""
        if "backend" in kw:
            self.backend = _normalise_backend(kw.pop("backend"))
        if "enabled" in kw:
            self.backend = _normalise_backend(kw.pop("enabled"))
        if "wrapper" in kw:
            self.backend = _normalise_backend(kw.pop("wrapper"))
        if kw:
            unknown = ", ".join(sorted(kw))
            raise ValueError(f"unknown API view config keys: {unknown}")

    @contextmanager
    def context(self, **kw: Any) -> Generator[APIViewConfig, None, None]:
        """Temporarily override API view settings."""
        snapshot = copy.copy(self.backend)
        try:
            if kw:
                self.configure(**kw)
            yield self
        finally:
            self.backend = snapshot

    def reset(self) -> None:
        """Reset to defaults."""
        self._init_defaults()

    def wrap_frame(self, data: Any, **metadata: Any) -> Any:
        """Apply the configured dataframe wrapper."""
        backend = self.backend
        if backend == "pandas":
            return data
        if callable(backend):
            return backend(data, **metadata)
        from .frame import _default_wrap_frame

        return _default_wrap_frame(data, **metadata)

    def enabled(self) -> bool:
        """Return whether API view wrapping is enabled."""
        return self.backend != "pandas"

    def summary(self) -> str:
        """Return a compact configuration summary."""
        backend = self.backend
        if callable(backend):
            name = getattr(backend, "__name__", backend.__class__.__name__)
        else:
            name = str(backend)
        return f"APIViewConfig(backend={name!r})"

    def __repr__(self) -> str:
        return self.summary()


PYCSAMT_API_VIEW = APIViewConfig()


def configure_api_view(**kw: Any) -> None:
    """Configure the global API view backend."""
    PYCSAMT_API_VIEW.configure(**kw)


def reset_api_view() -> None:
    """Reset the global API view backend."""
    PYCSAMT_API_VIEW.reset()
