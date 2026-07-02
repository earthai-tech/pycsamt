# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Map-view app launch helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ._core import MapData
from .loader import load_lines

_GUI_IMPORT_ERROR = (
    "The map-view platform requires the GUI dependencies "
    "(dash, dash-bootstrap-components). Install pycsamt with "
    "its app extras."
)


def launch_app(
    source: Any = None,
    *,
    host: str = "127.0.0.1",
    port: int = 8770,
    debug: bool = False,
    open_browser: bool = True,
    theme: str = "light",
    backend: str = "plotly",
    detect: str = "folder",
    recursive: bool = True,
    **load_kwargs: Any,
) -> None:
    """Launch the interactive map-view app.

    Parameters
    ----------
    source :
        Optional survey source.
        ``None`` starts the welcome screen.
        Otherwise accepts data, paths, sites, or mappings.
    host, port, debug, open_browser :
        Forwarded to :func:`pycsamt.app.mapview.launch`.
    theme, backend :
        Used when building a :class:`~pycsamt.map.MapView`.
    detect, recursive :
        Used for path/folder inputs.
    **load_kwargs :
        Extra options forwarded to ``MapView.from_folder`` or
        :func:`load_lines`.
    """
    launch_kwargs = {
        "host": host,
        "port": int(port),
        "debug": bool(debug),
        "open_browser": bool(open_browser),
    }
    if source is None:
        _launch_empty(**launch_kwargs)
        return

    view = _as_map_view(
        source,
        theme=theme,
        backend=backend,
        detect=detect,
        recursive=recursive,
        **load_kwargs,
    )
    view.launch(**launch_kwargs)


def launch_mapview(source: Any = None, **kwargs: Any) -> None:
    """Alias for :func:`launch_app`."""
    launch_app(source, **kwargs)


def open_app(source: Any = None, **kwargs: Any) -> None:
    """Backward-compatible alias for :func:`launch_app`."""
    launch_app(source, **kwargs)


def _launch_empty(**launch_kwargs: Any) -> None:
    try:
        from pycsamt.app.mapview import launch
    except ImportError as exc:  # pragma: no cover
        raise ImportError(_GUI_IMPORT_ERROR) from exc
    launch(**launch_kwargs)


def _as_map_view(
    source: Any,
    *,
    theme: str,
    backend: str,
    detect: str,
    recursive: bool,
    **load_kwargs: Any,
):
    from .view import MapView

    if isinstance(source, MapView):
        return source
    if isinstance(source, MapData):
        return MapView(source, theme=theme, backend=backend)
    if isinstance(source, (str, Path)):
        return MapView.from_folder(
            source,
            detect=detect,
            recursive=recursive,
            theme=theme,
            backend=backend,
            **load_kwargs,
        )
    return MapView(
        load_lines(source, **load_kwargs),
        theme=theme,
        backend=backend,
    )
