# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Optional plotting and geospatial backend loaders."""

from __future__ import annotations

from typing import Any


def require_plotly() -> Any:
    """Return :mod:`plotly.graph_objects`.

    Plotly is used only by the interactive map backend.
    It is loaded at the backend boundary with a clear error.
    """
    try:
        import plotly.graph_objects as go
    except ImportError as exc:
        raise ImportError(
            "Plotly is required for interactive pyCSAMT "
            "maps. "
            "Install pycsamt with its map dependencies."
        ) from exc
    return go


def require_plotly_subplots() -> Any:
    """Return Plotly's subplot factory."""
    try:
        from plotly.subplots import make_subplots
    except ImportError as exc:
        raise ImportError(
            "Plotly is required for multi-panel pyCSAMT "
            "maps. "
            "Install pycsamt with its map dependencies."
        ) from exc
    return make_subplots


def require_matplotlib_pyplot() -> Any:
    """Return :mod:`matplotlib.pyplot`."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError(
            "Matplotlib is required for static pyCSAMT maps. "
            "Install pycsamt with its map dependencies."
        ) from exc
    return plt


def require_pyproj_transformer() -> Any:
    """Return :class:`pyproj.Transformer`."""
    try:
        from pyproj import Transformer
    except ImportError as exc:
        raise ImportError(
            "pyproj is required for CRS reprojection in "
            "pycsamt.map."
        ) from exc
    return Transformer


def scipy_griddata() -> Any | None:
    """Return SciPy ``griddata`` when available.

    Contours can fall back to a small nearest-neighbour grid.
    Missing SciPy is not fatal for this helper.
    """
    try:
        from scipy.interpolate import griddata
    except ImportError:
        return None
    return griddata


def explain_plotly_image_error(
    _exc: Exception,
) -> ImportError:
    """Build a clear Plotly static image export error."""
    return ImportError(
        "Plotly static image export requires the kaleido "
        "package. "
        "Install kaleido or export to HTML instead."
    )
