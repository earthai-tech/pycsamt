# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Figure export helpers for :mod:`pycsamt.map`."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from ._backends import explain_plotly_image_error
from .config import ExportOptions

PLOTLY_IMAGE_KWARGS = {"scale", "width", "height"}


def write_html(
    fig: Any,
    path: str | Path,
    **kwargs: Any,
) -> Path:
    """Write a Plotly-like figure to HTML."""
    out = _prepare_output_path(path)
    if not hasattr(fig, "write_html"):
        raise TypeError(
            "HTML export requires a Plotly-like figure with "
            "a write_html method."
        )
    fig.write_html(str(out), **kwargs)
    return out


def write_image(
    fig: Any,
    path: str | Path,
    **kwargs: Any,
) -> Path:
    """Write a Plotly-like image."""
    out = _prepare_output_path(path)
    if not hasattr(fig, "write_image"):
        raise TypeError(
            "Static image export requires a Plotly-like "
            "figure "
            "with a write_image method, or use save_png for "
            "Matplotlib figures."
        )
    try:
        fig.write_image(str(out), **kwargs)
    except (ImportError, RuntimeError, ValueError) as exc:
        if _looks_like_kaleido_error(exc):
            raise explain_plotly_image_error(exc) from exc
        raise
    return out


def save_png(
    fig: Any,
    path: str | Path,
    **kwargs: Any,
) -> Path:
    """Save a Matplotlib-like or Plotly-like figure as PNG."""
    out = _prepare_output_path(path)
    if hasattr(fig, "savefig"):
        fig.savefig(
            str(out),
            format="png",
            **_matplotlib_kwargs(kwargs),
        )
        return out
    return write_image(fig, out, format="png", **kwargs)


def export_figure(fig: Any, options: ExportOptions) -> Path:
    """Export *fig* according to *options*."""
    path = Path(options.path)
    suffix = _export_suffix(path, options.format)
    if suffix in {".html", ".htm"}:
        return write_html(
            fig,
            path,
            include_plotlyjs=options.include_plotlyjs,
        )
    if suffix == ".png":
        kwargs: dict[str, Any] = {"scale": options.scale}
        if options.width is not None:
            kwargs["width"] = options.width
        if options.height is not None:
            kwargs["height"] = options.height
        return save_png(fig, path, **kwargs)
    if suffix == ".json":
        return write_json(fig, path)
    if suffix == ".dict":
        return write_dict(fig, path)
    if not suffix:
        msg = (
            "Export path must include an extension or "
            "ExportOptions.format must be set."
        )
        raise ValueError(msg)
    return write_image(fig, path)


def write_json(fig: Any, path: str | Path) -> Path:
    """Write a Plotly-like figure specification as JSON."""
    out = _prepare_output_path(path)
    if hasattr(fig, "write_json"):
        fig.write_json(str(out))
        return out
    if hasattr(fig, "to_json"):
        out.write_text(fig.to_json(), encoding="utf-8")
        return out
    out.write_text(
        json.dumps(figure_to_dict(fig), indent=2),
        encoding="utf-8",
    )
    return out


def write_dict(fig: Any, path: str | Path) -> Path:
    """Write a figure dictionary as JSON text."""
    out = _prepare_output_path(path)
    out.write_text(
        json.dumps(figure_to_dict(fig), indent=2),
        encoding="utf-8",
    )
    return out


def figure_to_dict(fig: Any) -> dict[str, Any]:
    """Return a serializable figure dictionary."""
    if hasattr(fig, "to_dict"):
        data = fig.to_dict()
        if isinstance(data, dict):
            return data
    if hasattr(fig, "to_plotly_json"):
        data = fig.to_plotly_json()
        if isinstance(data, dict):
            return data
    msg = "Figure dictionary export requires to_dict or to_plotly_json."
    raise TypeError(msg)


def _prepare_output_path(path: str | Path) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    return out


def _export_suffix(
    path: Path,
    fmt: str | None,
) -> str:
    if fmt:
        value = fmt.lower().strip()
        return value if value.startswith(".") else f".{value}"
    return path.suffix.lower()


def _matplotlib_kwargs(
    kwargs: dict[str, Any],
) -> dict[str, Any]:
    return {
        key: value
        for key, value in kwargs.items()
        if key not in PLOTLY_IMAGE_KWARGS
    }


def _looks_like_kaleido_error(exc: Exception) -> bool:
    msg = str(exc).lower()
    return "kaleido" in msg or "orca" in msg
