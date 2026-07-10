# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Figure export callbacks."""

from __future__ import annotations

import base64
import io

from dash import ALL, Input, Output, State, ctx, dcc
from dash.exceptions import PreventUpdate

from .._ids import IDs


def _b64_to_bytes(b64: str) -> bytes:
    return base64.b64decode(b64)


def _convert_png_to(
    b64_png: str,
    fmt: str,
) -> bytes:
    """
    Convert a base64 PNG to the requested format.
    Returns raw bytes for the target format.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    raw = base64.b64decode(b64_png)
    arr = np.frombuffer(raw, dtype=np.uint8)

    try:
        import cv2
        img = cv2.imdecode(arr, cv2.IMREAD_COLOR)
        img_rgb = cv2.cvtColor(
            img, cv2.COLOR_BGR2RGB
        )
    except ImportError:
        from PIL import Image
        img_rgb = np.array(
            Image.open(io.BytesIO(raw))
        )

    fig, ax = plt.subplots(
        figsize=(
            img_rgb.shape[1] / 100,
            img_rgb.shape[0] / 100,
        )
    )
    ax.imshow(img_rgb)
    ax.axis("off")
    fig.tight_layout(pad=0)

    buf = io.BytesIO()
    fmt_map = {
        "png": "png",
        "svg": "svg",
        "eps": "eps",
        "pdf": "pdf",
    }
    fig.savefig(
        buf,
        format=fmt_map.get(fmt, "png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)
    buf.seek(0)
    return buf.read()


_MIME = {
    "png": "image/png",
    "svg": "image/svg+xml",
    "eps": "application/postscript",
    "pdf": "application/pdf",
}


def register_preview(app) -> None:

    @app.callback(
        Output(IDs.EXPORT_DL, "data"),
        Input(
            {"type": "am-export", "fmt": ALL},
            "n_clicks",
        ),
        State(IDs.MODAL_FIG_KEY, "data"),
        State(IDs.STORE_FIGS, "data"),
        State(IDs.STORE_SETTINGS, "data"),
        prevent_initial_call=True,
    )
    def export_figure(
        n_clicks, fig_key, fig_store, settings
    ):
        if not any(n_clicks):
            raise PreventUpdate
        triggered = ctx.triggered_id
        if not triggered:
            raise PreventUpdate
        fmt = triggered.get("fmt", "png")
        if (
            not fig_key
            or not fig_store
            or fig_key not in fig_store
        ):
            raise PreventUpdate

        info = fig_store[fig_key]
        b64 = info.get("b64", "")
        title = info.get("title", "figure")
        safe = (
            title.replace(" ", "_")
            .replace("/", "_")
            .lower()
        )
        fname = f"{safe}.{fmt}"

        if fmt == "png":
            raw = base64.b64decode(b64)
            return dcc.send_bytes(
                raw, filename=fname
            )

        # convert for non-PNG formats
        raw = _convert_png_to(b64, fmt)
        return dcc.send_bytes(
            raw, filename=fname
        )
