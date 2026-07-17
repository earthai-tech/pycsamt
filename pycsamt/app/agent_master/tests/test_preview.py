# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.agent_master.callbacks.preview."""

from __future__ import annotations

import base64
import io

import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")


def _tiny_png_b64() -> str:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(1, 1))
    ax.plot([0, 1], [0, 1])
    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


class TestB64ToBytes:
    def test_decodes_base64_string(self):
        from pycsamt.app.agent_master.callbacks.preview import (
            _b64_to_bytes,
        )

        raw = _b64_to_bytes(base64.b64encode(b"hello").decode())
        assert raw == b"hello"


class TestConvertPngTo:
    def test_converts_to_svg(self):
        from pycsamt.app.agent_master.callbacks.preview import (
            _convert_png_to,
        )

        out = _convert_png_to(_tiny_png_b64(), "svg")
        assert out.strip().startswith(b"<?xml") or b"<svg" in out

    def test_converts_to_pdf(self):
        from pycsamt.app.agent_master.callbacks.preview import (
            _convert_png_to,
        )

        out = _convert_png_to(_tiny_png_b64(), "pdf")
        assert out.startswith(b"%PDF")

    def test_converts_to_eps(self):
        from pycsamt.app.agent_master.callbacks.preview import (
            _convert_png_to,
        )

        out = _convert_png_to(_tiny_png_b64(), "eps")
        assert b"%!" in out[:20]

    def test_unknown_format_defaults_to_png(self):
        from pycsamt.app.agent_master.callbacks.preview import (
            _convert_png_to,
        )

        out = _convert_png_to(_tiny_png_b64(), "bogus")
        assert out[:8] == b"\x89PNG\r\n\x1a\n"


class TestMimeMap:
    def test_covers_all_export_formats(self):
        from pycsamt.app.agent_master.callbacks.preview import _MIME

        assert _MIME["png"] == "image/png"
        assert _MIME["svg"] == "image/svg+xml"
        assert _MIME["eps"] == "application/postscript"
        assert _MIME["pdf"] == "application/pdf"


class TestRegisterPreview:
    def test_register_preview_is_callable(self):
        from pycsamt.app.agent_master.callbacks.preview import (
            register_preview,
        )

        assert callable(register_preview)

    def test_expected_output_wired(self):
        from pycsamt.app.agent_master._ids import IDs
        from pycsamt.app.agent_master.app import create_app

        app = create_app()
        cb_outputs = str(app.callback_map)
        assert IDs.EXPORT_DL in cb_outputs
