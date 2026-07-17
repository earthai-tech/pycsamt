# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.mapview.callbacks.inversion_import."""

from __future__ import annotations

import base64
import os

import pytest

pytest.importorskip("dash", reason="dash required")
pytest.importorskip("dash_bootstrap_components", reason="dbc required")


def _b64_payload(text: str) -> str:
    encoded = base64.b64encode(text.encode()).decode()
    return f"data:text/plain;base64,{encoded}"


class TestDecodeToTempdir:
    def test_writes_files_with_decoded_content(self):
        from pycsamt.app.mapview.callbacks.inversion_import import (
            _decode_to_tempdir,
        )

        tmpdir = _decode_to_tempdir(
            ["a.dat", "b.dat"],
            [_b64_payload("hello"), _b64_payload("world")],
        )
        try:
            with open(os.path.join(tmpdir, "a.dat")) as fh:
                assert fh.read() == "hello"
            with open(os.path.join(tmpdir, "b.dat")) as fh:
                assert fh.read() == "world"
        finally:
            import shutil

            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_sanitizes_nested_path_to_basename(self):
        from pycsamt.app.mapview.callbacks.inversion_import import (
            _decode_to_tempdir,
        )

        tmpdir = _decode_to_tempdir(
            ["sub/dir/c.dat"], [_b64_payload("nested")]
        )
        try:
            assert os.path.isfile(os.path.join(tmpdir, "c.dat"))
            assert not os.path.isdir(os.path.join(tmpdir, "sub"))
        finally:
            import shutil

            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_skips_unreadable_entries(self):
        from pycsamt.app.mapview.callbacks.inversion_import import (
            _decode_to_tempdir,
        )

        tmpdir = _decode_to_tempdir(
            ["good.dat", "bad.dat"],
            [_b64_payload("ok"), "not-a-valid-data-uri"],
        )
        try:
            assert os.path.isfile(os.path.join(tmpdir, "good.dat"))
            assert not os.path.exists(os.path.join(tmpdir, "bad.dat"))
        finally:
            import shutil

            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_empty_inputs_creates_empty_dir(self):
        from pycsamt.app.mapview.callbacks.inversion_import import (
            _decode_to_tempdir,
        )

        tmpdir = _decode_to_tempdir([], [])
        try:
            assert os.path.isdir(tmpdir)
            assert os.listdir(tmpdir) == []
        finally:
            import shutil

            shutil.rmtree(tmpdir, ignore_errors=True)


class TestRegisterInversionImport:
    def test_register_is_callable(self):
        from pycsamt.app.mapview.callbacks.inversion_import import (
            register_inversion_import,
        )

        assert callable(register_inversion_import)

    def test_expected_outputs_wired(self):
        from pycsamt.app.mapview._ids import IDs
        from pycsamt.app.mapview.app import create_app

        app = create_app()
        cb_outputs = str(app.callback_map)
        assert IDs.BTN_INV_CONFIRM in cb_outputs
        assert IDs.INV_STATUS in cb_outputs
