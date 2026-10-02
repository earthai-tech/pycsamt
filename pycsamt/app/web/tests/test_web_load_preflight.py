"""Tests for the web load-modal preflight helpers."""

from __future__ import annotations


def test_infer_line_counts_from_folder_paths():
    from pycsamt.app.web.callbacks.data import (
        _infer_line_counts,
    )

    names = [
        "WILLY_DATA/L18PLT/001.edi",
        "WILLY_DATA/L18PLT/002.edi",
        "WILLY_DATA/L22PLT/001.edi",
    ]

    assert _infer_line_counts(names) == {"L18PLT": 2, "L22PLT": 1}


def test_infer_line_counts_from_flat_files():
    from pycsamt.app.web.callbacks.data import (
        _infer_line_counts,
    )

    assert _infer_line_counts(["a.edi", "b.avg"]) == {"flat files": 2}


def test_ext_counts_groups_known_and_unknown_formats():
    from pycsamt.app.web.callbacks.data import _ext_counts

    counts = _ext_counts(["a.EDI", "b.avg", "c.J", "notes.txt"])

    assert counts == {"edi": 1, "avg": 1, "j": 1, "other": 1}


def test_upload_entries_keep_content_and_sanitize_editable_names():
    from pycsamt.app.web.callbacks.data import _upload_entries

    entries = _upload_entries(
        ["data:a", "data:b"],
        ["Line 1/a bad?.EDI", "b.avg"],
    )

    assert entries[0]["content"] == "data:a"
    assert entries[0]["original"] == "Line 1/a bad?.EDI"
    assert entries[0]["filename"] == "Line 1/a bad_.EDI"
    assert entries[1]["filename"] == "b.avg"
