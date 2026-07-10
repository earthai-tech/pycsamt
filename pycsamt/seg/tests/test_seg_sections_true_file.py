# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

from pathlib import Path

import pytest

from pycsamt.seg.base import EDIComponentBase
from pycsamt.seg.sections import (
    SECT_REGISTRY,
    iter_sections,
)


def _collect(path: Path) -> list[tuple[str, EDIComponentBase, int]]:
    return list(iter_sections(str(path)))


@pytest.mark.parametrize(
    "which, expect_any",
    [
        ("15125A_spe.edi", {">=SPECTRASECT"}),
        ("15125A_imp.edi", {">=MTSECT"}),
        # CSAMT may vary; accept any recognized section
        (
            "000CSA_csamt.edi",
            {
                ">=SPECTRASECT",
                ">=MTSECT",
                ">=TSERIESSECT",
                ">=OTHERSECT",
            },
        ),
    ],
)
def test_iter_sections_finds_expected_minimum(
    edi_path: Path, which: str, expect_any: set[str]
):
    p = edi_path / which
    if not p.exists():
        pytest.skip(f"Missing EDI: {p}")

    items = _collect(p)
    assert len(items) >= 1

    tags = {t for (t, _, _) in items}
    # For CSAMT we only require overlap
    assert tags & expect_any, (
        f"No expected tags in {tags} for {which}"
    )


@pytest.mark.parametrize(
    "which",
    ["15125A_spe.edi", "15125A_imp.edi", "000CSA_csamt.edi"],
)
def test_registry_dispatch_and_body_start_valid(
    edi_path: Path, which: str
):
    p = edi_path / which
    if not p.exists():
        pytest.skip(f"Missing EDI: {p}")

    txt = p.read_text(encoding="utf-8-sig", errors="replace")
    lines = txt.splitlines()
    n = len(lines)

    for tag, hdr, start in _collect(p):
        # tag present in registry
        assert tag in SECT_REGISTRY
        # header is instance of expected class
        expected = SECT_REGISTRY[tag]
        assert isinstance(hdr, expected)

        # body_start is in bounds and points to a tag
        assert 0 < start <= n
        # allow pointing to next section if no data blocks
        if start < n:
            assert lines[start].lstrip().startswith(">")

        # round-trip header write should start with a top tag
        write = getattr(hdr, "write", None)
        if callable(write):
            out = write()
            assert isinstance(out, list)
            assert len(out) >= 1
            assert out[0].lstrip().startswith(">=")


@pytest.mark.parametrize(
    "which, include, expect_tag",
    [
        ("15125A_spe.edi", [">=SPECTRASECT"], ">=SPECTRASECT"),
        ("15125A_imp.edi", [">=MTSECT"], ">=MTSECT"),
    ],
)
def test_iter_sections_include_filter(
    edi_path: Path, which: str, include: list[str], expect_tag: str
):
    p = edi_path / which
    if not p.exists():
        pytest.skip(f"Missing EDI: {p}")

    items = list(iter_sections(str(p), include=include))
    # at least one matching section should be yielded
    assert any(t == expect_tag for (t, _, _) in items)


@pytest.mark.parametrize(
    "which",
    ["15125A_spe.edi", "15125A_imp.edi", "000CSA_csamt.edi"],
)
def test_iter_sections_tags_in_file_order(
    edi_path: Path, which: str
):
    p = edi_path / which
    if not p.exists():
        pytest.skip(f"Missing EDI: {p}")

    items = _collect(p)
    # body_start indexes should be non-decreasing as we scan
    starts = [s for (_, _, s) in items]
    assert starts == sorted(starts)
