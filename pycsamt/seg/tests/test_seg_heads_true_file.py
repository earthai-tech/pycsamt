
from __future__ import annotations

from pathlib import Path

import pytest

from pycsamt.seg.heads import (
    Head,
    HeadMixin,
    Heads,
    Info,
    InfoMixin,
)


def _in(val, lo, hi):
    return (val is not None) and (lo < val < hi)


@pytest.mark.parametrize("which", [
    "15125A_imp.edi",
    "15125A_spe.edi",
    "000CSA_csamt.edi",
])
def test_head_parses_real_files(edi_path: Path, which: str):
    p = edi_path / which
    if not p.exists():
        pytest.skip(f"Missing EDI: {p}")
    h = Head.from_file(p)
    assert isinstance(h, Head)
    assert h.dataid is None or str(h.dataid).strip() != ""
    if h.lat is not None:
        assert _in(h.lat, -90.0, 90.0)
    if h.long is not None:
        assert _in(h.long, -180.0, 180.0)
    out = h.write()
    assert out and out[0].strip() == ">HEAD"


@pytest.mark.parametrize("which", [
    "15125A_imp.edi",
    "15125A_spe.edi",
    "000CSA_csamt.edi",
])
def test_info_parses_real_files(edi_path: Path, which: str):
    p = edi_path / which
    if not p.exists():
        pytest.skip(f"Missing EDI: {p}")
    info = Info.from_file(p)
    assert isinstance(info, Info)
    # Basic smoke on nested containers
    assert hasattr(info, "Source")
    assert hasattr(info, "Processing")
    out = info.write()
    assert out and out[0].strip() == ">INFO"


@pytest.mark.parametrize("which", [
    "15125A_imp.edi",
    "15125A_spe.edi",
    "000CSA_csamt.edi",
])
def test_heads_aggregator_order_and_write(
    edi_path: Path, which: str, tmp_path: Path
):
    p = edi_path / which
    if not p.exists():
        pytest.skip(f"Missing EDI: {p}")
    hs = Heads.from_file(p)
    assert isinstance(hs.head, Head)
    assert isinstance(hs.info, Info)
    lines = hs.write()
    text = "".join(lines)
    assert ">HEAD" in text and ">INFO" in text
    assert text.index(">HEAD") < text.index(">INFO")
    outp = tmp_path / f"roundtrip_{p.name}"
    outp.write_text(text, encoding="utf-8")
    assert outp.exists() and outp.stat().st_size > 0


@pytest.mark.parametrize("which", [
    "15125A_imp.edi",
    "15125A_spe.edi",
    "000CSA_csamt.edi",
])
def test_head_mixin_from_file(edi_path: Path, which: str):
    p = edi_path / which
    if not p.exists():
        pytest.skip(f"Missing EDI: {p}")

    class Host(HeadMixin):
        pass

    h = Host.from_file(p)
    assert isinstance(h, Head)
    assert h.dataid is None or str(h.dataid).strip() != ""


@pytest.mark.parametrize("which", [
    "15125A_imp.edi",
    "15125A_spe.edi",
    "000CSA_csamt.edi",
])
def test_info_mixin_from_file(edi_path: Path, which: str):
    p = edi_path / which
    if not p.exists():
        pytest.skip(f"Missing EDI: {p}")

    class Host(InfoMixin):
        pass

    info = Host.from_file(p)
    assert isinstance(info, Info)
    # Project may be absent in some files; just assert attribute exists
    assert hasattr(info.Source, "project")
