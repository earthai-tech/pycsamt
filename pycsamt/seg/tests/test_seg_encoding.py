"""Regression tests: EDI files with non-UTF-8 header text.

Chinese acquisition software (e.g. prMT) writes free-text header
fields such as ``ACQBY`` in GBK. The readers used to open files as
strict UTF-8 and failed with ``UnicodeDecodeError``; ``EDIBatch``
then silently skipped every file of a survey line.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest

from pycsamt.seg.edi import EDIFile
from pycsamt.seg.validation import read_edi_lines

ORG = "中南大学"


def _with_acqby(src: Path, dst: Path, codec: str) -> Path:
    """Copy *src* with ``ACQBY`` set to :data:`ORG` in *codec*."""
    text = src.read_text(encoding="utf-8", errors="replace")
    text, n = re.subn(
        r'(?m)^(\s*ACQBY\s*=).*$', rf'\1"{ORG}"', text, count=1
    )
    if n == 0:
        text = text.replace(">HEAD", f'>HEAD\n    ACQBY="{ORG}"', 1)
    dst.write_bytes(text.encode(codec))
    return dst


@pytest.mark.parametrize("codec", ["gbk", "utf-8", "utf-8-sig"])
def test_read_edi_lines_decodes(tmp_path: Path, codec: str) -> None:
    p = tmp_path / "x.edi"
    p.write_bytes(f'>HEAD\n ACQBY="{ORG}"\n>END\n'.encode(codec))
    lines = read_edi_lines(p)
    assert lines[0] == ">HEAD\n"  # BOM stripped for utf-8-sig
    assert lines[1].strip() == f'ACQBY="{ORG}"'


def test_read_edi_lines_latin1_never_fails(tmp_path: Path) -> None:
    p = tmp_path / "bad.edi"
    p.write_bytes(b">HEAD\n X=\xff\xfe\x81\n")
    assert len(read_edi_lines(p)) == 2


def test_read_edi_lines_custom_codecs_raise(tmp_path: Path) -> None:
    p = tmp_path / "bad.edi"
    p.write_bytes(b"\xff\xfe\x81")
    with pytest.raises(UnicodeDecodeError):
        read_edi_lines(p, encodings=("utf-8",))


def test_edifile_reads_gbk_header(
    edi_imp_file: Path, tmp_path: Path
) -> None:
    ref = EDIFile(edi_imp_file, verbose=0)
    gbk = _with_acqby(edi_imp_file, tmp_path / "gbk.edi", "gbk")
    with pytest.raises(UnicodeDecodeError):
        gbk.read_text(encoding="utf-8")
    ed = EDIFile(gbk, verbose=0)
    np.testing.assert_allclose(ed.Z.freq, ref.Z.freq)
    np.testing.assert_allclose(ed.Z.z, ref.Z.z)
    assert ed.get_section("head").acqby == ORG

    # Written back as UTF-8 without mojibake.
    ed.write(new_edifn="out.edi", savepath=str(tmp_path))
    out = (tmp_path / "out.edi").read_text(encoding="utf-8")
    assert ORG in out
    assert "�" not in out


def test_edibatch_loads_gbk_files(
    edi_imp_file: Path, tmp_path: Path
) -> None:
    from pycsamt.stratagem import EDIBatch

    for i in range(3):
        _with_acqby(edi_imp_file, tmp_path / f"S-{i * 30}.edi", "gbk")
    batch = EDIBatch(tmp_path).fit()
    assert len(batch) == 3
