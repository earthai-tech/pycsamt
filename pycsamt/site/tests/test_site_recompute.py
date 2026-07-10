from __future__ import annotations

import csv
from pathlib import Path

from pycsamt.seg.edi import EDIFile
from pycsamt.site import (
    EDIRecomputer,
    recompute_edi,
    recompute_edis,
)
from pycsamt.site.recompute import EDIRecomputeResult


def _copy_edi(src: Path, dst: Path) -> Path:
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(src.read_text(encoding="utf-8"), encoding="utf-8")
    return dst


def test_recompute_edi_returns_copy(simulated_edi: Path) -> None:
    edi = EDIFile(simulated_edi)
    out = recompute_edi(edi, rotate_angle=0.0, copy=True)

    assert out is not edi
    assert hasattr(out, "Z")
    assert out.Z.n_freq == edi.Z.n_freq


def test_recompute_line_folder_default_output(
    tmp_path: Path,
    simulated_edi: Path,
) -> None:
    root = tmp_path / "WILLY_DATA"
    line = root / "L18PLT"
    _copy_edi(simulated_edi, line / "L18A.edi")
    _copy_edi(simulated_edi, line / "L18B.edi")

    result = recompute_edis(
        line,
        template="{source_stem}.edi",
        overwrite=True,
        progress=False,
        verbose=0,
    )

    assert isinstance(result, EDIRecomputeResult)
    assert result.output_root == root / "recomputed_edis"
    assert len(result.records) == 2
    assert not result.failed
    assert (root / "recomputed_edis" / "L18PLT" / "L18A.edi").exists()
    assert (root / "recomputed_edis" / "L18PLT" / "L18B.edi").exists()
    assert (root / "recomputed_edis" / "manifest.csv").exists()


def test_recompute_survey_folder_preserves_lines(
    tmp_path: Path,
    simulated_edi: Path,
) -> None:
    root = tmp_path / "WILLY_DATA"
    for line in ("L18PLT", "L32PLT", "L22PLT"):
        _copy_edi(simulated_edi, root / line / f"{line}_A.edi")
        _copy_edi(simulated_edi, root / line / f"{line}_B.edi")

    result = EDIRecomputer(
        template="{source_stem}.edi",
        overwrite=True,
        verbose=0,
    ).run(root)

    assert len(result.records) == 6
    assert not result.failed
    for line in ("L18PLT", "L32PLT", "L22PLT"):
        outdir = root / "recomputed_edis" / line
        assert (outdir / f"{line}_A.edi").exists()
        assert (outdir / f"{line}_B.edi").exists()

    with (root / "recomputed_edis" / "manifest.csv").open(
        "r",
        encoding="utf-8",
        newline="",
    ) as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 6
    assert {row["line"] for row in rows} == {
        "L18PLT",
        "L32PLT",
        "L22PLT",
    }


def test_recompute_survey_folder_can_flatten_output(
    tmp_path: Path,
    simulated_edi: Path,
) -> None:
    root = tmp_path / "WILLY_DATA"
    _copy_edi(simulated_edi, root / "L18PLT" / "A.edi")
    _copy_edi(simulated_edi, root / "L32PLT" / "B.edi")

    result = recompute_edis(
        root,
        template="{line}_{source_stem}.edi",
        preserve_line_dirs=False,
        overwrite=True,
        verbose=0,
    )

    assert len(result.records) == 2
    assert (root / "recomputed_edis" / "L18PLT_A.edi").exists()
    assert (root / "recomputed_edis" / "L32PLT_B.edi").exists()
    assert not (root / "recomputed_edis" / "L18PLT").exists()


def test_recompute_collection_can_skip_writing(
    tmp_path: Path,
    simulated_edi: Path,
) -> None:
    p1 = _copy_edi(simulated_edi, tmp_path / "A.edi")
    p2 = _copy_edi(simulated_edi, tmp_path / "B.edi")
    edis = [EDIFile(p1), EDIFile(p2)]

    result = recompute_edis(edis, write=False, verbose=0)

    assert len(result.edis) == 2
    assert result.output_root is None
    assert result.paths == []
