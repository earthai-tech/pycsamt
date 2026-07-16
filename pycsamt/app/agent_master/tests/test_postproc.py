# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.agent_master.callbacks.postproc pure helpers."""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("dash", reason="dash required")


def _dup_edi(tmp_path: Path, src: Path, stem: str) -> Path:
    dst = tmp_path / f"{stem}.edi"
    dst.write_text(src.read_text(encoding="utf-8"), encoding="utf-8")
    return dst


def _sites(tmp_path: Path, simulated_edi: Path):
    from pycsamt.seg.edi import EDIFile
    from pycsamt.site.base import Site, Sites

    p1 = _dup_edi(tmp_path, simulated_edi, "T01")
    p2 = _dup_edi(tmp_path, simulated_edi, "T02")
    s1 = Site(EDIFile(p1))
    s2 = Site(EDIFile(p2))
    return Sites([s1.edi, s2.edi])


class TestExport:
    def test_writes_edi_files(self, tmp_path, simulated_edi):
        from pycsamt.app.agent_master.callbacks.postproc import _export

        sites = _sites(tmp_path, simulated_edi)
        dest = tmp_path / "out"
        paths = _export(sites, dest)
        assert len(paths) == 2
        assert all(Path(p).suffix == ".edi" for p in paths)


class TestCountValidImpedance:
    def test_counts_stations_with_z(self, tmp_path, simulated_edi):
        from pycsamt.app.agent_master.callbacks.postproc import (
            _count_valid_impedance,
        )

        sites = _sites(tmp_path, simulated_edi)
        assert _count_valid_impedance(sites) == 2


class TestValidateExport:
    def test_valid_export_reports_all_stations(self, tmp_path, simulated_edi):
        from pycsamt.app.agent_master.callbacks.postproc import (
            _export,
            _validate_export,
        )

        sites = _sites(tmp_path, simulated_edi)
        dest = tmp_path / "out2"
        _export(sites, dest)
        n_valid, n_total, detail = _validate_export(dest)
        assert n_valid == 2
        assert n_total == 2
        assert detail == ""

    def test_empty_folder_reports_zero_valid(self, tmp_path):
        from pycsamt.app.agent_master.callbacks.postproc import (
            _validate_export,
        )

        dest = tmp_path / "empty_out"
        dest.mkdir()
        n_valid, n_total, detail = _validate_export(dest)
        assert n_valid == 0
        assert n_total == 0
