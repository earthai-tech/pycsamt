# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations

from pathlib import Path

import pytest

from pycsamt.seg.meas import (
    DefineMeas,
    DefineMeasMixin,
    EMeasMixin,
    Emeasurement,
    HMeasMixin,
    Hmeasurement,
    MeasMixin,
)


_EDI_SAMPLES = [
    Path("MT") / "kap03lmt_edis" / "kap103.edi",
    Path("MT") / "SPECTRA" / "spectra01.edi",
    Path("CSAMT") / "csa000.edi",
]


@pytest.mark.parametrize(
    "which",
    _EDI_SAMPLES,
)
def test_define_meas_parses_real_files(edi_path: Path, which: str) -> None:
    p = edi_path / which
    if not p.exists():
        pytest.skip(f"Missing EDI: {p}")
    dm = DefineMeas.from_file(p)

    # at least one measurement defined
    total = len(dm.hmeas) + len(dm.emeas)
    assert total > 0

    # optional reference origin present or None
    if dm.reflat is not None:
        assert isinstance(dm.reflat, float)
    if dm.reflong is not None:
        assert isinstance(dm.reflong, float)


@pytest.mark.parametrize(
    "which",
    _EDI_SAMPLES,
)
def test_mixins_on_real_files(edi_path: Path, which: str) -> None:
    p = edi_path / which
    if not p.exists():
        pytest.skip(f"Missing EDI: {p}")

    class HostAll(MeasMixin):
        pass

    class HostD(DefineMeasMixin):
        pass

    class HostE(EMeasMixin):
        pass

    class HostH(HMeasMixin):
        pass

    dm = HostAll.from_file(p)
    assert isinstance(dm, DefineMeas)

    dm2 = HostD.from_file(p)
    assert isinstance(dm2, DefineMeas)

    e_list = HostE.from_file(p)
    h_list = HostH.from_file(p)
    assert isinstance(e_list, list)
    assert isinstance(h_list, list)

    # if lists are non-empty, check element types
    if e_list:
        assert isinstance(e_list[0], Emeasurement)
        assert e_list[0].to_line().startswith(">EMEAS ")
        d = e_list[0].to_dict()
        assert "id" in d and "chtype" in d

    if h_list:
        assert isinstance(h_list[0], Hmeasurement)
        assert h_list[0].to_line().startswith(">HMEAS ")
        d = h_list[0].to_dict()
        assert "id" in d and "chtype" in d
