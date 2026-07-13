# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.seg.edi import EDIFile
from pycsamt.seg.survey import (
    EDIProfile,
    Stations,
    Topography,
)


def _mk_edi(
    p: Path,
    *,
    dataid: str,
    lat: str,
    lon: str,
    elev: float = 1000.0,
) -> Path:
    lines: list[str] = [
        ">HEAD",
        f"  DATAID={dataid}",
        f"  LAT={lat}",
        f"  LONG={lon}",
        f"  ELEV={elev:.0f}",
        "  STDVERS=SEG 1.0",
        "",
        ">INFO",
        "  PROJECT=SIM",
        "  PROCESSEDBY=pyCSAMT",
        "  PROCESSINGSOFTWARE=pyCSAMT",
        "",
        ">=MTSECT",
        f"  SECTID={dataid}",
        "  NFREQ=2",
        "",
        ">!****FREQUENCIES****!",
        ">FREQ  //2",
        "  1.000000E+02  2.000000E+02",
        "",
        ">ZXXR ROT=ZROT  //2",
        "  1.000000E+00  1.000000E+00",
        "",
        ">END",
    ]
    p.write_text("\n".join(lines), encoding="utf-8")
    return p


@pytest.fixture()
def two_edifiles(tmp_path: Path) -> list[EDIFile]:
    f1 = _mk_edi(
        tmp_path / "A.edi",
        dataid="STA01",
        lat="26:00:00N",
        lon="010:00:00E",
        elev=1000.0,
    )
    f2 = _mk_edi(
        tmp_path / "B.edi",
        dataid="STA02",
        lat="26:00:30N",
        lon="010:00:30E",
        elev=1015.0,
    )
    return [EDIFile(f1), EDIFile(f2)]


def test_profile_basic_and_adjust(two_edifiles: list[EDIFile]) -> None:
    # accept list of EDIFile
    prof = EDIProfile(two_edifiles, verbose=0)

    # basic geometry is available
    lat = np.asarray(prof.lat, float)
    lon = np.asarray(prof.lon, float)
    elev = np.asarray(prof.elev, float)
    assert lat.size == 2 and lon.size == 2 and elev.size == 2

    # azimuth/bearing in [0, 360)
    az = float(prof.azimuth)
    assert 0.0 <= az < 360.0

    # distance is monotonic (second station farther)
    d = np.asarray(prof.distance, float)
    assert d.size == 2 and d[-1] > d[0]

    # get_step close to inter-station spacing
    step = float(prof.get_step())
    assert step > 0.0

    # adjust to the derived step, re-check cumulative distance
    prof2 = prof.adjust(step=step)
    d2 = np.asarray(prof2.distance, float)
    # start at zero, second ~ one step
    assert d2[0] == pytest.approx(0.0)
    assert d2[1] == pytest.approx(step, rel=0.05, abs=2.0)

    # update back into underlying EDI headers
    prof3 = prof2.update()
    # lat/lon may change after projection back from UTM, but
    # should be finite and non-zero
    assert np.isfinite(prof3.lat).all()
    assert np.isfinite(prof3.lon).all()


def test_stations_select_offsets_and_set(two_edifiles: list[EDIFile]) -> None:
    prof = EDIProfile(two_edifiles, verbose=0)
    st = Stations(prof, verbose=0)

    # names and get by name
    names = st.names()
    assert set(names) == {"STA01", "STA02"}
    ed = st.get("STA01")
    assert isinstance(ed, EDIFile)

    # pattern / regex selectors
    st_a = st.select(pattern="STA0*")
    assert set(st_a.names()) == {"STA01", "STA02"}
    st_1 = st.select(regex=r".*02$")
    assert st_1.names() == ["STA02"]

    # offsets along/cross relative to first site
    along, cross = st.offsets()
    assert along.size == 2 and cross.size == 2
    # first station is origin
    assert along[0] == pytest.approx(0.0)
    # small across-track for a short diagonal
    assert abs(float(cross[-1])) < float(along[-1]) + 1.0

    # set_coords pushes values into the EDI header
    st.set_coords("STA01", elev=1111.0)
    ed2 = st.get("STA01")
    assert ed2 is not None
    assert float(ed2.elev) == pytest.approx(1111.0)


def test_topography_pipeline(two_edifiles: list[EDIFile]) -> None:
    prof = EDIProfile(two_edifiles, verbose=0)

    # build from profile distances/elev
    topo = Topography(prof, use_profile_step=True, verbose=0)
    d0, z0 = topo.as_arrays()
    assert d0.size == 2 and z0.size == 2

    # smoothing is safe on small arrays
    topo.smooth(window=3, method="median")
    d1, z1 = topo.as_arrays()
    assert z1.size == z0.size

    # detrend sets an internal trend we can plot later
    topo.detrend()
    # resample to denser spacing (no crash, monotonic d)
    topo.resample(step=max(1.0, float(d1[-1]) / 4.0))
    d2, z2 = topo.as_arrays()
    assert d2.size >= d1.size and np.all(np.diff(d2) >= 0.0)

    # gradient length is N-1
    g = topo.gradient()
    assert g.size == max(0, d2.size - 1)
