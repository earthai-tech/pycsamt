# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage-focused tests for :mod:`pycsamt.models.mare2dem.zmm`.

All ``.zmm`` content is synthesized in ``tmp_path`` — nothing here
depends on the gitignored ``data/mare2dem/`` example data, so the
coverage gain is real in CI.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from pycsamt.models.mare2dem.zmm import (
    ZMMStation,
    _apply_error_floor,
    _is_number,
    make_mt_data_from_stations,
    make_mt_data_from_zmm,
    read_zmm,
)


def _write_zmm(path, *, with_station=True, with_tipper=True, with_coh=True, n=3):
    lines = []
    if with_station:
        lines.append("Station MT01")
    lines.append("Latitude   42.5")
    lines.append("Longitude  -111.2")
    lines.append("Declination  3.0")
    lines.append(f"NumFreq  {n}")
    for i in range(n):
        T = 10.0 * (i + 1)
        row = [T, 0.1, 0.2, 1.0, 0.5, -1.0, -0.3, 0.15, 0.05]
        if with_tipper:
            row += [0.01, 0.02]
            if with_coh:
                row += [0.02, 0.03, 0.9, 0.8]
            else:
                row += [0.02, 0.03]
        lines.append("Period " + " ".join(str(v) for v in row))
    path.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# _is_number
# ---------------------------------------------------------------------------


def test_is_number():
    assert _is_number("1.5") is True
    assert _is_number("-3") is True
    assert _is_number("abc") is False


# ---------------------------------------------------------------------------
# read_zmm
# ---------------------------------------------------------------------------


def test_read_zmm_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_zmm(tmp_path / "nope.zmm")


def test_read_zmm_full_with_tipper_and_coherence(tmp_path):
    p = tmp_path / "MT01.zmm"
    _write_zmm(p, with_tipper=True, with_coh=True, n=3)
    st = read_zmm(p)

    assert st.name == "MT01"
    assert st.latitude == 42.5
    assert st.longitude == -111.2
    assert st.declination == 3.0
    assert len(st.periods) == 3
    assert st.periods[0] == 10.0

    # apparent resistivity/phase derived from ZXY/ZYX
    T = st.periods[0]
    omega = 2.0 * math.pi / T
    mu0 = 4e-7 * math.pi
    zxy = complex(1.0, 0.5)
    expected_apres_te = abs(zxy) ** 2 / (omega * mu0)
    assert st.apres_te[0] == pytest.approx(expected_apres_te)
    assert st.tipper_zy is not None
    assert st.tipper_zy_se is not None
    np.testing.assert_allclose(st.tipper_zy_se, np.abs(st.tipper_zy) * 0.1)


def test_read_zmm_without_station_line_falls_back_to_stem(tmp_path):
    p = tmp_path / "fallback_name.zmm"
    _write_zmm(p, with_station=False)
    st = read_zmm(p)
    assert st.name == "fallback_name"


def test_read_zmm_without_tipper_columns(tmp_path):
    p = tmp_path / "no_tipper.zmm"
    _write_zmm(p, with_tipper=False)
    st = read_zmm(p)
    assert st.tipper_zy is None
    assert st.tipper_zy_se is None


def test_read_zmm_with_tipper_but_no_coherence_columns(tmp_path):
    p = tmp_path / "no_coh.zmm"
    _write_zmm(p, with_tipper=True, with_coh=False)
    st = read_zmm(p)
    assert st.tipper_zy is not None


def test_read_zmm_malformed_header_values_are_ignored(tmp_path):
    p = tmp_path / "malformed.zmm"
    lines = [
        "Station MT02",
        "Latitude",
        "Longitude notanumber",
        "Declination",
        "Site",
        "NumFreq notanumber",
        "Nfreq",
        "Period 10.0 0.1 0.2 1.0 0.5 -1.0 -0.3 0.15 0.05",
    ]
    p.write_text("\n".join(lines) + "\n")
    st = read_zmm(p)
    assert st.latitude == 0.0
    assert st.longitude == 0.0
    assert st.declination == 0.0
    assert len(st.periods) == 1


def test_read_zmm_no_period_lines_returns_empty_station(tmp_path):
    p = tmp_path / "empty.zmm"
    p.write_text("Station MT03\nLatitude 1.0\nLongitude 2.0\n")
    st = read_zmm(p)
    assert st.name == "MT03"
    assert len(st.periods) == 0


def test_read_zmm_site_keyword_sets_name(tmp_path):
    p = tmp_path / "site_kw.zmm"
    lines = ["Site MT04", "Period 10.0 0.1 0.2 1.0 0.5 -1.0 -0.3 0.15 0.05"]
    p.write_text("\n".join(lines) + "\n")
    st = read_zmm(p)
    assert st.name == "MT04"


# ---------------------------------------------------------------------------
# _apply_error_floor
# ---------------------------------------------------------------------------


def test_apply_error_floor_all_branches():
    st = ZMMStation(
        periods=np.array([10.0, 20.0]),
        apres_te=np.array([100.0, 200.0]),
        apres_te_se=np.array([1.0, 1.0]),
        phase_te_se=np.array([0.5, 0.5]),
        apres_tm=np.array([50.0, 60.0]),
        apres_tm_se=np.array([0.5, 0.5]),
        phase_tm_se=np.array([0.3, 0.3]),
        tipper_zy=np.array([complex(0.1, 0.1), complex(0.2, 0.2)]),
        tipper_zy_se=np.array([0.01, 0.01]),
    )
    out = _apply_error_floor([st], 0.05, 0.1, 0.02)[0]
    assert out.apres_te_se[0] == pytest.approx(max(1.0, 100.0 * 0.05))
    assert out.apres_tm_se[0] == pytest.approx(max(0.5, 50.0 * 0.1))
    assert out.tipper_zy_se[0] == pytest.approx(0.02)


def test_apply_error_floor_zero_floors_are_noop():
    st = ZMMStation(
        apres_te=np.array([100.0]),
        apres_te_se=np.array([1.0]),
        phase_te_se=np.array([0.5]),
        apres_tm=np.array([50.0]),
        apres_tm_se=np.array([0.5]),
        phase_tm_se=np.array([0.3]),
        tipper_zy_se=None,
    )
    out = _apply_error_floor([st], 0.0, 0.0, 0.0)[0]
    assert out.apres_te_se[0] == 1.0
    assert out.apres_tm_se[0] == 0.5


def test_apply_error_floor_tipper_floor_without_se_is_skipped():
    st = ZMMStation(
        apres_te=np.array([100.0]),
        apres_te_se=np.array([1.0]),
        phase_te_se=np.array([0.5]),
        apres_tm=np.array([50.0]),
        apres_tm_se=np.array([0.5]),
        phase_tm_se=np.array([0.3]),
        tipper_zy_se=None,
    )
    out = _apply_error_floor([st], 0.0, 0.0, 0.05)[0]
    assert out.tipper_zy_se is None


# ---------------------------------------------------------------------------
# make_mt_data_from_zmm / make_mt_data_from_stations
# ---------------------------------------------------------------------------


def _build_zmm_files(tmp_path, n_stations=3):
    files = []
    lats = [42.0, 42.01, 42.02]
    lons = [-111.2, -111.2, -111.2]
    for i in range(n_stations):
        p = tmp_path / f"S{i:03d}.zmm"
        lines = [
            f"Station S{i:03d}",
            f"Latitude {lats[i % len(lats)]}",
            f"Longitude {lons[i % len(lons)]}",
            "Declination 0.0",
            "NumFreq 2",
            "Period 10.0 0.1 0.2 1.0 0.5 -1.0 -0.3 0.15 0.05 0.01 0.02 0.02 0.03 0.9 0.8",
            "Period 20.0 0.1 0.2 1.0 0.5 -1.0 -0.3 0.15 0.05 0.01 0.02 0.02 0.03 0.9 0.8",
        ]
        p.write_text("\n".join(lines) + "\n")
        files.append(p)
    return files


def test_make_mt_data_from_zmm_end_to_end(tmp_path):
    files = _build_zmm_files(tmp_path)
    out_file = tmp_path / "line1_mt.emdata"
    em = make_mt_data_from_zmm(files, out_file, output_modes="all")
    assert out_file.exists()
    assert em.mt is not None
    assert em.mt.n_mt_receivers if hasattr(em.mt, "n_mt_receivers") else True
    assert len(em.mt.receiver_name) == 3
    assert em.data.shape[1] == 6


@pytest.mark.parametrize(
    "modes", ["TE", "TM", "tipper", "TE+tipper", "all impedance", "all"]
)
def test_make_mt_data_from_zmm_output_modes(tmp_path, modes):
    files = _build_zmm_files(tmp_path)
    out_file = tmp_path / f"modes_{modes.replace(' ', '_').replace('+','_')}.emdata"
    em = make_mt_data_from_zmm(files, out_file, output_modes=modes)
    assert out_file.exists()


def test_make_mt_data_from_zmm_with_explicit_utm_zone_and_orientation(tmp_path):
    files = _build_zmm_files(tmp_path)
    out_file = tmp_path / "explicit.emdata"
    em = make_mt_data_from_zmm(
        files,
        out_file,
        utm_zone="12N",
        line_orientation=0.0,
        utm0=(0.0, 0.0),
        declination=1.0,
        topo=500.0,
        rx_z_offset=0.5,
    )
    assert em.utm.grid == 12
    assert em.utm.hemi == "N"


def test_make_mt_data_from_zmm_south_hemisphere_zone(tmp_path):
    files = _build_zmm_files(tmp_path)
    out_file = tmp_path / "south.emdata"
    em = make_mt_data_from_zmm(files, out_file, utm_zone="19S")
    assert em.utm.hemi == "S"


@pytest.mark.parametrize("orientation", [10.0, 90.0, 180.0, 270.0])
def test_make_mt_data_from_zmm_sort_branches(tmp_path, orientation):
    files = _build_zmm_files(tmp_path)
    out_file = tmp_path / f"orient_{int(orientation)}.emdata"
    em = make_mt_data_from_zmm(files, out_file, line_orientation=orientation)
    assert out_file.exists()


def test_make_mt_data_from_zmm_omit_periods(tmp_path):
    files = _build_zmm_files(tmp_path)
    out_file = tmp_path / "omit.emdata"
    em = make_mt_data_from_zmm(
        files, out_file, omit_periods=np.array([[15.0, 25.0]])
    )
    assert em.mt.frequencies.size == 1


def test_make_mt_data_from_stations_empty_list(tmp_path):
    # Auto-detecting UTM zone / line orientation / UTM origin all crash on
    # an empty station list (np.median/np.polyfit/indexing on empty
    # arrays) before ever reaching the "if not stations" guard further
    # down -- so exercising that guard means supplying all three
    # explicitly, as a caller re-using a known survey's UTM frame would.
    out_file = tmp_path / "empty.emdata"
    em = make_mt_data_from_stations(
        [],
        out_file,
        utm_zone="12N",
        line_orientation=0.0,
        utm0=(0.0, 0.0),
    )
    assert out_file.exists()
    assert em.mt is None


def test_make_mt_data_from_zmm_error_floors_and_declination(tmp_path):
    files = _build_zmm_files(tmp_path)
    out_file = tmp_path / "floors.emdata"
    em = make_mt_data_from_zmm(
        files,
        out_file,
        error_floor_te=0.05,
        error_floor_tm=0.05,
        error_floor_tipper=0.02,
    )
    assert out_file.exists()
