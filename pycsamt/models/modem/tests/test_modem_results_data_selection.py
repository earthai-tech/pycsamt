# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""InversionResult observed-data ``.dat`` selection.

A ModEM ``-R`` "read and rewrite" run echoes the input data back out with
``GG_Lat``/``GG_Lon`` truncated to its fixed 3-decimal (``f9.3``) output
format -- ~100 m at mid latitudes -- which makes MapView profile lines
zig-zag. When the real input ``.dat`` is also in the folder,
``InversionResult`` must prefer it.
"""

from __future__ import annotations

import numpy as np

from pycsamt.models.modem.results import InversionResult

_HEADER = (
    "# synthetic\n"
    "# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error\n"
    "> Off_Diagonal_Impedance\n"
    "> exp(+i\\omega t)\n"
    "> [mV/km]/[nT]\n"
    "> 0.00\n"
    "> 0.000 0.000\n"
    "> 2 4\n"
)

# 4 stations, real 4-decimal lat/lon that straddle a 0.001-deg boundary.
_STATIONS = [
    ("S01", 32.1179, 119.1269, -300.0, 100.0),
    ("S02", 32.1188, 119.1265, -200.0, 101.0),
    ("S03", 32.1197, 119.1266, -100.0, 102.0),
    ("S04", 32.1206, 119.1264, 0.0, 103.0),
]
_PERIODS = (1.0e-4, 1.0e-3)


def _write_dat(path, *, decimals: int) -> None:
    lines = [_HEADER]
    for per in _PERIODS:
        for name, lat, lon, x, y in _STATIONS:
            for comp in ("ZXY", "ZYX"):
                lines.append(
                    f"  {per:.5E}  {name}  {lat:.{decimals}f}  {lon:.{decimals}f}  "
                    f"{x:.3f}  {y:.3f}  10.000  {comp}  "
                    f"1.000E+03  5.000E+02  1.000E+02\n"
                )
    path.write_text("".join(lines))


def test_prefers_full_precision_input_over_rewrite_echo(tmp_path):
    _write_dat(tmp_path / "ModEMData_v5.dat", decimals=4)
    _write_dat(tmp_path / "_rw_v5.dat", decimals=3)  # ModEM -R echo

    r = InversionResult(tmp_path, load_models=False, load_covariance=False)
    assert r.data_obs is not None
    lon = np.array([r.data_obs.site_lonlat[n][0] for n, *_ in _STATIONS])
    # the 4-decimal file keeps S02/S04 distinct from S01/S03
    assert not np.all(np.abs(lon - np.round(lon, 3)) < 1e-9)
    assert r.data_obs.site_lonlat["S01"][0] == 119.1269


def test_only_a_rewrite_echo_is_still_used(tmp_path):
    # No real input present -- the echo is all there is, so use it.
    _write_dat(tmp_path / "_rw_v5.dat", decimals=3)
    r = InversionResult(tmp_path, load_models=False, load_covariance=False)
    assert r.data_obs is not None
    assert len(r.data_obs.site_names) == 4
