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


def test_project_stem_with_short_number_is_not_an_iteration(tmp_path):
    """``BH_31.dat`` is a project name, not "iteration 31".

    ModEM writes its iteration counter with ``%03d`` (``_030``), so a
    stem ending in ``_<1-2 digits>`` must be treated as the observed
    data, not as a numbered response -- otherwise the observed and
    predicted files get swapped.
    """
    _write_dat(tmp_path / "BH_31.dat", decimals=4)  # real errors
    pred = (tmp_path / "BH_31.dat").read_text().replace(
        "1.000E+02\n", "1.000E+13\n"  # response file: masked errors
    )
    (tmp_path / "BH_31_NLCG_030.dat").write_text(pred)

    r = InversionResult(tmp_path, load_models=False, load_covariance=False)
    assert r.data_obs is not None and r.data_pred is not None
    obs_err = r.data_obs.blocks[0]["rows"][0][8]
    pred_err = r.data_pred.blocks[0]["rows"][0][8]
    assert obs_err < 1e9  # observed keeps its real error floor
    assert pred_err > 1e9  # predicted carries the masked sentinel
