# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

import numpy as np

from pycsamt.interp import HydroInterpreter, ResistivityModel
from pycsamt.inversion.results import InversionResult


def _model():
    rho = np.log10(np.array([
        [100.0, 12.0, 1500.0],
        [80.0, 25.0, 1200.0],
        [250.0, 70.0, 900.0],
        [1500.0, 1800.0, 2500.0],
    ]))
    return ResistivityModel.from_array(
        rho,
        x_centers=np.array([0.0, 500.0, 1000.0]),
        z_centers=np.array([25.0, 75.0, 150.0, 300.0]),
        station_x=np.array([0.0, 500.0, 1000.0]),
        station_names=["S1", "S2", "S3"],
        method="test",
    )


def test_hydro_interpreter_finds_aquifer_zones(tmp_path):
    hydro = HydroInterpreter(
        water_table_depth=20.0,
        aquifer_range=(30.0, 300.0),
        clay_max=20.0,
        min_zone_thickness=10.0,
    ).fit(_model())

    zones = hydro.aquifer_zones()
    assert zones
    assert hydro.unit_map.shape == (4, 3)
    assert hydro.unit_map[0, 1] == "clay"
    assert hydro.unit_map[0, 2] == "resistive basement"
    assert hydro.station_summary()[0]["station"] == "S1"

    cells = hydro.to_csv(tmp_path / "hydro_cells.csv")
    zone_file = hydro.zones_to_csv(tmp_path / "hydro_zones.csv")
    assert cells.exists()
    assert zone_file.exists()


def test_hydro_interpreter_accepts_inversion_result():
    inv_result = InversionResult(
        method="mt",
        dimension="2d",
        backend="builtin",
        model={
            "rho_2d": _model().rho_2d,
            "x_centers": _model().x_centers,
            "z_centers": _model().z_centers,
            "station_x": _model().station_x,
            "station_names": _model().station_names,
        },
    )
    hydro = HydroInterpreter(aquifer_range=(30.0, 300.0)).fit(inv_result)
    assert hydro.resistivity_model.n_x == 3
    assert len(hydro.logs) == 3
