# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Parity of the pyCSAMT adapters with the pcsf reference adapters.

pyCSAMT converts its own result objects (``InversionResult``), while
pcsf reads the native files directly. Both must produce the same
volume and geometry on pcsf's reference runs; the ModEM case guards
the north-to-south x order of ModEM files (``WS.inc`` writes
``i=Nx,1,-1``), which an earlier reader got wrong.
"""

from __future__ import annotations

import numpy as np
import pcsf
import pytest
from pcsf import datasets
from pcsf.adapters import modem_to_pcsf, occam2d_to_pcsf

import pycsamt.format as fmt
from pycsamt.format.adapters.modem3d import modem3d_to_pcsf
from pycsamt.format.adapters.occam2d import (
    occam2d_to_pcsf as occam2d_result_to_pcsf,
)


def test_format_objects_are_pcsf_objects():
    assert fmt.PCSFModel is pcsf.PCSFModel
    assert fmt.read_pcsf is pcsf.read_pcsf
    assert fmt.write_pcsm is pcsf.write_pcsm


def test_modem_adapter_matches_pcsf():
    from pycsamt.models.modem.results import InversionResult

    run = datasets.path("modem_tiny")
    ours = modem3d_to_pcsf(
        InversionResult(run), air_threshold_ohm_m=None
    )
    ref = modem_to_pcsf(
        run / "tiny_NLCG_002.rho",
        air_threshold_ohm_m=None,
        station_z_convention="elevation",
    )
    np.testing.assert_array_equal(ours.resistivity, ref.resistivity)
    for axis in ("x_nodes", "y_nodes", "z_nodes", "origin"):
        np.testing.assert_array_equal(
            getattr(ours.geometry, axis), getattr(ref.geometry, axis)
        )
    # The reference run varies along x, so a mirrored reading differs.
    assert not np.array_equal(
        ours.resistivity, ref.resistivity[:, :, ::-1]
    )


def test_occam2d_adapter_matches_pcsf():
    from pycsamt.models.occam2d.results import InversionResult

    run = datasets.path("occam2d_tiny")
    ours = occam2d_result_to_pcsf(InversionResult(run))
    ref = occam2d_to_pcsf(run)
    np.testing.assert_array_equal(
        ours.resistivity, ref.resistivity
    )
    np.testing.assert_array_equal(ours.geometry.x, ref.geometry.x)
    np.testing.assert_array_equal(ours.geometry.z, ref.geometry.z)


@pytest.mark.parametrize("name", datasets.NAMES)
def test_reference_runs_round_trip_through_pycsamt(name, tmp_path):
    from pcsf.adapters import mare2dem_to_pcsf

    run = datasets.path(name)
    if name == "mare2dem_tiny":
        pytest.importorskip("triangle")
        model = mare2dem_to_pcsf(run / "tiny.1.resistivity")
    elif name == "modem_tiny":
        model = modem_to_pcsf(
            run / "tiny_NLCG_002.rho",
            station_z_convention="elevation",
        )
    else:
        model = occam2d_to_pcsf(run)
    back = fmt.read_pcsf(fmt.write_pcsf(model, tmp_path / "m.pcsf"))
    np.testing.assert_array_equal(back.resistivity, model.resistivity)
