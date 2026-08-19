# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

"""Phase 9 release-level EMTF contracts and bundled synthetic fixtures."""

from __future__ import annotations

from pathlib import Path

import numpy as np

import pycsamt
from pycsamt.emtf import EMTF


_DATA = Path(pycsamt.__file__).resolve().parent / "data" / "emtf"


def test_release_fixtures_exist():
    expected = {
        "single_station_spectra.edi",
        "single_station_spectra.xml",
        "remote_reference_spectra.edi",
        "remote_reference_spectra.xml",
    }
    assert expected.issubset({path.name for path in _DATA.iterdir()})


def test_single_station_reference_pair_is_scientifically_consistent():
    from_edi = EMTF.from_edi_spectra(_DATA / "single_station_spectra.edi")
    from_xml = EMTF.from_xml(_DATA / "single_station_spectra.xml")

    np.testing.assert_allclose(from_edi.periods, from_xml.periods)
    np.testing.assert_allclose(from_edi.z, from_xml.z)
    np.testing.assert_allclose(from_edi.tipper, from_xml.tipper)
    np.testing.assert_allclose(
        from_edi.impedance.get_estimate("INVSIGCOV").data,
        from_xml.impedance.get_estimate("INVSIGCOV").data,
    )


def test_remote_reference_pair_is_scientifically_consistent():
    from_edi = EMTF.from_edi_spectra(_DATA / "remote_reference_spectra.edi")
    from_xml = EMTF.from_xml(_DATA / "remote_reference_spectra.xml")

    np.testing.assert_allclose(from_edi.periods, from_xml.periods)
    np.testing.assert_allclose(from_edi.z, from_xml.z)
    np.testing.assert_allclose(from_edi.tipper, from_xml.tipper)
    np.testing.assert_allclose(
        from_edi.tipper_tf.get_estimate("RESIDCOV").data,
        from_xml.tipper_tf.get_estimate("RESIDCOV").data,
    )
