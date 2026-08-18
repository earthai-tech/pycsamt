from __future__ import annotations

import numpy as np
import pytest

from pycsamt.airborne.mobilemt import (
    MOBILEMT_ADMITTANCE_TAG,
    MobileMTReferenceStation,
    build_mobilemt_emtf,
)
from pycsamt.metadata import LocationMeta, SiteMeta
from pycsamt.emtf import EMTF, EMTFEDIConversionError


def _doc():
    data = np.arange(12, dtype=float).reshape(2, 3, 2)
    data = data + 1j * (data + 0.25)
    ref = MobileMTReferenceStation(
        station_id="BASE01",
        site=SiteMeta(
            site_id="BASE01",
            location=LocationMeta(latitude=5.2, longitude=-3.1),
        ),
    )
    return build_mobilemt_emtf(
        data,
        frequency=[20.0, 100.0],
        variance=np.ones((2, 3, 2)),
        inverse_signal_covariance=np.tile(
            np.eye(2), (2, 1, 1)
        ).astype(complex),
        residual_covariance=np.tile(
            np.eye(3), (2, 1, 1)
        ).astype(complex),
        reference_station=ref,
    )


def test_mobilemt_custom_tf_roundtrips_through_emtf_xml():
    doc = _doc()
    xml = doc.to_xml()
    restored = EMTF.from_xml(xml)
    tf = restored.get_transfer_function(MOBILEMT_ADMITTANCE_TAG)
    assert tf is not None
    source = doc.get_transfer_function(MOBILEMT_ADMITTANCE_TAG)
    assert tf.shape == (2, 3, 2)
    assert tf.input_channels == ("Ex", "Ey")
    assert tf.output_channels == ("Hx", "Hy", "Hz")
    assert np.allclose(tf.data, source.data)
    assert np.allclose(
        tf.get_estimate("INVSIGCOV").data,
        source.get_estimate("INVSIGCOV").data,
    )
    assert np.allclose(
        tf.get_estimate("RESIDCOV").data,
        source.get_estimate("RESIDCOV").data,
    )
    notes = restored.metadata["notes"]["MobileMT"]
    assert notes["ReferenceStationId"] == "BASE01"
    assert notes["InputChannels"] == "Ex,Ey"
    assert notes["OutputChannels"] == "Hx,Hy,Hz"


def test_mobilemt_is_not_silently_coerced_to_edi():
    doc = _doc()
    with pytest.raises(EMTFEDIConversionError):
        doc.to_edi(on_loss="raise")

    with pytest.raises(EMTFEDIConversionError):
        doc.to_edi(on_loss="warn")
