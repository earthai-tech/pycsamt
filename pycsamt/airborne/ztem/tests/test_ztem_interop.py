from __future__ import annotations

import numpy as np
import pytest

from pycsamt.airborne.ztem import (
    ZTEMReferenceStation,
    build_ztem_emtf,
)
from pycsamt.emtf import EMTF, EMTFEDIConversionError
from pycsamt.metadata import LocationMeta, OrientationMeta, SiteMeta


def _doc():
    data = np.array(
        [
            [1.0 + 0.2j, 2.0 + 0.3j],
            [3.0 + 0.4j, 4.0 + 0.5j],
        ]
    )
    ref = ZTEMReferenceStation(
        station_id="BASE01",
        site=SiteMeta(
            site_id="BASE01",
            location=LocationMeta(latitude=5.2, longitude=-3.1),
        ),
    )
    return build_ztem_emtf(
        data,
        frequency=[30.0, 90.0],
        variance=np.ones((2, 1, 2)),
        inverse_signal_covariance=np.tile(
            np.eye(2), (2, 1, 1)
        ).astype(complex),
        residual_covariance=np.ones((2, 1, 1), dtype=complex),
        reference_station=ref,
        orientation=OrientationMeta(
            mode="orthogonal",
            angle_to_geographic_north=0.0,
        ),
    )


def test_ztem_standard_tipper_roundtrips_through_emtf_xml():
    doc = _doc()
    xml = doc.to_xml()
    assert "<T " in xml
    assert '<DataType name="T"' in xml
    restored = EMTF.from_xml(xml)
    tf = restored.tipper_tf
    source = doc.tipper_tf
    assert tf is not None
    assert tf.shape == (2, 1, 2)
    assert tf.input_channels == ("Hx", "Hy")
    assert tf.output_channels == ("Hz",)
    assert np.allclose(tf.data, source.data)
    assert np.allclose(
        tf.get_estimate("INVSIGCOV").data,
        source.get_estimate("INVSIGCOV").data,
    )
    assert np.allclose(
        tf.get_estimate("RESIDCOV").data,
        source.get_estimate("RESIDCOV").data,
    )
    notes = restored.metadata["notes"]["ZTEM"]
    assert notes["ReferenceStationId"] == "BASE01"
    assert notes["InputChannels"] == "Hx,Hy"
    assert notes["OutputChannels"] == "Hz"
    assert restored.processing.remote_reference.site == "BASE01"


def test_ztem_rotation_reuses_generic_tipper_engine():
    doc = _doc()
    rotated = doc.rotate(30.0)
    assert rotated.tipper_tf.shape == (2, 1, 2)
    assert rotated.orientation.angle_to_geographic_north == 30.0
    assert rotated.tipper_tf.get_estimate("INVSIGCOV").shape == (2, 2, 2)
    assert rotated.tipper_tf.get_estimate("RESIDCOV").shape == (2, 1, 1)

    restored = rotated.rotate(0.0)
    assert np.allclose(restored.tipper_tf.data, doc.tipper_tf.data)


def test_ztem_is_not_silently_coerced_to_standard_edi():
    doc = _doc()
    with pytest.raises(EMTFEDIConversionError):
        doc.to_edi(on_loss="raise")
    with pytest.raises(EMTFEDIConversionError):
        doc.to_edi(on_loss="warn")
