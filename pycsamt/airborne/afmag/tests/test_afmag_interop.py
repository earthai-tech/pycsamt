from __future__ import annotations

import numpy as np
import pytest

from pycsamt.airborne.afmag import (
    AFMAG_AP_TAG,
    AFMAG_TENSOR_TAG,
    AFMAG_TILT_TAG,
    AFMAGReferenceStation,
    build_airmt_emtf,
    build_original_afmag_emtf,
)
from pycsamt.emtf import EMTF, EMTFEDIConversionError
from pycsamt.metadata import LocationMeta, OrientationMeta, SiteMeta


def _tensor():
    return np.array(
        [
            [
                [1.0 + 0.2j, 0.1 + 0.05j],
                [0.2 - 0.1j, 1.1 + 0.1j],
                [0.3 + 0.02j, -0.2 + 0.03j],
            ],
            [
                [1.2 + 0.3j, 0.15 + 0.08j],
                [0.25 - 0.12j, 1.3 + 0.2j],
                [0.35 + 0.04j, -0.25 + 0.05j],
            ],
        ]
    )


def test_airmt_tensor_and_ap_roundtrip_through_emtf_xml():
    ref = AFMAGReferenceStation(
        station_id="BASE01",
        site=SiteMeta(
            site_id="BASE01",
            location=LocationMeta(latitude=5.2, longitude=-3.1),
        ),
    )
    doc = build_airmt_emtf(
        _tensor(),
        frequency=[24.0, 75.0],
        variance=np.ones((2, 3, 2)),
        inverse_signal_covariance=np.tile(
            np.eye(2), (2, 1, 1)
        ).astype(complex),
        residual_covariance=np.tile(
            np.eye(3), (2, 1, 1)
        ).astype(complex),
        reference_station=ref,
        orientation=OrientationMeta(
            mode="orthogonal",
            angle_to_geographic_north=0.0,
        ),
    )
    xml = doc.to_xml()
    assert '<DataType name="TI"' in xml
    assert '<DataType name="AP"' in xml

    restored = EMTF.from_xml(xml)
    tensor = restored.get_transfer_function(AFMAG_TENSOR_TAG)
    ap = restored.get_transfer_function(AFMAG_AP_TAG)
    assert tensor is not None
    assert ap is not None
    assert tensor.shape == (2, 3, 2)
    assert ap.shape == (2, 1, 1)
    assert np.allclose(
        tensor.data,
        doc.get_transfer_function(AFMAG_TENSOR_TAG).data,
    )
    assert np.allclose(
        ap.data,
        doc.get_transfer_function(AFMAG_AP_TAG).data,
    )
    assert restored.processing.remote_reference.site == "BASE01"
    assert restored.metadata["notes"]["AFMAG"]["Family"] == (
        "tensor_afmag_airmt"
    )


def test_original_afmag_scalar_roundtrips_through_emtf_xml():
    doc = build_original_afmag_emtf(
        [1.25, -2.5],
        frequency=[150.0, 510.0],
        variance=[0.1, 0.2],
    )
    xml = doc.to_xml()
    assert '<DataType name="AFTILT"' in xml
    restored = EMTF.from_xml(xml)
    tf = restored.get_transfer_function(AFMAG_TILT_TAG)
    assert tf is not None
    assert tf.shape == (2, 1, 1)
    assert tf.units == "deg"
    assert np.allclose(tf.data[:, 0, 0], [1.25, -2.5])
    assert np.allclose(
        tf.get_estimate("VAR").data[:, 0, 0],
        [0.1, 0.2],
    )


def test_afmag_products_are_not_silently_coerced_to_standard_edi():
    airmt = build_airmt_emtf(_tensor(), frequency=[24.0, 75.0])
    original = build_original_afmag_emtf(
        [1.0, 2.0],
        frequency=[150.0, 510.0],
    )
    for document in (airmt, original):
        with pytest.raises(EMTFEDIConversionError):
            document.to_edi(on_loss="raise")
        with pytest.raises(EMTFEDIConversionError):
            document.to_edi(on_loss="warn")
