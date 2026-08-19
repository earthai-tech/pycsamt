from __future__ import annotations

import numpy as np
import pytest

from pycsamt.core.base import TFBundle
from pycsamt.emtf import EMTF, TransferFunction
from pycsamt.metadata import (
    ChannelMeta,
    LocationMeta,
    OrientationMeta,
    ProcessingMeta,
    ProvenanceMeta,
    SiteLayout,
    SiteMeta,
    TransferFunctionQuality,
)


def test_emtf_accepts_format_neutral_metadata_without_seg_dependency():
    site = SiteMeta(
        project="Demo",
        survey="Survey-A",
        site_id="S001",
        name="Station 1",
        location=LocationMeta(
            latitude=5.2,
            longitude=-4.0,
            elevation=120.0,
        ),
    )
    layout = SiteLayout(
        input_channels=[
            ChannelMeta("Hx", "H", orientation=0.0),
            ChannelMeta("Hy", "H", orientation=90.0),
        ],
        output_channels=[
            ChannelMeta("Hz", "H"),
            ChannelMeta("Ex", "E", orientation=0.0),
            ChannelMeta("Ey", "E", orientation=90.0),
        ],
    )
    doc = EMTF(
        site=site,
        site_layout=layout,
        orientation=OrientationMeta(
            mode="orthogonal",
            angle_to_geographic_north=0.0,
        ),
        processing=ProcessingMeta(sign_convention="exp(+i w t)"),
        provenance=ProvenanceMeta(creating_application="pyCSAMT"),
        quality=TransferFunctionQuality(rating=4),
    )
    assert doc.site is site
    assert doc.station == "S001"
    assert doc.station_id == "S001"
    assert doc.lat == pytest.approx(5.2)
    assert doc.lon == pytest.approx(-4.0)
    assert doc.processing.sign_convention == "exp(+i ω t)"
    assert doc.site_layout.input_names == ("Hx", "Hy")


def test_legacy_bundle_remains_legacy_until_format_adapter_maps_metadata():
    bundle = TFBundle(
        freq=np.array([10.0, 1.0]),
        z=np.ones((2, 2, 2), dtype=complex),
        station="S002",
        station_id=2,
        lat=6.0,
        lon=-5.0,
        elev=80.0,
        attrs={"legacy": True},
    )
    doc = EMTF.from_bundle(bundle)
    # Phase 2 does not guess EMTF Site/Id versus Site/Name semantics from
    # legacy TFBundle fields. That mapping belongs to the EDI adapter.
    assert doc.site is None
    assert doc.station == "S002"
    assert doc.station_id == 2

    back = doc.to_bundle()
    assert back.station == "S002"
    assert back.station_id == 2
    assert back.lat == pytest.approx(6.0)
    assert back.attrs == {"legacy": True}


def test_duplicate_legacy_and_site_values_must_not_conflict():
    site = SiteMeta(
        name="S001",
        location=LocationMeta(latitude=5.0, longitude=-4.0),
    )
    EMTF(site=site, station="S001", lat=5.0, lon=-4.0)

    with pytest.raises(ValueError):
        EMTF(site=site, station="OTHER")
    with pytest.raises(ValueError):
        EMTF(site=site, lat=6.0)


def test_emtf_rejects_wrong_metadata_object_types():
    with pytest.raises(TypeError):
        EMTF(site={"name": "S001"})
    with pytest.raises(TypeError):
        EMTF(processing="processed")


def test_site_layout_does_not_constrain_arbitrary_transfer_function_matrix():
    layout = SiteLayout(
        input_channels=[
            ChannelMeta("Ex", "E"),
            ChannelMeta("Ey", "E"),
        ],
        output_channels=[
            ChannelMeta("Hx", "H"),
            ChannelMeta("Hy", "H"),
            ChannelMeta("Hz", "H"),
        ],
    )
    admittance = TransferFunction(
        name="admittance",
        data=np.zeros((4, 3, 2), dtype=complex),
        input_channels=layout.input_names,
        output_channels=layout.output_names,
    )
    doc = EMTF(site_layout=layout)
    doc.add_transfer_function(admittance)
    assert doc.get_transfer_function("admittance").shape == (4, 3, 2)
