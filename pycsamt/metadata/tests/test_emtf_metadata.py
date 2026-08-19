from __future__ import annotations

from datetime import datetime, timezone

import pytest

from pycsamt.metadata import (
    ChannelMeta,
    LocationMeta,
    OrientationMeta,
    Person,
    ProcessingMeta,
    ProvenanceMeta,
    QualityComment,
    RemoteReferenceMeta,
    SiteLayout,
    SiteMeta,
    Software,
    TransferFunctionQuality,
    normalize_sign_convention,
)


def test_location_meta_validation_and_coordinates():
    location = LocationMeta(
        latitude=45.532102,
        longitude=-69.864103,
        elevation=361.725,
        datum="WGS84",
        declination=15.3,
        declination_epoch=1995.0,
    )
    assert location.has_horizontal_coordinates
    assert location.elevation_units == "meters"

    with pytest.raises(ValueError):
        LocationMeta(latitude=95.0)
    with pytest.raises(ValueError):
        LocationMeta(longitude=181.0)


def test_site_meta_preserves_identity_and_time():
    start = datetime(2026, 1, 1, tzinfo=timezone.utc)
    end = datetime(2026, 1, 2, tzinfo=timezone.utc)
    site = SiteMeta(
        project="EarthScope",
        survey="Demo",
        year_collected=2026,
        country="USA",
        site_id="ID001",
        name="Example site",
        location=LocationMeta(latitude=40.0, longitude=-105.0),
        acquired_by="Field team",
        start=start,
        end=end,
        extra={"custom": "kept"},
    )
    assert site.preferred_name == "Example site"
    assert site.location.latitude == pytest.approx(40.0)
    assert site.extra["custom"] == "kept"

    with pytest.raises(ValueError):
        SiteMeta(year_collected=-1)
    with pytest.raises(ValueError):
        SiteMeta(start=end, end=start)


def test_orientation_separates_data_frame_from_site_layout():
    orientation = OrientationMeta(
        mode="orthogonal",
        angle_to_geographic_north=20.0,
        rotation_info="rotated after processing",
    )
    assert orientation.is_orthogonal
    assert not orientation.follows_site_layout

    site_layout = OrientationMeta(mode="site_layout")
    assert site_layout.mode == "sitelayout"
    assert site_layout.follows_site_layout

    with pytest.raises(ValueError):
        OrientationMeta(
            mode="sitelayout",
            angle_to_geographic_north=10.0,
        )


def test_channel_and_site_layout_support_general_input_output_geometry():
    hx = ChannelMeta(
        name="Hx",
        field_type="H",
        orientation=0.0,
        x=0.0,
        y=0.0,
        z=0.0,
    )
    hy = ChannelMeta(name="Hy", field_type="magnetic", orientation=90.0)
    ex = ChannelMeta(
        name="Ex",
        field_type="E",
        orientation=0.0,
        x=-50.0,
        y=0.0,
        x2=50.0,
        y2=0.0,
    )
    hz = ChannelMeta(name="Hz", field_type="H", orientation=0.0)
    layout = SiteLayout(
        input_channels=[hx, hy],
        output_channels=[hz, ex],
        input_units="m",
        output_units="m",
        input_reference="site",
        output_reference="site",
    )
    assert layout.input_names == ("Hx", "Hy")
    assert layout.output_names == ("Hz", "Ex")
    assert layout.n_input == 2
    assert layout.get_channel("ex", role="output") is ex
    assert ex.is_electric and ex.has_second_endpoint

    with pytest.raises(ValueError):
        SiteLayout(input_channels=[hx, hx])


def test_provenance_reuses_existing_person_metadata():
    creator = Person(
        name="Creator",
        email="creator@example.org",
        organization="Institute",
    )
    provenance = ProvenanceMeta(
        create_time="2026-08-15T12:00:00Z",
        creating_application="pyCSAMT 2.3",
        creator=creator,
        submitter=Person(name="Submitter"),
        extra={"archive_id": "A1"},
    )
    assert provenance.creator is creator
    assert provenance.extra["archive_id"] == "A1"

    with pytest.raises(TypeError):
        ProvenanceMeta(creator="not-a-person")


def test_processing_meta_keeps_unknown_sign_convention_instead_of_guessing():
    assert normalize_sign_convention("exp(+ i\\omega t)") == "exp(+i ω t)"
    assert normalize_sign_convention("exp(-i w t)") == "exp(-i ω t)"
    assert normalize_sign_convention("custom Fourier convention") == (
        "custom Fourier convention"
    )

    meta = ProcessingMeta(
        sign_convention="exp(+ i\\omega t)",
        processed_by="A. Processor",
        software=Software(name="EMTF", version="4.1"),
        remote_reference=RemoteReferenceMeta(
            reference_type="Robust Remote Reference",
            site="RR01",
        ),
        run_list=["run01", "run02"],
    )
    assert meta.sign_convention == "exp(+i ω t)"
    assert meta.remote_reference.site == "RR01"

    with pytest.raises(TypeError):
        ProcessingMeta(run_list="run01")


def test_transfer_function_quality_is_distinct_from_computed_data_quality():
    quality = TransferFunctionQuality(
        rating=5,
        good_from_period=8.0,
        good_to_period=20000.0,
        comments=[QualityComment("great TF", author="Analyst")],
        warning_flag=1,
        warnings=[QualityComment("cultural noise around 60 s")],
        extra={"source": "expert review"},
    )
    assert quality.rating_label == "great"
    assert quality.has_warning
    assert quality.comments[0].author == "Analyst"

    with pytest.raises(ValueError):
        TransferFunctionQuality(rating=6)
    with pytest.raises(ValueError):
        TransferFunctionQuality(
            good_from_period=100.0,
            good_to_period=10.0,
        )
