"""Borehole Builder page and integration contracts."""

from __future__ import annotations

from pycsamt.app.web.pages import borehole_builder
from pycsamt.format.borehole import document_from_builder, new_builder_draft


def test_builder_page_exposes_all_editors_and_accessible_feedback():
    rendered = str(borehole_builder.layout())

    for name in (
        "boreholes",
        "surveys",
        "intervals",
        "structures",
        "water",
        "construction",
        "samples",
        "assays",
    ):
        assert f"pcbh-builder-{name}" in rendered
    assert "pcbh-builder-validation" in rendered
    assert "role='alert'" in rendered or '"role":"alert"' in rendered
    assert "pcbh-builder-draft" in rendered
    assert "storage_type='local'" in rendered or '"storage_type":"local"' in rendered


def test_builder_callbacks_are_registered(web_app):
    callback_text = str(web_app.callback_map)
    assert "pcbh-builder-document" in callback_text
    assert "pcbh-builder-download" in callback_text
    assert "pcbh-builder-pcsf-download" in callback_text
    assert "pcbh-builder-csv-preview" in callback_text


def test_navigation_contains_builder_page(web_app):
    layout_text = str(web_app.layout)
    assert "nav-btn-borehole-builder" in layout_text
    assert "page-borehole-builder" in layout_text


def test_large_builder_draft_constructs_many_holes():
    draft = new_builder_draft()
    draft["project"]["crs_horizontal"] = "EPSG:32629"
    draft["boreholes"] = [
        {
            "id": f"BH-{index:04d}",
            "name": f"Hole {index}",
            "kind": "mining_exploration",
            "status": "planned",
            "x": 500000 + index,
            "y": 600000,
            "z": 100,
            "total_depth_md": 10,
        }
        for index in range(1_000)
    ]

    document = document_from_builder(draft)
    assert len(document.boreholes) == 1_000
