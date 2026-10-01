# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Deterministic request preconditions answered before any dispatch."""

from __future__ import annotations

import pytest

from pycsamt.app.agent_master._preconditions import (
    check_request,
    no_data_message,
)


def test_nonexistent_api_is_named_and_not_called():
    kind, message = check_request("Use pycsamt.run_workflow to process my survey.")
    assert kind == "answer"
    assert "`pycsamt.run_workflow` does not exist" in message


def test_existing_api_passes():
    assert check_request("How do I use pycsamt.emtools.ensure_sites?") is None


def test_missing_input_path_asks_for_a_source(tmp_path):
    missing = tmp_path / "not_here.edi"
    kind, message = check_request(f"Analyze {missing} and tell me its best station.")
    assert kind == "clarify" and "does not exist" in message
    real = tmp_path / "real.edi"
    real.write_text("")
    assert check_request(f"Analyze {real}.") is None


def test_urls_are_not_local_paths():
    assert check_request("See https://example.org/data/site.edi for format.") is None


@pytest.mark.parametrize("text", [
    "Apply these static-shift factors: NaN, zero, and minus two.",
    "Use factors of 1.2 and -0.5",
    "apply correction factors: inf",
])
def test_invalid_factors_are_rejected_without_correction(text):
    kind, message = check_request(text)
    assert kind == "answer"
    assert "finite and strictly positive" in message
    assert "No correction was applied" in message


@pytest.mark.parametrize("text", [
    "Apply these static-shift factors: 1.2, 0.8",
    "apply static shift correction with AMA",
    "what are static-shift factors?",
])
def test_valid_or_absent_factors_pass(text):
    assert check_request(text) is None


@pytest.mark.parametrize("text, angle, asks", [
    ("Rotate the data by 45 degrees.", None, True),
    ("Rotate the data by 45 degrees clockwise.", 45.0, False),
    ("rotate 30° counter-clockwise", -30.0, False),
    ("rotate by 20 degrees east of north", 20.0, False),
    ("rotate by -15 deg anticlockwise", -15.0, False),
    ("rotate to the strike", None, False),
    ("rotate the data to principal axes", None, False),
    ("Rotate the data.", None, False),
])
def test_rotation_request(text, angle, asks):
    from pycsamt.app.agent_master._preconditions import rotation_request

    result = rotation_request(text)
    assert result["angle"] == angle
    assert bool(result["question"]) is asks


def test_rotation_sign_matches_documented_convention():
    """Positive = clockwise from north: +90 deg turns the new x axis to east."""
    import numpy as np

    from pycsamt.seg.ops import rotate_tipper

    north = np.array([1.0, 0.0])
    np.testing.assert_allclose(rotate_tipper(north, 90.0), [0.0, -1.0], atol=1e-12)
    east = np.array([0.0, 1.0])
    np.testing.assert_allclose(rotate_tipper(east, 90.0), [1.0, 0.0], atol=1e-12)


class Registry:
    def lines(self):
        return ["L18PLT", "L22PLT"]


def test_no_data_messages_are_request_specific():
    unknown = no_data_message("qc", "Run QC on L99PLT.", "quality control", Registry())
    assert "Line L99PLT is not in the project registry" in unknown
    assert "L18PLT, L22PLT" in unknown
    known = no_data_message("qc", "Run QC on L22PLT.", "quality control", Registry())
    assert "Quality control needs station data" in known
    assert "installed ModEM executable" in no_data_message("modem", "Run ModEM", "ModEM", None)
    assert "invented data or settings" in no_data_message(
        "ai_inversion", "Invert my survey.", "1-D AI inversion", None)


def test_rotation_question_offers_one_click_replies():
    from pycsamt.app.agent_master._preconditions import rotation_request

    replies = [s["reply"] for s in rotation_request("Rotate the data by 45 degrees.")["suggestions"]]
    assert replies == ["Rotate the data by 45 degrees clockwise.",
                       "Rotate the data by 45 degrees counter-clockwise.",
                       "Rotate the data to the strike."]
    # each reply resolves without asking again
    assert all(rotation_request(r)["question"] is None for r in replies)
    assert rotation_request("rotate 30 deg clockwise")["suggestions"] == []


def test_unknown_line_offers_registered_lines():
    from pycsamt.app.agent_master._preconditions import line_suggestions

    got = line_suggestions("Run QC on L99PLT.", Registry())
    assert got == [{"label": "L18PLT", "reply": "Run QC on L18PLT."},
                   {"label": "L22PLT", "reply": "Run QC on L22PLT."}]
    assert line_suggestions("Run QC on L22PLT.", Registry()) == []
    assert line_suggestions("Run QC on L99PLT.", None) == []
