# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage tests for models/mare2dem/validation.py.

Pure filename/extension classification logic -- no bundled example data
required, so these tests always run (unlike most of the mare2dem suite,
which skips when ``data/mare2dem/`` is not present).
"""

from __future__ import annotations

import pytest

from pycsamt.models.mare2dem.validation import (
    Mare2DEMFileType,
    detect_file_type,
    is_emdata_file,
    is_log_file,
    is_response_file,
    is_resistivity_file,
    is_sensitivity_file,
    is_settings_file,
)


@pytest.mark.parametrize(
    "name, expected",
    [
        ("survey.emdata", Mare2DEMFileType.EMDATA),
        ("survey_MARE2DEM.emdata", Mare2DEMFileType.RESPONSE),
        ("Survey_MARE2DEM.EMDATA", Mare2DEMFileType.RESPONSE),
        ("iter00.resp", Mare2DEMFileType.RESPONSE),
        ("mare2dem.resistivity", Mare2DEMFileType.RESISTIVITY),
        ("mesh.poly", Mare2DEMFileType.POLY),
        ("mare2dem.settings", Mare2DEMFileType.SETTINGS),
        ("mare2dem.log", Mare2DEMFileType.LOG),
        ("mare2dem.logfile", Mare2DEMFileType.LOG),
        ("model.sensitivity", Mare2DEMFileType.SENSITIVITY),
        ("readme.txt", Mare2DEMFileType.UNKNOWN),
        ("no_extension", Mare2DEMFileType.UNKNOWN),
    ],
)
def test_detect_file_type(name, expected):
    assert detect_file_type(name) is expected


def test_detect_file_type_accepts_path_object(tmp_path):
    from pathlib import Path

    p = Path(tmp_path) / "demo.resistivity"
    assert detect_file_type(p) is Mare2DEMFileType.RESISTIVITY


def test_detect_file_type_case_insensitive_suffix_and_dir():
    assert detect_file_type("/some/DIR/demo.EMDATA") is Mare2DEMFileType.EMDATA
    assert (
        detect_file_type("/some/DIR/demo_mare2dem.EMDATA")
        is Mare2DEMFileType.RESPONSE
    )


def test_is_emdata_file():
    assert is_emdata_file("demo.emdata") is True
    assert is_emdata_file("demo_mare2dem.emdata") is False
    assert is_emdata_file("demo.resp") is False


def test_is_resistivity_file():
    assert is_resistivity_file("demo.resistivity") is True
    assert is_resistivity_file("demo.poly") is False


def test_is_settings_file():
    assert is_settings_file("mare2dem.settings") is True
    assert is_settings_file("mare2dem.log") is False


def test_is_log_file():
    assert is_log_file("demo.log") is True
    assert is_log_file("demo.logfile") is True
    assert is_log_file("demo.resistivity") is False


def test_is_response_file():
    assert is_response_file("demo.resp") is True
    assert is_response_file("demo_MARE2DEM.emdata") is True
    assert is_response_file("demo.emdata") is False


def test_is_sensitivity_file():
    assert is_sensitivity_file("demo.sensitivity") is True
    assert is_sensitivity_file("demo.resp") is False


def test_mare2dem_file_type_enum_members_unique():
    values = [member.value for member in Mare2DEMFileType]
    assert len(values) == len(set(values))
    assert Mare2DEMFileType.UNKNOWN in Mare2DEMFileType
