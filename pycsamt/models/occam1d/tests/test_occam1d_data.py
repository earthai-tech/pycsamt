"""Scientific and native-format tests for :mod:`occam1d.data`."""

import numpy as np
import pytest

from pycsamt.models.occam1d import Occam1DData


def _data(**overrides):
    values = {
        "frequency": [100.0, 10.0],
        "resistivity": [50.0, 80.0],
        "phase": [42.0, 45.0],
        "resistivity_error": [0.05, 0.10],
        "phase_error": [1.0, 2.0],
        "station": "S01",
    }
    values.update(overrides)
    return Occam1DData(**values)


def test_scientific_properties_and_diagnostics():
    data = _data()

    assert data.n_frequencies == 2
    assert data.n_resistivity == 2
    assert data.n_phase == 2
    assert data.n_data == 4
    assert data.frequency_bounds == (10.0, 100.0)
    np.testing.assert_allclose(data.period, [0.01, 0.1])
    assert len(data.to_records()) == 2
    assert data.diagnostics()["frequency_max_hz"] == 100.0
    assert "station='S01'" in data.summary()


def test_constructor_copies_input_arrays():
    frequency = np.array([100.0, 10.0])
    data = _data(frequency=frequency)
    frequency[0] = 1.0

    assert data.frequency[0] == 100.0


def test_orphan_errors_are_removed_with_warning():
    data = _data(
        resistivity=[50.0, np.nan],
        resistivity_error=[0.05, 0.20],
    )

    assert np.isnan(data.resistivity_error[1])
    assert len(data.warnings) == 1
    assert "Discarded 1" in data.warnings[0]


@pytest.mark.parametrize(
    "overrides, message",
    [
        ({"frequency": [10.0, 10.0]}, "unique"),
        ({"station": "  "}, "station"),
        ({"resistivity_error": [np.nan, 0.1]}, "requires"),
        ({"phase_error": [0.0, 1.0]}, "positive"),
        (
            {"resistivity": [np.nan, np.nan], "phase": [np.nan, np.nan]},
            "finite observation",
        ),
    ],
)
def test_invalid_scientific_states_are_rejected(overrides, message):
    with pytest.raises(ValueError, match=message):
        _data(**overrides)


def test_native_roundtrip_preserves_units_metadata_and_binding(tmp_path):
    original = _data(
        resistivity=[50.0, np.nan],
        resistivity_error=[0.05, np.nan],
        mode="xy",
    )

    path = original.write(tmp_path / "Data")
    restored = Occam1DData.read(path)

    assert original.is_bound
    assert restored.is_bound
    assert restored.path == path.resolve()
    assert restored.station == "S01"
    assert restored.mode == "xy"
    np.testing.assert_allclose(
        restored.resistivity, original.resistivity, equal_nan=True
    )
    np.testing.assert_allclose(
        restored.resistivity_error,
        original.resistivity_error,
        equal_nan=True,
    )


def test_reader_accepts_case_spacing_and_fortran_exponents(tmp_path):
    path = tmp_path / "Data"
    path.write_text(
        "Format: EMData_1.0\n"
        "! Station: A01; Mode: TM\n"
        " # FREQUENCIES : 1\n"
        "1.0D+02\n"
        "# Receivers: 1\n0 0 0 0 0 0\n"
        " # DATA : 2\n"
        "103 1 0 1 2.0D+00 2.17147241D-02\n"
        "104 1 0 1 4.5D+01 2.0D+00\n",
        encoding="utf8",
    )

    data = Occam1DData.read(path)

    assert data.frequency[0] == 100.0
    assert data.resistivity[0] == pytest.approx(100.0)
    assert data.resistivity_error[0] == pytest.approx(0.05)
    assert data.mode == "tm"


@pytest.mark.parametrize(
    "row, declared, message",
    [
        ("103 2 0 1 2 0.1\n", 1, "outside"),
        ("103 1 0 1 nan 0.1\n", 1, "Non-finite"),
        ("103 1 0 1 1000 0.1\n", 1, "overflows"),
        ("103 1 0 1 2\n", 1, "six fields"),
        (
            "103 1 0 1 2 0.1\n103 1 0 1 2 0.1\n",
            2,
            "Duplicate",
        ),
        ("103 1 0 1 2 0.1\n", 2, "Declared"),
    ],
)
def test_reader_rejects_corrupt_native_rows(tmp_path, row, declared, message):
    path = tmp_path / "Data"
    path.write_text(
        "Format: EMData_1.1\n"
        "# Frequencies: 1\n100\n"
        "# Receivers: 1\n0 0 0 0 0 0\n"
        f"# Data: {declared}\n"
        "! Type Freq# Tx# Rx# Data Std_Error\n"
        f"{row}",
        encoding="utf8",
    )

    with pytest.raises(ValueError, match=message):
        Occam1DData.read(path)


def test_unknown_native_type_is_reported_as_warning(tmp_path):
    path = tmp_path / "Data"
    path.write_text(
        "Format: EMData_1.1\n"
        "# Frequencies: 1\n100\n"
        "# Receivers: 1\n0 0 0 0 0 0\n"
        "# Data: 2\n"
        "103 1 0 1 2 0.1\n"
        "999 1 0 1 3 0.2\n",
        encoding="utf8",
    )

    data = Occam1DData.read(path)

    assert data.n_data == 1
    assert "unsupported" in data.warnings[0]
