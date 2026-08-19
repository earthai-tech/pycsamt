"""Scientific, dialect, and alignment tests for Occam1D responses."""

import math

import numpy as np
import pytest

from pycsamt.models.occam1d import Occam1DData, Occam1DResponse


def _response():
    return Occam1DResponse(
        [
            [1, 1, 103, 0.1, 2.0, 2.1, -1.0],
            [1, 1, 104, 2.0, 45.0, 44.0, 0.5],
            [1, 2, 103, 0.2, 3.0, 2.8, 1.0],
        ]
    )


def _data():
    return Occam1DData(
        [100.0, 10.0],
        [100.0, 1000.0],
        [45.0, 50.0],
        [0.1 * math.log(10.0), 0.2 * math.log(10.0)],
        [2.0, 2.0],
    )


def test_response_properties_statistics_and_diagnostics():
    response = _response()

    assert response.n_data == 3
    assert response.n_residuals == 3
    assert response.type_codes.tolist() == [103, 104]
    assert response.rms == pytest.approx(math.sqrt(2.25 / 3))
    assert response.mean_residual == pytest.approx(1 / 6)
    assert response.mean_absolute_residual == pytest.approx(2.5 / 3)
    assert response.chi_square == pytest.approx(2.25)
    assert response.misfit_per_frequency()[1] == pytest.approx(
        math.sqrt(0.625)
    )
    assert response.misfit_per_type()[103] == pytest.approx(1.0)
    assert response.diagnostics()["n_data"] == 3
    assert "types=[103, 104]" in response.summary()


def test_constructor_copies_rows_and_empty_response_is_valid():
    rows = np.asarray([[1, 1, 103, 0.1, 2.0, 2.1, -1.0]])
    response = Occam1DResponse(rows)
    rows[0, 4] = 9.0

    assert response.observed[0] == 2.0
    empty = Occam1DResponse()
    assert empty.rows.shape == (0, 7)
    assert np.isnan(empty.rms)


@pytest.mark.parametrize(
    "rows, message",
    [
        ([[1, 1, 103]], "shape"),
        ([[0, 1, 103, 1, 1, 1, 0]], "positive integers"),
        ([[1, 1.5, 103, 1, 1, 1, 0]], "integer-valued"),
        ([[1, 1, 103, -1, 1, 1, 0]], "non-negative"),
        ([[1, 1, 103, 1, 1, np.inf, 0]], "infinite"),
        (
            [
                [1, 1, 103, 1, 1, 1, 0],
                [1, 1, 103, 1, 1, 1, 0],
            ],
            "unique",
        ),
    ],
)
def test_invalid_response_tables_are_rejected(rows, message):
    with pytest.raises(ValueError, match=message):
        Occam1DResponse(rows)


def test_canonical_reader_supports_fortran_exponents_and_headers(tmp_path):
    path = tmp_path / "RESP_01.resp"
    path.write_text(
        "response table\n"
        "1 1 103 1.0D-01 2.0D+00 2.1D+00 -1.0D+00\n"
        "1 1 104 2 45 44 5.0D-01\n",
        encoding="utf8",
    )

    response = Occam1DResponse.read(path)

    assert response.is_bound
    assert response.path == path.resolve()
    assert response.n_data == 2
    assert response.rms == pytest.approx(math.sqrt(0.625))


def test_emdata_seven_column_dialect_is_not_misclassified(tmp_path):
    path = tmp_path / "response.data"
    path.write_text(
        "103 1 0 1 2.0 0.1 2.1\n"
        "104 1 0 1 45.0 2.0 44.0\n",
        encoding="utf8",
    )

    response = Occam1DResponse.read(path)

    assert response.site_index.tolist() == [1, 1]
    assert response.frequency_index.tolist() == [1, 1]
    assert response.type_code.tolist() == [103, 104]
    np.testing.assert_allclose(response.residuals, [-1.0, 0.5])


def test_strict_and_permissive_malformed_row_handling(tmp_path):
    path = tmp_path / "response.resp"
    path.write_text(
        "1 1 103 bad 2 2.1 -1\n"
        "1 1 104 2 45 44 0.5\n",
        encoding="utf8",
    )

    with pytest.raises(ValueError, match="bad"):
        Occam1DResponse.read(path)

    response = Occam1DResponse.read(path, strict=False)
    assert response.n_data == 1
    assert "Skipped 1" in response.warnings[0]


def test_physical_conversion_errors_records_selection_and_csv(tmp_path):
    response = _response()
    observed, modeled = response.physical_values()

    assert observed.tolist() == [100.0, 45.0, 1000.0]
    assert modeled[0] == pytest.approx(10**2.1)
    assert modeled[1] == 44.0
    errors = response.physical_errors()
    assert errors[0] == pytest.approx(100 * math.log(10) * 0.1)
    assert errors[1] == 2.0
    rho = response.select(type_code=103)
    assert rho.n_data == 2
    assert len(response.to_records(physical=True)) == 3
    path = response.to_csv(tmp_path / "response.csv")
    assert len(path.read_text(encoding="utf8").splitlines()) == 4


def test_computed_residuals_exclude_zero_or_missing_errors():
    response = Occam1DResponse(
        [
            [1, 1, 103, 0.0, 2.0, 2.1, np.nan],
            [1, 1, 104, np.nan, 45.0, 44.0, np.nan],
        ]
    )

    assert np.all(np.isnan(response.computed_residuals()))
    assert response.n_residuals == 0
    assert len(response.warnings) == 1


def test_unknown_type_is_preserved_with_warning():
    response = Occam1DResponse([[1, 1, 999, 1, 2, 3, -1]])

    assert response.type_codes.tolist() == [999]
    assert "999" in response.warnings[0]
    observed, modeled = response.physical_values()
    assert observed[0] == 2
    assert modeled[0] == 3


def test_response_validates_against_sounding():
    response = _response()
    response.validate_against(_data())

    bad = response.rows.copy()
    bad[0, 4] = 2.2
    with pytest.raises(ValueError, match="row 1"):
        Occam1DResponse(bad).validate_against(_data())
    bad = response.rows.copy()
    bad[0, 1] = 3
    with pytest.raises(ValueError, match="outside"):
        Occam1DResponse(bad).validate_against(_data())
    with pytest.raises(TypeError, match="Occam1DData"):
        response.validate_against(object())


def test_physical_resistivity_overflow_is_reported():
    response = Occam1DResponse([[1, 1, 103, 1, 1000, 2, 0]])

    with pytest.raises(ValueError, match="exceeds"):
        response.physical_values()
