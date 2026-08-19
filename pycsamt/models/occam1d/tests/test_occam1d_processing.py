"""Tests for site-to-Occam1D scientific processing."""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from pycsamt.models.occam1d import Occam1DConfig
from pycsamt.models.occam1d.processing import (
    extract_sounding,
    normalise_sites,
)


def _impedance_site(name="S01"):
    z = np.zeros((3, 2, 2), dtype=complex)
    z[:, 0, 1] = [1 + 1j, 2 + 2j, 3 + 3j]
    z[:, 1, 0] = -z[:, 0, 1]
    return SimpleNamespace(
        name=name,
        freq=np.array([1.0, 100.0, 10.0]),
        z=z,
        z_err=np.full((3, 2, 2), 0.01),
    )


def test_component_uses_canonical_pycsamt_impedance_units():
    site = _impedance_site()

    data = extract_sounding(site, Occam1DConfig(mode="xy"))

    order = [1, 2, 0]
    expected = 0.2 * np.abs(site.z[:, 0, 1]) ** 2 / site.freq
    np.testing.assert_allclose(data.resistivity, expected[order])
    np.testing.assert_allclose(data.phase, 45.0)
    assert data.metadata["resistivity_error_source"] == (
        "propagated impedance error"
    )
    assert data.metadata["frequency_order"] == "descending"


def test_determinant_response_and_error_propagation():
    site = _impedance_site()

    data = extract_sounding(site, Occam1DConfig(mode="determinant"))

    order = [1, 2, 0]
    determinant = np.sqrt(
        -site.z[:, 0, 1] * site.z[:, 1, 0]
    )
    expected = 0.2 * np.abs(determinant) ** 2 / site.freq
    np.testing.assert_allclose(data.resistivity, expected[order])
    assert np.all(data.resistivity_error >= 0.05)
    assert np.all(data.phase_error >= 2.0)


def test_determinant_phase_uses_conventional_mt_branch():
    site = _impedance_site()
    site.z[:, 0, 1] = [1 - 1j, 2 - 2j, 3 - 3j]
    site.z[:, 1, 0] = -site.z[:, 0, 1]

    data = extract_sounding(site, Occam1DConfig(mode="determinant"))

    np.testing.assert_allclose(data.phase, 45.0)
    assert np.all((data.phase >= 0.0) & (data.phase <= 90.0))


def test_direct_physical_vectors_take_precedence_and_errors_are_relative():
    site = SimpleNamespace(
        name="P01",
        freq=[10.0, 1.0],
        rho=[100.0, 200.0],
        phase=[40.0, 50.0],
        rho_err=[10.0, 40.0],
        phase_err=[1.0, 3.0],
    )
    config = Occam1DConfig(
        mode="xy", error_floor_rho=0.05, error_floor_phase=2.0
    )

    data = extract_sounding(site, config)

    np.testing.assert_allclose(data.resistivity, [100.0, 200.0])
    np.testing.assert_allclose(data.resistivity_error, [0.1, 0.2])
    np.testing.assert_allclose(data.phase_error, [2.0, 3.0])
    assert "supplied absolute" in data.metadata["resistivity_error_source"]


def test_tensor_physical_values_override_derived_impedance():
    site = _impedance_site()
    site.rho = np.full((3, 2, 2), 999.0)
    site.phase = np.full((3, 2, 2), 33.0)

    data = extract_sounding(site, Occam1DConfig(mode="xy"))

    np.testing.assert_allclose(data.resistivity, 999.0)
    np.testing.assert_allclose(data.phase, 33.0)


def test_tm_phase_is_mapped_to_conventional_branch():
    site = _impedance_site()

    data = extract_sounding(site, Occam1DConfig(mode="tm"))

    np.testing.assert_allclose(data.phase, 45.0)


def test_normalise_sites_loads_a_single_edi_file_path():
    repository = Path(__file__).resolve().parents[4]
    source = repository / "data" / "3edis" / "new_E1_1.edi"

    sites = normalise_sites(source)
    data = extract_sounding(sites[0], Occam1DConfig(mode="determinant"))

    assert len(sites) == 1
    assert data.station == "E1_1"
    assert data.n_frequencies == 53
    assert np.all((data.phase >= 0.0) & (data.phase <= 90.0))


def test_frequency_and_observation_filtering_and_provenance():
    site = SimpleNamespace(
        name="F01",
        freq=[1.0, 10.0, 100.0, np.nan],
        rho=[100.0, np.nan, 300.0, 400.0],
        phase=[40.0, np.nan, 50.0, 60.0],
    )
    config = Occam1DConfig(mode="xy", freq_min=5.0, freq_max=200.0)

    data = extract_sounding(site, config)

    np.testing.assert_allclose(data.frequency, [100.0])
    assert data.metadata["source_frequency_count"] == 4
    assert data.metadata["selected_frequency_count"] == 1


def test_missing_uncertainties_receive_floors_only_for_observations():
    site = SimpleNamespace(
        freq=[10.0, 1.0],
        rho=[100.0, np.nan],
        phase=[np.nan, 45.0],
    )

    data = extract_sounding(site, Occam1DConfig(mode="xy"))

    np.testing.assert_allclose(data.resistivity_error[0], 0.05)
    assert np.isnan(data.resistivity_error[1])
    assert np.isnan(data.phase_error[0])
    np.testing.assert_allclose(data.phase_error[1], 2.0)


@pytest.mark.parametrize(
    "site, message",
    [
        (SimpleNamespace(rho=[1.0]), "frequency vector"),
        (
            SimpleNamespace(freq=[1.0], rho=[1.0], phase=[45.0]),
            "Determinant mode requires",
        ),
        (
            SimpleNamespace(
                freq=[1.0, 2.0],
                z=np.zeros((3, 2, 2), dtype=complex),
            ),
            "got",
        ),
        (
            SimpleNamespace(
                freq=[1.0],
                z_err=np.ones((1, 2, 2)),
                rho=[1.0],
                phase=[45.0],
            ),
            "without an impedance",
        ),
    ],
)
def test_invalid_sources_raise_descriptive_errors(site, message):
    with pytest.raises(ValueError, match=message):
        extract_sounding(site)


def test_duplicate_selected_frequencies_are_rejected():
    site = SimpleNamespace(
        freq=[10.0, 10.0],
        rho=[100.0, 200.0],
        phase=[40.0, 50.0],
    )

    with pytest.raises(ValueError, match="unique"):
        extract_sounding(site, Occam1DConfig(mode="xy"))


def test_public_configuration_and_mode_are_validated():
    site = _impedance_site()

    with pytest.raises(TypeError, match="config"):
        extract_sounding(site, config={})
    with pytest.raises(TypeError, match="mode"):
        extract_sounding(site, mode=1)
    with pytest.raises(ValueError, match="mode must be"):
        extract_sounding(site, mode="invalid")


def test_normalise_sites_supports_single_mapping_and_iterable():
    first = _impedance_site("A")
    second = _impedance_site("B")

    assert normalise_sites(first) == [first]
    assert normalise_sites({"A": first, "B": second}) == [first, second]
    assert normalise_sites(item for item in (first, second)) == [first, second]
    with pytest.raises(TypeError, match="None"):
        normalise_sites(None)
    with pytest.raises(TypeError, match="source must"):
        normalise_sites(42)
