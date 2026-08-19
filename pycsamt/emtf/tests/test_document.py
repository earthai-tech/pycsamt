from __future__ import annotations

import numpy as np
import pytest

from pycsamt.core.base import TFBundle
from pycsamt.emtf import EMTF, StatisticalEstimate, TransferFunction


def _bundle() -> TFBundle:
    n = 3
    freq = np.array([100.0, 10.0, 1.0])
    z = np.arange(n * 4).reshape(n, 2, 2).astype(float) + 1j
    z_err = np.full((n, 2, 2), 0.1)
    tip = np.zeros((n, 2), dtype=complex)
    tip[:, 0] = 0.1 + 0.2j
    tip_err = np.full((n, 2), 0.02)
    rho = np.full((n, 2, 2), 100.0)
    phase = np.full((n, 2, 2), 45.0)
    return TFBundle(
        freq=freq,
        z=z,
        z_err=z_err,
        tipper=tip,
        tipper_err=tip_err,
        rho=rho,
        phase=phase,
        station="S001",
        station_id=1,
        lat=5.2,
        lon=-4.0,
        elev=120.0,
        azimuth=10.0,
        attrs={"source": "legacy"},
    )


def test_emtf_can_hold_multiple_transfer_functions():
    periods = np.array([0.1, 1.0])
    z = TransferFunction(
        name="Z",
        data=np.zeros((2, 2, 2), dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Ex", "Ey"),
        periods=periods,
    )
    t = TransferFunction(
        name="T",
        data=np.zeros((2, 1, 2), dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Hz",),
        periods=periods,
    )
    doc = EMTF(periods=periods)
    doc.add_transfer_function(z).add_transfer_function(t)
    assert doc.impedance is z
    assert doc.tipper_tf is t
    assert set(doc.tags) == {"impedance", "tipper"}


def test_document_rejects_inconsistent_period_grids():
    doc = EMTF(periods=[1.0, 2.0])
    tf = TransferFunction(
        name="tipper",
        data=np.zeros((2, 1, 2), dtype=complex),
        input_channels=("Hx", "Hy"),
        output_channels=("Hz",),
        periods=[1.0, 3.0],
    )
    with pytest.raises(ValueError):
        doc.add_transfer_function(tf)


def test_legacy_bundle_to_emtf_preserves_core_arrays_and_metadata():
    source = _bundle()
    doc = EMTF.from_bundle(source)

    assert np.allclose(doc.frequency, source.freq)
    assert np.allclose(doc.z, source.z)
    assert np.allclose(doc.z_err, source.z_err)
    assert doc.tipper.shape == (3, 1, 2)
    assert doc.tipper_err.shape == (3, 1, 2)
    assert doc.station == "S001"
    assert doc.lat == pytest.approx(5.2)
    assert doc.attrs["source"] == "legacy"

    # Existing numerical containers remain the compatibility surface.
    assert doc.Z is not None
    assert doc.Z.z.shape == (3, 2, 2)
    assert doc.Tip is not None
    assert doc.Tip.tipper.shape == (3, 1, 2)


def test_emtf_bundle_roundtrip_preserves_legacy_subset():
    source = _bundle()
    back = EMTF.from_bundle(source).to_bundle()
    assert np.allclose(back.freq, source.freq)
    assert np.allclose(back.z, source.z)
    assert np.allclose(back.z_err, source.z_err)
    assert np.allclose(back.tipper, source.tipper)
    assert np.allclose(back.tipper_err, source.tipper_err)
    assert np.allclose(back.rho, source.rho)
    assert np.allclose(back.phase, source.phase)
    assert back.station == source.station
    assert back.station_id == source.station_id


def test_generalized_bundle_payload_is_retained_in_memory():
    tf = TransferFunction(
        name="admittance",
        data=np.zeros((2, 3, 2), complex),
        input_channels=("Ex", "Ey"),
        output_channels=("Hx", "Hy", "Hz"),
    )
    bundle = TFBundle(transfer_functions={"admittance": tf})
    assert bundle.is_empty() is False

    doc = EMTF.from_bundle(bundle)
    assert doc.get_transfer_function("admittance").shape == (2, 3, 2)
    back = doc.to_bundle()
    assert isinstance(back.transfer_functions["admittance"], TransferFunction)


def test_standard_error_is_not_relabelled_as_emtf_variance():
    doc = EMTF.from_bundle(_bundle())
    estimate = doc.impedance.get_estimate("standard_error")
    assert isinstance(estimate, StatisticalEstimate)
    assert estimate.name == "STD"
    assert doc.impedance.get_estimate("variance") is None


def test_bundle_roundtrip_does_not_invent_derived_arrays():
    source = TFBundle(
        freq=np.array([10.0]),
        z=np.ones((1, 2, 2), dtype=complex),
    )
    back = EMTF.from_bundle(source).to_bundle()
    assert back.rho is None
    assert back.phase is None
