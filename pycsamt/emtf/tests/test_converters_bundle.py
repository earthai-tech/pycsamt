from __future__ import annotations

import numpy as np
import pytest

from pycsamt.core.base import TFBundle
from pycsamt.emtf import EMTF
from pycsamt.emtf.converters.bundle import bundle_to_emtf, emtf_to_bundle


def _sample_bundle() -> TFBundle:
    freq = np.array([100.0, 10.0])
    z = np.stack(
        [np.eye(2, dtype=complex) * (i + 1) for i in range(2)]
    )
    return TFBundle(
        freq=freq,
        z=z,
        station="S01",
        lat=5.1,
        lon=-3.2,
    )


def test_bundle_to_emtf_delegates_to_from_bundle():
    bundle = _sample_bundle()
    doc = bundle_to_emtf(bundle)
    assert isinstance(doc, EMTF)
    assert doc.station == "S01"
    np.testing.assert_allclose(doc.z, bundle.z)


def test_emtf_to_bundle_delegates_to_to_bundle():
    bundle = _sample_bundle()
    doc = bundle_to_emtf(bundle)
    back = emtf_to_bundle(doc)
    assert isinstance(back, TFBundle)
    assert back.station == "S01"
    np.testing.assert_allclose(back.z, bundle.z)


def test_emtf_to_bundle_rejects_non_emtf_input():
    with pytest.raises(TypeError, match="pycsamt.emtf.EMTF"):
        emtf_to_bundle(_sample_bundle())
