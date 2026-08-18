from __future__ import annotations

import pycsamt


def test_airborne_is_lazy_root_subpackage():
    assert hasattr(pycsamt, "airborne")
    assert pycsamt.airborne.AirborneEMDataset.__name__ == "AirborneEMDataset"
