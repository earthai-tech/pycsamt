from __future__ import annotations

from pycsamt.core.base import TFBundle


def test_old_positional_tfbundle_signature_is_preserved():
    attrs = {"legacy": True}
    bundle = TFBundle(
        [1.0],
        [[[1.0]]],
        [[[0.1]]],
        [[0.0, 0.0]],
        [[0.1, 0.1]],
        [100.0],
        [45.0],
        "S01",
        1,
        5.0,
        -4.0,
        100.0,
        0.0,
        attrs,
    )
    assert bundle.station == "S01"
    assert bundle.attrs is attrs
    assert bundle.transfer_functions == {}
    assert bundle.estimates == {}


def test_legacy_empty_logic_is_unchanged():
    assert TFBundle().is_empty() is True
    assert TFBundle(z=[[1]]).is_empty() is False
    assert TFBundle(rho=[1.0]).is_empty() is True
    assert TFBundle(rho=[1.0], phase=[2.0]).is_empty() is False
