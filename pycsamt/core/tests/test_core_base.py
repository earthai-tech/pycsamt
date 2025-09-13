# -*- coding: utf-8 -*-
from __future__ import annotations

import pytest

from pycsamt.core import base as b
from pycsamt.core import config as regmod


def test_public_api_all():
    expected = {
        "TFBundle",
        "SupportsToBundle",
        "SupportsFromBundle",
        "ensure_station",
        "pick_adapter_key",
        "to_edi",
    }
    assert set(b.__all__) == expected


def test_tfbundle_is_empty_logic():
    # empty when neither Z nor (rho, phase) present
    x = b.TFBundle()
    assert x.is_empty() is True

    # Z present -> not empty
    y = b.TFBundle(z=[[1]])
    assert y.is_empty() is False

    # rho without phase -> still empty
    z = b.TFBundle(rho=[1.0])
    assert z.is_empty() is True

    # both rho and phase -> not empty
    w = b.TFBundle(rho=[1.0], phase=[2.0])
    assert w.is_empty() is False


def test_ensure_station_with_policy_override():
    pol = regmod.StationNamePolicy(prefix="STA", pad=4)
    name = b.ensure_station("  S#01  ", None, policy=pol)
    assert name == "S01"  # cleaned, not synthesized

    synthesized = b.ensure_station(None, 7, policy=pol)
    assert synthesized == "STA0007"


def test_pick_adapter_key_inference():
    # Craft dummy classes with module names to trigger inference
    Avg = type("Avg", (), {})
    Avg.__module__ = "pycsamt.zonge.avg"
    Jf = type("Jf", (), {})
    Jf.__module__ = "pycsamt.jones.j"
    Edi = type("Edi", (), {})
    Edi.__module__ = "pycsamt.seg.edi"

    assert b.pick_adapter_key(Avg()) == "avg"
    assert b.pick_adapter_key(Jf()) == "j"
    assert b.pick_adapter_key(Edi()) == "edi"

    # Hint should override detection
    assert b.pick_adapter_key(object(), hint="avg") == "avg"


def test_to_edi_errors_when_no_key_or_adapter():
    class Unknown:  # no recognizable module name
        pass

    with pytest.raises(RuntimeError):
        b.to_edi(Unknown())

    # Explicit key but no adapter registered
    with pytest.raises(RuntimeError):
        b.to_edi(object(), key="avg")


def test_to_edi_with_registered_adapter(monkeypatch):
    # Register a simple adapter under the key inferred by pick_adapter_key
    calls = {}

    def dummy_avg_adapter(src, **kw):
        calls["called"] = True
        calls["kw"] = kw
        return {"ok": True, "src": src}

    regmod.register_adapter("avg", dummy_avg_adapter)

    Avg = type("Avg", (), {})
    Avg.__module__ = "pycsamt.zonge.avg"

    out = b.to_edi(Avg(), some="arg")
    assert out["ok"] is True
    assert calls.get("called") is True
    assert calls.get("kw", {}).get("some") == "arg"
