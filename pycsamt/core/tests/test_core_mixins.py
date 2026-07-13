from __future__ import annotations

import numpy as np
import pytest

from pycsamt.core import config as cfg
from pycsamt.core import mixins as m


def test_public_api_all():
    expected = {
        "bundle_from_edi",
        "BundleMixin",
        "BundleContainerMixin",
    }
    assert set(m.__all__) == expected


class _EDIStub:
    def __init__(self, n: int = 3) -> None:
        self.freq = np.linspace(1, 10, n)
        self.z = np.random.rand(n, 2, 2) + 1j * np.random.rand(n, 2, 2)
        self.z_err = np.full((n, 2, 2), 0.1)
        self.resistivity = np.full((n, 2, 2), 100.0)
        self.phase = np.full((n, 2, 2), 45.0)
        # tipper with shape (n, 1, 2) to test normalization
        self.tipper = np.zeros((n, 1, 2), dtype=complex)
        self.tipper_err = np.full((n, 2), 0.05)
        self.station = " S#01 "
        self.lat = 1.23
        self.lon = 4.56
        self.elev = 789.0
        self.azimuth = 0.0


def test_bundle_from_edi_extracts_and_normalizes():
    edi = _EDIStub(n=4)
    b = m.bundle_from_edi(edi)
    assert b.freq is not None and len(b.freq) == 4
    assert b.z is not None and b.z.shape == (4, 2, 2)
    assert b.z_err is not None and b.z_err.shape == (4, 2, 2)
    # tipper normalized to (n, 2)
    assert b.tipper is not None and b.tipper.shape == (4, 2)
    assert b.station.strip(" ") == "S#01"
    assert b.lat == pytest.approx(1.23)
    assert b.lon == pytest.approx(4.56)
    assert b.elev == pytest.approx(789.0)


class Host(m.BundleMixin):
    def __init__(self, bundle: m.TFBundle | None = None) -> None:
        self._bundle = bundle or m.TFBundle()

    def to_bundle(self) -> m.TFBundle:
        return self._bundle

    @classmethod
    def from_bundle(cls, bundle: m.TFBundle):
        return cls(bundle)


def test_bundle_mixin_ensure_station_name():
    name = Host.ensure_station_name("  K-01  ", None)
    assert name == "K-01"


def test_from_edi_single_and_collection():
    # single
    h = Host.from_edi(_EDIStub(2))
    assert isinstance(h, Host)
    assert h._bundle.freq is not None and len(h._bundle.freq) == 2

    # collection
    coll = [_EDIStub(1), _EDIStub(3)]
    out = Host.from_edi(coll)
    assert isinstance(out, list) and len(out) == 2
    assert isinstance(out[0], Host)


def test_to_edi_calls_adapter(monkeypatch):
    calls = {"n": 0}

    def adapter(obj, **kw):
        calls["n"] += 1
        # return the bundle to check dispatch
        return obj.to_bundle()

    cfg.register_adapter("mixins_test", adapter)

    h1 = Host(m.TFBundle(station="A"))
    out = h1.to_edi(key="mixins_test")

    assert calls["n"] == 1
    assert isinstance(out, m.TFBundle)
    assert out.station == "A"


class Item:
    def __init__(self, name: str):
        self._b = m.TFBundle(station=name)

    def to_bundle(self) -> m.TFBundle:
        return self._b


class Container(m.BundleContainerMixin):
    def __init__(self):
        self._items = {"A": Item("A"), "B": Item("B")}

    def items(self):  # mapping-like
        return list(self._items.items())


def test_container_iter_bundles_and_to_edi_collection():
    calls = {"n": 0}

    def adapter(obj, **kw):
        calls["n"] += 1
        return obj.to_bundle()

    cfg.register_adapter("mixins_coll", adapter)

    c = Container()
    bundles = list(c.iter_bundles())
    assert len(bundles) == 2
    assert {b.station for b in bundles} == {"A", "B"}

    edis = c.to_edi_collection(key="mixins_coll")
    assert len(edis) == 2
    assert calls["n"] == 2
    assert {e.station for e in edis} == {"A", "B"}
