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


def test_get_returns_none_when_nothing_matches():
    class Empty:
        pass

    b = m.bundle_from_edi(Empty())
    assert b.freq is None
    assert b.z is None
    assert b.station is None


def test_bundle_from_edi_leaves_non_3d_tipper_untouched():
    class E:
        tipper = np.zeros((5, 2), dtype=complex)  # already (n, 2)

    b = m.bundle_from_edi(E())
    assert b.tipper.shape == (5, 2)


def test_bundle_from_edi_tipper_3d_non_1x2_left_untouched():
    class E:
        tipper = np.zeros((5, 2, 2), dtype=complex)  # ndim==3 but not (1, 2)

    b = m.bundle_from_edi(E())
    assert b.tipper.shape == (5, 2, 2)


def test_bundle_from_edi_tipper_shape_probe_exception_is_swallowed():
    class _BadTipper:
        ndim = 3  # lies about ndim, has no .shape

    class E:
        tipper = _BadTipper()

    b = m.bundle_from_edi(E())
    assert isinstance(b.tipper, _BadTipper)


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


def test_looks_collection_heuristics():
    assert m._looks_collection([1, 2]) is True
    assert m._looks_collection((1, 2)) is True
    assert m._looks_collection({1, 2}) is True
    assert m._looks_collection("a string") is False
    assert m._looks_collection(b"bytes") is False
    # generic iterable without a `z` attribute -> collection
    assert m._looks_collection(iter([1, 2])) is True
    # has __iter__ but also a `z` attribute -> treated as a single item,
    # not a collection
    class IterableWithZ:
        z = 1

        def __iter__(self):
            return iter([])

    assert m._looks_collection(IterableWithZ()) is False
    # neither iterable nor list-like -> False
    assert m._looks_collection(object()) is False


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


def test_container_items_skips_entries_without_to_bundle():
    class Container2(m.BundleContainerMixin):
        def __init__(self):
            self._items = {"A": Item("A"), "B": object()}

        def items(self):
            return list(self._items.items())

    c = Container2()
    bundles = list(c.iter_bundles())
    assert len(bundles) == 1
    assert bundles[0].station == "A"


class ListContainer(m.BundleContainerMixin):
    """No items() method: falls back to plain __iter__."""

    def __init__(self, items):
        self._items = items

    def __iter__(self):
        return iter(self._items)


def test_container_iter_bundles_fallback_to_plain_iteration():
    c = ListContainer([Item("X"), object(), Item("Y")])
    bundles = list(c.iter_bundles())
    assert {b.station for b in bundles} == {"X", "Y"}


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
