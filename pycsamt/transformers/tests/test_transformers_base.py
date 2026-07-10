from __future__ import annotations

import numpy as np

from pycsamt.core.base import TFBundle
from pycsamt.core.config import config_context
from pycsamt.transformers import _base as t


def test_public_api_all():
    assert set(t.__all__) == {"TransformerMixin"}


class DummyTransformer(t.TransformerMixin):
    def __init__(self):
        self.calls: list[str] = []
        self._bundle_in: TFBundle | None = None

    def extract(self, source):
        self.calls.append("extract")
        self._bundle_in = source
        return source

    def emit_edi(self, bundle):
        self.calls.append("emit")
        # return a small dict so we can inspect
        return {"bundle": bundle, "ok": True}

    def post_emit(self, edi_obj, source, bundle):
        self.calls.append("post")
        return edi_obj

    def compute_res_from_z(self, b: TFBundle) -> TFBundle:
        self.calls.append("res_from_z")
        if b.z is not None and b.freq is not None:
            n = len(b.freq)
            b.rho = np.full((n, 2, 2), 100.0)
            b.phase = np.full((n, 2, 2), 45.0)
        return b

    def compute_z_from_res(self, b: TFBundle) -> TFBundle:
        self.calls.append("z_from_res")
        if b.rho is not None and b.phase is not None:
            n = len(b.rho)
            b.z = np.ones((n, 2, 2), complex)
        return b


def test_finalize_orders_and_dedups_and_station():
    tr = DummyTransformer()
    # craft frequencies unsorted with one near-duplicate
    freq = np.array([1.0, 10.0, 10.0 * (1 + 5e-10), 5.0])
    z = np.arange(freq.size * 4).reshape(freq.size, 2, 2)
    b = TFBundle(freq=freq, z=z, station=None, station_id=7)

    # default is desc sorting; near-duplicate gets dropped
    out = tr._finalize(b)
    assert out.station == "S007"
    assert np.all(np.diff(out.freq) <= 0)
    # dedup removed one element → size becomes 3
    assert out.freq.size == 3
    assert out.z.shape[0] == 3

    # now enforce ascending order
    with config_context(freq_order="asc"):
        out2 = tr._finalize(TFBundle(freq=freq, z=z))
        assert np.all(np.diff(out2.freq) >= 0)


def test_fill_missing_hooks_res_from_z_and_back():
    tr = DummyTransformer()
    f = np.array([8.0, 4.0, 2.0])
    z = np.ones((3, 2, 2), complex)

    # compute_res_from_z path
    with config_context(
        compute_res_from_z=True,
        compute_z_from_res=True,
    ):
        b1 = TFBundle(freq=f, z=z)
        o1 = tr._finalize(b1)
        assert "res_from_z" in tr.calls
        assert o1.rho is not None and o1.phase is not None

    # compute_z_from_res path
    tr.calls.clear()
    rp_rho = np.full((3, 2, 2), 50.0)
    rp_phi = np.full((3, 2, 2), 30.0)
    with config_context(
        compute_res_from_z=False,
        compute_z_from_res=True,
    ):
        b2 = TFBundle(freq=f, rho=rp_rho, phase=rp_phi)
        o2 = tr._finalize(b2)
        assert "z_from_res" in tr.calls
        assert o2.z is not None


def test_transform_calls_in_order():
    tr = DummyTransformer()
    f = np.array([3.0, 2.0, 1.0])
    b = TFBundle(freq=f, z=np.zeros((3, 2, 2)))

    out = tr.transform(b, name="X", station_id=1)
    # assert tr.calls[:3] == ["extract", "emit", "post"]
    i_ext = tr.calls.index("extract")
    i_emit = tr.calls.index("emit")
    i_post = tr.calls.index("post")
    assert i_ext < i_emit < i_post

    assert out["ok"] is True
    assert out["bundle"].station == "X"


def test_finalize_with_no_freq_is_robust():
    tr = DummyTransformer()
    b = TFBundle(z=np.ones((2, 2, 2)))
    out = tr._finalize(b, name="A")
    assert out.station == "A"
    # no freq → ordering/dedup are no-ops
    assert out.freq is None
