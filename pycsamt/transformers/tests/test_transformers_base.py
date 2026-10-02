from __future__ import annotations

import numpy as np
import pytest

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


def test_default_post_emit_and_compute_hooks_are_noop():
    mixin = t.TransformerMixin()
    sentinel = object()
    assert mixin.post_emit(sentinel, None, None) is sentinel

    b = TFBundle(freq=[1.0])
    assert mixin.compute_res_from_z(b) is b
    assert mixin.compute_z_from_res(b) is b


def test_order_and_dedup_freq_numpy_branch_touches_all_optional_fields():
    freq = np.array([1.0, 3.0, 2.0])
    n = freq.size
    z = np.arange(n * 4).reshape(n, 2, 2).astype(complex)
    z_err = np.ones((n, 2, 2))
    tipper = np.ones((n, 1, 2), complex)
    tipper_err = np.ones((n, 1, 2))
    rho = np.arange(n * 4).reshape(n, 2, 2).astype(float)
    phase = np.arange(n * 4).reshape(n, 2, 2).astype(float)
    b = TFBundle(
        freq=freq,
        z=z,
        z_err=z_err,
        tipper=tipper,
        tipper_err=tipper_err,
        rho=rho,
        phase=phase,
    )
    mixin = t.TransformerMixin()
    ordered = mixin._order_freq(b)
    # default freq_order is desc
    assert np.all(np.diff(ordered.freq) <= 0)
    assert ordered.freq.tolist() == [3.0, 2.0, 1.0]
    assert ordered.z.shape[0] == n
    assert ordered.z_err.shape[0] == n
    assert ordered.tipper.shape[0] == n
    assert ordered.tipper_err.shape[0] == n
    assert ordered.rho.shape[0] == n
    assert ordered.phase.shape[0] == n

    deduped = mixin._dedup_freq(ordered)
    # no near-duplicates present, nothing dropped, all optional arrays kept
    assert deduped.freq.size == n
    assert deduped.z_err.shape[0] == n
    assert deduped.tipper.shape[0] == n
    assert deduped.tipper_err.shape[0] == n
    assert deduped.rho.shape[0] == n
    assert deduped.phase.shape[0] == n


def test_order_freq_list_mode_fallback_without_numpy(monkeypatch):
    monkeypatch.setattr(t, "np", None)
    b = TFBundle(
        freq=[1.0, 3.0, 2.0],
        z=[10, 30, 20],
        z_err=[0.1, 0.3, 0.2],
        tipper=[100, 300, 200],
        tipper_err=[1, 3, 2],
        rho=[1000, 3000, 2000],
        phase=[45, 46, 47],
    )
    out = t.TransformerMixin()._order_freq(b)
    assert out.freq == [3.0, 2.0, 1.0]
    assert out.z == [30, 20, 10]
    assert out.z_err == [0.3, 0.2, 0.1]
    assert out.tipper == [300, 200, 100]
    assert out.tipper_err == [3, 2, 1]
    assert out.rho == [3000, 2000, 1000]
    assert out.phase == [46, 47, 45]


def test_dedup_freq_list_mode_fallback_without_numpy(monkeypatch):
    monkeypatch.setattr(t, "np", None)
    # 0.1 and 0.1 + 5e-10 are only distinguishable if the "denom" floor is
    # 1.0 (as in the numpy branch); this is a regression check for a bug
    # where the sub-1.0 relative-tolerance denominator diverged from the
    # numpy code path (see fix in _dedup_freq's list-mode fallback).
    b = TFBundle(
        freq=[0.1, 0.1 + 5e-10, 5.0],
        z=[1, 2, 3],
        z_err=[0.1, 0.2, 0.3],
        tipper=[1, 2, 3],
        tipper_err=[1, 2, 3],
        rho=[1, 2, 3],
        phase=[1, 2, 3],
    )
    out = t.TransformerMixin()._dedup_freq(b)
    assert out.freq == [0.1, 5.0]
    assert out.z == [1, 3]
    assert out.z_err == [0.1, 0.3]
    assert out.tipper == [1, 3]
    assert out.tipper_err == [1, 3]
    assert out.rho == [1, 3]
    assert out.phase == [1, 3]


def test_ensure_head_info_definemeas_mtsect_helpers():
    from pycsamt.seg.edi import EDIFile

    mixin = t.TransformerMixin()
    ed = EDIFile(verbose=0)

    head1 = mixin._ensure_head(ed, station="S01", empty=1e32)
    assert head1.dataid == "S01"
    assert head1.stdvers == "SEG 1.0"
    assert head1.progvers == "PYCSAMT"
    assert head1.empty == 1e32
    assert head1.lat == 0.0 and head1.long == 0.0 and head1.elev == 0.0

    # second call finds the existing section and does not recreate it
    head2 = mixin._ensure_head(ed, station="OTHER", empty=1e32)
    assert head2 is head1
    assert head2.dataid == "S01"

    info1 = mixin._ensure_info(ed, survey_id="SURV1")
    assert ed.get_section("info") is info1
    info2 = mixin._ensure_info(ed, survey_id="OTHER")
    assert info2 is info1

    head1.lat, head1.long, head1.elev = 10.0, 20.0, 300.0
    dm1 = mixin._ensure_definemeas(ed, units="M", reftype="CART")
    assert dm1.units == "M" and dm1.reftype == "CART"
    assert dm1.reflat == pytest.approx(10.0)
    assert dm1.reflong == pytest.approx(20.0)
    assert dm1.refelev == pytest.approx(300.0)
    dm2 = mixin._ensure_definemeas(ed)
    assert dm2 is dm1

    mt1 = mixin._ensure_mtsect(ed, sectid="S01", nfreq=5)
    assert mt1.sectid == "S01" and mt1.nfreq == 5
    mt2 = mixin._ensure_mtsect(ed, sectid="S02", nfreq=8)
    assert mt2 is mt1
    assert mt1.sectid == "S02" and mt1.nfreq == 8


def test_ensure_definemeas_without_existing_head_section():
    from pycsamt.seg.edi import EDIFile

    ed = EDIFile(verbose=0)
    dm = t.TransformerMixin()._ensure_definemeas(ed)
    assert dm.units == "M"
    assert dm.reftype == "CART"
