"""Small branch tests that push near-threshold emtools modules over 80%."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from pycsamt.emtools import csumt, impedance


def test_csumt_schedule_empty_and_frequency_plot():
    assert csumt.frequency_schedule([-2.0, 0.0], 100.0).size == 0

    fig, ax = plt.subplots()
    out = csumt.plot_frequency_schedule(
        [1.0, 20.0, 1000.0],
        300.0,
        f_min=1.0e4,
        f_max=6.0e5,
        ax=ax,
        title="schedule",
    )
    assert out is ax
    assert ax.get_yscale() == "log"
    assert ax.get_title() == "schedule"
    assert len(ax.collections) >= 1
    plt.close(fig)


def test_csumt_unwrap_falls_back_to_original_object():
    obj = object()
    assert csumt._unwrap(obj) is obj


def test_impedance_empty_site_placeholders(monkeypatch):
    monkeypatch.setattr(impedance, "ensure_sites", lambda *a, **k: object())
    monkeypatch.setattr(impedance, "_iter_items", lambda obj: iter(()))

    polar = impedance.plot_phasor_wheel(object())
    assert polar.name == "polar"
    assert polar.texts[0].get_text() == "no sites"

    heat = impedance.plot_offdiag_antisym_residual(object())
    assert heat.texts[0].get_text() == "no data"
    plt.close(polar.figure)
    plt.close(heat.figure)


def test_impedance_z_block_legacy_signature(monkeypatch):
    calls = []

    def old_get_z(ed, *args, **kwargs):
        calls.append(kwargs)
        if kwargs:
            raise TypeError("legacy signature")
        return ("Z", np.zeros((1, 2, 2), complex), np.ones(1))

    monkeypatch.setattr(impedance, "_get_z_block", old_get_z)
    assert impedance._zblk_flex(object(), need_err=True)[0] == "Z"
    assert calls == [{"with_errors": True}, {}]
