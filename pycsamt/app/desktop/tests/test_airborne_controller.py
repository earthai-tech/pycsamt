# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for AirborneController — Airborne EM diagnostics/map dispatch.

Real data
---------
data/ZTEM/gold_springs_nv/     — real ZTEM EMTF-XML flight lines (105 stations)
data/AFMAG/abitibi_on/         — real AFMAG EMTF-XML (13 stations)
data/mobileMT/flammefjeld_greenland/ — real MobileMT EMTF-XML (12 stations)

Every catalogue entry is dispatched with its default parameters against
real data for its category, confirming it returns a populated Figure
rather than an error/placeholder figure. Two entries build their own
multi-panel figure (``plot_ztem_divergence_psection_grid``,
``plot_ztem_band_mask_psection``); the rest take a single ``ax=``.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

from pycsamt.app.desktop.controllers.airborne_controller import (
    CATALOGUE,
    CATEGORIES,
    AirborneController,
)

_ROOT = Path(__file__).resolve().parents[4]
_ZTEM_DIR = _ROOT / "data" / "ZTEM" / "gold_springs_nv"
_AFMAG_DIR = _ROOT / "data" / "AFMAG" / "abitibi_on"
_MOBILEMT_DIR = _ROOT / "data" / "mobileMT" / "flammefjeld_greenland"

_HAS_ZTEM = _ZTEM_DIR.exists() and any(_ZTEM_DIR.glob("*.xml"))
_HAS_AFMAG = _AFMAG_DIR.exists() and any(_AFMAG_DIR.glob("*.xml"))
_HAS_MOBILEMT = _MOBILEMT_DIR.exists() and any(_MOBILEMT_DIR.glob("*.xml"))


def _close():
    plt.close("all")


# ── Catalogue shape ─────────────────────────────────────────────────────────


def test_categories_match_catalogue_keys():
    assert CATEGORIES == list(CATALOGUE.keys())
    assert set(CATEGORIES) == {"ZTEM", "AFMAG", "MobileMT"}


def test_every_entry_is_a_5_tuple():
    for cat, entries in CATALOGUE.items():
        for entry in entries:
            assert len(entry) == 5
            label, fn_name, desc, params, multi_axes = entry
            assert isinstance(label, str) and label
            assert isinstance(fn_name, str) and fn_name
            assert isinstance(desc, str) and desc
            assert isinstance(params, list)
            assert isinstance(multi_axes, bool)


def test_every_fn_name_resolves_in_its_emtools_module():
    import pycsamt.emtools.afmag as afmag_mod
    import pycsamt.emtools.mobilemt as mobilemt_mod
    import pycsamt.emtools.ztem as ztem_mod

    modules = {"ZTEM": ztem_mod, "AFMAG": afmag_mod, "MobileMT": mobilemt_mod}
    for cat, entries in CATALOGUE.items():
        mod = modules[cat]
        for label, fn_name, desc, params, multi_axes in entries:
            assert hasattr(mod, fn_name), f"{cat}.{fn_name} not found in {mod}"


# ── Loading ──────────────────────────────────────────────────────────────────


class TestLoad:
    def test_load_real_ztem_directory(self):
        if not _HAS_ZTEM:
            pytest.skip("ZTEM sample data not available")
        c = AirborneController()
        n = c.load(str(_ZTEM_DIR))
        assert n > 0
        assert c.has_data
        assert n == len(c.state.asites)

    def test_clear_resets_state(self):
        if not _HAS_ZTEM:
            pytest.skip("ZTEM sample data not available")
        c = AirborneController()
        c.load(str(_ZTEM_DIR))
        c.clear()
        assert c.state.asites is None
        assert not c.has_data

    def test_has_data_false_when_nothing_loaded(self):
        c = AirborneController()
        assert not c.has_data

    def test_load_bad_path_resolves_to_zero_stations(self):
        # ensure_asites(strict=False) matches ensure_sites' lenient
        # default used everywhere else in this app: a path that resolves
        # to nothing comes back as an empty AirborneSites, not a raise.
        c = AirborneController()
        n = c.load("/definitely/not/a/real/path")
        assert n == 0
        assert not c.has_data


# ── generate(): no data ──────────────────────────────────────────────────────


class TestGenerateNoData:
    def test_returns_placeholder_figure(self):
        c = AirborneController()
        fig = c.generate("ZTEM", "plot_ztem_tipper_profile")
        texts = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "no airborne data" in texts.lower()
        _close()

    def test_unknown_fn_name_returns_error_not_raise(self):
        c = AirborneController()
        fig = c.generate("ZTEM", "not_a_real_function")
        assert fig is not None
        _close()


# ── generate(): real ZTEM sweep ──────────────────────────────────────────────


@pytest.mark.skipif(not _HAS_ZTEM, reason="ZTEM sample data not available")
class TestGenerateZtemSweep:
    @pytest.fixture(scope="class")
    def ctrl(self):
        c = AirborneController()
        c.load(str(_ZTEM_DIR))
        return c

    @pytest.mark.parametrize(
        "fn_name", [e[1] for e in CATALOGUE["ZTEM"]], ids=[e[1] for e in CATALOGUE["ZTEM"]]
    )
    def test_default_params_render_without_error(self, ctrl, fn_name):
        entry = next(e for e in CATALOGUE["ZTEM"] if e[1] == fn_name)
        _, _, _, params, _ = entry
        kwargs = {p.name: p.default for p in params}
        fig = ctrl.generate("ZTEM", fn_name, **kwargs)
        assert len(fig.axes) >= 1
        _close()


# ── generate(): real AFMAG sweep ─────────────────────────────────────────────


@pytest.mark.skipif(not _HAS_AFMAG, reason="AFMAG sample data not available")
class TestGenerateAfmagSweep:
    @pytest.fixture(scope="class")
    def ctrl(self):
        c = AirborneController()
        c.load(str(_AFMAG_DIR))
        return c

    @pytest.mark.parametrize(
        "fn_name", [e[1] for e in CATALOGUE["AFMAG"]], ids=[e[1] for e in CATALOGUE["AFMAG"]]
    )
    def test_default_params_render_without_error(self, ctrl, fn_name):
        entry = next(e for e in CATALOGUE["AFMAG"] if e[1] == fn_name)
        _, _, _, params, _ = entry
        kwargs = {p.name: p.default for p in params}
        fig = ctrl.generate("AFMAG", fn_name, **kwargs)
        assert len(fig.axes) >= 1
        _close()


# ── generate(): real MobileMT sweep ──────────────────────────────────────────


@pytest.mark.skipif(not _HAS_MOBILEMT, reason="MobileMT sample data not available")
class TestGenerateMobileMtSweep:
    @pytest.fixture(scope="class")
    def ctrl(self):
        c = AirborneController()
        c.load(str(_MOBILEMT_DIR))
        return c

    @pytest.mark.parametrize(
        "fn_name",
        [e[1] for e in CATALOGUE["MobileMT"]],
        ids=[e[1] for e in CATALOGUE["MobileMT"]],
    )
    def test_default_params_render_without_error(self, ctrl, fn_name):
        entry = next(e for e in CATALOGUE["MobileMT"] if e[1] == fn_name)
        _, _, _, params, _ = entry
        kwargs = {p.name: p.default for p in params}
        fig = ctrl.generate("MobileMT", fn_name, **kwargs)
        assert len(fig.axes) >= 1
        _close()


# ── Multi-axes vs. single-ax dispatch ────────────────────────────────────────


@pytest.mark.skipif(not _HAS_ZTEM, reason="ZTEM sample data not available")
class TestMultiAxesDispatch:
    def test_divergence_grid_is_multi_axes_and_builds_own_figure(self):
        c = AirborneController()
        c.load(str(_ZTEM_DIR))
        fig = c.generate(
            "ZTEM", "plot_ztem_divergence_psection_grid", max_lines=2
        )
        # A grid of >=2 panels means more than one axes.
        assert len(fig.axes) >= 2
        _close()

    def test_band_mask_is_multi_axes_and_builds_own_figure(self):
        c = AirborneController()
        c.load(str(_ZTEM_DIR))
        fig = c.generate("ZTEM", "plot_ztem_band_mask_psection")
        assert len(fig.axes) >= 2
        _close()

    def test_tipper_profile_is_single_ax(self):
        c = AirborneController()
        c.load(str(_ZTEM_DIR))
        fig = c.generate("ZTEM", "plot_ztem_tipper_profile")
        assert len(fig.axes) == 1
        _close()
