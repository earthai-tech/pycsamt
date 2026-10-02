# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.desktop.controllers.correction_views -- the Qt-free
Before/After comparison renderer used by CorrectionWindow.

Synthetic stations with analytically known responses are used so the
numbers can be checked exactly: a pure static shift by factor k must give
ρ_after/ρ_before == k at every period and Δφ == 0.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from pycsamt.app.desktop.controllers.correction_controller import _LIGHT
from pycsamt.app.desktop.controllers.correction_views import (
    ALL_STATIONS,
    PlotUnavailable,
    extract_responses,
    figure_blank_reason,
    render_curves,
    render_section,
)


class _Z:
    def __init__(self, z, freq):
        self.z = z
        self.freq = freq


class _Station:
    def __init__(self, name, z, freq):
        self.station = name
        self.Z = _Z(z, freq)


def _survey(n=4, nf=12, scale=None, rotate_phase_deg=0.0, nan_yx=False):
    """n stations of a 100 ohm-m half-space-like response (phi = 45 deg)."""
    freq = np.logspace(3, -1, nf)
    T = 1.0 / freq
    out = []
    for i in range(n):
        rho = 100.0 * (1 + i)  # distinct per station
        amp = np.sqrt(rho / (0.2 * T))  # rho = 0.2 T |Z|^2
        k = 1.0 if scale is None else scale[i]
        ph = np.radians(45.0 + rotate_phase_deg)
        z = np.zeros((nf, 2, 2), dtype=complex)
        z[:, 0, 1] = np.sqrt(k) * amp * np.exp(1j * ph)
        z[:, 1, 0] = -np.sqrt(k) * amp * np.exp(1j * ph)  # 3rd quadrant
        if nan_yx:
            z[:, 1, 0] = np.nan
        out.append(_Station(f"S{i:02d}", z, freq))
    return out


@pytest.fixture
def fig():
    f = plt.figure(figsize=(9, 6), layout="constrained")
    yield f
    plt.close(f)


class TestExtract:
    def test_rho_and_phase_follow_field_unit_convention(self):
        r = extract_responses(_survey(n=1))["S00"]
        np.testing.assert_allclose(r.rho["xy"], 100.0, rtol=1e-9)
        np.testing.assert_allclose(r.phi["xy"], 45.0, atol=1e-9)
        # YX phase is folded from the 3rd quadrant into the 1st
        np.testing.assert_allclose(r.phi["yx"], 45.0, atol=1e-9)
        assert np.all(np.diff(r.T) > 0)  # sorted by period

    def test_duplicate_station_names_kept_apart(self):
        s = _survey(n=2)
        s[1].station = s[0].station
        assert len(extract_responses(s)) == 2

    def test_none_is_empty(self):
        assert extract_responses(None) == {}


class TestCurves:
    @pytest.mark.parametrize("mode", ["Before / After", "Overlay"])
    def test_modes_draw_real_data(self, fig, mode):
        render_curves(fig, _survey(), _survey(scale=[2, 2, 2, 2]), mode=mode,
                      theme=_LIGHT)
        assert figure_blank_reason(fig) is None

    def test_before_after_shares_y_scale_across_states(self, fig):
        render_curves(fig, _survey(), _survey(scale=[10] * 4),
                      mode="Before / After", quantities=("rho",), theme=_LIGHT)
        ax_b, ax_a = [ax for ax in fig.axes if ax.get_label() != "<colorbar>"][:2]
        assert ax_b.get_ylim() == ax_a.get_ylim()

    def test_static_shift_diff_is_exact_ratio_and_zero_phase(self, fig):
        k = [0.5, 2.0, 4.0, 1.0]
        render_curves(fig, _survey(), _survey(scale=k), mode="Diff",
                      quantities=("rho", "phi"), components=("xy",),
                      theme=_LIGHT)
        ax_rho, ax_phi = fig.axes[:2]
        ratios = sorted(
            float(np.median(ln.get_ydata()))
            for ln in ax_rho.get_lines() if len(ln.get_ydata()) > 2
        )
        np.testing.assert_allclose(ratios, sorted(k), rtol=1e-9)
        for ln in ax_phi.get_lines():
            if len(ln.get_ydata()) > 2:
                np.testing.assert_allclose(ln.get_ydata(), 0.0, atol=1e-9)
        assert any("unchanged" in t.get_text() for t in ax_phi.texts)

    def test_phase_change_is_visible_in_diff(self, fig):
        render_curves(fig, _survey(), _survey(rotate_phase_deg=10.0),
                      mode="Diff", quantities=("phi",), theme=_LIGHT)
        ys = [ln.get_ydata() for ln in fig.axes[0].get_lines()
              if len(ln.get_ydata()) > 2]
        np.testing.assert_allclose(np.concatenate(ys), 10.0, atol=1e-9)

    def test_diff_of_identical_data_is_explained(self, fig):
        with pytest.raises(PlotUnavailable, match="identical"):
            render_curves(fig, _survey(), _survey(), mode="Diff", theme=_LIGHT)

    def test_single_station_selection(self, fig):
        render_curves(fig, _survey(), _survey(scale=[3] * 4), mode="Overlay",
                      station="S02", quantities=("rho",),
                      components=("xy", "yx"), theme=_LIGHT)
        assert len(fig.axes[0].get_lines()) == 4  # 2 comps x 2 states

    def test_unknown_station_is_explained(self, fig):
        with pytest.raises(PlotUnavailable, match="not in this dataset"):
            render_curves(fig, _survey(), _survey(), station="NOPE",
                          theme=_LIGHT)

    def test_all_nan_component_is_explained(self, fig):
        s = _survey(nan_yx=True)
        with pytest.raises(PlotUnavailable, match="No valid YX"):
            render_curves(fig, s, s, components=("yx",), theme=_LIGHT)

    def test_all_stations_mode_has_station_key(self, fig):
        render_curves(fig, _survey(), _survey(scale=[2] * 4), theme=_LIGHT)
        assert any(ax.get_ylabel() == "Station" for ax in fig.axes)

    def test_no_data_is_explained(self, fig):
        with pytest.raises(PlotUnavailable, match="No data"):
            render_curves(fig, None, None, theme=_LIGHT)


class TestSection:
    @pytest.mark.parametrize("mode", ["Before / After", "Overlay", "Diff"])
    def test_modes_draw_real_data(self, fig, mode):
        render_section(fig, _survey(), _survey(scale=[2, 1, 1, 0.5]),
                       mode=mode, quantities=("rho", "phi"),
                       affected_stations=["S01"], station="S02", theme=_LIGHT)
        assert figure_blank_reason(fig) is None

    def test_before_after_uses_one_colour_scale(self, fig):
        render_section(fig, _survey(), _survey(scale=[10] * 4),
                       mode="Before / After", quantities=("rho",), theme=_LIGHT)
        from matplotlib.collections import QuadMesh

        meshes = [c for ax in fig.axes if ax.get_label() != "<colorbar>"
                  for c in ax.collections if isinstance(c, QuadMesh)]
        clims = {tuple(np.round(m.get_clim(), 9)) for m in meshes}
        assert len(meshes) == 2 and len(clims) == 1

    def test_diff_is_centred_on_zero(self, fig):
        render_section(fig, _survey(), _survey(scale=[2, 1, 1, 0.5]),
                       mode="Diff", quantities=("rho",), theme=_LIGHT)
        (mesh,) = [c for c in fig.axes[0].collections if hasattr(c, "get_clim")]
        lo, hi = mesh.get_clim()
        assert lo == pytest.approx(-hi)

    def test_one_station_is_explained(self, fig):
        s = _survey(n=1)
        with pytest.raises(PlotUnavailable, match="at least two stations"):
            render_section(fig, s, s, theme=_LIGHT)


class TestBlankReason:
    def test_text_only_figure_is_blank(self, fig):
        fig.add_subplot(111).text(0.5, 0.5, "No valid Z data")
        assert figure_blank_reason(fig) == "No valid Z data"

    def test_figure_with_line_is_not_blank(self, fig):
        fig.add_subplot(111).plot([1, 2], [1, 2])
        assert figure_blank_reason(fig) is None

    def test_empty_figure_is_blank(self, fig):
        assert figure_blank_reason(fig)


def test_all_stations_sentinel_means_no_selection(fig):
    render_curves(fig, _survey(), _survey(scale=[2] * 4), station=ALL_STATIONS,
                  quantities=("rho",), mode="Overlay", theme=_LIGHT)
    assert len(fig.axes[0].get_lines()) == 8  # 4 stations x 2 states
