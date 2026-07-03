#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
examples/demo/demo_emtools.py
==============================

pycsamt.emtools — visualisation gallery
-----------------------------------------
Generates one PNG per emtools plotting function using the bundled WILLY_DATA
AMT survey (5 profile lines, 128 stations).  All figures land in
``examples/demo/figures/``.

Run from the project root::

    python examples/demo/demo_emtools.py

Requirements
------------
pycsamt (installed or on PYTHONPATH), matplotlib, numpy, scipy.
"""
from __future__ import annotations

import os
import sys
import warnings

import matplotlib
matplotlib.use("Agg")           # headless — no display required
import matplotlib.pyplot as plt
import matplotlib.cm as mcm
import numpy as np

warnings.filterwarnings("ignore")

# ── resolve paths relative to this file ────────────────────────────────────
_DEMO_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT      = os.path.dirname(os.path.dirname(_DEMO_DIR))
sys.path.insert(0, _ROOT)

DATA_ROOT   = os.path.join(_ROOT, "data", "AMT", "WILLY_DATA")
OUT_DIR     = os.path.join(_DEMO_DIR, "figures", "emtools")
TIPPER_DIR  = os.path.join(_ROOT, "data", "AMT", "TIPPER")
TIPPER_OUT  = os.path.join(_DEMO_DIR, "figures", "tipper")
os.makedirs(OUT_DIR,    exist_ok=True)
os.makedirs(TIPPER_OUT, exist_ok=True)

# ── emtools imports ─────────────────────────────────────────────────────────
from pycsamt.emtools._core       import ensure_sites, _iter_items, _name
from pycsamt.emtools.inspect     import plot_coverage
from pycsamt.emtools.qc          import (
    plot_confidence_band_summary,
    plot_confidence_profile,
    plot_frequency_confidence_psection,
    plot_station_confidence_dashboard,
)
from pycsamt.emtools.frequency   import (
    edit_frequencies_by_confidence,
    plot_coverage_quality_heatmap,
    plot_apparent_depth_psection,
    plot_frequency_edit_decisions,
    plot_frequency_edit_summary,
)
from pycsamt.emtools.remove_noise import (
    apply_emap_filter,
    confidence_gated_emap_filter,
    emap_filter_report,
    plot_emap_filter_profile,
    plot_emap_filter_psection,
)
from pycsamt.emtools.impedance   import plot_phasor_wheel, plot_determinant_track
from pycsamt.emtools.tensor      import (
    plot_phase_tensor_psection,
    plot_phase_tensor_summary,
    plot_phase_tensor_rose,
    plot_phase_tensor_map,
    plot_phase_tensor_strip,
    plot_phase_tensor_strip_grid,
    plot_theta_vs_period,
    plot_dimensionality_psection,
    plot_skew_ellipt_density,
    plot_theta_rose_grid,
)
from pycsamt.emtools.skew        import (
    plot_skew_traffic_psection,
    plot_skew_percentile_ribbon,
)
from pycsamt.emtools.dimensionality import plot_dim_confidence_grid
from pycsamt.emtools.anisotropy     import plot_anisotropy
from pycsamt.emtools.strike         import (
    plot_strike_rose,
    plot_strike_ribbon,
    plot_strike_profile,
    plot_strike_mapsticks,
    plot_strike_analysis,
)
from pycsamt.emtools.tf             import (
    plot_tipper_hodograms,
    plot_tipper_polar,
    plot_induction_arrows,
    plot_induction_map,
    plot_induction_section,
    plot_induction_convention,
    plot_induction_rose,
)
from pycsamt.emtools.plot           import plot_raw_sites_1d
from pycsamt.api.control            import PYCSAMT_CONTROL

DPI = 150
N_FIGURES = 0

# ── helper ──────────────────────────────────────────────────────────────────
def save(obj, name: str, title: str = "") -> None:
    """Accept a Figure *or* an Axes, then save it to the gallery."""
    global N_FIGURES
    fig = obj if isinstance(obj, plt.Figure) else obj.get_figure()
    path = os.path.join(OUT_DIR, name)
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    N_FIGURES += 1
    print(f"  ✔  {name}")


# ── load data ────────────────────────────────────────────────────────────────
print("━" * 62)
print("  pycsamt.emtools  —  demo figure gallery")
print("━" * 62)
print("\nLoading EDI data …")

L22   = ensure_sites(os.path.join(DATA_ROOT, "L22PLT"))       # 25 stations, AMT
L22_PATH = os.path.join(DATA_ROOT, "L22PLT")

# All 5 lines (128 stations) with explicit line groups for multi-line plots
S_all = ensure_sites(DATA_ROOT, recursive=True)
all_names = [_name(ed, i) for i, ed in enumerate(_iter_items(S_all))]
groups_by_line = {
    "L18": sorted(n for n in all_names if n.startswith("18-")),
    "L22": sorted(n for n in all_names if n.startswith("22-")),
    "L26": sorted(n for n in all_names if n.startswith("26-")),
    "L30": sorted(n for n in all_names if n.startswith("30-")),
    "L34": sorted(n for n in all_names if n.startswith("34-")),
}
print(f"  L22PLT : {len(groups_by_line['L22'])} stations")
print(f"  ALL    : {len(all_names)} stations across {len(groups_by_line)} lines\n")

# ────────────────────────────────────────────────────────────────────────────
# Section 1 — Frequency coverage & data quality
# ────────────────────────────────────────────────────────────────────────────
print("── Section 1 — Coverage & data quality ──")

save(
    plot_coverage(L22, figsize=(12, 3.8)),
    "fig01_coverage.png",
    "Fig 01 — Frequency coverage by station  (L22PLT)",
)

save(
    plot_coverage_quality_heatmap(L22, figsize=(12, 4.2)),
    "fig02_coverage_heatmap.png",
    "Fig 02 — Coverage quality heatmap  (L22PLT)",
)

save(
    plot_apparent_depth_psection(L22, figsize=(12, 4.8)),
    "fig03_apparent_depth.png",
    "Fig 03 — Apparent-depth pseudo-section  (L22PLT)",
)

save(
    plot_confidence_profile(
        L22,
        method="composite",
        figsize=(12, 4.2),
        ci_hi=0.90,
        ci_lo=0.50,
    ),
    "fig03b_station_confidence_composite.png",
    "Fig 03b — Composite station-confidence profile  (L22PLT)",
)

save(
    plot_frequency_confidence_psection(
        L22,
        method="composite",
        figsize=(12, 4.8),
        ci_hi=0.90,
        ci_lo=0.50,
    ),
    "fig03c_frequency_confidence_section.png",
    "Fig 03c — Frequency-confidence pseudo-section  (L22PLT)",
)

save(
    plot_confidence_band_summary(
        L22,
        method="composite",
        figsize=(8.5, 4.0),
        ci_hi=0.90,
        ci_lo=0.50,
    ),
    "fig03d_confidence_period_summary.png",
    "Fig 03d — Period-band confidence summary  (L22PLT)",
)

save(
    plot_station_confidence_dashboard(
        L22,
        station="22-14BF",
        method="composite",
        figsize=(10.5, 6.0),
        ci_hi=0.90,
        ci_lo=0.50,
    ),
    "fig03e_station_confidence_dashboard.png",
    "Fig 03e — Frequency-confidence dashboard  (22-14BF)",
)

# Frequency confidence editing demo.  Separate loads keep the unedited
# baseline independent from the edited source objects.
L22_edit_before = ensure_sites(L22_PATH)
L22_edit_source = ensure_sites(L22_PATH)
frequency_edit = edit_frequencies_by_confidence(
    L22_edit_source,
    before_sites=L22_edit_before,
    mode="recover",
    method="composite",
    ci_hi=0.90,
    ci_lo=0.50,
    reject="drop",
)
frequency_edit.report.to_csv(
    os.path.join(OUT_DIR, "table03f_frequency_edit_report.csv"),
    index=False,
)
frequency_edit.decisions.to_csv(
    os.path.join(OUT_DIR, "table03g_frequency_edit_decisions.csv"),
    index=False,
)
print("  ✔  table03f_frequency_edit_report.csv")
print("  ✔  table03g_frequency_edit_decisions.csv")

save(
    plot_frequency_edit_summary(
        L22_edit_before,
        frequency_edit.sites,
        method="composite",
        ci_hi=0.90,
        ci_lo=0.50,
        figsize=(10, 4.0),
    ),
    "fig03f_frequency_edit_summary.png",
    "Fig 03f — Frequency confidence edit summary  (L22PLT)",
)

save(
    plot_frequency_edit_decisions(
        L22_edit_before,
        frequency_edit.sites,
        method="composite",
        ci_hi=0.90,
        ci_lo=0.50,
        figsize=(12, 4.8),
    ),
    "fig03g_frequency_edit_decisions.png",
    "Fig 03g — Frequency confidence edit decisions  (L22PLT)",
)

L22_emap_before = ensure_sites(L22_PATH)
L22_emap_source = ensure_sites(L22_PATH)
L22_flma = apply_emap_filter(
    L22_emap_source,
    method="flma",
    window=5,
    component="xy",
)
emap_report = emap_filter_report(
    L22_emap_before,
    L22_flma,
    component="xy",
    period_s=0.01,
)
emap_report.to_csv(
    os.path.join(OUT_DIR, "table03h_emap_flma_report.csv"),
    index=False,
)
print("  ✔  table03h_emap_flma_report.csv")

save(
    plot_emap_filter_profile(
        L22_emap_before,
        L22_flma,
        method="flma",
        component="xy",
        period_s=0.01,
        figsize=(10.5, 4.0),
    ),
    "fig03h_emap_flma_profile.png",
    "Fig 03h — EMAP FLMA station profile  (L22PLT)",
)

save(
    plot_emap_filter_psection(
        L22_emap_before,
        L22_flma,
        method="flma",
        component="xy",
        figsize=(11.5, 8.2),
    ),
    "fig03i_emap_flma_psection.png",
    "Fig 03i — EMAP FLMA pseudo-section comparison  (L22PLT)",
)

L22_gated_before = ensure_sites(L22_PATH)
L22_gated_source = ensure_sites(L22_PATH)
gated_emap = confidence_gated_emap_filter(
    L22_gated_source,
    before_sites=L22_gated_before,
    method="flma",
    confidence_method="composite",
    component="xy",
    window=5,
    ci_hi=0.90,
    ci_lo=0.50,
)
gated_emap.report.to_csv(
    os.path.join(OUT_DIR, "table03j_emap_gated_report.csv"),
    index=False,
)
gated_emap.decisions.to_csv(
    os.path.join(OUT_DIR, "table03k_emap_gated_decisions.csv"),
    index=False,
)
print("  ✔  table03j_emap_gated_report.csv")
print("  ✔  table03k_emap_gated_decisions.csv")

save(
    plot_emap_filter_psection(
        L22_gated_before,
        gated_emap.sites,
        method="confidence-gated FLMA",
        component="xy",
        figsize=(11.5, 8.2),
    ),
    "fig03j_emap_gated_flma_psection.png",
    "Fig 03j — Confidence-gated EMAP FLMA pseudo-section  (L22PLT)",
)

# ────────────────────────────────────────────────────────────────────────────
# Section 2 — Apparent resistivity & phase (direct from Z tensor)
# ────────────────────────────────────────────────────────────────────────────
print("\n── Section 2 — Apparent resistivity & phase ──")

raw_panel_stations = ["22-1BF", "22-14BF", "22-24BF"]

with PYCSAMT_CONTROL.context(
    phase__range=(-180.0, 180.0),
    rho__view="log10",
    x__view="log10_period",
):
    save(
        plot_raw_sites_1d(
            L22,
            stations=raw_panel_stations,
            components=("xx", "xy", "yx", "yy"),
            raw=True,
            ncols_groups=3,
            figsize_scale=(4.2, 3.2),
            show_component_legend=True,
        ),
        "fig04_raw_1d_station_panels.png",
        "Fig 04 — Raw 1-D rho/phase panels  (L22PLT)",
    )

with PYCSAMT_CONTROL.context(
    phase__range=(-90.0, 90.0),
    rho__view="log10",
    x__view="period",
):
    save(
        plot_raw_sites_1d(
            L22,
            stations=raw_panel_stations,
            components=("xy", "yx"),
            raw=False,
            ncols_groups=3,
            figsize_scale=(3.2, 3.0),
            show_component_legend=True,
        ),
        "fig04b_component_style_1d_panels.png",
        "Fig 04b — Component-coloured 1-D rho/phase panels  (L22PLT)",
    )

# Build ρa / φ overlay directly from Site.rho / Site.phase / Site.freq
# (all 25 L22PLT stations, XY component)
sites_list = list(_iter_items(L22))
cmap_st    = mcm.get_cmap("tab20", len(sites_list))

fig, (ax_r, ax_p) = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
for k, ed in enumerate(sites_list):
    nm  = _name(ed, k)
    per = 1.0 / ed.freq
    col = cmap_st(k)
    ax_r.loglog(per, ed.rho[:, 0, 1], ".", ms=2.5, color=col, label=nm)
    ax_p.semilogx(per, ed.phase[:, 0, 1], ".", ms=2.5, color=col)

ax_r.set_ylabel(r"$\rho_a^{XY}$  (Ω·m)")
ax_r.set_title("Fig 04 — Apparent resistivity & phase (XY component, all stations)",
               fontsize=9, pad=6)
ax_r.legend(ncol=5, fontsize=5.5, loc="upper right",
            markerscale=3, framealpha=0.7)
ax_p.set_ylabel("Phase XY  (°)")
ax_p.set_xlabel("Period  (s)")
ax_p.set_ylim(0, 90)
fig.tight_layout()
save(fig, "fig04_rhoa_phi.png")

# ────────────────────────────────────────────────────────────────────────────
# Section 3 — Tipper & impedance diagnostics
# ────────────────────────────────────────────────────────────────────────────
print("\n── Section 3 — Impedance diagnostics & station map ──")

# Note: WILLY_DATA EDI files contain no Hz/tipper component, so
# plot_tipper_hodograms is replaced with plot_strike_mapsticks which
# places per-station strike sticks on a lon/lat map.
save(
    plot_strike_mapsticks(L22, method="consensus", figsize=(6, 6)),
    "fig05_strike_mapsticks.png",
    "Fig 05 — Geoelectric strike map-sticks  (L22PLT)",
)

save(
    plot_phasor_wheel(L22, figsize=(10, 5)),
    "fig06_phasor_wheel.png",
    "Fig 06 — Impedance phasor wheel  (L22PLT)",
)

save(
    plot_determinant_track(L22, figsize=(10, 4.5)),
    "fig07_determinant_track.png",
    "Fig 07 — Determinant impedance track  (L22PLT)",
)

# ────────────────────────────────────────────────────────────────────────────
# Section 4 — Phase tensor analysis
# ────────────────────────────────────────────────────────────────────────────
print("\n── Section 4 — Phase tensor analysis ──")

save(
    plot_phase_tensor_psection(L22, figsize=(12, 5.5)),
    "fig08_pt_psection.png",
    "Fig 08 — Phase-tensor ellipse pseudo-section  (L22PLT)",
)

save(
    plot_phase_tensor_summary(L22, figsize=(13, 10)),
    "fig09_pt_summary.png",
    "Fig 09 — Phase-tensor summary: φmax | φmin | β  (L22PLT)",
)

save(
    plot_phase_tensor_rose(
        L22,
        bins             = 36,
        bar_style        = "gradient",
        cmap             = "plasma",
        outer_ring_lw    = 2.5,
        outer_ring_color = "0.12",
        n_rings          = 4,
        ring_label_angle = 22.5,
        ring_label_fontsize = 7.5,
        ring_label_fmt   = "{:.0f}",
        spoke_every      = 45.0,
        compass_labels   = "NESW",
        compass_fontsize = 9.0,
        compass_fontweight = "bold",
        show_mean        = True,
        mean_color       = "crimson",
        mean_lw          = 2.2,
        show_secondary   = True,
        secondary_ls     = "--",
        show_annotation  = True,
        show_n           = True,
        annotation_fontsize = 8.5,
        figsize          = (6, 6),
    ),
    "fig10_pt_rose.png",
    "Fig 10 — Phase-tensor θ rose  (L22PLT, all periods)",
)

save(
    plot_theta_vs_period(L22, figsize=(9, 4.2)),
    "fig11_theta_vs_period.png",
    "Fig 11 — Strike angle θ vs period  (L22PLT, all stations)",
)

# 11b — full L18PLT pseudo-section (28 stations): the field-data case that
# motivated the min_aspect floor + reserved-legend-strip fixes (real EDI
# data has many near-degenerate φ_min cells that used to render invisible).
L18 = ensure_sites(os.path.join(DATA_ROOT, "L18PLT"))
save(
    plot_phase_tensor_psection(L18, figsize=(14, 5.5)),
    "fig11b_pt_psection_L18.png",
    "Fig 11b — Phase-tensor ellipse pseudo-section  (L18PLT, 28 stations)",
)

# 11c — geographic phase-tensor map at a representative short period
save(
    plot_phase_tensor_map(L22, period=0.01, figsize=(8, 7)),
    "fig11c_pt_map.png",
    "Fig 11c — Phase-tensor map at T=0.01 s  (L22PLT)",
)

# 11d — single-station ellipse strip (period on x, one row)
save(
    plot_phase_tensor_strip(L22, station="22-14BF", figsize=(7, 1.6)),
    "fig11d_pt_strip_single.png",
    "Fig 11d — Phase-tensor ellipse strip  (22-14BF)",
)

# 11e — multi-profile ellipse-strip grid: 4 representative (evenly spaced)
# stations per line, across all 5 WILLY_DATA profiles, one shared colorbar
def _pick_representative(names, k=4):
    names = sorted(names)
    if len(names) <= k:
        return names
    idx = sorted(set(np.linspace(0, len(names) - 1, k).round().astype(int)))
    return [names[i] for i in idx]

strip_profiles = {
    f"Profile {ln}": _pick_representative(names)
    for ln, names in groups_by_line.items()
}
save(
    plot_phase_tensor_strip_grid(
        S_all,
        profiles=strip_profiles,
        c_by="skew", cmap="RdBu_r",
        panel_size=(3.4, 1.05),
        suptitle="Fig 11e — Phase-tensor ellipse strips by profile  (WILLY_DATA)",
    ),
    "fig11e_pt_strip_grid.png",
)

# ────────────────────────────────────────────────────────────────────────────
# Section 5 — Dimensionality classification
# ────────────────────────────────────────────────────────────────────────────
print("\n── Section 5 — Dimensionality ──")

save(
    plot_dimensionality_psection(L22, figsize=(12, 4.2)),
    "fig12_dim_psection.png",
    "Fig 12 — Dimensionality pseudo-section: 1D / 2D / 3D  (L22PLT)",
)

save(
    plot_dim_confidence_grid(L22, figsize=(12, 4.5)),
    "fig13_dim_confidence.png",
    "Fig 13 — Dimensionality confidence grid  (L22PLT)",
)

# ────────────────────────────────────────────────────────────────────────────
# Section 6 — Skewness analysis
# ────────────────────────────────────────────────────────────────────────────
print("\n── Section 6 — Skewness analysis ──")

save(
    plot_skew_traffic_psection(L22, figsize=(12, 4.2)),
    "fig14_skew_traffic.png",
    "Fig 14 — Skewness traffic-light pseudo-section  (L22PLT)",
)

save(
    plot_skew_percentile_ribbon(L22, figsize=(9, 3.8)),
    "fig15_skew_ribbon.png",
    "Fig 15 — Skewness percentile ribbon  (L22PLT)",
)

save(
    plot_skew_ellipt_density(L22, figsize=(7, 5.5)),
    "fig16_skew_density.png",
    "Fig 16 — Skewness β vs ellipticity λ density  (L22PLT)",
)

# ────────────────────────────────────────────────────────────────────────────
# Section 7 — Anisotropy
# ────────────────────────────────────────────────────────────────────────────
print("\n── Section 7 — Anisotropy ──")

save(
    plot_anisotropy(L22, figsize=(11, 4.5)),
    "fig17_anisotropy.png",
    "Fig 17 — Anisotropy ratio pseudo-section  (L22PLT)",
)

# ────────────────────────────────────────────────────────────────────────────
# Section 8 — Geoelectric strike analysis
# ────────────────────────────────────────────────────────────────────────────
print("\n── Section 8 — Geoelectric strike analysis ──")

# 18  Strike rose — all 5 profile lines side-by-side
save(
    plot_strike_rose(
        S_all,
        groups           = groups_by_line,
        method           = "consensus",
        bins             = 24,
        bar_style        = "gradient",
        cmap             = "YlOrRd",
        outer_ring_lw    = 2.5,
        outer_ring_color = "0.12",
        n_rings          = 3,
        spoke_every      = 45.0,
        compass_labels   = "NESW",
        compass_fontsize = 7.5,
        mean_color       = "crimson",
        mean_lw          = 2.2,
        show_secondary   = True,
        secondary_ls     = "--",
        annotation_fontsize = 8.0,
        show_n_stations  = True,
        subplot_size     = 3.0,
        n_cols           = 5,
        suptitle         = "Fig 18 — Geoelectric strike rose diagrams  "
                           "(WILLY_DATA — 5 profile lines, 128 stations)",
        suptitle_fontsize = 10,
    ),
    "fig18_strike_rose_multiline.png",
)

# 19  Strike rose — frequency-band decomposition on L22PLT
save(
    plot_strike_rose(
        L22,
        bar_style        = "bands",
        freq_bands       = [(1e-4, 1e-2), (1e-2, 1e0)],
        band_labels      = ["Short period  (0.1–10 ms)", "Long period  (10 ms–1 s)"],
        bins             = 24,
        outer_ring_lw    = 2.5,
        outer_ring_color = "0.12",
        spoke_every      = 45.0,
        compass_labels   = "NESW",
        annotation_fontsize = 8.5,
        show_n_stations  = True,
        subplot_size     = 4.2,
        suptitle         = "Fig 19 — Strike rose: short vs long period  (L22PLT)",
        suptitle_fontsize = 10,
    ),
    "fig19_strike_rose_bands.png",
)

# 20  Strike colour ribbon — strike vs log-period for all stations
save(
    plot_strike_ribbon(L22, method="sweep", figsize=(12, 3.8)),
    "fig20_strike_ribbon.png",
    "Fig 20 — Strike colour ribbon  (L22PLT, sweep method)",
)

# 21  Strike profile — median ± IQR along the station profile
save(
    plot_strike_profile(L22, method="consensus", figsize=(10, 3.8)),
    "fig21_strike_profile.png",
    "Fig 21 — Geoelectric strike profile  (L22PLT)",
)

# 22  Phase-tensor θ rose grid — one rose per frequency decade
save(
    plot_theta_rose_grid(L22),
    "fig22_theta_rose_grid.png",
    "Fig 22 — Phase-tensor θ rose grid by frequency decade  (L22PLT)",
)

# ────────────────────────────────────────────────────────────────────────────
# Section 9 — Tipper / induction-arrow analysis  (HBH03_IMP.edi)
# ────────────────────────────────────────────────────────────────────────────
print("\n── Section 9 — Tipper & induction-arrow analysis ──")

TIPPER_EDI = os.path.join(TIPPER_DIR, "HBH03_IMP.edi")
S_tip = ensure_sites(TIPPER_EDI)
N_TIP = 0

def save_tip(obj, name: str, title: str = "") -> None:
    global N_TIP
    fig = obj if isinstance(obj, plt.Figure) else obj.get_figure()
    path = os.path.join(TIPPER_OUT, name)
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    N_TIP += 1
    print(f"  ✔  tipper/{name}")


# T01 — Tipper hodograms: real and imaginary vectors per period band
save_tip(
    plot_tipper_hodograms(
        S_tip,
        n_bands  = 4,
        normalize = False,
        figsize  = (10, 8),
    ),
    "fig_t01_tipper_hodograms.png",
    "Fig T01 — Tipper hodograms by period band  (HBH03)",
)

# T02 — Tipper polar: magnitude vs azimuth for each frequency
save_tip(
    plot_tipper_polar(
        S_tip,
        component = "real",
        cmap      = "plasma",
        figsize   = (6, 6),
    ),
    "fig_t02_tipper_polar_real.png",
    "Fig T02 — Tipper polar magnitude (real)  (HBH03)",
)

save_tip(
    plot_tipper_polar(
        S_tip,
        component = "imag",
        cmap      = "viridis",
        figsize   = (6, 6),
    ),
    "fig_t02b_tipper_polar_imag.png",
    "Fig T02b — Tipper polar magnitude (imaginary)  (HBH03)",
)

# T03 — Induction arrows section at representative periods
save_tip(
    plot_induction_arrows(
        S_tip,
        periods  = [1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1],
        scale    = 2.0,
        figsize  = (11, 4.5),
    ),
    "fig_t03_induction_arrows.png",
    "Fig T03 — Induction arrows section  (HBH03)",
)

# T04 — Induction map at a single representative period
save_tip(
    plot_induction_map(
        S_tip,
        period  = 0.01,
        scale   = 0.5,
        figsize = (7, 5),
    ),
    "fig_t04_induction_map.png",
    "Fig T04 — Induction arrow map  (HBH03, T=0.01 s)",
)

# T05 — Tipper magnitude section vs period
save_tip(
    plot_induction_section(
        S_tip,
        component = "real",
        n_periods = 20,
        figsize   = (9, 4.5),
    ),
    "fig_t05_induction_section.png",
    "Fig T05 — Tipper section (real component)  (HBH03)",
)

# T06 — Induction convention diagram (returns 2×2 axes array)
_conv_axs = plot_induction_convention(
    S_tip,
    period  = 0.01,
    figsize = (10, 8),
)
save_tip(
    _conv_axs.flat[0].get_figure(),
    "fig_t06_induction_convention.png",
    "Fig T06 — Induction arrow convention  (HBH03, T=0.01 s)",
)

# T07 — Rose diagram of induction arrow azimuths (real component)
save_tip(
    plot_induction_rose(
        S_tip,
        component = "real",
        nbins     = 36,
        figsize   = (5.5, 5.5),
    ),
    "fig_t07_induction_rose_real.png",
    "Fig T07 — Induction rose (real)  (HBH03, all periods)",
)

# T08 — Rose diagram of induction arrow azimuths (imaginary component)
save_tip(
    plot_induction_rose(
        S_tip,
        component = "imag",
        nbins     = 36,
        figsize   = (5.5, 5.5),
    ),
    "fig_t08_induction_rose_imag.png",
    "Fig T08 — Induction rose (imaginary)  (HBH03, all periods)",
)

# T09 — Combined Strike(Z) / PT Azimuth / Tipper Strike analysis
save_tip(
    plot_strike_analysis(
        S_tip,
        style   = "pycsamt",
        method  = "sweep",
        bins    = 36,
        suptitle = "Strike Analysis — HBH03  (Z | PT Azimuth | Tipper)",
        subplot_size = 4.0,
    ),
    "fig_t09_strike_analysis.png",
    "Fig T09 — Strike analysis rose (Z / PT / Tipper)  (HBH03)",
)

# T10 — Strike analysis with short-period band only
save_tip(
    plot_strike_analysis(
        S_tip,
        style   = "pycsamt",
        method  = "sweep",
        band    = (1e-4, 1e-2),
        bins    = 36,
        suptitle = "Strike Analysis — HBH03  (T = 0.1–10 ms)",
        subplot_size = 4.0,
    ),
    "fig_t10_strike_analysis_shortband.png",
    "Fig T10 — Strike analysis: short-period band  (HBH03)",
)

print(f"\n  {N_TIP} tipper figures saved → {TIPPER_OUT}")

# ────────────────────────────────────────────────────────────────────────────
print("\n" + "━" * 62)
print(f"  {N_FIGURES} emtools figures saved →  {OUT_DIR}")
print(f"  {N_TIP} tipper  figures saved →  {TIPPER_OUT}")
print("━" * 62)
