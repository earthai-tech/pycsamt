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

DATA_ROOT  = os.path.join(_ROOT, "data", "AMT", "WILLY_DATA")
OUT_DIR    = os.path.join(_DEMO_DIR, "figures")
os.makedirs(OUT_DIR, exist_ok=True)

# ── emtools imports ─────────────────────────────────────────────────────────
from pycsamt.emtools._core       import ensure_sites, _iter_items, _name
from pycsamt.emtools.inspect     import plot_coverage
from pycsamt.emtools.frequency   import (
    plot_coverage_quality_heatmap,
    plot_apparent_depth_psection,
)
from pycsamt.emtools.impedance   import plot_phasor_wheel, plot_determinant_track
from pycsamt.emtools.tensor      import (
    plot_phase_tensor_psection,
    plot_phase_tensor_summary,
    plot_phase_tensor_rose,
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
)

DPI = 150

# ── helper ──────────────────────────────────────────────────────────────────
def save(obj, name: str, title: str = "") -> None:
    """Accept a Figure *or* an Axes, optionally add a suptitle, then save."""
    fig = obj if isinstance(obj, plt.Figure) else obj.get_figure()
    if title:
        fig.suptitle(title, y=1.02, fontsize=10, fontweight="normal")
    path = os.path.join(OUT_DIR, name)
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✔  {name}")


# ── load data ────────────────────────────────────────────────────────────────
print("━" * 62)
print("  pycsamt.emtools  —  demo figure gallery")
print("━" * 62)
print("\nLoading EDI data …")

L22   = ensure_sites(os.path.join(DATA_ROOT, "L22PLT"))       # 25 stations, AMT

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

# ────────────────────────────────────────────────────────────────────────────
# Section 2 — Apparent resistivity & phase (direct from Z tensor)
# ────────────────────────────────────────────────────────────────────────────
print("\n── Section 2 — Apparent resistivity & phase ──")

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
print("\n" + "━" * 62)
print(f"  22 figures saved →  {OUT_DIR}")
print("━" * 62)
