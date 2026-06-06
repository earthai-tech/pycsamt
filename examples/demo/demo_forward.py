#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
examples/demo/demo_forward.py
==============================

pycsamt.forward — 1-D and 2-D forward modelling gallery
---------------------------------------------------------
Generates one PNG per forward plotting function, covering:

* 1-D MT model depth profiles and soundings
* Multi-model overlays with geology priors
* 2-D finite-difference MT forward models and pseudo-sections
* Lateral response profiles
* Style variants (default, publication, dark)

All figures land in ``examples/demo/figures/``.

Run from the project root::

    python examples/demo/demo_forward.py

Requirements
------------
pycsamt (installed or on PYTHONPATH), matplotlib, numpy, scipy.
"""
from __future__ import annotations

import os
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

warnings.filterwarnings("ignore")

# ── resolve paths ─────────────────────────────────────────────────────────────
_DEMO_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT      = os.path.dirname(os.path.dirname(_DEMO_DIR))
sys.path.insert(0, _ROOT)

OUT_DIR = os.path.join(_DEMO_DIR, "figures", "forward")
os.makedirs(OUT_DIR, exist_ok=True)

# ── forward imports ──────────────────────────────────────────────────────────
from pycsamt.forward import (
    MT1DForward, CSAMT1DForward, TEM1DForward,
    LayeredModel, GEOLOGY_PRIORS,
    Grid2D, MT2DForward,
    Grid3D, MT3DForward,
    plot_response_1d,
    plot_model_1d,
    plot_response_and_model_1d,
    plot_model_2d,
    plot_pseudosection_2d,
    plot_response_profiles,
    plot_model_3d,
    plot_response_map_3d,
    plot_response_section_3d,
    plot_tensor_components_3d,
)
from pycsamt.api.style   import use_style, reset_style, configure_style
from pycsamt.api.control import PYCSAMT_CONTROL, configure_control

DPI        = 150
N_FIGURES  = 0

# ── save helper ───────────────────────────────────────────────────────────────
def save(obj, name: str, title: str = "") -> None:
    """Accept a Figure, Axes, or ndarray of Axes; optionally add suptitle."""
    global N_FIGURES
    if isinstance(obj, plt.Figure):
        fig = obj
    elif isinstance(obj, np.ndarray):
        fig = obj.flat[0].get_figure()
    else:
        fig = obj.get_figure()
    if title:
        fig.suptitle(title, y=1.02, fontsize=10, fontweight="normal")
    path = os.path.join(OUT_DIR, name)
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    N_FIGURES += 1
    print(f"  ✔  {name}")


# ── shared frequency grids ────────────────────────────────────────────────────
FREQS_MT    = np.logspace(-3, 4, 35)      # broad MT/AMT band
FREQS_CSAMT = np.logspace(1, 4, 25)       # CSAMT: 10 Hz – 10 kHz
TIMES_TEM   = np.logspace(-6, -2, 30)     # TEM gate times

# ═══════════════════════════════════════════════════════════════════════════════
print("━" * 64)
print("  pycsamt.forward  —  forward modelling gallery")
print("━" * 64)

# ─────────────────────────────────────────────────────────────────────────────
# Section 1 — 1-D model depth profiles
# ─────────────────────────────────────────────────────────────────────────────
print("\n── Section 1 — 1-D model depth profiles ──")

# Representative layered models
M_SEDIMENTARY = LayeredModel(
    [1_000., 20., 5., 300.],
    [200., 600., 1_500.],
    name="sedimentary",
)
M_CRYSTALLINE = LayeredModel(
    [800., 8_000., 600.],
    [2_000., 15_000.],
    name="crystalline",
)
M_GEOTHERMAL = LayeredModel(
    [500., 8., 250., 3_000.],
    [100., 400., 2_500.],
    name="geothermal",
)
M_HALFSPACE = LayeredModel([100.], [], name="halfspace")
M_CONDUCTIVE = LayeredModel(
    [200., 5., 400., 100.],
    [150., 500., 2_000.],
    name="conductive-layer",
)

# Fig 01 — single model depth profile (sedimentary)
save(
    plot_model_1d(M_SEDIMENTARY, title="Sedimentary model"),
    "fig_fwd01_model_1d_single.png",
    "Fig 01 — 1-D resistivity–depth profile  (sedimentary)",
)

# Fig 02 — five geology-inspired models overlaid
save(
    plot_model_1d(
        [M_SEDIMENTARY, M_CRYSTALLINE, M_GEOTHERMAL, M_CONDUCTIVE, M_HALFSPACE],
        labels=["sedimentary", "crystalline", "geothermal",
                "conductive-layer", "halfspace"],
        figsize=(4.5, 6),
    ),
    "fig_fwd02_models_overlay.png",
    "Fig 02 — 1-D model overlay  (5 geological scenarios)",
)

# Fig 03 — random models from GEOLOGY_PRIORS
rng_models = []
rng_labels = []
for scenario in ("sedimentary", "crystalline", "geothermal", "marine", "permafrost"):
    m = LayeredModel.from_geology(scenario, seed=42)
    m.name = scenario
    rng_models.append(m)
    rng_labels.append(scenario)

save(
    plot_model_1d(rng_models, labels=rng_labels, figsize=(4.5, 6)),
    "fig_fwd03_models_geology_priors.png",
    "Fig 03 — GEOLOGY_PRIORS random models  (one sample each, seed=42)",
)

# ─────────────────────────────────────────────────────────────────────────────
# Section 2 — 1-D MT/CSAMT soundings
# ─────────────────────────────────────────────────────────────────────────────
print("\n── Section 2 — 1-D soundings ──")

# Compute MT responses
R_SED  = MT1DForward(FREQS_MT).run(M_SEDIMENTARY)
R_CRYS = MT1DForward(FREQS_MT).run(M_CRYSTALLINE)
R_GEO  = MT1DForward(FREQS_MT).run(M_GEOTHERMAL)
R_COND = MT1DForward(FREQS_MT).run(M_CONDUCTIVE)
R_HS   = MT1DForward(FREQS_MT).run(M_HALFSPACE)

# Fig 04 — classic 2-panel sounding (sedimentary)
save(
    plot_response_1d(R_SED, title="Sedimentary MT sounding"),
    "fig_fwd04_response_1d_mt_sed.png",
    "Fig 04 — MT1D apparent resistivity & phase  (sedimentary model)",
)

# Fig 05 — geothermal sounding (strong phase suppression)
save(
    plot_response_1d(R_GEO, title="Geothermal MT sounding"),
    "fig_fwd05_response_1d_mt_geo.png",
    "Fig 05 — MT1D apparent resistivity & phase  (geothermal model)",
)

# Fig 06 — 3-panel composite: model + ρ_a + phase  (validate & save)
save(
    plot_response_and_model_1d(
        R_COND, M_CONDUCTIVE,
        title="Conductive-layer model — validate & save",
        figsize=(11, 4.5),
    ),
    "fig_fwd06_validate_1d_conductive.png",
    "Fig 06 — 1-D forward validate view  (conductive-layer model)",
)

# Fig 07 — 3-panel validate for geothermal
save(
    plot_response_and_model_1d(
        R_GEO, M_GEOTHERMAL,
        title="Geothermal model — validate & save",
        figsize=(11, 4.5),
    ),
    "fig_fwd07_validate_1d_geothermal.png",
    "Fig 07 — 1-D forward validate view  (geothermal model)",
)

# Fig 08 — multi-response overlay with MultilineStyle
# Build a multi-response comparison plot manually
fig, axs = plt.subplots(2, 1, figsize=(9, 6), sharex=True,
                         constrained_layout=True)
from pycsamt.api.style import PYCSAMT_STYLE
responses = [R_SED, R_CRYS, R_GEO, R_COND, R_HS]
labels    = ["sedimentary", "crystalline", "geothermal",
             "conductive-layer", "halfspace"]
colors    = PYCSAMT_STYLE.multiline.colors(5)
x         = PYCSAMT_CONTROL.x.transform(FREQS_MT)
for k, (resp, lab) in enumerate(zip(responses, labels)):
    kw = dict(color=colors[k], lw=1.6, alpha=0.88, label=lab)
    axs[0].plot(x, np.log10(resp.rho_a), **kw)
    axs[1].plot(x, resp.phase, **kw)
for ax in axs:
    ax.grid(True, which="both", ls=":", lw=0.4, color="0.75")
    ax.set_axisbelow(True)
axs[0].set_ylabel(r"$\log_{10}\rho_a$  ($\Omega\cdot$m)", fontsize=9)
axs[1].set_ylabel(r"Phase ($^\circ$)", fontsize=9)
axs[1].set_xlabel(PYCSAMT_CONTROL.x.label(), fontsize=9)
axs[0].legend(fontsize=8, ncol=2, framealpha=0.8)
save(fig, "fig_fwd08_response_comparison_mt.png",
     "Fig 08 — MT1D sounding comparison  (5 geological scenarios)")

# Fig 09 — CSAMT sounding (uses higher frequencies)
R_SED_CSAMT = CSAMT1DForward(FREQS_CSAMT).run(M_SEDIMENTARY)
save(
    plot_response_and_model_1d(
        R_SED_CSAMT, M_SEDIMENTARY,
        title="CSAMT sounding — sedimentary model",
        figsize=(11, 4.5),
    ),
    "fig_fwd09_validate_1d_csamt.png",
    "Fig 09 — CSAMT1D validate view  (sedimentary model, 10 Hz–10 kHz)",
)

# Fig 10 — TEM step-off response (manual, no rho_a)
R_TEM = TEM1DForward(TIMES_TEM, loop_radius=50.0).run(M_CONDUCTIVE)
fig_tem, ax_tem = plt.subplots(figsize=(7, 4), constrained_layout=True)
ax_tem.loglog(R_TEM.times, np.abs(R_TEM.dBz_dt),
              color=PYCSAMT_STYLE.mt.te.color, lw=1.6, marker="o",
              ms=3.5, mfc="white", mew=1.0, alpha=0.9, label="TEM step-off")
ax_tem.set_xlabel("Time  (s)", fontsize=9)
ax_tem.set_ylabel(r"$|d\mathbf{B}_z/dt|$  (T/s)", fontsize=9)
ax_tem.set_title("TEM1D step-off dBz/dt response", fontsize=9, pad=6)
ax_tem.legend(fontsize=8)
ax_tem.grid(True, which="both", ls=":", lw=0.4, color="0.75")
ax_tem.set_axisbelow(True)
save(fig_tem, "fig_fwd10_response_1d_tem.png",
     "Fig 10 — TEM1D step-off response  (conductive-layer model, 50 m loop)")

# ─────────────────────────────────────────────────────────────────────────────
# Section 3 — 2-D MT forward models and pseudo-sections
# ─────────────────────────────────────────────────────────────────────────────
print("\n── Section 3 — 2-D MT forward modelling ──")

FREQS_2D = np.logspace(-2, 2, 18)   # 0.01 – 100 Hz (18 freq points)

# ── Build two 2-D models ──
# Model A: conductive body (simulates a fault zone / graphite schist)
GRID_FAULT = Grid2D.with_anomaly(
    bg_rho       = 500.0,
    anomaly_rho  = 3.0,
    anomaly_bounds = (2_500., 5_500., 300., 1_800.),
    nx=45, nz=32,
    x_max=9_000., z_max=5_000.,
    n_stations=12, n_pad=8,
    name="fault-zone conductor",
)

# Model B: resistive body (simulates a salt dome / igneous intrusion)
GRID_SALT = Grid2D.with_anomaly(
    bg_rho       = 50.0,
    anomaly_rho  = 5_000.0,
    anomaly_bounds = (3_000., 6_000., 500., 3_000.),
    nx=45, nz=32,
    x_max=9_000., z_max=5_000.,
    n_stations=12, n_pad=8,
    name="resistive intrusion",
)

print("  Running 2-D FD solver — fault-zone conductor …")
RESP_FAULT = MT2DForward(FREQS_2D, GRID_FAULT, verbose=False).run()

print("  Running 2-D FD solver — resistive intrusion …")
RESP_SALT  = MT2DForward(FREQS_2D, GRID_SALT,  verbose=False).run()

# Fig 11 — 2-D model view: fault-zone conductor
save(
    plot_model_2d(GRID_FAULT, figsize=(11, 4)),
    "fig_fwd11_model_2d_fault.png",
    "Fig 11 — 2-D resistivity model  (fault-zone conductor, 3 Ω·m on 500 Ω·m)",
)

# Fig 12 — 2-D model view: resistive intrusion
save(
    plot_model_2d(GRID_SALT, figsize=(11, 4)),
    "fig_fwd12_model_2d_salt.png",
    "Fig 12 — 2-D resistivity model  (resistive intrusion, 5 000 Ω·m on 50 Ω·m)",
)

# Fig 13 — TE pseudo-section: fault conductor
save(
    plot_pseudosection_2d(RESP_FAULT, mode="te", quantity="rho_a",
                          figsize=(11, 5)),
    "fig_fwd13_psection_te_rho_fault.png",
    "Fig 13 — TE apparent-resistivity pseudo-section  (fault-zone conductor)",
)

# Fig 14 — TM pseudo-section: fault conductor
save(
    plot_pseudosection_2d(RESP_FAULT, mode="tm", quantity="rho_a",
                          figsize=(11, 5)),
    "fig_fwd14_psection_tm_rho_fault.png",
    "Fig 14 — TM apparent-resistivity pseudo-section  (fault-zone conductor)",
)

# Fig 15 — TE phase pseudo-section: fault conductor
save(
    plot_pseudosection_2d(RESP_FAULT, mode="te", quantity="phase",
                          figsize=(11, 5)),
    "fig_fwd15_psection_te_phase_fault.png",
    "Fig 15 — TE phase pseudo-section  (fault-zone conductor)",
)

# Fig 16 — TE pseudo-section with contour overlay: resistive intrusion
save(
    plot_pseudosection_2d(RESP_SALT, mode="te", quantity="rho_a",
                          n_contours=8, figsize=(11, 5)),
    "fig_fwd16_psection_te_rho_salt_contour.png",
    "Fig 16 — TE apparent-resistivity pseudo-section  (resistive intrusion + contours)",
)

# Fig 17 — lateral profiles: TE at 5 selected periods
save(
    plot_response_profiles(RESP_FAULT, mode="te", quantity="rho_a",
                           n_freqs_shown=5, figsize=(9, 4)),
    "fig_fwd17_profiles_te_rho_fault.png",
    "Fig 17 — Lateral TE ρ_a profiles at 5 periods  (fault-zone conductor)",
)

# Fig 18 — lateral profiles: TM at same periods
save(
    plot_response_profiles(RESP_FAULT, mode="tm", quantity="rho_a",
                           n_freqs_shown=5, figsize=(9, 4)),
    "fig_fwd18_profiles_tm_rho_fault.png",
    "Fig 18 — Lateral TM ρ_a profiles at 5 periods  (fault-zone conductor)",
)

# Fig 19 — combined 2-row: model + TE + TM pseudo-section (fault)
fig19, axs19 = plt.subplots(3, 1, figsize=(12, 13),
                              constrained_layout=True)
plot_model_2d(GRID_FAULT,    ax=axs19[0], show_stations=True, title="Resistivity model")
plot_pseudosection_2d(RESP_FAULT, ax=axs19[1], mode="te", quantity="rho_a",
                      show_stations=True,
                      title=r"TE — $\log_{10}\rho_a$")
plot_pseudosection_2d(RESP_FAULT, ax=axs19[2], mode="tm", quantity="rho_a",
                      show_stations=True,
                      title=r"TM — $\log_{10}\rho_a$")
save(fig19, "fig_fwd19_validate_2d_fault.png",
     "Fig 19 — 2-D forward validate view  (fault-zone conductor)")

# Fig 20 — combined 3-row for resistive intrusion
fig20, axs20 = plt.subplots(3, 1, figsize=(12, 13),
                              constrained_layout=True)
plot_model_2d(GRID_SALT,    ax=axs20[0], show_stations=True, title="Resistivity model")
plot_pseudosection_2d(RESP_SALT, ax=axs20[1], mode="te", quantity="rho_a",
                      show_stations=True,
                      title=r"TE — $\log_{10}\rho_a$")
plot_pseudosection_2d(RESP_SALT, ax=axs20[2], mode="tm", quantity="rho_a",
                      show_stations=True,
                      title=r"TM — $\log_{10}\rho_a$")
save(fig20, "fig_fwd20_validate_2d_salt.png",
     "Fig 20 — 2-D forward validate view  (resistive intrusion)")

# ─────────────────────────────────────────────────────────────────────────────
# Section 4 — Style variants (same fault model, three looks)
# ─────────────────────────────────────────────────────────────────────────────
print("\n── Section 4 — Style variants ──")

# Fig 21 — publication style (grayscale / high-contrast)
use_style("publication")
fig21, axs21 = plt.subplots(1, 3, figsize=(13, 4.5), constrained_layout=True)
plot_model_1d(M_SEDIMENTARY, ax=axs21[0],
              title="Depth profile")
plot_response_1d(R_SED, axes=axs21[1:3])
save(fig21, "fig_fwd21_style_publication.png",
     "Fig 21 — Publication style  (sedimentary model, grayscale)")
reset_style()

# Fig 22 — dark style
use_style("dark")
fig22, axs22 = plt.subplots(1, 3, figsize=(13, 4.5), constrained_layout=True,
                             facecolor="#1a1a2e")
for ax in axs22:
    ax.set_facecolor("#1a1a2e")
plot_model_1d(M_CONDUCTIVE, ax=axs22[0],
              title="Depth profile")
plot_response_1d(R_COND, axes=axs22[1:3])
save(fig22, "fig_fwd22_style_dark.png",
     "Fig 22 — Dark style  (conductive-layer model)")
reset_style()

# Fig 23 — period axis (instead of log10-period)
configure_control(x__view="period")
save(
    plot_response_and_model_1d(R_GEO, M_GEOTHERMAL,
                               title="Period axis (linear scale)",
                               figsize=(11, 4.5)),
    "fig_fwd23_axis_period.png",
    "Fig 23 — PYCSAMT_CONTROL: period axis  (geothermal model)",
)
configure_control(x__view="log10_period")   # restore

# ─────────────────────────────────────────────────────────────────────────────
# Section 5 — 3-D quasi-3D forward gallery
# ─────────────────────────────────────────────────────────────────────────────
print("\n── Section 5 — 3-D forward gallery ──")

FREQS_3D = np.logspace(0, 3, 12)   # 1 Hz – 1 kHz

# ── Build two 3-D models ─────────────────────────────────────────────────────
# Halfspace: validation baseline
G3_HALF = Grid3D.halfspace(
    rho=100.0, nx=8, ny=8, nz=6,
    x_max=4_000.0, y_max=4_000.0, z_max=2_000.0,
    n_pad=3, nx_stations=4, ny_stations=4,
    name="halfspace-3D",
)

# 3-D conductive block embedded in a resistive background
G3_BLOCK = Grid3D.block_anomaly(
    bg_rho=500.0, anomaly_rho=5.0,
    bounds=(1_200.0, 2_800.0, 1_200.0, 2_800.0, 300.0, 1_200.0),
    nx=8, ny=8, nz=6,
    x_max=4_000.0, y_max=4_000.0, z_max=2_000.0,
    n_pad=3, nx_stations=4, ny_stations=4,
    name="block-anomaly-3D",
)

# ── Fig 24 — 3-D halfspace model slices ──────────────────────────────────────
save(
    plot_model_3d(G3_HALF,
                  title="3-D halfspace model  (100 Ω·m)"),
    "fig_fwd24_model_3d_halfspace.png",
    "Fig 24 — 3-D model slices  (uniform halfspace)",
)

# ── Fig 25 — 3-D block-anomaly model slices ───────────────────────────────────
save(
    plot_model_3d(G3_BLOCK,
                  title="3-D block-anomaly model  (500/5 Ω·m)"),
    "fig_fwd25_model_3d_block.png",
    "Fig 25 — 3-D model slices  (conductive block anomaly)",
)

# ── Run quasi-3D forward ──────────────────────────────────────────────────────
print("  running quasi-3D forward (halfspace)…")
R3_HALF  = MT3DForward(FREQS_3D, G3_HALF,  verbose=False).run()
print("  running quasi-3D forward (block anomaly)…")
R3_BLOCK = MT3DForward(FREQS_3D, G3_BLOCK, verbose=False).run()

# ── Fig 26 — map-view ρ_a  (halfspace, Z_xy) ─────────────────────────────────
save(
    plot_response_map_3d(R3_HALF, freq_idx=0, component="xy",
                         title=r"Halfspace — map  $\rho_a$  [Z_XY]  (T = 1 s)"),
    "fig_fwd26_map_3d_halfspace.png",
    r"Fig 26 — quasi-3D map view  $\rho_a$  (halfspace)",
)

# ── Fig 27 — map-view ρ_a  (block anomaly, Z_xy) ─────────────────────────────
save(
    plot_response_map_3d(R3_BLOCK, freq_idx=4, component="xy",
                         title=r"Block anomaly — map  $\rho_a$  [Z_XY]"),
    "fig_fwd27_map_3d_block_rho.png",
    r"Fig 27 — quasi-3D map view  $\rho_a$  (conductive block)",
)

# ── Fig 28 — map-view phase  (block anomaly, Z_xy) ────────────────────────────
save(
    plot_response_map_3d(R3_BLOCK, freq_idx=4, component="xy",
                         quantity="phase",
                         title=r"Block anomaly — map phase  [Z_XY]"),
    "fig_fwd28_map_3d_block_phase.png",
    r"Fig 28 — quasi-3D map view  phase  (conductive block)",
)

# ── Fig 29 — pseudo-section  (halfspace, Z_xy) ────────────────────────────────
save(
    plot_response_section_3d(R3_HALF, component="xy",
                              title=r"Halfspace — 3-D pseudo-section  $\rho_a$  [Z_XY]"),
    "fig_fwd29_psection_3d_halfspace.png",
    r"Fig 29 — quasi-3D pseudo-section  $\rho_a$  (halfspace)",
)

# ── Fig 30 — pseudo-section  (block anomaly, Z_xy) ────────────────────────────
save(
    plot_response_section_3d(R3_BLOCK, component="xy", n_contours=4,
                              title=r"Block anomaly — 3-D pseudo-section  $\rho_a$  [Z_XY]"),
    "fig_fwd30_psection_3d_block.png",
    r"Fig 30 — quasi-3D pseudo-section  $\rho_a$  (conductive block)",
)

# ── Fig 31 — pseudo-section phase  (block anomaly, Z_yx) ──────────────────────
save(
    plot_response_section_3d(R3_BLOCK, component="yx", quantity="phase",
                              title=r"Block anomaly — pseudo-section phase  [Z_YX]"),
    "fig_fwd31_psection_3d_block_phase.png",
    r"Fig 31 — quasi-3D pseudo-section  phase  Z_YX  (block anomaly)",
)

# ── Fig 32 — full tensor 2×2 map panel  (block anomaly) ──────────────────────
save(
    plot_tensor_components_3d(R3_BLOCK, freq_idx=4,
                               title=r"Block anomaly — full impedance tensor  $\rho_a$"),
    "fig_fwd32_tensor_3d_block.png",
    r"Fig 32 — quasi-3D full tensor components  $\rho_a$  (block anomaly)",
)

# ── Fig 33 — full tensor phase panel  (block anomaly) ────────────────────────
save(
    plot_tensor_components_3d(R3_BLOCK, freq_idx=4, quantity="phase",
                               title="Block anomaly — full impedance tensor  phase"),
    "fig_fwd33_tensor_3d_block_phase.png",
    "Fig 33 — quasi-3D full tensor components  phase  (block anomaly)",
)

# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "━" * 64)
print(f"  {N_FIGURES} figures saved →  {OUT_DIR}")
print("━" * 64)
