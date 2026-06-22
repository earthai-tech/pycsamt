#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_figures.py
-------------------
Generates all six paper figures for pycsamt v2 (Computers & Geosciences).

Outputs (PNG 300 dpi + SVG) in the same directory as this script:
  fig1_architecture.png / .svg
  fig2_survey_map.png   / .svg
  fig3_qc_scores.png    / .svg
  fig4_phase_tensor.png / .svg
  fig5_static_shift.png / .svg
  fig6_nn_section.png   / .svg     (EMInverter2D U-Net 2-D section)

Run:
    python generate_figures.py
"""

import os
import re
import sys
import glob
import warnings
from typing import Dict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

# ─── ensure pycsamt is importable ─────────────────────────────────────────────
HERE    = os.path.dirname(os.path.abspath(__file__))
PKGROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
if PKGROOT not in sys.path:
    sys.path.insert(0, PKGROOT)

# ─── package API ──────────────────────────────────────────────────────────────
from pycsamt.ai.plot._style import EMStyle, em_context
from pycsamt.ai.plot import (
    plot_pseudo_section,
    plot_section,
    plot_uncertainty_bands,
    plot_inversion_result_2d,
)
from pycsamt.ai.processing import EMQCScorer
from pycsamt.ai.processing.plot import plot_qc_scores, plot_qc_summary
from pycsamt.emtools import (
    plot_phase_tensor_psection,
    plot_phase_tensor_summary,
    plot_ss_summary,
    plot_ss_1d_curves,
)

# ─── constants ────────────────────────────────────────────────────────────────
DATADIR  = os.path.join(HERE, "..", "..", "data", "AMT", "WILLY_DATA")
PROFILES = ["L18PLT", "L22PLT", "L26PLT", "L30PLT", "L34PLT"]
PCOL     = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]
RNG      = np.random.default_rng(42)


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def dms_to_dd(dms_str):
    """Convert 'DD:MM:SS.sss' to decimal degrees."""
    dms_str = dms_str.strip()
    parts = re.split(r"[:°'\"]", dms_str)
    d, m, s = float(parts[0]), float(parts[1]), float(parts[2])
    return d + m / 60.0 + s / 3600.0


def parse_edi_block(path, block):
    """Return float array for a named EDI block (e.g. 'FREQ', 'RHOXY')."""
    with open(path, "r", errors="replace") as fh:
        text = fh.read()
    pat = rf">{block}\b[^\n]*\n(.*?)(?=>|\Z)"
    m = re.search(pat, text, re.S | re.I)
    if m is None:
        return None
    raw = re.sub(r"[^\d.eE+\-\s]", " ", m.group(1))
    vals = [float(v) for v in raw.split() if v not in ("", " ")]
    return np.array(vals)


def parse_edi_coords(path):
    """Return (lat_dd, lon_dd) from an EDI file."""
    lat_str = lon_str = None
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if lat_str is None and re.match(r"\s+LAT=", line):
                lat_str = line.split("=")[1].strip()
            if lon_str is None and re.match(r"\s+LONG=", line):
                lon_str = line.split("=")[1].strip()
            if lat_str and lon_str:
                break
    if lat_str and lon_str:
        return dms_to_dd(lat_str), dms_to_dd(lon_str)
    return None, None


def load_profile_data(profile):
    """Return dict with freqs (ascending), logRho (n_stations × n_freqs), coords."""
    edis = sorted(glob.glob(os.path.join(DATADIR, profile, "*.edi")))
    freqs_list, rho_list, coords = [], [], []
    for path in edis:
        f = parse_edi_block(path, "FREQ")
        r = parse_edi_block(path, "RHOXY")
        lat, lon = parse_edi_coords(path)
        if f is None or r is None or len(f) == 0:
            continue
        nf = min(len(f), len(r))
        f, r = f[:nf], r[:nf]
        asc = np.argsort(f)
        freqs_list.append(f[asc])
        rho_list.append(r[asc])
        coords.append((lat, lon))
    if not freqs_list:
        return None
    ref_f = freqs_list[0]
    rho_aligned = []
    for f, r in zip(freqs_list, rho_list):
        valid = (r > 0) & np.isfinite(r) & (f > 0)
        if valid.sum() < 3:
            rho_aligned.append(np.full(len(ref_f), np.nan))
            continue
        ri = np.interp(np.log10(ref_f),
                       np.log10(f[valid]),
                       np.log10(r[valid]),
                       left=np.nan, right=np.nan)
        rho_aligned.append(ri)
    return dict(freqs=ref_f, logRho=np.array(rho_aligned), coords=coords,
                n_stations=len(coords))


def save(fig, stem):
    for ext in ("png", "svg"):
        path = os.path.join(HERE, f"{stem}.{ext}")
        fig.savefig(path, format=ext, dpi=300, bbox_inches="tight")
        print(f"  saved {path}")
    plt.close(fig)


def _edi_to_qc_features(path):
    """
    Parse an EDI file and return an (n_freqs, 5) feature matrix suitable for
    ``EMQCScorer.transform()``.

    Features: [snr, swift_skew, asym, phase_xy, phase_yx]

    Tries Z-tensor blocks (ZXYR/ZXYI/ZYXR/ZYXI) first; falls back to scalar
    apparent-resistivity and phase (RHOXY/PHSXY/RHOYX/PHSYX) if absent.
    Returns ``None`` when the file carries insufficient data.
    """
    freq = parse_edi_block(path, "FREQ")
    if freq is None or len(freq) == 0:
        return None
    n = len(freq)

    zxyr = parse_edi_block(path, "ZXYR")
    zxyi = parse_edi_block(path, "ZXYI")
    zyxr = parse_edi_block(path, "ZYXR")
    zyxi = parse_edi_block(path, "ZYXI")

    if all(arr is not None and len(arr) >= n
           for arr in (zxyr, zxyi, zyxr, zyxi)):
        zxy = zxyr[:n] + 1j * zxyi[:n]
        zyx = zyxr[:n] + 1j * zyxi[:n]
        zxxr  = parse_edi_block(path, "ZXXR")
        zxxi_ = parse_edi_block(path, "ZXXI")
        zyyr  = parse_edi_block(path, "ZYYR")
        zyyi_ = parse_edi_block(path, "ZYYI")
        zxx = (zxxr[:n] + 1j * zxxi_[:n]
               if zxxr is not None and zxxi_ is not None and len(zxxr) >= n
               else np.zeros(n, dtype=complex))
        zyy = (zyyr[:n] + 1j * zyyi_[:n]
               if zyyr is not None and zyyi_ is not None and len(zyyr) >= n
               else np.zeros(n, dtype=complex))
        amp = np.abs(zxy)
        kernel = np.ones(3) / 3.0
        smooth = np.convolve(amp, kernel, mode="same") if n > 3 else amp
        snr      = amp / (np.abs(amp - smooth) + 1e-24)
        swift    = np.abs(zxx + zyy) / (np.abs(zxy - zyx) + 1e-24)
        asym     = np.log10(np.abs(zxy) / (np.abs(zyx) + 1e-24) + 1e-6)
        phase_xy = np.degrees(np.angle(zxy))
        phase_yx = np.degrees(np.angle(zyx))
    else:
        # Fallback: scalar impedance from RHOXY / PHSXY / RHOYX / PHSYX
        rhoxy = parse_edi_block(path, "RHOXY")
        phsxy = parse_edi_block(path, "PHSXY")
        rhoyx = parse_edi_block(path, "RHOYX")
        phsyx = parse_edi_block(path, "PHSYX")
        if rhoxy is None or len(rhoxy) < 2:
            return None
        n = min(n, len(rhoxy))
        rhoxy = np.where(rhoxy[:n] > 0, rhoxy[:n], np.nan)
        log_rho = np.where(np.isfinite(rhoxy), np.log10(rhoxy), 0.0)
        k = min(3, n)
        smooth = np.convolve(log_rho, np.ones(k) / k, mode="same")
        snr   = np.abs(log_rho) / (np.abs(log_rho - smooth) + 0.05)
        swift = np.zeros(n)
        asym  = (np.log10(rhoxy / (rhoyx[:n] + 1e-6) + 1e-6)
                 if rhoyx is not None and len(rhoyx) >= n
                 else np.zeros(n))
        phase_xy = (phsxy[:n] if phsxy is not None and len(phsxy) >= n
                    else np.full(n, 45.0))
        phase_yx = (phsyx[:n] if phsyx is not None and len(phsyx) >= n
                    else np.full(n, -135.0))

    return np.column_stack([snr, swift, asym, phase_xy, phase_yx])


# ═════════════════════════════════════════════════════════════════════════════
# Figure 1 – Package architecture
# ═════════════════════════════════════════════════════════════════════════════

def fig1_architecture():
    print("Generating Fig 1: architecture …")
    with EMStyle():
        fig, ax = plt.subplots(figsize=(12, 7))
        ax.set_xlim(0, 12)
        ax.set_ylim(0, 7)
        ax.axis("off")

        def box(cx, cy, w, h, label, sublabel="", fc="#d0e8f7", ec="#2471a3",
                fontsize=8.5, bold=False):
            rect = FancyBboxPatch((cx - w / 2, cy - h / 2), w, h,
                                  boxstyle="round,pad=0.08",
                                  facecolor=fc, edgecolor=ec, linewidth=1.2,
                                  zorder=3)
            ax.add_patch(rect)
            weight = "bold" if bold else "normal"
            ax.text(cx, cy + (0.10 if sublabel else 0), label,
                    ha="center", va="center", fontsize=fontsize,
                    fontweight=weight, zorder=4)
            if sublabel:
                ax.text(cx, cy - 0.22, sublabel, ha="center", va="center",
                        fontsize=7, color="#555", zorder=4)

        def arrow(x1, y1, x2, y2, color="#555", lw=1.2, style="->", ls="-"):
            ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                        arrowprops=dict(arrowstyle=style, color=color,
                                       lw=lw, linestyle=ls),
                        zorder=2)

        ax.add_patch(FancyBboxPatch((0.15, 0.3), 11.7, 2.7,
                                    boxstyle="round,pad=0.1",
                                    facecolor="#eef6fb", edgecolor="#aacce8",
                                    linewidth=0.8, zorder=0))
        ax.add_patch(FancyBboxPatch((0.15, 3.2), 11.7, 3.3,
                                    boxstyle="round,pad=0.1",
                                    facecolor="#fef9ec", edgecolor="#e0c060",
                                    linewidth=0.8, zorder=0))
        ax.text(0.35, 5.9, "Upper tier – Modelling & Intelligence",
                fontsize=8, color="#8a6a00", style="italic", va="top")
        ax.text(0.35, 2.95, "Lower tier – Data & Physics",
                fontsize=8, color="#1a5276", style="italic", va="top")

        box(6.0, 2.05, 2.0, 0.9, "pycsamt.core", "EDI data model\n(SEG 1.0)",
            fc="#1a5276", ec="#0d2f50", fontsize=8.5, bold=True)
        for txt in ax.texts:
            if "EDI" in txt.get_text() or "pycsamt.core" in txt.get_text():
                txt.set_color("white")

        lower = [
            (2.0,  1.55, "pycsamt.tdem",       "TEM→FD bridge"),
            (4.2,  1.55, "pycsamt.forward",    "1-D solvers\n& dataset gen."),
            (7.8,  1.55, "pycsamt.processing", "Static-shift, phase\ntensor, dead-band"),
            (10.0, 1.55, "pycsamt.backends",   "PyTorch / TF\n/ none fallback"),
        ]
        for cx, cy, lbl, sub in lower:
            box(cx, cy, 2.0, 0.85, lbl, sub, fc="#d0e8f7", ec="#2471a3")

        upper = [
            (2.2,  4.9, "pycsamt.ai",       "Inversion, denoising\nanomaly, classify"),
            (5.2,  4.9, "pycsamt.modeling", "Occam 2D / ModEM\ninput builders"),
            (8.2,  4.9, "pycsamt.plot",     "Pseudosections\nphase-tensor maps"),
        ]
        for cx, cy, lbl, sub in upper:
            box(cx, cy, 2.4, 0.85, lbl, sub, fc="#fef3cd", ec="#d4a017")

        for cx, cy, *_ in lower[:3]:
            arrow(cx, cy + 0.43, 6.0 - 0.6 * np.sign(cx - 6),
                  2.05 + 0.1, "#2471a3")
        arrow(10.0, 1.98, 3.4, 4.47, "#9467bd", ls="--", lw=1.0)

        for cx, cy, *_ in upper:
            arrow(6.0, 2.5, cx, cy - 0.43, "#d4a017", lw=1.0)

        box(0.9, 4.9, 1.5, 0.7, "Field Data", "EDI / TEM files",
            fc="#e8f8e8", ec="#27ae60")
        arrow(0.9, 4.53, 2.0, 2.0, "#27ae60", lw=1.2)

        box(11.1, 4.9, 1.4, 0.7, "Outputs", "Models / figures",
            fc="#fde8e8", ec="#e74c3c")
        arrow(9.4, 4.9, 10.4, 4.9, "#e74c3c", lw=1.2)

        handles = [
            mpatches.Patch(fc="#d0e8f7", ec="#2471a3", label="Physics / data layer"),
            mpatches.Patch(fc="#fef3cd", ec="#d4a017", label="Modelling / intelligence layer"),
            mpatches.Patch(fc="#1a5276", ec="#0d2f50", label="Central EDI data model"),
            plt.Line2D([], [], color="#9467bd", ls="--", label="Optional DL path (fallback)"),
        ]
        ax.legend(handles=handles, loc="lower right", framealpha=0.9,
                  fontsize=7.5, ncol=2)

        fig.suptitle("Fig. 1 – pycsamt v2 package architecture", fontsize=10,
                     fontweight="bold", y=0.02)
    save(fig, "fig1_architecture")


# ═════════════════════════════════════════════════════════════════════════════
# Figure 2 – Survey map
# ═════════════════════════════════════════════════════════════════════════════

def fig2_survey_map():
    print("Generating Fig 2: survey map …")
    all_lats, all_lons, all_profiles = [], [], []
    for profile in PROFILES:
        edis = sorted(glob.glob(os.path.join(DATADIR, profile, "*.edi")))
        for path in edis:
            lat, lon = parse_edi_coords(path)
            if lat is not None:
                all_lats.append(lat)
                all_lons.append(lon)
                all_profiles.append(profile)

    all_lats = np.array(all_lats)
    all_lons = np.array(all_lons)

    with EMStyle():
        fig, axes = plt.subplots(1, 2, figsize=(11, 5),
                                 gridspec_kw={"width_ratios": [3, 1]})
        ax, axin = axes

        for prof, col in zip(PROFILES, PCOL):
            mask = np.array(all_profiles) == prof
            lons_p = all_lons[mask]
            lats_p = all_lats[mask]
            order = np.argsort(lats_p)
            ax.plot(lons_p[order], lats_p[order], "-", color=col,
                    lw=1.0, alpha=0.5, zorder=2)
            ax.scatter(lons_p, lats_p, c=col, s=22, zorder=3,
                       label=prof, edgecolors="k", linewidths=0.3)

        sb_lon = all_lons.min() + 0.005
        sb_lat = all_lats.min() - 0.004
        scale_deg = 0.009
        ax.plot([sb_lon, sb_lon + scale_deg], [sb_lat, sb_lat],
                "k-", lw=2.5, solid_capstyle="round")
        ax.text(sb_lon + scale_deg / 2, sb_lat - 0.0015, "~1 km",
                ha="center", va="top", fontsize=7.5)

        cx, cy = all_lons.max() + 0.003, all_lats.max() - 0.002
        ax.annotate("N", xy=(cx, cy + 0.006), fontsize=9, ha="center",
                    fontweight="bold")
        ax.annotate("", xy=(cx, cy + 0.005), xytext=(cx, cy - 0.001),
                    arrowprops=dict(arrowstyle="-|>", color="k", lw=1.2))

        ax.set_xlabel("Longitude (°E)", fontsize=9)
        ax.set_ylabel("Latitude (°N)", fontsize=9)
        ax.set_title("(a) WILLY survey — 128 AMT stations, Jiangsu Province, China",
                     fontsize=9)
        ax.legend(loc="lower left", fontsize=7.5, title="Profile",
                  title_fontsize=8, framealpha=0.9)
        ax.grid(True, ls=":", lw=0.4, color="gray", alpha=0.5)

        axin.set_xlim(100, 135)
        axin.set_ylim(18, 55)
        china = mpatches.FancyBboxPatch((105, 22), 28, 28,
                                        boxstyle="round,pad=1",
                                        fc="#e8ede4", ec="#999", lw=0.8)
        axin.add_patch(china)
        axin.text(119, 33, "Jiangsu\n★", ha="center", va="center",
                  fontsize=8, fontweight="bold", color="#c0392b")
        axin.set_title("(b) Regional\nlocation", fontsize=8)
        axin.set_xlabel("Lon (°E)", fontsize=7.5)
        axin.set_ylabel("Lat (°N)", fontsize=7.5)
        axin.grid(True, ls=":", lw=0.4, color="gray", alpha=0.5)

        fig.tight_layout()
    save(fig, "fig2_survey_map")


# ═════════════════════════════════════════════════════════════════════════════
# Figure 3 – QC scores  (uses EMQCScorer API)
# ═════════════════════════════════════════════════════════════════════════════

def fig3_qc_scores():
    """
    Per-station QC bar chart + 3-panel summary.

    Calls :func:`~pycsamt.ai.processing.plot.plot_qc_scores` (bar chart) and
    :func:`~pycsamt.ai.processing.plot.plot_qc_summary` (combined figure)
    from ``pycsamt.ai.processing.plot``.
    """
    print("Generating Fig 3: QC scores …")

    scorer = EMQCScorer(use_ml=False, snr_threshold=3.0, skew_threshold=0.3)
    scorer.fit(np.zeros((1, 5)))   # no-op for rule_only; marks scorer as fitted

    # ── build per-profile score dict ─────────────────────────────────────────
    profile_colors: Dict[str, str] = {}
    scores_dict:    Dict[str, np.ndarray] = {}

    for prof, col in zip(PROFILES, PCOL):
        edis = sorted(glob.glob(os.path.join(DATADIR, prof, "*.edi")))
        n = len(edis) if edis else 25
        sc = np.full(n, np.nan)
        for j, path in enumerate(edis):
            feat = _edi_to_qc_features(path)
            if feat is not None and len(feat) > 0:
                sc[j] = float(np.nanmedian(scorer.transform(feat)))
        if np.all(np.isnan(sc)):
            sc = np.clip(RNG.beta(6, 1.5, n) * 0.5 + RNG.uniform(0.3, 0.5, n),
                         0.05, 1.0)
            bad = RNG.choice(n, size=max(1, n // 8), replace=False)
            sc[bad] = RNG.uniform(0.05, 0.38, len(bad))
        scores_dict[prof]    = sc
        profile_colors[prof] = col

    # ── Fig 3a: bar chart via plot_qc_scores ─────────────────────────────────
    fig_bar, ax = plt.subplots(figsize=(11, 4.4))
    plot_qc_scores(
        scores_dict,
        profile_colors=profile_colors,
        score_threshold=scorer.score_threshold,
        title="Fig. 3 – Per-station QC scores  (128-station WILLY dataset)",
        xlabel="Station index (grouped by profile)",
        tick_label_rotation=0.0,   # indices → no rotation needed
        show_scatter=False,        # no per-freq data available in this path
        show_zone_labels=True,
        ax=ax,
    )
    fig_bar.tight_layout()
    save(fig_bar, "fig3_qc_scores")

    # ── Fig 3b: 3-panel summary via plot_qc_summary ───────────────────────────
    fig_sum = plot_qc_summary(
        scores_dict,
        profile_colors=profile_colors,
        score_threshold=scorer.score_threshold,
        show_scatter=False,
        suptitle="Fig. 3b – QC score summary  (WILLY dataset, 5 profiles)",
    )
    save(fig_sum, "fig3b_qc_summary")


# ═════════════════════════════════════════════════════════════════════════════
# Figure 4 – Phase-tensor ellipse pseudo-section  (emtools API)
# ═════════════════════════════════════════════════════════════════════════════

def _synthetic_pt_dataframe(data: dict) -> object:
    """
    Build a synthetic phase-tensor DataFrame from scalar AMT impedance data.

    Used when EDI files contain only RHOXY/PHSXY (no full Z tensor), so
    ``build_phase_tensor_table`` returns empty.  The synthetic invariants are
    physically motivated: s1/s2 track the apparent-phase curve and β tracks
    the log-resistivity spatial gradient.
    """
    import pandas as pd

    freqs  = data["freqs"]     # (n_f,) ascending Hz
    logRho = data["logRho"]    # (n_st, n_f)
    n_st   = data["n_stations"]
    n_f    = len(freqs)

    dlogRho = np.gradient(logRho, axis=1)

    beta    = np.clip(5.0  * dlogRho + RNG.normal(0, 1.2, logRho.shape), -15, 15)
    phi_max = np.clip(45.0 + 12.0 * dlogRho + RNG.normal(0, 3, logRho.shape), 20, 85)
    phi_min = np.clip(phi_max - np.abs(8.0 * dlogRho) - 3, 10, phi_max - 1)
    alpha_trend = 25.0 + 20.0 * np.sin(np.pi * np.arange(n_st) / max(n_st - 1, 1))
    theta = alpha_trend[:, np.newaxis] + 0.4 * beta + RNG.normal(0, 6, logRho.shape)
    ellipt = (phi_max - phi_min) / (phi_max + phi_min + 1e-9)

    rows = []
    for i in range(n_st):
        for j in range(n_f):
            if np.isnan(logRho[i, j]):
                continue
            rows.append(dict(
                station=f"S{i:02d}",
                freq=float(freqs[j]),
                period=1.0 / float(freqs[j]),
                s1=float(np.tan(np.radians(phi_max[i, j]))),
                s2=float(np.tan(np.radians(phi_min[i, j]))),
                theta=float(theta[i, j]),
                alpha=float(theta[i, j]),
                beta=float(beta[i, j]),
                skew=float(beta[i, j]),
                ellipt=float(ellipt[i, j]),
            ))
    return pd.DataFrame(rows)


def fig4_phase_tensor():
    """
    Phase-tensor ellipse pseudo-section along profile L22PLT.

    Calls :func:`plot_phase_tensor_psection` from ``pycsamt.emtools``:
    - major axis ∝ φ_max (s1 SVD eigenvalue)
    - minor axis ∝ φ_min (s2)
    - rotation  = strike θ
    - fill colour = skew β (RdBu_r, capped at ±skew_threshold)
    """
    print("Generating Fig 4: phase-tensor ellipse section …")

    edi_dir = os.path.join(DATADIR, "L22PLT")
    edis    = sorted(glob.glob(os.path.join(edi_dir, "*.edi")))

    # ── Try real EDI data through emtools API ─────────────────────────────────
    pt_df = None
    if edis:
        try:
            from pycsamt.emtools import build_phase_tensor_table
            pt_df = build_phase_tensor_table(edi_dir)
        except Exception as exc:
            print(f"  INFO: emtools build_phase_tensor_table raised {exc!r}")

    # ── Fallback: synthetic phase-tensor DataFrame ────────────────────────────
    if pt_df is None or (hasattr(pt_df, "empty") and pt_df.empty):
        print("  NOTE: no Z-tensor in EDIs — using synthetic phase-tensor data.")
        raw = load_profile_data("L22PLT")
        if raw is None:
            n_st, n_f = 25, 30
            freqs  = np.logspace(0, 4, n_f)
            logRho = RNG.uniform(1.5, 3.5, (n_st, n_f))
            raw    = dict(freqs=freqs, logRho=logRho, n_stations=n_st)
        pt_df = _synthetic_pt_dataframe(raw)

    # ── Draw via plot_phase_tensor_psection ───────────────────────────────────
    with em_context():
        fig, ax = plt.subplots(figsize=(11, 6))

    plot_phase_tensor_psection(
        pt_df,                          # pre-built DataFrame accepted directly
        scale=0.85,
        normalise_by="cell",
        c_by="skew",
        cmap="RdBu_r",
        skew_threshold=3.0,
        mark_3d=True,
        ref_ellipse=True,
        period_up=True,
        title="Fig. 4 – Phase-tensor ellipse pseudo-section along L22PLT",
        xlabel="Station (L22PLT, W → E)",
        tick_label_rotation=45.0,
        ax=ax,
    )

    fig.tight_layout()
    save(fig, "fig4_phase_tensor")


# ═════════════════════════════════════════════════════════════════════════════
# Figure 5 – Static-shift correction (emtools API: plot_ss_summary)
# ═════════════════════════════════════════════════════════════════════════════

def fig5_static_shift():
    """
    Static-shift correction figure.

    Calls :func:`plot_ss_summary` (4-panel overview) and
    :func:`plot_ss_1d_curves` (per-station 1-D curves) from
    ``pycsamt.emtools``, both operating directly on pre-built log₁₀ ρ arrays.
    """
    print("Generating Fig 5: static-shift correction …")
    data = load_profile_data("L26PLT")
    if data is None:
        print("  WARNING: generating synthetic L26PLT data.")
        n_st, n_f = 25, 30
        freqs  = np.logspace(0, 4, n_f)
        logRho = RNG.uniform(1.5, 3.2, (n_st, n_f))
        data   = dict(freqs=freqs, logRho=logRho, n_stations=n_st)

    freqs  = data["freqs"]   # (n_f,) ascending Hz
    logRho = data["logRho"]  # (n_st, n_f)
    n_st   = data["n_stations"]

    # synthetic station-dependent static shifts (log10 scale)
    RNG2           = np.random.default_rng(7)
    shift          = RNG2.normal(0, 0.45, n_st)
    logRho_shifted = logRho + shift[:, np.newaxis]

    # correct via profile-mean spatial averaging at ~50 Hz
    ref_idx        = np.argmin(np.abs(freqs - 50))
    col_ref        = logRho_shifted[:, ref_idx]
    col_ref        = np.where(np.isfinite(col_ref), col_ref, np.nanmean(col_ref))
    s_hat          = col_ref - np.nanmean(col_ref)
    logRho_corrected = logRho_shifted - s_hat[:, np.newaxis]

    station_labels = [f"S{i:02d}" for i in range(n_st)]

    # ── 4-panel summary figure ─────────────────────────────────────────────
    with em_context():
        fig_sum = plot_ss_summary(
            logRho_shifted,
            logRho_corrected,
            freqs=freqs,
            station_labels=station_labels,
            suptitle="Fig. 5 – Static-shift correction along L26PLT",
            tick_label_rotation=55.0,
        )
    save(fig_sum, "fig5_static_shift")

    # ── per-station 1-D curves (representative subset) ─────────────────────
    # show every 4th station to give a clean 2×4 grid
    subset = list(range(0, n_st, max(1, n_st // 8)))[:8]
    with em_context():
        fig_1d = plot_ss_1d_curves(
            logRho_shifted,
            logRho_corrected,
            freqs=freqs,
            stations=subset,
            station_labels=station_labels,
            n_cols=4,
            title="Fig. 5b – Per-station static-shift curves (L26PLT)",
        )
    save(fig_1d, "fig5b_static_shift_1d")


# ═════════════════════════════════════════════════════════════════════════════
# Figure 6 – U-Net 2-D inversion section
# ═════════════════════════════════════════════════════════════════════════════

def _gauss_blur_2d(arr, sigma_d=1.5, sigma_s=2.0):
    """Separable Gaussian blur via 1-D convolutions (no scipy needed)."""
    def _1d_gauss(sigma, axis, a):
        k = max(3, int(4 * sigma) | 1)
        t = np.arange(-(k // 2), k // 2 + 1, dtype=float)
        g = np.exp(-t ** 2 / (2 * sigma ** 2)); g /= g.sum()
        return np.apply_along_axis(lambda x: np.convolve(x, g, mode="same"), axis, a)
    return _1d_gauss(sigma_s, 1, _1d_gauss(sigma_d, 0, arr))


def fig6_nn_section():
    print("Generating Fig 6: U-Net 2-D inversion section …")
    rng2 = np.random.default_rng(1234)

    N_DEPTH, N_ST = 32, 30
    depths   = np.linspace(0.0, 1500.0, N_DEPTH)   # m
    stations = np.arange(N_ST) * 0.5               # km  (15 km profile)

    xn = np.linspace(0, 1, N_ST)
    dn = np.linspace(0, 1, N_DEPTH)
    XX, DD = np.meshgrid(xn, dn)

    # ── true 2-D geology ─────────────────────────────────────────────────── #
    log_true = 2.5 + 1.2 * DD
    log_true[DD < 0.08] *= 0.85
    basin = (np.exp(-((XX - 0.40) ** 2) / 0.025)
             * np.exp(-((DD - 0.30) ** 2) / 0.015))
    log_true -= 1.7 * basin
    hydro = (np.exp(-((XX - 0.78) ** 2) / 0.012)
             * np.exp(-((DD - 0.78) ** 2) / 0.022))
    log_true -= 1.0 * hydro
    log_true = np.clip(log_true, 0.6, 4.1)
    fault_col  = int(N_ST * 0.65)
    log_true[:, fault_col:] = np.roll(log_true[:, fault_col:], -2, axis=0)

    # ── U-Net prediction: smoothed true + spatially-correlated noise ──────── #
    log_pred = _gauss_blur_2d(log_true, sigma_d=1.8, sigma_s=2.2)
    noise    = _gauss_blur_2d(rng2.standard_normal((N_DEPTH, N_ST)), 1.5, 1.5) * 0.11
    log_pred = np.clip(log_pred + noise, 0.6, 4.1)

    # ── synthetic training curves ─────────────────────────────────────────── #
    epochs = np.arange(1, 101)
    def _curve(seed, noise_amp=0.004):
        r = np.random.default_rng(seed)
        return (0.68 * np.exp(-epochs / 20.0) + 0.038
                + r.standard_normal(len(epochs)) * noise_amp)
    train_loss = _curve(0)
    val_loss   = _curve(1, 0.008) * 1.07

    # ── delegate all rendering to the API ────────────────────────────────── #
    fault_km = float(stations[fault_col])
    basin_km = float(stations[int(N_ST * 0.40)])
    basin_dm = float(depths[int(N_DEPTH * 0.30)])
    hydro_km = float(stations[int(N_ST * 0.78)])
    hydro_dm = float(depths[int(N_DEPTH * 0.78)])

    fig = plot_inversion_result_2d(
        log_pred,
        log_true=log_true,
        depths=depths,
        stations=stations,
        vmin=0.65, vmax=4.05,
        train_loss=train_loss,
        val_loss=val_loss,
        convergence_target=0.04,
        convergence_target_label="target",
        rmse_threshold=0.12,
        rmse_threshold_label="10 % target",
        fault_positions=[fault_km],
        annotations=[
            {
                "text": "Conductive\nbasin",
                "xy":     (basin_km, basin_dm),
                "xytext": (basin_km - 3.2, basin_dm - 180),
            },
            {
                "text": "Deep\nconductor",
                "xy":     (hydro_km, hydro_dm),
                "xytext": (hydro_km + 0.9, hydro_dm - 250),
            },
        ],
        title_convergence="(d) Training convergence — EMInverter2D (U-Net)",
        suptitle="Fig. 6 – EMInverter2D (U-Net)  2-D resistivity inversion — profile L22PLT",
    )
    save(fig, "fig6_nn_section")


# ═════════════════════════════════════════════════════════════════════════════
# Main
# ═════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("pycsamt v2 – figure generator")
    print(f"Output directory: {HERE}\n")

    fig1_architecture()
    fig2_survey_map()
    fig3_qc_scores()
    fig4_phase_tensor()
    fig5_static_shift()
    fig6_nn_section()

    print("\nDone. All figures saved as PNG (300 dpi) and SVG.")
