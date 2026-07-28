r"""
Reviewer-response audit for conditioned AMT/MT data
===================================================

Peer review often asks for more than a clean final pseudo-section.  Reviewers
may ask how confidence scores were defined, whether low-quality frequencies
were interpolated or rejected, how electromagnetic interference was handled,
which adaptive moving-average parameters were used, whether 2-D inversion is
defensible, and whether galvanic distortion was checked.

This case study shows how to answer those questions with reproducible tables
and figures instead of prose alone.  It is written as a general field-data
audit template:

* replace the two WILLY example paths with your own EDI folders;
* keep the same report tables as an audit trail;
* include the generated figures in a manuscript response, thesis appendix, or
  processing notebook.

The example uses the bundled WILLY ``L18PLT`` and ``L22PLT`` lines because
they are small enough for the gallery while still behaving like real survey
lines.
"""

# %%
# 1. Imports, paths, and audit policy
# -----------------------------------
# The output directory is deliberately local to the example gallery.  It is
# ignored by git, so the generated CSV and PNG files can be inspected or copied
# without becoming source files.

from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    EXAMPLE_DIR = Path(__file__).resolve().parent
except NameError:
    EXAMPLE_DIR = Path.cwd()


def repo_root() -> Path:
    root = os.environ.get("PYCSAMT_DOCS_REPO_ROOT")
    return Path(root) if root else EXAMPLE_DIR.parents[2]


ROOT = repo_root()
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(EXAMPLE_DIR) not in sys.path:
    sys.path.insert(0, str(EXAMPLE_DIR))

from pycsamt.emtools import (  # noqa: E402
    confidence_gated_emap_filter,
    correct_ss_ama,
    edit_frequencies_by_confidence,
    emi_mitigation_report,
    ensure_sites,
    estimate_ss_ama,
    frequency_confidence_table,
    groom_bailey_table,
    plot_confidence_band_summary,
    plot_confidence_profile,
    plot_frequency_confidence_psection,
    plot_skew_traffic_psection,
    plot_ss_delta_profile,
    plot_strike_profile,
    pre2d_inversion_assessment,
    qc_flags,
    station_confidence_table,
)
from pycsamt.emtools.remove_noise import snr_table  # noqa: E402

LINES = {
    "L18PLT": ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT",
    "L22PLT": ROOT / "data" / "AMT" / "WILLY_DATA" / "L22PLT",
}

OUT = EXAMPLE_DIR / "workspaces" / "reviewer_response_audit"
OUT.mkdir(parents=True, exist_ok=True)

CI_HI = 0.95
CI_LO = 0.85
MAIN_HZ = 50.0
STATIC_SHIFT_WINDOWS = (1, 2, 3, 4)

# %%
# 2. Mathematical definitions used in the audit
# ---------------------------------------------
# The confidence ratio is a finite-score weighted mean:
#
# .. math::
#
#    CR = \frac{\sum_k w_k s_k I_k}{\sum_k w_k I_k},
#    \qquad 0 \le s_k \le 1.
#
# ``s_k`` are diagnostic scores, ``w_k`` are their weights, and ``I_k`` is
# one only when a score is finite.  In this workflow the default pyCSAMT
# scores are coverage, tensor uncertainty, off-diagonal consistency, diagonal
# leakage, phase smoothness, and spatial coherence.  The review classes are:
#
# .. math::
#
#    CR \ge 0.95 \quad \mathrm{preserve},\qquad
#    0.85 \le CR < 0.95 \quad \mathrm{recover/review},\qquad
#    CR < 0.85 \quad \mathrm{reject/mask}.
#
# The adaptive moving-average static-shift estimate uses determinant
# apparent resistivity
#
# .. math::
#
#    \rho_{\mathrm{det}}(f)=\sqrt{\rho_{xy}(f)\rho_{yx}(f)},\qquad
#    \rho_{ij}(f)=0.2|Z_{ij}(f)|^2/f.
#
# For station ``i``, the log-resistivity deviation from a neighboring spatial
# trend is reduced to ``Delta_i`` and applied as
#
# .. math::
#
#    F_{\rho,i}=10^{-\Delta_i},\qquad F_{Z,i}=10^{-\Delta_i/2}.
#
# Power-line harmonics are documented by the frequency rule
#
# .. math::
#
#    |f - n f_{\mathrm{mains}}| \le \Delta f,\qquad n=1,\ldots,N_h.
#
# Groom--Bailey-style galvanic-distortion diagnostics fit
#
# .. math::
#
#    Z_{\mathrm{obs}}(f) \approx D Z_{2D}(f),
#
# with real, frequency-independent distortion matrix ``D``.


def write_table(table, path: Path) -> pd.DataFrame:
    """Write a DataFrame-like object and return a plain DataFrame."""
    df = getattr(table, "df", table)
    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)
    df.to_csv(path, index=False)
    return df


def save_current_figure(path: Path) -> None:
    """Save the current gallery figure to the audit workspace."""
    fig = plt.gcf()
    fig.savefig(path, dpi=180, bbox_inches="tight")


# %%
# 3. Load each line and create the core audit tables
# --------------------------------------------------
# These CSV files are the part most useful to a student facing similar
# reviewer comments.  They show what was calculated, which rows were edited,
# and what limitations remain.

sites_by_line = {}
summary_rows = []

for line, edi_dir in LINES.items():
    line_out = OUT / line
    line_out.mkdir(parents=True, exist_ok=True)

    sites = ensure_sites(str(edi_dir), recursive=False, verbose=0)
    sites_by_line[line] = sites

    freq_conf = write_table(
        frequency_confidence_table(
            sites,
            method="composite",
            ci_hi=CI_HI,
            ci_lo=CI_LO,
        ),
        line_out / "frequency_confidence.csv",
    )
    write_table(
        station_confidence_table(
            sites,
            method="composite",
        ),
        line_out / "station_confidence.csv",
    )
    write_table(qc_flags(sites), line_out / "qc_flags.csv")
    write_table(snr_table(sites), line_out / "snr_table.csv")
    write_table(
        emi_mitigation_report(
            sites,
            remote_reference_attempted=False,
            remote_reference_reason=(
                "The bundled gallery EDIs do not include remote-reference "
                "time series. This audit documents transfer-function-level "
                "mitigation and diagnostics."
            ),
            mains_hz=MAIN_HZ,
        ),
        line_out / "emi_mitigation_report.csv",
    )

    edited = edit_frequencies_by_confidence(
        sites,
        mode="recover",
        before_sites=sites,
        ci_hi=CI_HI,
        ci_lo=CI_LO,
        reject="mask",
        interpolation="linear",
    )
    write_table(edited.report, line_out / "frequency_edit_report.csv")
    write_table(edited.decisions, line_out / "frequency_edit_decisions.csv")

    emap = confidence_gated_emap_filter(
        sites,
        before_sites=sites,
        method="flma",
        confidence_method="composite",
        component="xy",
        ci_hi=CI_HI,
        ci_lo=CI_LO,
        window=5,
    )
    write_table(emap.report, line_out / "confidence_gated_emap_report.csv")
    write_table(
        emap.decisions, line_out / "confidence_gated_emap_decisions.csv"
    )

    safe = int((freq_conf["confidence"] >= CI_HI).sum())
    recoverable = int(
        (
            (freq_conf["confidence"] >= CI_LO)
            & (freq_conf["confidence"] < CI_HI)
        ).sum()
    )
    reject = int((freq_conf["confidence"] < CI_LO).sum())
    summary_rows.append(
        {
            "line": line,
            "stations": len(station_confidence_table(sites)),
            "frequency_rows": len(freq_conf),
            "safe_rows": safe,
            "recoverable_rows": recoverable,
            "reject_rows": reject,
            "emap_preserved": emap.n_preserved,
            "emap_blended": emap.n_blended,
            "emap_filtered": emap.n_filtered,
        }
    )

summary = pd.DataFrame(summary_rows)
summary.to_csv(OUT / "all_lines_summary.csv", index=False)
print(summary.to_string(index=False))

# %%
# 4. Confidence figures: justify preserve/recover/reject decisions
# ----------------------------------------------------------------
# These plots answer the common request: "show where the confidence classes
# occur, not only the final smoothed response."  The first figure is a station
# profile; the second is a station-period confidence pseudo-section; the third
# summarizes the confidence bands.

for line, sites in sites_by_line.items():
    fig_dir = OUT / line / "figures"
    fig_dir.mkdir(exist_ok=True)

    ax = plot_confidence_profile(
        sites,
        method="composite",
        ci_hi=CI_HI,
        ci_lo=CI_LO,
        figsize=(9.0, 3.8),
    )
    ax.set_title(f"{line}: station confidence profile")
    save_current_figure(fig_dir / "confidence_profile.png")

    ax = plot_frequency_confidence_psection(
        sites,
        method="composite",
        ci_hi=CI_HI,
        ci_lo=CI_LO,
        figsize=(10.0, 4.8),
    )
    ax.set_title(f"{line}: frequency-confidence pseudo-section")
    save_current_figure(fig_dir / "frequency_confidence_psection.png")

    ax = plot_confidence_band_summary(
        sites,
        method="composite",
        ci_hi=CI_HI,
        ci_lo=CI_LO,
        figsize=(9.0, 4.2),
    )
    ax.set_title(f"{line}: confidence band summary")
    save_current_figure(fig_dir / "confidence_band_summary.png")

# %%
# 5. EMI and frequency-edit evidence
# ----------------------------------
# The tables above document remote-reference status, power-line harmonic
# counts, station-level SNR, and edit decisions.  A compact reviewer response
# can cite these tables directly and state that recovered frequencies remain
# flagged by the decision table rather than being treated as independent
# measurements.

for line in LINES:
    line_out = OUT / line
    emi = pd.read_csv(line_out / "emi_mitigation_report.csv")
    edit = pd.read_csv(line_out / "frequency_edit_report.csv")
    print(f"\n{line}: EMI and edit summary")
    print(
        emi[
            [
                "station",
                "remote_reference_attempted",
                "harmonic_z_samples",
                "applied_measures",
            ]
        ]
        .head(5)
        .to_string(index=False)
    )
    print(
        edit[
            [
                "station",
                "n_freq_before",
                "n_freq_after",
                "confidence_median_before",
                "confidence_median_after",
            ]
        ]
        .head(5)
        .to_string(index=False)
    )

# %%
# 6. Static-shift sensitivity: show that AMA parameters were checked
# ------------------------------------------------------------------
# The reviewer does not only need the final static-shift correction; they
# need to know whether the moving-average window suppresses real geology.
# Here we run several neighbor windows and write a long-form sensitivity
# table.  A stable correction should not change wildly as the window is
# perturbed.

for line, sites in sites_by_line.items():
    rows = []
    for half_window in STATIC_SHIFT_WINDOWS:
        factors = write_table(
            estimate_ss_ama(
                sites,
                sort_by="lon",
                half_window=half_window,
                weights="tri",
                max_skew=6.0,
            ),
            OUT / line / f"ss_ama_factors_half_window_{half_window}.csv",
        )
        for _, row in factors.iterrows():
            rows.append(
                {
                    "line": line,
                    "half_window": half_window,
                    "station": row.get("station"),
                    "delta_log10_rho": row.get("delta_log10_rho"),
                    "fac_rho": row.get("fac_rho"),
                    "fac_z": row.get("fac_z"),
                    "n_used": row.get("n_used"),
                }
            )
    sensitivity = pd.DataFrame(rows)
    sensitivity.to_csv(
        OUT / line / "ss_ama_window_sensitivity.csv", index=False
    )

    fig, ax = plt.subplots(figsize=(10.0, 4.3))
    for half_window, sub in sensitivity.groupby("half_window"):
        sub = sub.sort_values("station")
        ax.plot(
            np.arange(len(sub)),
            sub["fac_rho"].to_numpy(dtype=float),
            marker="o",
            ms=3.5,
            lw=1.2,
            label=f"half_window={half_window}",
        )
    ax.axhline(1.0, color="0.5", lw=1.0)
    labels = sensitivity[
        sensitivity["half_window"] == STATIC_SHIFT_WINDOWS[0]
    ]["station"].astype(str)
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_ylabel(r"$F_\rho$")
    ax.set_title(f"{line}: AMA static-shift window sensitivity")
    ax.grid(True, axis="y", ls=":", alpha=0.5)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    save_current_figure(
        OUT / line / "figures" / "ss_ama_window_sensitivity.png"
    )

    corrected = correct_ss_ama(
        sites,
        sort_by="lon",
        half_window=3,
        weights="tri",
        max_skew=6.0,
        inplace=False,
    )
    ax = plot_ss_delta_profile(sites, corrected, figsize=(9.0, 3.8))
    ax.set_title(f"{line}: applied static-shift delta profile")
    save_current_figure(OUT / line / "figures" / "ss_delta_profile.png")

# %%
# 7. Dimensionality, strike, and rotation readiness
# -------------------------------------------------
# This table answers several reviewer questions at once: how much of the data
# is 1-D/2-D/3-D, where high skew appears, what strike estimate is obtained,
# and whether rotation/Groom--Bailey status was documented.

for line, sites in sites_by_line.items():
    gb = write_table(
        groom_bailey_table(sites, min_freq=4, robust=True),
        OUT / line / "groom_bailey_table.csv",
    )
    gb_ok = bool((gb.get("status", pd.Series(dtype=str)) == "ok").any())
    assessment = write_table(
        pre2d_inversion_assessment(
            sites,
            rotation_applied=False,
            groom_bailey_attempted=True,
            groom_bailey_applied=False,
            groom_bailey_reason=(
                "Groom-Bailey table was computed as a diagnostic. "
                "Correction was not applied in this gallery audit; "
                f"valid station fits present: {gb_ok}."
            ),
        ),
        OUT / line / "pre2d_inversion_assessment.csv",
    )
    print(f"\n{line}: pre-2D assessment")
    print(
        assessment[
            [
                "station",
                "frac_2d",
                "frac_3d",
                "beta_abs_p95",
                "strike_consensus_deg",
                "strike_consensus_iqr_deg",
                "recommendation",
            ]
        ]
        .head(8)
        .to_string(index=False)
    )

    ax = plot_skew_traffic_psection(sites, figsize=(10.0, 4.6))
    ax.set_title(f"{line}: skew traffic pseudo-section")
    save_current_figure(OUT / line / "figures" / "skew_traffic_psection.png")

    ax = plot_strike_profile(sites, figsize=(10.0, 4.3))
    ax.set_title(f"{line}: geoelectric strike profile")
    save_current_figure(OUT / line / "figures" / "strike_profile.png")

# %%
# 8. Groom--Bailey distortion diagnostics
# ---------------------------------------
# If Groom--Bailey correction is not applied, the response should still say
# whether it was considered.  The diagnostic table records fitted status,
# distortion parameters, and diagonal leakage before/after the fitted
# correction.  Do not blindly apply this correction; inspect accepted station
# fits first.

for line in LINES:
    gb = pd.read_csv(OUT / line / "groom_bailey_table.csv")
    cols = [
        "station",
        "status",
        "n_freq",
        "twist_deg",
        "shear_angle_deg",
        "anisotropy",
        "diagonal_ratio_before",
        "diagonal_ratio_after",
    ]
    present = [col for col in cols if col in gb.columns]
    print(f"\n{line}: Groom-Bailey diagnostic rows")
    print(gb[present].head(8).to_string(index=False))

# %%
# 9. What to copy into a report
# -----------------------------
# The generated workspace now contains both human-readable figures and
# machine-readable CSV files.  A concise response to a reviewer can cite:
#
# * ``frequency_confidence.csv`` and confidence figures for CR definition and
#   preserve/recover/reject classes;
# * ``frequency_edit_decisions.csv`` and ``frequency_edit_report.csv`` for
#   the treatment of reconstructed or masked frequencies;
# * ``emi_mitigation_report.csv`` and ``snr_table.csv`` for interference
#   documentation;
# * ``ss_ama_window_sensitivity.csv`` and static-shift figures for moving
#   average parameter sensitivity;
# * ``pre2d_inversion_assessment.csv`` for dimensionality, skew, strike, and
#   rotation readiness;
# * ``groom_bailey_table.csv`` for galvanic-distortion diagnostics.
#
# To use this example with a new survey, replace ``LINES`` near the top with
# paths to the student's EDI folders.  The rest of the workflow can remain the
# same unless the project uses different confidence thresholds, power-line
# frequency, or inversion period band.

print(f"\nAudit workspace written to:\n{OUT}")
