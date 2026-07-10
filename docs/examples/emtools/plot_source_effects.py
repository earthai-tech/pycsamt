r"""
CSAMT source overprint and near-field effects (:mod:`pycsamt.emtools.source_effects`)
======================================================================================

:mod:`pycsamt.emtools.source_effects` quantifies how the artificial
CSAMT transmitter contaminates a sounding, using two independent
formulas from two different papers: the Yan & Fu (2004) /
Da et al. (2016) ground-wave/surface-wave amplitude ratio
:math:`\beta_{Ey}` (a source-overprint index, threshold 3%), and the
Wang & Lin (2023) skin-depth field-zone classification with a
near-field correction built on the equatorial HED transfer function
:math:`F(p)`. Unlike natural-source MT, CSAMT needs the
source-receiver offset :math:`r` to evaluate either formula, and no
real offset is bundled with these docs (a standard EDI file does not
carry it). This example reuses the same honest choice already made in
the ``fieldzone`` example: a representative 2 km offset for **L18PLT**
(``data/AMT/WILLY_DATA/``), chosen for illustration and stated plainly
rather than read from survey metadata.
"""

# %%
# 1. The core concept: how :math:`\beta_{Ey}` depends on offset
# --------------------------------------------------------------------
# :func:`~pycsamt.emtools.source_effects.overprint_beta` is pure
# arithmetic — no sites needed. At fixed resistivity, sweeping
# frequency at three representative offsets (the same 500 m / 2 km /
# 8 km trio used in the ``fieldzone`` example) shows the expected
# physical trend: a closer transmitter pushes the "safe" (low
# :math:`\beta`) frequency band down.

import matplotlib.pyplot as plt
import numpy as np

from pycsamt.emtools import (
    BETA_THRESH_PCT,
    correct_near_field,
    detect_source_overprint,
    normalize_response,
    overprint_beta,
    plot_normalized_response,
    plot_overprint_section,
    source_overprint_table,
)

freq_sweep = np.logspace(-1, 3, 60)
fig, ax = plt.subplots(figsize=(7.5, 4.2))
for offset, color in zip([500.0, 2000.0, 8000.0], ["#d62728", "#ff7f0e", "#2ca02c"]):
    beta = overprint_beta(rho=300.0, freq=freq_sweep, offset=offset)
    ax.loglog(freq_sweep, beta, color=color, label=f"r={offset:g} m")
    above = freq_sweep[beta > BETA_THRESH_PCT]
    if above.size:
        print(f"offset={offset:g} m: beta>{BETA_THRESH_PCT:g}% for f up to {above.max():.3g} Hz")
ax.axhline(BETA_THRESH_PCT, color="0.3", ls="--", lw=1, label=f"threshold {BETA_THRESH_PCT:g}%")
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel(r"$\beta_{Ey}$ (%)")
ax.set_title(r"$\beta_{Ey}$ vs. frequency at $\rho$=300 $\Omega\cdot$m")
ax.grid(True, which="both", alpha=0.25)
ax.legend(fontsize=8)

# %%
# **Reading this output/figure.** At the same 300 :math:`\Omega\cdot`m
# half-space, the frequency below which :math:`\beta_{Ey}` exceeds 3%
# shrinks as the offset grows: up to 1000 Hz at 500 m, only 392 Hz at
# 2 km, and just 27.6 Hz at 8 km. A closer transmitter contaminates a
# *wider* frequency band — exactly the intuition behind treating offset
# as the single most important survey-design choice, echoing the
# ``fieldzone`` example's finding for the skin-depth field zones.

# %%
# 2. Per-frequency overprint on a real sounding
# ------------------------------------------------------
# :func:`~pycsamt.emtools.source_effects.detect_source_overprint`
# applies the same formula, station by station, using each station's
# own real observed :math:`\rho_a` and frequency instead of a fixed
# half-space value.

import warnings  # noqa: E402

from _datasets import load_survey  # noqa: E402

survey = load_survey("amt_l18plt")

with warnings.catch_warnings():
    warnings.simplefilter("ignore")  # no source_offset metadata in real EDIs
    detail = detect_source_overprint(survey, source_offset=2000.0)

print(detail["beta_pct"].describe())
print(f"flagged (beta > {BETA_THRESH_PCT:g}%): "
      f"{int(detail['overprint_flag'].sum())}/{len(detail)} station-frequency rows")

# %%
# **Reading this output.** At the assumed 2 km offset, 944 of 1484
# rows (64%) are flagged — a real, if severe, consequence of treating
# an ordinary CSAMT-band AMT-style line as if it had this specific
# transmitter geometry: :func:`overprint_beta` in section 1 already
# showed 2 km keeps :math:`\beta` under 3% only above ~392 Hz, and this
# line's real resistivities push that crossover frequency around even
# further station to station. This is exactly why the offset has to be
# a real, measured survey parameter in practice — it is not a detail
# that can be guessed after the fact from the impedance alone.

# %%
# 3. Per-station summary and the da2016 slope criterion
# ------------------------------------------------------------------
# :func:`~pycsamt.emtools.source_effects.source_overprint_table`
# collapses each station to :math:`\beta_{max}`/:math:`\beta_{mean}`
# plus a low-/high-frequency log-log slope comparison (da2016): a
# strongly negative ``slope_delta`` signals a resistivity contrast
# beneath the source. The default ``f_split=1.0`` Hz sits right at this
# survey's own lower frequency limit (~1.008 Hz), leaving no
# low-frequency samples at all — raising ``f_split`` to 50 Hz actually
# splits the band in two.

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    table_default_split = source_overprint_table(survey, source_offset=2000.0)
    table = source_overprint_table(survey, source_offset=2000.0, f_split=50.0)

print("f_split=1.0 Hz (default): all lf_slope NaN?",
      table_default_split["lf_slope"].isna().all())
print(table[["station", "beta_max_pct", "overprint_frac", "slope_delta"]]
      .sort_values("slope_delta").head(3))
print(table[["station", "beta_max_pct", "overprint_frac", "slope_delta"]]
      .sort_values("slope_delta").tail(3))

# %%
# **Reading this output.** With a sensible split, ``18-020A`` has the
# most negative ``slope_delta`` (-1.13) — its low-frequency
# :math:`\rho_a` trend rises far more steeply than its high-frequency
# one, the da2016 signature of a contrast beneath the source. Every
# station's ``beta_max_pct`` sits at essentially the same ceiling
# (~50%, the formula's own numerical cap at very small offset-to-skin
# -depth ratios), so ``overprint_frac`` (the *fraction* of frequencies
# flagged) is the more discriminating per-station number here, not the
# maximum.

# %%
# 4. The overprint pseudo-section
# --------------------------------------
# :func:`~pycsamt.emtools.source_effects.plot_overprint_section` shows
# :math:`\beta_{Ey}` for every station and period at once, with a
# dashed white contour at the 3% threshold.

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    plot_overprint_section(survey, source_offset=2000.0)

# %%
# **Reading this figure.** The 3% threshold contour (dashed white)
# sweeps through most of the middle of the period range at nearly every
# station — consistent with section 2's 64% flagged fraction — with a
# visibly darker, more saturated patch around ``18-021B`` at short
# periods where :math:`\beta_{Ey}` climbs highest. The same station
# reappears in section 6 below from a completely different formula.

# %%
# 5. Wang & Lin (2023): normalized response and field zone
# ------------------------------------------------------------------
# :func:`~pycsamt.emtools.source_effects.normalize_response` divides
# observed :math:`\rho_a` by a reference half-space value and subtracts
# a reference phase, while independently classifying each
# (station, frequency) point into near/transition/far zones using the
# Wang & Lin skin-depth thresholds (0.5\ :math:`\delta`, 4\ :math:`\delta`)
# — the same zone definitions used in the ``fieldzone`` example.

norm = normalize_response(survey, rho_ref=300.0, source_offset=2000.0)
print(norm["zone"].value_counts())

# %%
# **Reading this output.** At the same assumed 2 km offset, the 1484
# station-frequency points split roughly evenly: 603 far, 468
# transition, 413 near — a genuinely mixed survey where no single zone
# dominates, the same kind of offset-sensitive mixed outcome the
# ``fieldzone`` example found when comparing assumed offsets directly.

# %%
# 6. Two-panel normalized pseudo-section
# ----------------------------------------------
# :func:`~pycsamt.emtools.source_effects.plot_normalized_response`
# plots :math:`\rho_n=\rho_{obs}/\rho_{ref}` and
# :math:`\varphi_{diff}=\varphi_{obs}-\varphi_{ref}` side by side.

fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.0))
plot_normalized_response(survey, rho_ref=300.0, source_offset=2000.0, axes=axes)
fig.tight_layout()

# %%
# **Reading this figure.** The left panel's darkest red patch sits at
# short periods around ``18-021B``/``18-022U`` — the same station
# ``plot_overprint_section`` (section 4) flagged from the unrelated
# :math:`\beta_{Ey}` formula. The right (phase) panel splits roughly in
# two along the line: stations up to about ``18-017U`` run mostly red
# (phase above the 45-degree reference) while later stations show more
# blue — a spatial pattern in phase behaviour that neither the
# amplitude-only :math:`\beta_{Ey}` view nor ``rho_n`` alone reveals.

# %%
# 7. Near-field correction
# -------------------------------
# :func:`~pycsamt.emtools.source_effects.correct_near_field` divides
# :math:`Z` by the complex near-field factor
# :math:`F(p)=1-3/p^2+3/p^3`. In the far field :math:`F\to 1` and the
# correction does nothing; in the near field it can be large.

from pycsamt.emtools import ensure_sites  # noqa: E402
from pycsamt.emtools._core import (  # noqa: E402
    _get_z_block,
    _iter_items,
    _name,
)


def rho_xy(sites, name="18-016A"):
    s = ensure_sites(sites, recursive=False)
    for i, ed in enumerate(_iter_items(s)):
        if _name(ed, i) == name:
            _, z, fr = _get_z_block(ed)
            return 0.2 * np.abs(z[:, 0, 1]) ** 2 / fr, fr
    raise KeyError(name)


before, fr0 = rho_xy(survey, "18-016A")
corrected = correct_near_field(survey, source_offset=2000.0)
after, _ = rho_xy(corrected, "18-016A")
print(f"18-016A: max |change| in log10(rho) after near-field correction = "
      f"{np.nanmax(np.abs(np.log10(before) - np.log10(after))):.2f}")

# %%
# **Reading this output.** A swing of more than 6 decades in
# :math:`\log_{10}\rho_a` — the correction can shrink the apparent
# resistivity by a factor of roughly 2.6 million at the frequencies
# where this station sits deepest in the near field for a 2 km offset.
# That is not a bug: :math:`F(p)` genuinely diverges as the near-field
# geometric term dominates, and a correction this large is itself a
# useful diagnostic — it flags exactly how far the raw, uncorrected
# sounding departed from the plane-wave assumption at those
# frequencies, rather than a small, reassuring nudge.

# %%
# 8. Advanced: do the two independent formulas agree?
# ------------------------------------------------------------
# Section 2's :math:`\beta_{Ey}` overprint flag (Yan & Fu 2004) and
# section 5's skin-depth zone (Wang & Lin 2023) come from two unrelated
# papers and formulas. Merging them point by point checks whether they
# actually agree on which frequencies are source-contaminated.

merged = detail.merge(
    norm[["station", "freq_hz", "zone"]], on=["station", "freq_hz"],
)
agreement = merged.groupby("zone")["overprint_flag"].mean()
print(agreement)

# %%
# **Reading this output.** Agreement is essentially total for the
# unambiguous cases: 100% of "near" and 100% of "transition" points are
# also flagged by :math:`\beta_{Ey}`, while only 10.4% of "far" points
# are. Two formulas built on different physical arguments — one an
# amplitude-ratio criterion, the other a skin-depth threshold — reach
# the same practical conclusion about which frequencies are usable at
# this assumed offset, which is a genuine, useful cross-check before
# trusting either one alone on a real survey where the offset is known
# rather than assumed.
