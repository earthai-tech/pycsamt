r"""
CSAMT field-zone classification (:mod:`pycsamt.emtools.fieldzone`)
======================================================================

:mod:`pycsamt.emtools.fieldzone` answers a question specific to
controlled-source AMT: at a given frequency and source-receiver offset
*r*, is the measurement in the plane-wave (far) zone the standard MT
apparent-resistivity formula assumes, or has the receiver drifted into
the near/transition zone where that assumption breaks down? Both
:func:`~pycsamt.emtools.fieldzone.classify_field_zones` and
:func:`~pycsamt.emtools.fieldzone.near_field_factor` reduce to the same
dimensionless parameter,

.. math::

    |k \cdot r| = \frac{r}{\delta_B}, \qquad
    \delta_B = 356\sqrt{\rho_a / f}

(Chen & Yan 2005) — the source-receiver distance measured in Bostick
skin depths.

This is a genuinely CSAMT-specific concept: it needs a real
source-receiver offset, which a natural-source AMT survey like
**L18PLT** (``data/AMT/WILLY_DATA/``) does not record (there is no
transmitter). As in the ``csumt`` example, this page applies the same
physics to L18PLT's real, CSAMT-band (1 Hz-10.4 kHz) apparent
resistivity using a few representative assumed offsets — a legitimate
way to see how the method behaves, while being upfront that the exact
offset numbers are chosen for illustration, not read from survey
metadata.
"""

# %%
# 1. The |k.r| parameter itself
# ------------------------------------
# Before any real data: for a fixed apparent resistivity, |k.r| depends
# only on frequency and the chosen offset. Larger offsets or higher
# frequencies push a sounding toward the far field faster.

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pycsamt.emtools.fieldzone import (
    classify_field_zones,
    near_field_factor,
    plot_field_zones,
)

BOSTICK_C = 356.0
RHO_DEMO = 300.0
freq = np.logspace(0, 4, 300)  # 1 Hz-10 kHz, matches L18PLT's band
delta_b = BOSTICK_C * np.sqrt(RHO_DEMO / freq)

fig, ax = plt.subplots(figsize=(7, 5))
for offset, color in zip([500.0, 2000.0, 8000.0], ["#d62728", "#ff7f0e", "#2ca02c"]):
    ax.loglog(freq, offset / delta_b, color=color, label=f"r={offset:g} m")
ax.axhspan(3.0, 1e6, color="#2ca02c", alpha=0.08)
ax.axhspan(0.3, 3.0, color="#ff7f0e", alpha=0.10)
ax.axhspan(1e-6, 0.3, color="#d62728", alpha=0.08)
ax.axhline(3.0, color="0.3", ls="--", lw=0.8)
ax.axhline(0.3, color="0.3", ls="--", lw=0.8)
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel(r"$|k \cdot r|$")
ax.set_title(fr"$|k \cdot r|$ vs. frequency  ($\rho_a$={RHO_DEMO:g} $\Omega\cdot$m)")
ax.legend(fontsize=8)
ax.grid(True, which="both", alpha=0.3)

# %%
# **Reading this figure.** At a fixed 300 Ω·m, an 8 km offset only
# leaves the far-field band (green) below about 5 Hz — it is in
# far field for almost the entire displayed range. A 500 m offset is
# far worse off: it drops out of far field already below 1.4 kHz and
# reaches near-field territory (red) below about 14 Hz. The offset you
# assume changes which part of your own recorded band you should trust
# at face value.

# %%
# 2. One real station's |k.r| curve
# -------------------------------------------
# :func:`~pycsamt.emtools.fieldzone.classify_field_zones` computes |k.r|
# from each station's actual measured apparent resistivity rather than
# an assumed constant — using a representative 2 km offset.

from _datasets import load_survey  # noqa: E402

survey = load_survey("amt_l18plt")
OFFSET_DEMO = 2000.0
zones = classify_field_zones(survey, OFFSET_DEMO)
station = "18-001A"
d = zones[zones["station"] == station].sort_values("period_s")

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.axhspan(3.0, 1e4, color="#2ca02c", alpha=0.08)
ax.axhspan(0.3, 3.0, color="#ff7f0e", alpha=0.10)
ax.axhspan(1e-4, 0.3, color="#d62728", alpha=0.08)
ax.loglog(d["period_s"], d["kr"], "o-", ms=3, color="0.2")
ax.set_xlabel("Period (s)")
ax.set_ylabel(r"$|k \cdot r|$")
ax.set_title(f"{station} — measured |k.r|  (r={OFFSET_DEMO:g} m)")

# %%
# **Reading this figure.** |k.r| falls from about 65 at the shortest
# period, dipping as low as 0.08 around 0.24 s before settling around
# 0.1-0.2 at the longest periods — this station's own sounding crosses
# all three zones within its recorded band, entirely because apparent
# resistivity and frequency both change with period, not because the
# assumed offset changes.

# %%
# 3. Cross-checking against the near-field correction factor
# -----------------------------------------------------------------------
# :func:`~pycsamt.emtools.fieldzone.near_field_factor` computes a
# continuous bias factor |F(p)| from the full complex wavenumber,
# independent of the |k.r| threshold rule. If the threshold-based
# zoning is doing its job, |F| should sit close to 1 in the far zone
# and depart sharply as the near zone is entered.

nff = near_field_factor(survey, OFFSET_DEMO)
merged = zones.merge(nff[["station", "freq_hz", "nf_factor"]], on=["station", "freq_hz"])
colors = {"far": "#2ca02c", "transition": "#ff7f0e", "near": "#d62728"}

fig, ax = plt.subplots(figsize=(6.5, 5))
for zone in ("far", "transition", "near"):
    m = merged["zone"] == zone
    ax.scatter(merged.loc[m, "kr"], merged.loc[m, "nf_factor"],
               s=10, alpha=0.5, color=colors[zone], label=zone)
ax.axhline(1.0, color="0.3", ls=":", lw=1)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$|k \cdot r|$")
ax.set_ylabel(r"near-field factor $|F(p)|$")
ax.legend(fontsize=8)
ax.set_title("L18PLT — near-field factor vs. |k.r|, all stations")

# %%
# **Reading this figure.** The two independent computations agree
# closely: |F| averages 0.99 in the far zone (essentially unbiased),
# climbs to a mean of about 13 in transition, and explodes to a mean
# above 900 (up to nearly 17,000 at the most extreme point) in the near
# zone — a smooth, monotonic blow-up as |k.r| falls below 1, exactly
# where the plane-wave approximation is expected to fail. That the
# threshold rule and the continuous formula tell the same story is a
# useful internal consistency check before trusting either on its own.

# %%
# 4. The module's pseudo-section
# --------------------------------------
# :func:`~pycsamt.emtools.fieldzone.plot_field_zones` is the headline
# view: every station's zone across the whole recorded band, with
# dashed white |k.r| contours.

plot_field_zones(survey, OFFSET_DEMO)

# %%
# **Reading this figure.** Every station shows the same broad pattern
# as the single-station curve in section 2 — green (far) at the bottom
# (shortest periods), warming through orange to red toward the top
# (longest periods) — because every station in this line shares a
# comparable resistivity range and the same assumed 2 km offset.

# %%
# 5. Advanced: how much does the assumed offset matter?
# -----------------------------------------------------------------
# Section 1 showed this qualitatively; here it is quantified across the
# whole survey for three representative offsets.

offsets = [500.0, 2000.0, 8000.0]
fracs = []
for off in offsets:
    z = classify_field_zones(survey, off)
    vc = z["zone"].value_counts(normalize=True)
    fracs.append({"offset": off, **{k: vc.get(k, 0.0) for k in ("far", "transition", "near")}})
frac_df = pd.DataFrame(fracs).set_index("offset")
print(frac_df.round(3))

fig, ax = plt.subplots(figsize=(6.5, 4.5))
x = np.arange(len(offsets))
bottom = np.zeros(len(offsets))
for zone in ("far", "transition", "near"):
    vals = frac_df[zone].to_numpy()
    ax.bar(x, vals, bottom=bottom, color=colors[zone], label=zone, width=0.6)
    bottom += vals
ax.set_xticks(x, [f"{o:g} m" for o in offsets])
ax.set_ylabel("fraction of (station, frequency) pairs")
ax.set_xlabel("assumed source offset")
ax.legend(fontsize=8, loc="upper right")
ax.set_title("L18PLT — zone mix vs. assumed offset")

# %%
# **Reading this figure.** The near-field fraction drops from 35% at
# 500 m to essentially 0% at 8 km, while the far-field fraction rises
# from 26% to 71% over the same range — the *same* recorded data reads
# as mostly contaminated or mostly clean depending entirely on an
# offset number that, for a real CSAMT survey, must come from the
# actual transmitter-receiver geometry. Getting it wrong doesn't just
# shift a percentage: it silently changes which frequencies you should
# have discarded before inversion.

# %%
# 6. Advanced: near vs. far offset, side by side
# ----------------------------------------------------------
# The same pseudo-section at the two extreme offsets from section 5,
# sharing one figure via ``ax``.

fig, (axa, axb) = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
plot_field_zones(survey, 500.0, ax=axa)
axa.set_title("r = 500 m")
plot_field_zones(survey, 8000.0, ax=axb)
axb.set_title("r = 8000 m")
fig.tight_layout()

# %%
# **Reading this figure.** The 500 m panel is dominated by orange/red
# down through the middle of the band; the 8 km panel is nearly all
# green until the very longest periods. Same stations, same
# resistivities — the offset alone decides which one a practitioner
# would trust.

# %%
# 7. Advanced: what the near-field bias would do to a sounding curve
# ------------------------------------------------------------------------------
# The module's own docstring describes ``nf_factor`` as the bias factor
# on apparent resistivity in the near field. Dividing the measured
# curve by :math:`|F(p)|^2` gives a rough sense of how far off an
# uncorrected near-field reading could be — illustrative here, not a
# substitute for checking the exact convention in Chen & Yan (2005)
# before applying it to a real inversion.

st_merged = merged[merged["station"] == station].sort_values("period_s")
corrected = st_merged["rho_a_ohmm"] / st_merged["nf_factor"] ** 2

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.loglog(st_merged["period_s"], st_merged["rho_a_ohmm"], "o-", ms=3,
          color="0.3", label="measured (uncorrected)")
ax.loglog(st_merged["period_s"], corrected, "s--", ms=3,
          color="#d62728", label=r"/ $|F|^2$ (illustrative)")
ax.set_xlabel("Period (s)")
ax.set_ylabel(r"$\rho_a$ ($\Omega\cdot$m)")
ax.legend(fontsize=8)
ax.set_title(f"{station} — measured vs. near-field-corrected ρ_a")

# %%
# **Reading this figure.** The two curves overlap almost exactly at
# short period (the far zone, where |F| ≈ 1) and diverge sharply at the
# longest periods, where section 2 already showed this station sits
# deep in the near-field zone — precisely the frequencies where a plain
# plane-wave inversion would be misled by this station's own data.
