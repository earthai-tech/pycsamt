r"""
Geoelectric strike estimation and visualization (:mod:`pycsamt.emtools.strike`)
================================================================================

:mod:`pycsamt.emtools.strike` estimates the geoelectric strike
direction — the preferred 2-D structural axis inferred from the MT
impedance tensor — three independent ways (impedance-tensor rotation
sweep, phase-tensor azimuth, and a weighted consensus blend), applies
that estimate to rotate data onto strike, and renders it through five
plot styles: a per-frequency ribbon, along-line profile, geographic
map-sticks, and single/multi-line rose diagrams, finishing with a
three-panel diagnostic comparable to MTPy's ``StrikeAnalysis`` plot.
This example uses **L18PLT** (``data/AMT/WILLY_DATA/``) throughout,
and brings in its sibling line **L22PLT** for the multi-line rose
comparison in section 6.
"""

# %%
# 1. Three ways to estimate strike
# ------------------------------------
# :func:`~pycsamt.emtools.strike.estimate_strike_sweep` rotates the
# impedance tensor in 1-degree steps and picks the angle that best
# diagonalises (or best off-diagonalises, depending on ``metric``) it
# per frequency, then reports the per-station median.
# :func:`~pycsamt.emtools.strike.estimate_strike_phase_tensor` instead
# reads the phase-tensor skew angle directly.
# :func:`~pycsamt.emtools.strike.estimate_strike_consensus` blends the
# two.

import numpy as np
from _datasets import dataset_path, load_survey

from pycsamt.emtools import (
    estimate_strike_consensus,
    estimate_strike_phase_tensor,
    estimate_strike_sweep,
    plot_strike_analysis,
    plot_strike_mapsticks,
    plot_strike_profile,
    plot_strike_ribbon,
    plot_strike_rose,
    plot_strike_rose_by_line,
    rotate_to_strike,
    strike_curve_sweep,
)

survey = load_survey("amt_l18plt")

t_sweep = estimate_strike_sweep(survey)
t_pt = estimate_strike_phase_tensor(survey)
t_cons = estimate_strike_consensus(survey)

print(t_sweep.head())
print(
    f"sweep     ang range: {t_sweep['ang'].min():.1f} to "
    f"{t_sweep['ang'].max():.1f} deg  (iqr median {t_sweep['iqr'].median():.1f})"
)
print(
    f"phase-tns ang range: {t_pt['ang'].min():.1f} to "
    f"{t_pt['ang'].max():.1f} deg  (iqr median {t_pt['iqr'].median():.1f})"
)
print(
    f"consensus ang range: {t_cons['ang'].min():.1f} to "
    f"{t_cons['ang'].max():.1f} deg  (iqr median {t_cons['iqr'].median():.1f})"
)

# %%
# **Reading this output.** Every station's *reported* angle
# (``ang``) already falls inside a sensible -90 to +90 degree window,
# but the *stability* behind that number differs enormously between
# methods: the sweep method's per-station IQR (its scatter across
# frequency) has a median of 104.5 degrees — more than half the entire
# axial range — while the phase-tensor method's median IQR is only
# 35.4 degrees. The rotation sweep is far noisier station-to-station
# on this dataset than the phase-tensor azimuth.

# %%
# 2. The correlation trap: why raw ``corrcoef`` misleads for strike
# ------------------------------------------------------------------------
# Strike is axial (0-180 degrees, with 180 degrees equivalent to 0):
# a naive Pearson correlation between two methods' raw angles ignores
# that wraparound and can report nonsense.

m = t_sweep.merge(t_pt, on="station", suffixes=("_sw", "_pt"))
naive_corr = np.corrcoef(m["ang_sw"], m["ang_pt"])[0, 1]
print(f"naive corrcoef(sweep, pt): {naive_corr:.3f}")

axial_diff = ((m["ang_sw"] - m["ang_pt"] + 90) % 180) - 90
print(f"median |axial diff| sweep vs pt: {axial_diff.abs().median():.1f} deg")
print(f"mean   |axial diff| sweep vs pt: {axial_diff.abs().mean():.1f} deg")
print(axial_diff.describe())

# %%
# **Reading this output.** The naive correlation comes back a weak
# 0.155 — which would suggest the two methods barely agree — but that
# understates things, because it treats angles near +90 and -90 as
# maximally different when, axially, they are almost the same
# direction. Computing the actual axial difference station-by-station
# tells a better story: the median absolute difference is only 19.0
# degrees and the interquartile range of that difference is roughly
# -12 to +19 degrees — the two methods mostly agree to within about 20
# degrees, consistent with the well-known result that phase-tensor
# azimuth is more robust to 3-D/noise effects than a raw
# impedance-tensor rotation sweep (matching the IQR gap found in
# section 1), not that
# they disagree.

# %%
# 3. Rotating data onto strike
# ---------------------------------
# :func:`~pycsamt.emtools.strike.rotate_to_strike` rotates every
# station's impedance tensor by its own estimated strike angle.
# Fixed along the way: the function located each station's angle
# correctly but then rotated the wrong object — a freshly-wrapped
# ``Sites`` collection rather than the underlying EDI item — so
# :func:`~pycsamt.site.edit.rotate` could never find the ``.Z``
# section it needs and the whole operation was a silent no-op; every
# call to ``rotate_to_strike`` returned data byte-identical to its
# input. Now fixed, rotating the underlying EDI item directly.

before = estimate_strike_consensus(survey)
rotated = rotate_to_strike(survey, method="consensus")
after = estimate_strike_consensus(rotated)

print(f"mean |strike| before rotation: {before['ang'].abs().mean():.1f} deg")
print(f"mean |strike| after  rotation: {after['ang'].abs().mean():.1f} deg")

# %%
# **Reading this output.** Re-estimating the consensus strike on the
# rotated data pulls the mean ``|angle|`` from 36.6 down to 24.6 degrees —
# a real, if imperfect, reduction (imperfect because the consensus
# estimate mixes two methods that do not each go exactly to zero under
# a single rotation, and because strike estimation itself is noisy).
# Before the fix, this number never moved at all.

# %%
# 4. Per-frequency strike: the curve and the ribbon
# -------------------------------------------------------
# :func:`~pycsamt.emtools.strike.strike_curve_sweep` returns one
# smoothed sweep angle per station *per frequency* — the same
# rotation-sweep idea as section 1, but without collapsing across
# frequency first. :func:`~pycsamt.emtools.strike.plot_strike_ribbon`
# renders it as a station x period image: hue encodes strike angle,
# saturation encodes local stability (white = high local variance).

curve = strike_curve_sweep(survey)
print(curve.shape, list(curve.columns))
print(
    f"stations x frequencies: {curve['station'].nunique()} x "
    f"{curve.groupby('station').size().iloc[0]}"
)

ax = plot_strike_ribbon(survey, method="sweep")

# %%
# 5. Rose diagrams — one survey, then a per-line comparison
# ------------------------------------------------------------------
# :func:`~pycsamt.emtools.strike.plot_strike_rose` groups stations by
# profile line automatically when no ``groups`` are given, using a
# station-name prefix heuristic (letters, then digits — e.g.
# ``"E1S01"`` groups to ``"E1"``). L18PLT's station names follow the
# opposite convention (``"18-001A"`` — digits first), so every station
# becomes its own singleton "group"; ``plot_strike_rose`` falls back
# to a single ``"All"`` rose rather than showing nothing.

fig = plot_strike_rose(survey)
ax_all = fig.get_axes()[0]
print("single-rose annotation:", [t.get_text() for t in ax_all.texts])

# %%
# **Reading this output.** The annotation reads ``143.1 deg, n=28`` —
# every one of the 28 stations folded into one axial histogram. That
# is a reasonable fallback, but it cannot compare *between* lines.
# :func:`~pycsamt.emtools.strike.plot_strike_rose_by_line` has no such
# fallback (it shows "no groups" instead), so a real multi-line
# comparison needs an explicit ``groups`` mapping — which is also the
# natural way to use it, since the function is built to compare
# profile lines against each other. Bringing in the sibling line
# **L22PLT**:

p18, p22 = dataset_path("amt_l18plt"), dataset_path("amt_l22plt")
files18 = sorted(p18.glob("*.edi"))
files22 = sorted(p22.glob("*.edi"))
combined = files18 + files22
groups = {
    "L18PLT": [f.stem for f in files18],
    "L22PLT": [f.stem for f in files22],
}
fig2 = plot_strike_rose_by_line(combined, groups=groups)
for ax in fig2.get_axes():
    print(ax.get_title(), [t.get_text() for t in ax.texts])

# %%
# **Reading this output.** L18PLT (28 stations, 143.1 degrees) and
# L22PLT (25 stations, 144.1 degrees) come back barely a degree apart
# — a geologically sensible check for two nearby lines from the same
# survey.

# %%
# 6. Frequency-band decomposition
# ------------------------------------
# ``bar_style="bands"`` stacks separate histograms per period band in
# one rose, showing whether shallow and deep structure share a strike.

fig3 = plot_strike_rose(
    survey,
    bar_style="bands",
    freq_bands=[(0.001, 0.01), (0.01, 1.0)],
    band_labels=["short period", "long period"],
)
print("bands figure axes:", len(fig3.get_axes()))

# %%
# 7. Geographic view: map-sticks
# ------------------------------------
# :func:`~pycsamt.emtools.strike.plot_strike_mapsticks` draws one
# short line segment per station at its real (lon, lat), oriented
# along its estimated strike — darker/more opaque where the estimate
# is more stable (lower IQR).

ax = plot_strike_mapsticks(survey)
print(f"map extent: lon {ax.get_xlim()}, lat {ax.get_ylim()}")

# %%
# 8. Along-line profile
# --------------------------
# :func:`~pycsamt.emtools.strike.plot_strike_profile` plots strike
# (with an IQR ribbon) against station order. Fixed along the way: the
# same recurring bug found elsewhere in this gallery (see :doc:`/user_guide/emtools/ss`,
# :doc:`/user_guide/emtools/dimensionality`) — ``sort_by="lon"``/``"lat"``/``"auto"``
# checked only flat ``.lon``/``.lat`` attributes that real ``Site``
# objects do not have (coordinates live in ``.coords``), so every
# station order silently collapsed to alphabetical-by-name. Now fixed.

ax_name = plot_strike_profile(survey, sort_by="name")
order_name = [t.get_text() for t in ax_name.get_xticklabels()]
ax_lon = plot_strike_profile(survey, sort_by="lon")
order_lon = [t.get_text() for t in ax_lon.get_xticklabels()]
print("name order:", order_name[:6])
print("lon  order:", order_lon[:6])
print("orders differ now?", order_name != order_lon)

# %%
# **Reading this output.** ``sort_by="lon"`` (and the default
# ``"auto"``) now genuinely differs from ``sort_by="name"`` — the true
# west-to-east station order for this near-north-south line
# (``18-020A, 18-024U, 18-022U, ...``) bears little resemblance to
# alphabetical order, because longitude barely varies along a line
# that runs mostly north-south (the same axis-choice caveat already
# documented in :doc:`/user_guide/emtools/ss`).

# %%
# 9. Combined strike/PT-azimuth diagnostic
# ----------------------------------------
# :func:`~pycsamt.emtools.strike.plot_strike_analysis` puts Strike (Z)
# and PT Azimuth side by side in one figure, adding a third Tipper
# Strike panel only when *survey* actually carries a vertical-field
# channel somewhere — rather than reserving a panel that would only
# ever read "no data".

fig4 = plot_strike_analysis(survey)
for ax in fig4.get_axes():
    print(ax.get_title(), [t.get_text() for t in ax.texts])

# %%
# **Reading this output.** Only two panels are drawn — L18PLT, like the
# other AMT lines in ``data/AMT/WILLY_DATA/``, has no vertical-field
# (tipper) channel at all (the same reason :doc:`/user_guide/emtools/tf`'s
# induction-vector example uses **KAP03** instead). Strike (Z) and PT
# Azimuth still agree reasonably well here (140.8 vs 147.4 degrees),
# consistent with the axial-difference analysis in section 2.
