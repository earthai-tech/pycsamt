"""
Frequency editing and EMAP denoising
=====================================

Once the confidence diagnostics (see
:doc:`plot_coverage_and_confidence`) show *where* the data is weak, the
next step is to act on it: recover or drop individual frequencies by
confidence, and suppress spatially incoherent noise with an
electromagnetic-array-profiling (EMAP) filter. This example runs both
edits on line **L22PLT**, keeping an untouched baseline alongside each
edited copy so the before/after figures are honest comparisons.

The editing routines live in :mod:`pycsamt.emtools.frequency`
(confidence-driven frequency edits) and
:mod:`pycsamt.emtools.remove_noise` (EMAP filtering and its
confidence-gated variant).
"""

# sphinx_gallery_thumbnail_number = 5

# %%
# Confidence-driven frequency editing
# -----------------------------------
# :func:`~pycsamt.emtools.frequency.edit_frequencies_by_confidence`
# scores every ``(station, frequency)`` cell and, in ``"recover"`` mode,
# keeps high-confidence cells, attempts to recover marginal ones, and
# drops the rest. We load two independent copies of L22PLT — one frozen
# baseline, one to edit — so nothing leaks between them.

import logging
from contextlib import contextmanager

from _datasets import load_sites

from pycsamt.emtools.frequency import (
    edit_frequencies_by_confidence,
    plot_frequency_edit_decisions,
    plot_frequency_edit_summary,
)


@contextmanager
def quiet_recompute_logs():
    """Silence the impedance object's expected recompute diagnostics.

    Recover-mode editing masks rejected cells; setting the (now
    non-finite) Z then logs a "cannot recompute rho/phi" error per masked
    cell. Those are internal diagnostics, not failures — scope a global
    ``logging.disable`` around the edit so the example output stays
    readable without affecting any other gallery script.
    """
    logging.disable(logging.ERROR)
    try:
        yield
    finally:
        logging.disable(logging.NOTSET)


from pycsamt.emtools.remove_noise import (
    apply_emap_filter,
    confidence_gated_emap_filter,
    emap_filter_report,
    plot_emap_filter_profile,
    plot_emap_filter_psection,
)

before = load_sites("amt_l22plt")
source = load_sites("amt_l22plt")

with quiet_recompute_logs():
    edit = edit_frequencies_by_confidence(
        source,
        before_sites=before,
        mode="recover",
        method="composite",
        ci_hi=0.90,
        ci_lo=0.50,
        reject="drop",
    )

# The edit returns a per-cell decision table and a per-station report —
# inspect them just as you would on your own data.
print("report columns:", list(edit.report.columns))
print(edit.report.head().to_string(index=False))
print("\ndecision table columns:", list(edit.decisions.columns))
print(edit.decisions.head().to_string(index=False))

# %%
# 1. Edit summary
# ---------------
# :func:`~pycsamt.emtools.frequency.plot_frequency_edit_summary`
# contrasts the confidence distribution before and after the edit, so
# you can confirm the operation improved the line as a whole.

plot_frequency_edit_summary(
    before,
    edit.sites,
    method="composite",
    ci_hi=0.90,
    ci_lo=0.50,
    figsize=(10, 4.0),
)

# %%
# 2. Per-cell edit decisions
# --------------------------
# :func:`~pycsamt.emtools.frequency.plot_frequency_edit_decisions` maps
# the actual keep / recover / drop decision onto the station-by-period
# grid — a precise audit of what the edit changed and where.

plot_frequency_edit_decisions(
    before,
    edit.sites,
    method="composite",
    ci_hi=0.90,
    ci_lo=0.50,
    figsize=(12, 4.8),
)

# %%
# EMAP FLMA denoising
# -------------------
# :func:`~pycsamt.emtools.remove_noise.apply_emap_filter` applies a
# first-order moving-average (FLMA) EMAP filter along the profile,
# averaging each station's response with its neighbours to suppress
# static and spatially random noise on the ``xy`` component.

emap_before = load_sites("amt_l22plt")
emap_source = load_sites("amt_l22plt")
flma = apply_emap_filter(emap_source, method="flma", window=5, component="xy")

report = emap_filter_report(
    emap_before,
    flma,
    component="xy",
    period_s=0.01,
)
print(report.head().to_string(index=False))

# %%
# 3. Filter effect on one station profile
# ---------------------------------------
# :func:`~pycsamt.emtools.remove_noise.plot_emap_filter_profile` overlays
# the raw and filtered apparent resistivity along the line at a fixed
# period, showing how much smoothing the window applies.

plot_emap_filter_profile(
    emap_before,
    flma,
    method="flma",
    component="xy",
    period_s=0.01,
    figsize=(10.5, 4.0),
)

# %%
# 4. Before / after pseudo-section
# --------------------------------
# :func:`~pycsamt.emtools.remove_noise.plot_emap_filter_psection` gives
# the full-band picture: raw and FLMA-filtered pseudo-sections
# side-by-side over every station and period.

plot_emap_filter_psection(
    emap_before,
    flma,
    method="flma",
    component="xy",
    figsize=(11.5, 8.2),
)

# %%
# 5. Confidence-gated EMAP
# ------------------------
# :func:`~pycsamt.emtools.remove_noise.confidence_gated_emap_filter`
# combines both ideas: it only lets the FLMA filter modify low-confidence
# cells, leaving trustworthy data untouched. This is the recommended
# denoiser — it cleans the noisy corners of the band without smearing the
# parts that were already good.

gated_before = load_sites("amt_l22plt")
gated_source = load_sites("amt_l22plt")
with quiet_recompute_logs():
    gated = confidence_gated_emap_filter(
        gated_source,
        before_sites=gated_before,
        method="flma",
        confidence_method="composite",
        component="xy",
        window=5,
        ci_hi=0.90,
        ci_lo=0.50,
    )
print("gated decision table columns:", list(gated.decisions.columns))
print(gated.decisions.head().to_string(index=False))

plot_emap_filter_psection(
    gated_before,
    gated.sites,
    method="confidence-gated FLMA",
    component="xy",
    figsize=(11.5, 8.2),
)
