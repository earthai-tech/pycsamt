r"""
Uncertainty quantification
==========================

Every hydro-geophysical number carries the uncertainty of the petrophysical
assumptions behind it. :class:`~pycsamt.interp.MonteCarloHydro` propagates
priors on the rock-physics parameters (pore-water resistivity, cementation
exponent, ...) through the whole workflow, producing P10/P50/P90 maps and
profiles instead of single values. This example quantifies the uncertainty
on the synthetic section's aquifer properties.
"""

# %%
# Monte-Carlo propagation
# -----------------------
# :class:`~pycsamt.interp.UncertaintyBounds` sets the prior ranges; each
# draw re-runs the hydro model, and ``run`` aggregates the ensemble into an
# :class:`~pycsamt.interp.UncertaintyResult` with percentile maps and
# coefficient-of-variation fields.

from _interp_data import demo_model

from pycsamt.interp import (
    MonteCarloHydro,
    PetrophysicalConfig,
    UncertaintyBounds,
)

# Use the uncertainty section (2nd figure) as the card thumbnail.
# sphinx_gallery_thumbnail_number = 2

rm = demo_model()
cfg = PetrophysicalConfig(rho_w=20.0, porosity_prior=0.25)
bounds = UncertaintyBounds(rho_w_range=(10.0, 60.0), m_range=(1.6, 2.2))
unc = MonteCarloHydro(rm, cfg, bounds, n_samples=200).run()

import numpy as np

print("worst-case CV of K:", round(float(np.nanmax(unc.cv_K)), 2))
print(
    "water-table P10-P90 spread (m), first 5:",
    np.round((unc.p90_wt - unc.p10_wt)[:5], 1),
)

# %%
# Uncertainty section
# -------------------
# :class:`~pycsamt.interp.plot.PlotUncertaintySection` pairs the P50 estimate
# (top) with its spread (bottom) for a chosen quantity — so a confident,
# high-``K`` aquifer cell is visibly distinct from a high-``K`` cell that is
# merely a lucky draw.

from pycsamt.interp.plot import (
    PlotUncertaintyHistogram,
    PlotUncertaintyProfile,
    PlotUncertaintySection,
)

PlotUncertaintySection(unc, quantity="K").plot()

# %%
# Water-table and transmissivity with error bands
# -----------------------------------------------
# :class:`~pycsamt.interp.plot.PlotUncertaintyProfile` is the water-table /
# transmissivity profile from :doc:`plot_4_hydro_geophysics`, now with
# P10-P90 envelopes — the honest version of a well-siting figure.

PlotUncertaintyProfile(unc, reference_depth=20.0).plot()

# %%
# Distribution of a target property
# ---------------------------------
# :class:`~pycsamt.interp.plot.PlotUncertaintyHistogram` shows the full
# ensemble distribution of a scalar target (e.g. peak transmissivity),
# making the shape of the uncertainty — skew, multi-modality — explicit
# rather than collapsing it to a single error bar.

PlotUncertaintyHistogram(unc).plot()

# %%
# **Reading it.** Uncertainty is largest where the aquifer is thin or its
# resistivity sits near a unit boundary, and smallest in the thick clean
# sand — so report calibrated P10-P90 ranges, not point values, whenever the
# interpretation feeds a drilling or management decision. This closes the
# interpretation workflow: model -> lithology -> hydrogeology -> properties ->
# calibration -> monitoring -> uncertainty.
