r"""
Hydro-geophysics: aquifer properties
====================================

:class:`~pycsamt.interp.EMHydroModel` is the quantitative step: it runs a
petrophysical model (Archie by default) over the resistivity section to
estimate **porosity, water saturation, hydraulic conductivity, the water
table, and transmissivity** — the numbers a hydrogeologist actually needs.
The result is an :class:`~pycsamt.interp.EMHydroResult`, which the
:mod:`pycsamt.interp.plot` classes turn into finished figures.
"""

# %%
# Run the model
# -------------
# A :class:`~pycsamt.interp.PetrophysicalConfig` sets the rock physics
# (pore-water resistivity, porosity prior, cementation). ``fit`` returns the
# result object with one map/profile per property.

from _interp_data import demo_model

from pycsamt.interp import EMHydroModel, PetrophysicalConfig

# Use the hydraulic-conductivity section (2nd figure) as the thumbnail.
# sphinx_gallery_thumbnail_number = 2

rm = demo_model()
cfg = PetrophysicalConfig(rho_w=20.0, porosity_prior=0.25)
result = EMHydroModel(rm, cfg, method_tag="AMT").fit()

import numpy as np

print(
    "porosity range:",
    np.round([result.porosity.min(), result.porosity.max()], 2),
)
print(
    "K range (m/s):",
    np.format_float_scientific(np.nanmin(result.hydraulic_K), 1),
    "-",
    np.format_float_scientific(np.nanmax(result.hydraulic_K), 1),
)
print(
    "water table (m), first 5 stations:", np.round(result.water_table[:5], 1)
)

# %%
# Hydraulic-conductivity section
# ------------------------------
# :class:`~pycsamt.interp.plot.PlotHydroSection` images any quantitative map.
# Hydraulic conductivity ``K`` is the headline aquifer property — high in the
# clean sand aquifer, low in the clay.

from pycsamt.interp.plot import (
    PlotAquiferCharacterization,
    PlotHydroSection,
    PlotWaterTableProfile,
)

PlotHydroSection(result, quantity="K").plot()

# %%
# Water-saturation section
# ------------------------
# The same view for saturation ``Sw`` separates the saturated aquifer from
# the drained overburden above the water table.

PlotHydroSection(result, quantity="saturation").plot()

# %%
# Water table and transmissivity
# ------------------------------
# :class:`~pycsamt.interp.plot.PlotWaterTableProfile` reduces the section to
# the two most-requested profiles along the line: the water-table depth and
# the aquifer transmissivity (:math:`T`, on a log scale) — the deliverable
# for well-siting.

PlotWaterTableProfile(result, reference_depth=20.0).plot()

# %%
# Dar-Zarrouk aquifer characterization
# ------------------------------------
# :class:`~pycsamt.interp.plot.PlotAquiferCharacterization` stacks the
# classic Dar-Zarrouk parameters — transverse resistance ``TR`` and
# longitudinal conductance ``S`` — with the water table and transmissivity,
# a compact one-figure aquifer summary.

PlotAquiferCharacterization(result).plot()

# %%
# **Reading it.** ``K`` and ``T`` both peak where the sand aquifer is
# thickest and cleanest, and collapse over the clay — so the best well
# targets are the high-``T`` segments of the profile. The
# :doc:`petrophysics <plot_5_petrophysics>` example shows the rock-physics
# model behind these numbers, and :doc:`uncertainty <plot_8_uncertainty>`
# puts error bars on them.
