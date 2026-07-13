"""
Profiles and pseudo-sections
============================

Between the flat station map and a full 3-D volume sit the **per-line**
views: a profile of resistivity/phase along each survey line, and a
pseudo-section that images one line as station (x) against period (a proxy
for depth, y). These are the interpretive workhorses — the first place
lateral and vertical structure show up together.

All figures are interactive Plotly panels; hover for values, zoom into a
band of interest.
"""

# %%
# Load the survey
# ---------------

import os

from pycsamt.map import MapView

# sphinx_gallery_thumbnail_path = '_static/map_thumbs/plot_2_profiles_and_sections.png'

DATA = os.path.join(
    os.environ.get("PYCSAMT_DOCS_REPO_ROOT", "."), "data", "AMT", "WILLY_DATA"
)
mv = MapView.from_folder(DATA, recursive=True)
print(f"{mv.n_stations} stations, lines: {mv.lines}")

# %%
# Line profiles
# -------------
# :meth:`MapView.profile <pycsamt.map.MapView.profile>` draws apparent
# resistivity and phase for every station, grouped by line — the quickest
# read on how each line behaves across the frequency band.

fig = mv.profile()
fig.update_layout(height=560)
fig

# %%
# Resistivity pseudo-section
# --------------------------
# :meth:`MapView.pseudosection <pycsamt.map.MapView.pseudosection>` images
# apparent resistivity as station vs period. High frequencies (top) sense
# shallow structure, low frequencies (bottom) sense deep — so the vertical
# axis reads as increasing depth.

fig = mv.pseudosection(quantity="rho")
fig.update_layout(height=520)
fig

# %%
# Phase pseudo-section
# --------------------
# The phase pseudo-section of the same line complements resistivity: phase
# tracks the *rate of change* of resistivity with depth, so a phase high
# often flags a conductor the resistivity image only hints at.

fig = mv.pseudosection(quantity="phase")
fig.update_layout(height=520)
fig

# %%
# Distance along the profile
# --------------------------
# Setting ``x_axis="distance"`` spaces stations by their true separation
# (metres) instead of by index — important when station spacing is uneven,
# so lateral features are not distorted.

fig = mv.pseudosection(quantity="rho", x_axis="distance")
fig.update_layout(height=520)
fig

# %%
# **Next.** A pseudo-section flattens one line. The
# :doc:`fence diagram <plot_3_fence_3d>` hangs every line's section in a
# single 3-D scene at its true map position.
