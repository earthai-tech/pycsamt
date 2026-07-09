"""
Station maps and overlays
=========================

The starting point for every :mod:`pycsamt.map` view is a **station map** —
a plan view of the survey with each site plotted at its true coordinates.
The high-level :class:`~pycsamt.map.MapView` façade loads a whole survey
once and renders it many ways; this example builds the map up from bare
locations to a contoured apparent-resistivity surface.

Every figure below is a **live Plotly map**: drag to pan, scroll to zoom,
and hover a station to read its values.
"""

# %%
# Load the survey
# ---------------
# :meth:`MapView.from_folder <pycsamt.map.MapView.from_folder>` reads every
# EDI line under a folder into one geo-referenced view. WILLY_DATA has five
# lines and 128 stations.

import os

from pycsamt.map import MapView

# Card thumbnail: a pre-rendered static PNG (interactive figures below have
# no raster the gallery can thumbnail).
# sphinx_gallery_thumbnail_path = '_static/map_thumbs/plot_1_station_maps.png'

DATA = os.path.join(os.environ.get("PYCSAMT_DOCS_REPO_ROOT", "."),
                    "data", "AMT", "WILLY_DATA")
mv = MapView.from_folder(DATA, recursive=True)
print(f"{mv.n_stations} stations across {len(mv.lines)} lines: {mv.lines}")
print("geo-referenced:", mv.has_geo)

# %%
# Station locations
# -----------------
# The simplest overlay, ``index``, colours each station by its position
# along the line — enough to see the survey geometry. Hover any point for
# its ID and coordinates.

fig = mv.station(overlay="index")
fig.update_layout(height=560)
fig

# %%
# Coloured by apparent resistivity
# --------------------------------
# Switch the overlay to ``rho`` and pick a frequency: now each station is
# coloured by its apparent resistivity at ~100 Hz. Lateral resistivity
# variation across the survey becomes visible at a glance.

fig = mv.station(overlay="rho", component="xy", frequency=100.0)
fig.update_layout(height=560)
fig

# %%
# Add an interpolated contour surface
# -----------------------------------
# ``show_contours=True`` interpolates the station values onto a grid and
# draws filled contours beneath the markers — a continuous resistivity map
# rather than discrete points.

fig = mv.station(overlay="rho", component="xy", frequency=100.0,
                 show_contours=True, contour_levels=14)
fig.update_layout(height=600)
fig

# %%
# Phase instead of resistivity
# ----------------------------
# The same map keyed on impedance ``phase`` highlights different structure —
# phase responds to *gradients* in resistivity, so it often sharpens
# boundaries the resistivity map smooths over.

fig = mv.station(overlay="phase", component="xy", frequency=100.0,
                 show_contours=True)
fig.update_layout(height=600)
fig

# %%
# Frequency is depth: compare two bands
# -------------------------------------
# Higher frequencies sense shallower structure. Rendering the same survey at
# ~1000 Hz (shallow) shows a different pattern than the ~100 Hz map above —
# the basis for reading a station map as a depth-dependent slice.

fig = mv.station(overlay="rho", component="xy", frequency=1000.0,
                 show_contours=True, title="Apparent resistivity ~1000 Hz (shallow)")
fig.update_layout(height=600)
fig

# %%
# **Next.** These plan views collapse all depth into one frequency. The
# following examples open up the third dimension — first as per-line
# :doc:`profiles and pseudo-sections <plot_2_profiles_and_sections>`, then as
# full 3-D :doc:`fence diagrams <plot_3_fence_3d>` and volumes.
