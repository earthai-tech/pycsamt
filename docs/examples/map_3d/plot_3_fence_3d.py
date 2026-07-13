"""
3-D fence diagrams
==================

The **fence diagram** is the natural first step into 3-D: each survey line
is drawn as a vertical cross-section ("fence panel") hanging at its true map
position, so you see all five WILLY_DATA lines together in one rotatable
scene. It is the most faithful 3-D view of raw survey data — no volumetric
interpolation between lines, just the real per-line pseudo-sections placed
in space.

Drag to orbit the scene, scroll to zoom, and hover any panel for its
resistivity and depth.
"""

# %%
# Load the survey
# ---------------

import os

from pycsamt.map import MapView

# sphinx_gallery_thumbnail_path = '_static/map_thumbs/plot_3_fence_3d.png'

DATA = os.path.join(
    os.environ.get("PYCSAMT_DOCS_REPO_ROOT", "."), "data", "AMT", "WILLY_DATA"
)
mv = MapView.from_folder(DATA, recursive=True)
print(f"{mv.n_stations} stations across {len(mv.lines)} lines")

# %%
# The fence
# ---------
# :meth:`MapView.map3d <pycsamt.map.MapView.map3d>` with ``mode="fence"``
# builds one panel per line. Colour is apparent resistivity; the vertical
# axis is pseudo-depth from the period-to-depth transform.

fig = mv.map3d(mode="fence")
fig.update_layout(height=640, scene_aspectmode="cube")
fig

# %%
# Add the stations
# ----------------
# ``show_stations=True`` drops a marker at every station on top of the
# panels, tying the 3-D fence back to the plan-view survey geometry.

fig = mv.map3d(mode="fence", show_stations=True, station_size=3)
fig.update_layout(height=640, scene_aspectmode="cube")
fig

# %%
# Phase fence
# -----------
# The same fence keyed on ``phase`` rather than resistivity — useful for
# cross-checking conductor positions between the two quantities in 3-D.

fig = mv.map3d(mode="fence", quantity="phase")
fig.update_layout(height=640, scene_aspectmode="cube")
fig

# %%
# Constrain the colour range
# --------------------------
# Passing an explicit ``rho_range`` fixes the colour scale across a survey,
# so several fences (or several surveys) are directly comparable rather than
# each auto-scaling to its own extremes.

fig = mv.map3d(mode="fence", rho_range=(1.0, 1000.0))
fig.update_layout(height=640, scene_aspectmode="cube")
fig

# %%
# **Next.** A fence shows the measured lines only. The
# :doc:`volume, block, and depth views <plot_4_volume_block_depth>`
# interpolate *between* the lines to fill a continuous 3-D resistivity body.
