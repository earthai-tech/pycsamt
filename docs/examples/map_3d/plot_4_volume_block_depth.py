"""
3-D blocks, depth slices, and surfaces
======================================

Where the fence diagram shows only the measured lines, the volumetric modes
interpolate *between* lines to fill a continuous 3-D resistivity body, then
render it three ways: a solid **block**, a stack of horizontal **depth**
slices, and an iso-**surface**. Each is the same underlying volume seen
through a different lens, and each answers a different question.

All scenes are interactive — orbit, zoom, and hover to probe the volume.
"""

# %%
# Load the survey
# ---------------

import os

from pycsamt.map import MapView

# sphinx_gallery_thumbnail_path = '_static/map_thumbs/plot_4_volume_block_depth.png'

DATA = os.path.join(
    os.environ.get("PYCSAMT_DOCS_REPO_ROOT", "."), "data", "AMT", "WILLY_DATA"
)
mv = MapView.from_folder(DATA, recursive=True)
print(f"{mv.n_stations} stations across {len(mv.lines)} lines")

# %%
# The resistivity block
# ---------------------
# ``mode="block"`` renders the interpolated volume as a semi-transparent
# Plotly ``Volume`` — the most immediate "solid earth" view. Opacity lets
# you see conductive/resistive bodies suspended inside.

fig = mv.map3d(mode="block")
fig.update_layout(height=640, scene_aspectmode="cube")
fig

# %%
# Depth slices
# ------------
# ``mode="depth"`` cuts the volume into a stack of horizontal
# constant-depth maps. Reading top to bottom is like peeling the survey
# downward one resistivity map at a time.

fig = mv.map3d(mode="depth", n_slices=6)
fig.update_layout(height=640, scene_aspectmode="cube")
fig

# %%
# Iso-surface
# -----------
# ``mode="surface"`` draws iso-resistivity surfaces — the 3-D contour of a
# chosen resistivity value, which isolates the shape of a conductor or a
# resistive basement.

fig = mv.map3d(mode="surface")
fig.update_layout(height=640, scene_aspectmode="cube")
fig

# %%
# The block, with stations
# ------------------------
# Overlaying the station markers ties the interpolated body back to where
# the data actually constrain it — a reminder that resistivity *between*
# lines is inferred, not measured.

fig = mv.map3d(mode="block", show_stations=True, station_size=3)
fig.update_layout(height=640, scene_aspectmode="cube")
fig

# %%
# **Next.** So far the volume floats at nominal depth. The
# :doc:`topography example <plot_5_topography_3d>` drapes it over real
# surface elevation.
