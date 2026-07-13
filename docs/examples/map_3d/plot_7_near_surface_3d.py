"""
Near-surface 3-D view with topography and stations
===================================================

For interpretation, the **top ~1 km** is usually what matters — it is where
AMT resolution is highest and where targets sit. This example builds a
focused 3-D view of just that near-surface window, draped over real
topography and annotated with the survey stations, so structure is read
directly against the ground surface and the sites that constrain it.

Every scene is an interactive Plotly view: orbit to look under the terrain,
hover for resistivity and depth.
"""

# %%
# Load the survey
# ---------------

import os

from pycsamt.map import MapView

# sphinx_gallery_thumbnail_path = '_static/map_thumbs/plot_7_near_surface_3d.png'

DATA = os.path.join(
    os.environ.get("PYCSAMT_DOCS_REPO_ROOT", "."), "data", "AMT", "WILLY_DATA"
)
mv = MapView.from_folder(DATA, recursive=True)
print(f"{mv.n_stations} stations across {len(mv.lines)} lines")

# %%
# The near-surface block
# ----------------------
# ``depth_range=(0, 1000)`` restricts the interpolated volume to the top
# kilometre. With topography on (the default) and ``show_stations=True``,
# the resistivity body hangs beneath the true ground surface with a marker
# at every site.

fig = mv.map3d(
    mode="block", depth_range=(0, 1000), show_stations=True, station_size=4
)
fig.update_layout(height=660, scene_aspectmode="cube")
fig

# %%
# Near-surface fence with stations
# --------------------------------
# The same 1 km window as a fence: each line's shallow cross-section, draped
# on topography, with its stations. This ties the near-surface resistivity
# straight back to the acquisition geometry.

fig = mv.map3d(
    mode="fence", depth_range=(0, 1000), show_stations=True, station_size=4
)
fig.update_layout(height=660, scene_aspectmode="cube")
fig

# %%
# Shallow depth slices, labelled
# ------------------------------
# Depth slices through the top kilometre give a peel-down sequence of maps.
# Adding ``station_labels=True`` writes the site IDs, turning the scene into
# a self-contained near-surface figure.

fig = mv.map3d(
    mode="depth",
    depth_range=(0, 1000),
    n_slices=5,
    show_stations=True,
    station_size=3,
    station_labels=True,
)
fig.update_layout(height=660, scene_aspectmode="cube")
fig

# %%
# A cleaner target view
# ---------------------
# Dropping the terrain surface (``show_terrain=False``) while keeping the
# topographic *shift* leaves an unobstructed view of the shallow body with
# its stations — often the clearest angle for presenting a near-surface
# anomaly.

fig = mv.map3d(
    mode="block",
    depth_range=(0, 1000),
    show_terrain=False,
    show_stations=True,
    station_size=4,
    opacity=0.7,
)
fig.update_layout(height=660, scene_aspectmode="cube")
fig

# %%
# **See also.** To isolate the conductive or resistive parts of this volume,
# continue to :doc:`plot_8_resistivity_bodies`.
