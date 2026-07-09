"""
Advanced 3-D compositions
=========================

The final step: fine control over how the volume is rendered and viewed.
Every 3-D mode accepts the full :class:`~pycsamt.map.VolumeMapOptions`
vocabulary — colour map and range, cell-filtering iso-window, slice and
surface density, opacity, station overlays, colour clipping, vertical
exaggeration, and camera azimuth. This example composes them into
presentation-quality scenes.

.. admonition:: Depth is capped to 2 km for visibility
   :class: important

   The interpolated pseudo-depth volume reaches roughly **30 km**, so the
   shallow structure that actually matters renders as a vanishingly thin
   sliver at the top of a tall, empty box. **Every scene on this page caps
   the view to the top 2 km** with ``depth_range=(0, 2000)`` (set once as
   ``DEPTH_CAP`` below). Change that single value — e.g. ``(0, 1000)`` for a
   1 km near-surface view, or ``(0, 5000)`` to see deeper — to reframe every
   figure. Pairing a tight ``depth_range`` with ``aspectmode="cube"`` is the
   single most effective way to get a legible 3-D resistivity scene.

The scenes are live; use them to find the viewpoint that best tells the
story, then export with :meth:`MapView.export`.
"""

# %%
# Load the survey and set the depth cap
# -------------------------------------
# ``DEPTH_CAP`` is threaded through every call below, so the whole page
# re-frames from one place.

import os

from pycsamt.map import MapView

# sphinx_gallery_thumbnail_path = '_static/map_thumbs/plot_6_advanced_compositions.png'

DEPTH_CAP = (0.0, 2000.0)          # top 2 km — see the admonition above

DATA = os.path.join(os.environ.get("PYCSAMT_DOCS_REPO_ROOT", "."),
                    "data", "AMT", "WILLY_DATA")
mv = MapView.from_folder(DATA, recursive=True)
print(f"{mv.n_stations} stations across {len(mv.lines)} lines")
print(f"depth cap: {DEPTH_CAP[1]:.0f} m")

# %%
# 1. A custom-styled block
# ------------------------
# Override the colour map, fix the resistivity range for reproducibility,
# and raise the opacity to firm up the body. ``log_color=True`` spreads the
# scale over resistivity's many decades. Capped to the top 2 km, the shallow
# structure fills the frame instead of hiding at the top of a 30 km box.

fig = mv.map3d(mode="block", depth_range=DEPTH_CAP, cmap="Turbo",
               rho_range=(1.0, 5000.0), log_color=True, opacity=0.55,
               show_stations=True, station_size=3)
fig.update_layout(height=660, scene_aspectmode="cube")
fig

# %%
# 2. Sharpen contrast with colour clipping
# ----------------------------------------
# ``rho_range`` chooses *which cells render*; ``value_range`` instead clips
# the *colour scale* (``cmin``/``cmax`` in log10 ρ). Narrowing it saturates
# the extremes and throws the mid-range resistivity contrast into relief —
# useful when one very resistive or very conductive cell would otherwise
# flatten the whole colour bar.

fig = mv.map3d(mode="block", depth_range=DEPTH_CAP, cmap="RdYlBu_r",
               value_range=(1.5, 3.3), opacity=0.6, show_stations=True,
               station_size=3)
fig.update_layout(height=660, scene_aspectmode="cube")
fig

# %%
# 3. A dense depth stack
# ----------------------
# More slices give a finer sense of the vertical resistivity gradient — ten
# constant-depth maps through the top 2 km, each lightly transparent so
# deeper slices show through.

fig = mv.map3d(mode="depth", depth_range=DEPTH_CAP, n_slices=10,
               opacity=0.8, cmap="RdYlBu_r")
fig.update_layout(height=660, scene_aspectmode="cube")
fig

# %%
# 4. Focused iso-surfaces
# -----------------------
# Narrowing ``rho_range`` and raising ``surface_count`` isolates a specific
# resistivity band — here the transition into the conductive zone — and
# renders it as several nested shells with crisp boundaries.

fig = mv.map3d(mode="surface", depth_range=DEPTH_CAP, surface_count=14,
               rho_range=(1.0, 800.0), opacity=0.45, show_stations=True,
               station_size=3)
fig.update_layout(height=660, scene_aspectmode="cube")
fig

# %%
# 5. A phase volume
# -----------------
# Every mode also renders ``quantity="phase"``. Because phase responds to
# resistivity *gradients* rather than absolute values, a phase block often
# outlines boundaries — the top of a conductor, a fault contact — more
# sharply than the resistivity block does.

fig = mv.map3d(mode="block", depth_range=DEPTH_CAP, quantity="phase",
               cmap="Viridis", opacity=0.6, show_stations=True,
               station_size=3)
fig.update_layout(height=660, scene_aspectmode="cube")
fig

# %%
# 6. Camera and exaggeration
# --------------------------
# ``azimuth`` rotates the scene and ``aspectmode="cube"`` exaggerates the
# (shallow) depth axis so thin structure stays legible. With station labels
# this is the fence figure to export for a report.

fig = mv.map3d(mode="fence", depth_range=DEPTH_CAP, azimuth=45.0,
               aspectmode="cube", show_stations=True, station_labels=True)
fig.update_layout(height=680)
fig

# %%
# 7. A presentation composite
# ---------------------------
# Everything together: the top 2 km, draped on topography, with a firm
# colour block, station markers, and a chosen camera — the kind of scene
# you would drop straight into a report or slide.

fig = mv.map3d(mode="block", depth_range=DEPTH_CAP, cmap="Turbo",
               log_color=True, opacity=0.65, topography=True,
               show_terrain=True, terrain_opacity=0.5,
               show_stations=True, station_size=3, azimuth=25.0)
fig.update_layout(height=700, scene_aspectmode="cube")
fig

# %%
# Exporting
# ---------
# Any of these figures round-trips to a standalone interactive HTML file (or
# a static image) with one call — the same object you see above::
#
#     mv.export("survey_volume.html", view="map3d", mode="block",
#               depth_range=(0, 2000))
#
# **That completes the tour** — from a flat station map through profiles,
# fences, volumes, topography, near-surface and body-isolation views, to
# fully styled 3-D compositions. See the :doc:`Map tools guide
# </map/index>` for the complete API and the Dash workbench
# (:func:`pycsamt.map.launch_mapview`) for a live, no-code version.
