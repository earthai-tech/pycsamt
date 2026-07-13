"""
3-D maps with topography
========================

Real surveys are not flat. :mod:`pycsamt.map` can drape the 3-D views over
**surface elevation**, so cross-sections and volumes hang beneath the true
terrain rather than a flat datum — essential for shallow targets where a few
tens of metres of relief matter. WILLY_DATA carries per-station elevations
(37–224 m), so topography here is real, not synthetic.

The scenes are interactive: orbit to see the terrain surface and the survey
draped beneath it.
"""

# %%
# Load the survey
# ---------------

import os

import numpy as np

from pycsamt.map import MapView

# sphinx_gallery_thumbnail_path = '_static/map_thumbs/plot_5_topography_3d.png'

DATA = os.path.join(
    os.environ.get("PYCSAMT_DOCS_REPO_ROOT", "."), "data", "AMT", "WILLY_DATA"
)
mv = MapView.from_folder(DATA, recursive=True)
elev = mv.table()["Elevation"].to_numpy(dtype=float)
print(f"elevation range: {np.nanmin(elev):.0f}–{np.nanmax(elev):.0f} m")

# %%
# Baseline: a flat datum
# ----------------------
# Turning topography *off* pins every panel to a flat z=0 datum — the
# reference to compare against. This is the pure cross-section geometry,
# ignoring surface relief.

fig = mv.map3d(mode="fence", topography=False, show_terrain=False)
fig.update_layout(height=640, scene_aspectmode="cube")
fig

# %%
# Draped over real terrain
# ------------------------
# Topography is *on by default*: each panel is shifted to sit beneath its
# station's true elevation, and ``show_terrain`` adds the interpolated ground
# surface. Compare with the flat fence above — the panels now follow the
# relief.

fig = mv.map3d(mode="fence")  # topography=True is the default
fig.update_layout(height=660, scene_aspectmode="cube")
fig

# %%
# Volume beneath the terrain
# --------------------------
# The block view honours topography too — the resistivity body hangs under
# the real ground surface, so depths read from the surface down. A lighter
# terrain lets the volume show through.

fig = mv.map3d(
    mode="block", terrain_opacity=0.5, show_stations=True, station_size=3
)
fig.update_layout(height=660, scene_aspectmode="cube")
fig

# %%
# Supplying your own elevations
# -----------------------------
# When EDI headers lack elevation (or you have a better DEM), attach one with
# :meth:`MapView.with_elevations <pycsamt.map.MapView.with_elevations>` —
# a ``{station_id: elevation_m}`` mapping. :meth:`MapView.fetch_elevations`
# can also pull a DEM online. Here we exaggerate the relief to show the
# effect.

boosted = {
    s: float(200 + 3 * (e - np.nanmean(elev)))
    for s, e in zip(mv.stations, elev)
}
fig = mv.with_elevations(boosted).map3d(mode="fence")
fig.update_layout(height=660, scene_aspectmode="cube")
fig

# %%
# **Next.** The :doc:`advanced compositions <plot_6_advanced_compositions>`
# example tunes colour, slicing, iso-values, and the camera for
# presentation-quality figures.
