"""
Isolating resistive, conductive, and background bodies
======================================================

A single resistivity volume mixes conductors, resistors, and the background
medium into one cloud. The most useful interpretive trick is to **peel them
apart**: render only the cells in a chosen resistivity window, so a
conductor or a resistive block stands alone in space. ``rho_range`` does
exactly this — it sets the volume's ``isomin``/``isomax``, hiding everything
outside the band.

Here the WILLY_DATA volume (ρ ≈ 4 – 10⁵ Ω·m) is split into conductive
(< 100 Ω·m), background (100 – 1000 Ω·m), and resistive (> 1000 Ω·m) bodies,
all draped on topography. Every scene is interactive.
"""

# %%
# Load the survey
# ---------------

import os

from pycsamt.map import MapView

# sphinx_gallery_thumbnail_path = '_static/map_thumbs/plot_8_resistivity_bodies.png'

DATA = os.path.join(
    os.environ.get("PYCSAMT_DOCS_REPO_ROOT", "."), "data", "AMT", "WILLY_DATA"
)
mv = MapView.from_folder(DATA, recursive=True)
print(f"{mv.n_stations} stations across {len(mv.lines)} lines")

# %%
# The 2-D slices through the volume
# ---------------------------------
# Before isolating bodies, the ``depth`` mode shows the volume as a stack of
# horizontal 2-D maps — a reference for where the conductors and resistors
# sit at each level.

fig = mv.map3d(mode="depth", n_slices=6)
fig.update_layout(height=640, scene_aspectmode="cube")
fig

# %%
# The full resistivity block
# --------------------------
# The complete volume, draped on topography — every resistivity class
# present at once. The examples below carve this apart.

fig = mv.map3d(mode="block", opacity=0.5, show_stations=True, station_size=3)
fig.update_layout(height=640, scene_aspectmode="cube")
fig

# %%
# Conductive block only
# ---------------------
# ``rho_range=(1, 100)`` keeps only cells below 100 Ω·m — the conductors
# (clay, alteration, fluids) isolated from everything else, so their shape
# and depth are unobstructed.

fig = mv.map3d(
    mode="block",
    rho_range=(1.0, 100.0),
    opacity=0.75,
    show_stations=True,
    station_size=3,
)
fig.update_layout(height=640, scene_aspectmode="cube")
fig

# %%
# Resistive block only
# --------------------
# The complementary view: ``rho_range=(1000, 100000)`` shows only the
# resistive bodies (fresh basement, intrusions), which here form the bulk of
# the deeper survey.

fig = mv.map3d(
    mode="block",
    rho_range=(1000.0, 100_000.0),
    opacity=0.55,
    show_stations=True,
    station_size=3,
)
fig.update_layout(height=640, scene_aspectmode="cube")
fig

# %%
# Background medium only
# ----------------------
# Isolating the intermediate band ``rho_range=(100, 1000)`` leaves the
# "host" medium — everything that is neither a marked conductor nor a strong
# resistor. Comparing the three isolations shows how the survey volume
# partitions.

fig = mv.map3d(
    mode="block",
    rho_range=(100.0, 1000.0),
    opacity=0.5,
    show_stations=True,
    station_size=3,
)
fig.update_layout(height=640, scene_aspectmode="cube")
fig

# %%
# Iso-surface of a single class
# -----------------------------
# For a crisp boundary rather than a semi-transparent cloud, ``mode="surface"``
# with the same ``rho_range`` draws the iso-surface enclosing the conductive
# body — the cleanest way to communicate a target's geometry.

fig = mv.map3d(
    mode="surface",
    rho_range=(1.0, 100.0),
    surface_count=6,
    opacity=0.6,
    show_stations=True,
    station_size=3,
)
fig.update_layout(height=640, scene_aspectmode="cube")
fig

# %%
# **Takeaway.** ``rho_range`` turns one volume into a set of
# resistivity-class views — conductor, background, resistor — each draped on
# topography and referenced to the stations. Combine with
# :doc:`plot_7_near_surface_3d` to isolate a class *within* the near-surface
# window.
