r"""
The resistivity model
=====================

Everything in :mod:`pycsamt.interp` starts from a
:class:`~pycsamt.interp.ResistivityModel` — a 2-D section of
:math:`\log_{10}\rho` on a ``(depth, distance)`` grid, plus station
positions and metadata. It usually comes from an inversion
(:meth:`ResistivityModel.from_occam2d
<pycsamt.interp._base.ResistivityModel.from_occam2d>`), but you can build one
directly from arrays, which is what this whole gallery does so every
deliverable traces back to a known model.

This first example builds the shared synthetic section and reads it as
1-D resistivity-depth profiles.
"""

# %%
# Build a model from arrays
# -------------------------
# :meth:`ResistivityModel.from_array
# <pycsamt.interp._base.ResistivityModel.from_array>` takes a
# ``(n_z, n_x)`` array of :math:`\log_{10}\rho` with cell-centre depth and
# distance axes. Here the section is a classic hydrogeological sequence: a
# resistive dry overburden, a conductive sand aquifer, a very conductive
# clay aquitard, and a resistive basement.

import numpy as np
import matplotlib.pyplot as plt

from _interp_data import demo_model

# Use the resistivity section (2nd figure) as the card thumbnail.
# sphinx_gallery_thumbnail_number = 2

rm = demo_model()
print(rm)
print("grid (n_z, n_x):", rm.rho_2d.shape,
      " depth 0-{:.0f} m,  distance 0-{:.0f} m".format(
          rm.z_centers.max(), rm.x_centers.max()))

# %%
# The section itself
# ------------------
# A quick look at the model with an EM-style colour scale: cool colours are
# conductive (water, clay), warm colours resistive (dry ground, basement).

fig, ax = plt.subplots(figsize=(10, 4.2), constrained_layout=True)
im = ax.pcolormesh(rm.x_centers, rm.z_centers, rm.rho_2d,
                   cmap="Spectral", shading="auto")
ax.plot(rm.station_x, np.full_like(rm.station_x, rm.z_centers.min()),
        "kv", ms=6, clip_on=False)
ax.invert_yaxis()
ax.set_xlabel("distance (m)")
ax.set_ylabel("depth (m)")
ax.set_title("Synthetic resistivity section (station markers on top)")
fig.colorbar(im, ax=ax, label=r"$\log_{10}\rho$  ($\Omega\cdot$m)")

# %%
# 1-D resistivity-depth profiles
# ------------------------------
# :class:`~pycsamt.interp.plot.PlotResistivityDepthProfile` extracts the
# sounding beneath a station and plots resistivity against depth — the most
# direct way to pick layer boundaries. Passing a ``ResistivityModel`` shows
# the raw log; later examples pass an ``EMHydroResult`` to add zone shading.

from pycsamt.interp.plot import PlotResistivityDepthProfile

PlotResistivityDepthProfile(rm).plot()

# %%
# **Reading it.** The profile drops sharply into the aquifer (~40 Ohm-m near
# 20 m), bottoms out in the clay (~12 Ohm-m), then rises steeply into the
# resistive basement (>1000 Ohm-m). Those four segments are exactly what the
# :doc:`lithology <plot_2_lithology>` and :doc:`hydro <plot_3_hydro_interpretation>`
# examples turn into named units and aquifer zones.
