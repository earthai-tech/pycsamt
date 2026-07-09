r"""
Plot styles and axis control
============================

Every :mod:`pycsamt.forward` figure is drawn through the library's shared
styling layer, so the *same* forward result can be re-rendered for
different audiences without touching the plotting code. This short example
shows the built-in ``publication`` and ``dark`` styles applied to a 1-D
sounding, and how :obj:`~pycsamt.api.control.PYCSAMT_CONTROL` switches the
frequency axis between log-period and linear period.

.. note::

   :func:`~pycsamt.api.style.use_style` and
   :func:`~pycsamt.api.control.configure_control` set **global** state.
   Always pair them with :func:`~pycsamt.api.style.reset_style` (and
   restore the control) so later figures are unaffected — done here with a
   ``try/finally`` block.
"""

# %%
# Setup
# -----
# One sedimentary and one conductive-layer model, each with its MT
# response, are enough to demonstrate every style.

import numpy as np
import matplotlib.pyplot as plt

from pycsamt.forward import (
    LayeredModel, MT1DForward, plot_model_1d, plot_response_1d,
    plot_response_and_model_1d,
)
from pycsamt.api.style import use_style, reset_style
from pycsamt.api.control import configure_control

# Use the dark-style figure (2nd) as the section thumbnail.
# sphinx_gallery_thumbnail_number = 2

M_SEDIMENTARY = LayeredModel([1_000., 20., 5., 300.], [200., 600., 1_500.],
                             name="sedimentary")
M_CONDUCTIVE = LayeredModel([200., 5., 400., 100.], [150., 500., 2_000.],
                            name="conductive-layer")
M_GEOTHERMAL = LayeredModel([500., 8., 250., 3_000.], [100., 400., 2_500.],
                            name="geothermal")

FREQS_MT = np.logspace(-3, 4, 35)
R_SED = MT1DForward(FREQS_MT).run(M_SEDIMENTARY)
R_COND = MT1DForward(FREQS_MT).run(M_CONDUCTIVE)
R_GEO = MT1DForward(FREQS_MT).run(M_GEOTHERMAL)

# %%
# 1. Publication style
# --------------------
# :func:`~pycsamt.api.style.use_style` swaps in a high-contrast, print-ready
# palette. We compose the depth profile and the two response panels into a
# single row using the ``ax=``/``axes=`` arguments.

use_style("publication")
try:
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.5), constrained_layout=True)
    plot_model_1d(M_SEDIMENTARY, ax=axs[0], title="Depth profile")
    plot_response_1d(R_SED, axes=axs[1:3])
    fig.suptitle("Publication style - sedimentary model", y=1.03, fontsize=11)
finally:
    reset_style()

# %%
# 2. Dark style
# -------------
# The ``dark`` style targets slides and screens. The figure/axes
# backgrounds are set explicitly so the whole canvas — not just the plot
# area — picks up the dark theme.

use_style("dark")
try:
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.5), constrained_layout=True,
                            facecolor="#1a1a2e")
    for ax in axs:
        ax.set_facecolor("#1a1a2e")
    plot_model_1d(M_CONDUCTIVE, ax=axs[0], title="Depth profile")
    plot_response_1d(R_COND, axes=axs[1:3])
    fig.suptitle("Dark style - conductive-layer model", y=1.03,
                 fontsize=11, color="white")
finally:
    reset_style()

# %%
# 3. Switching the frequency axis
# -------------------------------
# :func:`~pycsamt.api.control.configure_control` changes how every figure
# maps frequency to the x-axis. Here we switch from the default
# ``log10_period`` to a linear ``period`` axis, then restore it so the
# setting does not leak into other examples.

configure_control(x__view="period")
try:
    fig = plot_response_and_model_1d(
        R_GEO, M_GEOTHERMAL,
        title="Period axis (linear scale) - geothermal model",
        figsize=(11, 4.5),
    )
finally:
    configure_control(x__view="log10_period")   # restore default
