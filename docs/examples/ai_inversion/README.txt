.. _ai_inversion_gallery:

AI inversion
------------

Deep-learning inversion with :mod:`pycsamt.ai` — training neural networks
to map electromagnetic responses **directly** to subsurface resistivity,
a fast surrogate for the classical iterative inverters. The gallery builds
from the ground up, 1-D to 3-D and beyond:

* **Train a 1-D inverter** — generate a synthetic training set with the
  forward solver, fit a CNN, and validate predicted layer models;
* **Compare architectures** — CNN vs ResNet vs FCN on the same data, as a
  grid of convergence, residual, and per-layer error panels;
* **Uncertainty and calibration** — deep ensembles, Monte-Carlo dropout,
  and conformal prediction intervals with coverage diagnostics;
* **2-D section inversion** — a U-Net that predicts a full resistivity
  cross-section from a pseudo-section of responses;
* **3-D graph-network inversion** — a graph convolutional network that
  inverts a whole multi-line survey at once, using station geometry as
  graph context.

Every example trains a *real* network on synthetic data during the
documentation build, kept small so each page runs in seconds. On your own
data you would use larger training sets and more epochs, but the API and
the validation figures are exactly the same. See the :ref:`AI-inversion
user guide <user_guide_ai_inversion>` for the narrative reference.

.. note::

   These examples require a deep-learning backend (PyTorch or TensorFlow).
   Install one with ``pip install pycsamt[torch]``.
