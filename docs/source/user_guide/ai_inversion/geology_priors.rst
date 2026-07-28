.. _ai_inversion_geology_priors:

Correlated geological priors
============================

.. admonition:: Status
   :class: note

   Implemented. Statistical recovery of the requested correlation
   ranges is verified by unit tests, but the M2 gate's broader
   diagnostic gallery (histograms, variograms, class/feature coverage,
   depth occupancy across many realizations) is not yet assembled as
   a page or report.

Training a network on tiled 1-D models teaches it a per-station
relationship, not a spatial one: nothing in that training data ever
shows the network what a laterally continuous resistivity structure
looks like.  :mod:`pycsamt.ai.geology` exists to generate complete
2-D sections and 3-D volumes directly, so response simulation and
training happen on genuine geology instead of on independent columns
stitched together after the fact.

The package builds on a shared
:class:`~pycsamt.ai.geology.GeologyGrid` (canonical ``(z, x)`` or
``(z, y, x)`` cell-centre grid) and composes several generators on
top of it: anisotropic correlated Gaussian fields
(:class:`~pycsamt.ai.geology.CorrelatedField`, via
:func:`~pycsamt.ai.geology.generate_gaussian_field`), layered
geology with dipping interfaces
(:func:`~pycsamt.ai.geology.generate_layered_geology`), ellipsoidal
lenses and bodies (:func:`~pycsamt.ai.geology.insert_lenses`), and
topographic surfaces
(:class:`~pycsamt.ai.geology.TopographicSurface`).  Every generator
takes an explicit random seed and exposes its generation provenance,
so a training dataset can always be traced back to exactly the
configuration that produced it.

.. code-block:: pycon

   >>> import numpy as np
   >>> from pycsamt.ai.geology import GeologyGrid, GaussianCorrelation
   >>> from pycsamt.ai.geology import generate_gaussian_field
   >>> grid = GeologyGrid.regular_2d(nx=32, nz=16, dx_m=100, dz_m=50)
   >>> model = GaussianCorrelation(length_x_m=500, length_z_m=100)
   >>> field = generate_gaussian_field(grid, model, seed=12)
   >>> field.values.shape
   (16, 32)
   >>> round(float(np.mean(field.values)), 6)
   0.0
   >>> round(float(np.std(field.values)), 6)
   1.0

:func:`~pycsamt.ai.geology.directional_variogram` gives an empirical
semivariogram along one axis, used to check that a generated field's
realized correlation length actually matches the requested
:class:`~pycsamt.ai.geology.GaussianCorrelation`, rather than trusting
the generator's parameters at face value.

Planned expansion
-----------------

A gallery of 2-D and 3-D realizations spanning the generator's
parameter space, empirical variogram fits confirming recovered
correlation lengths against requested ones, and a layered-geology
plus lens/fault composition example built up step by step.
