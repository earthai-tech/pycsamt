.. _emtools_gradient_imaging:

pycsamt.emtools.gradient_imaging — Gradient-Based CSAMT Pseudo-Sections
============================================================================

:mod:`pycsamt.emtools.gradient_imaging` implements the frequency,
spatial, and joint frequency-spatial gradient apparent-resistivity
formulations of Zhang, Farquharson & Liu (2021, *Geophysical
Prospecting*). Instead of plotting :math:`\rho_a` itself, it plots how
:math:`\rho_a` *changes* — along the line, with frequency/depth, or
both at once — which the paper argues sharpens boundary delineation
and suppresses spurious background variation that a plain
:math:`\rho_a` pseudo-section carries even over a uniform half-space.

Functions
---------

- :func:`~pycsamt.emtools.gradient_imaging.rho_spatial_gradient`
- :func:`~pycsamt.emtools.gradient_imaging.rho_frequency_gradient`
- :func:`~pycsamt.emtools.gradient_imaging.rho_joint_gradient`
- :func:`~pycsamt.emtools.gradient_imaging.plot_gradient_section`

Full signatures and parameter documentation live in the
:doc:`API reference <../api/emtools>`.

Worked Example
--------------

Built on **L18PLT** (``data/AMT/WILLY_DATA/``), the same CSAMT-band
line used in the ``anisotropy``/``csumt``/``diag``/``fieldzone``
examples. It moves from one station pair's spatial gradient and one
station's frequency gradient, to the joint gradient combining both, to
the module's three pseudo-sections side by side (which agree on where
the two main hot spots sit), and finishes with three checks that go
beyond just calling the plotting function: whether the joint
gradient's claimed background suppression actually shows up in this
real line's numbers (yes — about half the spatial gradient's standard
deviation), how much the choice of impedance component changes the
result (a lot — ``yx`` alone is roughly 7x noisier than the default
geometric-mean ``det``), and a same-survey comparison against
neighbouring line L22PLT.

.. include:: ../examples/emtools/plot_gradient_imaging.rst
   :start-after: .. _sphx_glr_examples_emtools_plot_gradient_imaging.py:
   :end-before: .. _sphx_glr_download_examples_emtools_plot_gradient_imaging.py:
