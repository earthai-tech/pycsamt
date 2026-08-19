.. _models:

Classical model integrations
=============================

The :mod:`pycsamt.models` package connects pyCSAMT to three established
electromagnetic modelling and inversion ecosystems: Occam2D, ModEM, and
MARE2DEM. These integrations keep native files and external executables
visible while providing Python tools to build inputs, validate projects, run
solvers when explicitly requested, load results, and create review figures.

This section is the classical engine layer of the
:doc:`../inversion/index` guide. It is separate from
:doc:`../forward/index`, which focuses on predicting synthetic responses, and
from :doc:`../ai_inversion/index`, which covers learned inversion methods.

.. admonition:: External solvers
   :class: important

   pyCSAMT provides integration code and project tooling. Availability of an
   external solver, its executable license, compilation requirements, and
   computational resources remain the user's responsibility. External
   execution is explicit rather than automatic.

.. toctree::
   :maxdepth: 3
   :class: pycsamt-guide-toc

   overview
   choosing_backend
   configuration_and_io
   compilation
   occam2d
   occam1d
   modem
   mare2dem
