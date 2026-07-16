.. _forward:

Forward Modelling
==================

:mod:`pycsamt.forward` contains the in-process forward modelling layer of
pyCSAMT v2. It provides resistivity model containers, 1-D/2-D/quasi-3-D
electromagnetic solvers, synthetic dataset generators, noise models, and
diagnostic plotting utilities.

Forward modelling answers a controlled physical question:

   If this earth model were true, what response would the survey record?

That question sits underneath every serious inversion workflow. Inversion
searches for a model that explains observed data; forward modelling computes
the data predicted by a known model. A good inversion setup therefore depends
on forward modelling for sensitivity tests, synthetic recovery experiments,
training data, uncertainty design, and quality control.

.. toctree::
   :maxdepth: 3
   :class: pycsamt-guide-toc

   overview
   concepts
   configuration
   solvers_and_grids
   synthetic_datasets
   plotting
   forward_to_inversion
