.. _user_guide_inversion:

Inversion
=========

Inversion turns processed electromagnetic observations into a resistivity
model. In pyCSAMT this guide is the decision point between three related but
different paths: the backend-neutral :mod:`pycsamt.inversion` API for
built-in and adapter workflows, the classical engine integrations in
:mod:`pycsamt.models` for Occam2D, ModEM, and MARE2DEM native projects, and
the learned workflows in :mod:`pycsamt.ai`, documented separately as
:doc:`../ai_inversion/index`.

Forward modelling is not an inversion path. It predicts responses from a
known model and is documented in :doc:`../forward/index`; it is often used to
design, test, or validate inversion choices. The common Python interface
itself lives at :doc:`../../api/inversion`.

.. toctree::
   :maxdepth: 3
   :class: pycsamt-guide-toc

   overview
   ../models/index
   ../ai_inversion/index
