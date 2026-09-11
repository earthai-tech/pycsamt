.. _emtools_gallery:

EM Tools Guide
==============

``pycsamt.emtools`` is the science-facing toolbox for MT, AMT, and CSAMT
processing. Use this page as a task map: start with loading and inspection,
then move through quality control, frequency editing, noise/static-shift
conditioning, tensor/strike analysis, source diagnostics, survey design, and
publication plots. Each page opens a narrative guide with examples, figures,
and links to the public API. The same workflows are also available as
runnable gallery examples in :ref:`EM tools examples <emtools_examples>`. For
the complete callable reference, see :doc:`../../api/emtools`.

Every rule-based tool here has a trained, machine-learning counterpart in
:mod:`pycsamt.ai.processing`, documented separately as
:doc:`../ai_processing/index`: an :term:`isolation forest` complementing
:doc:`qc`'s confidence ratio, an :term:`autoencoder` complementing
:doc:`remove_noise`'s filters, and a
:term:`multi-layer perceptron` complementing
:doc:`dimensionality`'s threshold rule.

.. toctree::
   :maxdepth: 3
   :class: pycsamt-guide-toc

   inspect
   qc
   frequency
   remove_noise
   ss
   gb
   tensor
   strike
   dimensionality
   skew
   anisotropy
   impedance
   tf
   afmag
   ztem
   mobilemt
   spectra
   source_effects
   source_array
   fieldzone
   csumt
   gradient_imaging
   lcurve
   plot
   advanced
   diag
   ../ai_processing/index
