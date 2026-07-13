.. _emtools_gallery:

EM Tools Guide
==============

``pycsamt.emtools`` is the science-facing toolbox for MT, AMT, and CSAMT
processing.  Use this page as a task map: start with loading and inspection,
then move through quality control, frequency editing, noise/static-shift
conditioning, tensor/strike analysis, source diagnostics, and publication
plots.

Each card opens a narrative page with examples, figures, and links to the
public API.  The same workflows are also available as runnable gallery
examples in :ref:`EM tools examples <emtools_examples>`.  For the complete
callable reference, see :doc:`../api/emtools`.

Start Here
----------

.. grid:: 1 1 2 3
   :gutter: 3

   .. grid-item-card:: Inspect EDI Surveys
      :link: inspect
      :link-type: doc
      :img-top: ../../_static/icons/loading.svg
      :class-card: pycsamt-card sd-text-center

      Load EDI folders, summarize stations, check missing sections, and draw
      quick response plots before processing.

   .. grid-item-card:: Quality Control
      :link: qc
      :link-type: doc
      :img-top: ../../_static/icons/diagnostic.svg
      :class-card: pycsamt-card sd-text-center

      Build confidence tables, QC flags, SNR summaries, and confidence
      pseudo-sections for audit-ready data review.

   .. grid-item-card:: Frequency Editing
      :link: frequency
      :link-type: doc
      :img-top: ../../_static/icons/edit.svg
      :class-card: pycsamt-card sd-text-center

      Select bands, align grids, drop duplicates, and recover or mask
      low-confidence frequencies with decision tables.

Condition And Correct
---------------------

.. grid:: 1 1 2 3
   :gutter: 3

   .. grid-item-card:: Noise Removal
      :link: remove_noise
      :link-type: doc
      :img-top: ../../_static/icons/processing.svg
      :class-card: pycsamt-card sd-text-center

      Suppress power-line harmonics, smooth response trends, apply EMAP
      filters, and document mitigation choices.

   .. grid-item-card:: Static Shift
      :link: ss
      :link-type: doc
      :img-top: ../../_static/icons/static-shift.svg
      :class-card: pycsamt-card sd-text-center

      Estimate AMA, LOESS, bilateral, or reference-median factors and compare
      corrected curves across a profile.

   .. grid-item-card:: Galvanic Distortion
      :link: gb
      :link-type: doc
      :img-top: ../../_static/icons/distorsion.svg
      :class-card: pycsamt-card sd-text-center

      Fit Groom--Bailey-style distortion parameters and decide whether a
      correction is defensible.

Tensor, Strike, And Dimensionality
----------------------------------

.. grid:: 1 1 2 3
   :gutter: 3

   .. grid-item-card:: Phase Tensor Tools
      :link: tensor
      :link-type: doc
      :img-top: ../../_static/icons/amt-csamt-mt-electromag.svg
      :class-card: pycsamt-card sd-text-center

      Build phase-tensor tables, maps, pseudo-sections, rose grids, and
      tensor operations.

   .. grid-item-card:: Strike Analysis
      :link: strike
      :link-type: doc
      :img-top: ../../_static/icons/location.svg
      :class-card: pycsamt-card sd-text-center

      Estimate geoelectric strike, compare methods, plot stability, and
      rotate tensors for 2-D workflows.

   .. grid-item-card:: Dimensionality
      :link: dimensionality
      :link-type: doc
      :img-top: ../../_static/icons/selection.svg
      :class-card: pycsamt-card sd-text-center

      Classify 1-D/2-D/3-D behavior, build pre-2-D inversion assessments, and
      identify localized 3-D effects.

   .. grid-item-card:: Skew
      :link: skew
      :link-type: doc
      :img-top: ../../_static/icons/polar-uncertainty.svg
      :class-card: pycsamt-card sd-text-center

      Compute Bahr skewness, select low-skew bands, and visualize skew
      traffic across station and period.

   .. grid-item-card:: Anisotropy
      :link: anisotropy
      :link-type: doc
      :img-top: ../../_static/icons/anisotropie.svg
      :class-card: pycsamt-card sd-text-center

      Measure apparent anisotropy and Swift-style skew/strike indicators for
      interpretation screening.

   .. grid-item-card:: Impedance Diagnostics
      :link: impedance
      :link-type: doc
      :img-top: ../../_static/icons/impedance.svg
      :class-card: pycsamt-card sd-text-center

      Inspect phasors, determinant tracks, and off-diagonal antisymmetry.

Source, Spectra, And Transfer Functions
---------------------------------------

.. grid:: 1 1 2 3
   :gutter: 3

   .. grid-item-card:: Transfer Functions
      :link: tf
      :link-type: doc
      :img-top: ../../_static/icons/transfer-tipper.svg
      :class-card: pycsamt-card sd-text-center

      Plot tippers, induction arrows, hodograms, rose diagrams, and
      multi-period induction maps.

   .. grid-item-card:: Spectra
      :link: spectra
      :link-type: doc
      :img-top: ../../_static/icons/spectra.svg
      :class-card: pycsamt-card sd-text-center

      Work with spectra-derived coherence, PSD, SNR, and transfer-function
      diagnostics.

   .. grid-item-card:: Source Effects
      :link: source_effects
      :link-type: doc
      :img-top: ../../_static/icons/source-control.svg
      :class-card: pycsamt-card sd-text-center

      Detect source overprint, normalize responses, and evaluate near-field
      correction behavior.

   .. grid-item-card:: Source Array Design
      :link: source_array
      :link-type: doc
      :img-top: ../../_static/icons/power.svg
      :class-card: pycsamt-card sd-text-center

      Explore array factors, beam steering, PAS patterns, and SNR gain for
      controlled-source layouts.

   .. grid-item-card:: Field Zones
      :link: fieldzone
      :link-type: doc
      :img-top: ../../_static/icons/map-and-profile.svg
      :class-card: pycsamt-card sd-text-center

      Classify near-, transition-, and far-field behavior along a profile.

Survey Design And Inversion Readiness
-------------------------------------

.. grid:: 1 1 2 3
   :gutter: 3

   .. grid-item-card:: CS/AMT Survey Design
      :link: csumt
      :link-type: doc
      :img-top: ../../_static/icons/controls.svg
      :class-card: pycsamt-card sd-text-center

      Estimate Bostick depth, frequency schedules, depth coverage, and
      vertical resolution.

   .. grid-item-card:: Gradient Imaging
      :link: gradient_imaging
      :link-type: doc
      :img-top: ../../_static/icons/gradient.svg
      :class-card: pycsamt-card sd-text-center

      Use spatial, frequency, and joint apparent-resistivity gradients for
      rapid structure screening.

   .. grid-item-card:: L-Curve
      :link: lcurve
      :link-type: doc
      :img-top: ../../_static/icons/bounds-vs-ridge-summary.svg
      :class-card: pycsamt-card sd-text-center

      Summarize regularization trade-offs and choose stable inversion
      parameters.

Plotting And Composite Diagnostics
----------------------------------

.. grid:: 1 1 2 3
   :gutter: 3

   .. grid-item-card:: Site Plots
      :link: plot
      :link-type: doc
      :img-top: ../../_static/icons/gallery-icon-images.svg
      :class-card: pycsamt-card sd-text-center

      Plot raw sites, response/tipper panels, survey comparisons, and
      observed-versus-predicted fits.

   .. grid-item-card:: Advanced Figures
      :link: advanced
      :link-type: doc
      :img-top: ../../_static/icons/advanced.svg
      :class-card: pycsamt-card sd-text-center

      Build composite diagnostics, Mohr circles, distortion radar plots,
      fingerprint plots, and multi-panel publication figures.

   .. grid-item-card:: Diagnostic Coverage
      :link: diag
      :link-type: doc
      :img-top: ../../_static/icons/external-validation.svg
      :class-card: pycsamt-card sd-text-center

      Evaluate coverage scores, uncertainty intervals, and polar error
      diagnostics against expected response behavior.

.. toctree::
   :hidden:
   :maxdepth: 1

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
