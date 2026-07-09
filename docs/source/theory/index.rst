.. _theory:

Scientific Background
=====================

This section explains the core electromagnetic and inversion concepts
behind pyCSAMT workflows. The library is practical by design, but its
outputs are only meaningful when the physical assumptions are
understood — these pages cover the response functions, field-method
differences, distortion effects, inversion ideas, and time-domain
concepts that appear throughout the rest of the documentation.

Use this section when you need to understand why a workflow asks for a
particular component, error floor, dimensionality, correction, or
diagnostic plot.

Concept Pages
-------------

.. grid:: 1 1 2 3
   :gutter: 3
   :class-container: cta-tiles

   .. grid-item-card:: CSAMT, AMT, and MT overview
      :link: csamt_amt_mt_overview
      :link-type: doc
      :img-top: ../_static/icons/amt-csamt-mt-electromag.svg
      :class-card: pycsamt-card sd-text-center

      How CSAMT, AMT, MT, CSEM, and TDEM relate: source types,
      frequency-domain assumptions, and the survey-design consequences.

   .. grid-item-card:: Impedance tensor
      :link: impedance_tensor
      :link-type: doc
      :img-top: ../_static/icons/impedance.svg
      :class-card: pycsamt-card sd-text-center

      Field coupling, tensor components, apparent resistivity and
      phase, dimensionality indicators, and rotation — the language of
      EDI data.

   .. grid-item-card:: Static shift
      :link: static_shift
      :link-type: doc
      :img-top: ../_static/icons/static-shift.svg
      :class-card: pycsamt-card sd-text-center

      Near-surface galvanic distortion: symptoms, correction
      strategies, and the uncertainty that remains after correction.

   .. grid-item-card:: Inversion concepts
      :link: inversion_concepts
      :link-type: doc
      :img-top: ../_static/icons/inversion-concepts.svg
      :class-card: pycsamt-card sd-text-center

      Forward models, misfit and RMS, regularization, roughness, and
      how to read inversion diagnostics before trusting a model.

   .. grid-item-card:: TDEM basics
      :link: tdem_basics
      :link-type: doc
      :img-top: ../_static/icons/tdem.svg
      :class-card: pycsamt-card sd-text-center

      Transient EM diffusion, time gates, transmitter/receiver
      geometry, and how TDEM data align with frequency-domain results.

Recommended Reading Paths
-------------------------

For MT, AMT, or CSAMT impedance workflows:

#. Start with :doc:`csamt_amt_mt_overview`.
#. Read :doc:`impedance_tensor`.
#. Read :doc:`static_shift` if apparent resistivity curves are shifted between
   nearby stations or if near-surface heterogeneity is expected.
#. Read :doc:`inversion_concepts` before running Occam2D, ModEM, or MARE2DEM.

For TDEM workflows:

#. Start with :doc:`tdem_basics`.
#. Read :doc:`inversion_concepts` to understand how transient data become an
   inversion data vector.
#. Read :doc:`csamt_amt_mt_overview` if TDEM products are being compared with
   frequency-domain EM results.

For model-backend decisions:

#. Read :doc:`csamt_amt_mt_overview` to understand method assumptions.
#. Read :doc:`inversion_concepts` to understand regularization and
   dimensionality.
#. Continue to :doc:`../user_guide/models/choosing_backend`.

Core Quantities
---------------

Many pyCSAMT workflows revolve around a small set of physical quantities.

.. list-table::
   :header-rows: 1
   :widths: 24 36 40

   * - Quantity
     - Typical notation
     - Why it matters
   * - Resistivity
     - :math:`\rho`
     - The main interpreted property in most EM inversions.
   * - Conductivity
     - :math:`\sigma = 1 / \rho`
     - Controls EM diffusion and current concentration.
   * - Angular frequency
     - :math:`\omega = 2 \pi f`
     - Links harmonic fields, impedance, skin depth, and phase.
   * - Impedance tensor
     - :math:`\mathbf{Z}`
     - Relates horizontal electric and magnetic fields in MT, AMT, and CSAMT
       style workflows.
   * - Apparent resistivity
     - :math:`\rho_a`
     - A frequency-dependent response estimate, not a direct layer
       resistivity.
   * - Phase
     - :math:`\phi`
     - Measures phase lag between field components and helps identify
       conductive or resistive structure.
   * - RMS misfit
     - :math:`RMS`
     - Summarizes data fit relative to assigned uncertainties.

Relationship To Other Sections
------------------------------

Theory pages connect directly to the applied documentation:

* :doc:`../getting_started/data_formats` explains how field data enter
  pyCSAMT.
* :doc:`../user_guide/models/index` explains how physical assumptions map to
  external modelling backends.
* :doc:`../user_guide/pipeline/index` explains how repeated processing and
  inversion workflows are organized.
* :doc:`../tutorials/index` gives worked examples that apply these concepts.

.. toctree::
   :maxdepth: 1
   :hidden:

   csamt_amt_mt_overview
   impedance_tensor
   static_shift
   inversion_concepts
   tdem_basics
