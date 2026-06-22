.. _theory:

Scientific Background
=====================

This section explains core electromagnetic and inversion concepts used by
pyCSAMT workflows.

pyCSAMT is a practical modelling and interpretation library, but the outputs
are only meaningful when the physical assumptions are understood. The pages in
this section explain the response functions, field-method differences,
distortion effects, inversion ideas, and time-domain concepts that appear
throughout the rest of the documentation.

Use this section when you need to understand why a workflow asks for a
particular component, error floor, dimensionality, correction, or diagnostic
plot.

Concept Map
-----------

.. list-table::
   :header-rows: 1
   :widths: 26 44 30

   * - Topic
     - What it explains
     - Read before
   * - :doc:`csamt_amt_mt_overview`
     - The relationship between CSAMT, AMT, MT, CSEM, and TDEM; source types;
       frequency-domain assumptions; and survey-design consequences.
     - Choosing acquisition assumptions, model backends, or processing paths.
   * - :doc:`impedance_tensor`
     - Electric and magnetic field coupling, impedance tensor components,
       apparent resistivity, phase, dimensionality indicators, and rotation.
     - Reading EDI data, selecting TE/TM or full tensor components, plotting
       response fits.
   * - :doc:`static_shift`
     - Near-surface galvanic distortion, static-shift symptoms, correction
       strategies, and uncertainty that remains after correction.
     - Correcting apparent resistivity curves or interpreting shallow
       resistivity contrasts.
   * - :doc:`inversion_concepts`
     - Forward models, data misfit, RMS, regularization, model roughness,
       uncertainty floors, and how to read inversion diagnostics.
     - Occam2D, ModEM, MARE2DEM, and pipeline inversion workflows.
   * - :doc:`tdem_basics`
     - Transient EM diffusion, time gates, transmitter/receiver geometry,
       TDEM data fields, and frequency-domain compatibility.
     - Working with TDEM readers, transforms, or TDEM inversion inputs.

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
#. Continue to :doc:`../models/choosing_backend`.

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
* :doc:`../models/index` explains how physical assumptions map to external
  modelling backends.
* :doc:`../pipeline/index` explains how repeated processing and inversion
  workflows are organized.
* :doc:`../tutorials/index` gives worked examples that apply these concepts.

Contents
--------

.. toctree::
   :maxdepth: 1

   csamt_amt_mt_overview
   impedance_tensor
   static_shift
   inversion_concepts
   tdem_basics
