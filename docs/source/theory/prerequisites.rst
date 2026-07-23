.. _theory_prerequisites:

Prerequisites
=============

Before working through the response functions, distortion effects, and
inversion ideas in the rest of this section, it helps to know which page
answers which question and to share a common notation for the handful of
quantities that recur on nearly every page. This page is that map: a
reading order depending on the field method in use, and the small
vocabulary of symbols assumed everywhere else.

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

Many pyCSAMT workflows revolve around a small set of physical quantities,
each with a notation kept consistent across every theory page:

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

These are notation, not numbers -- the actual constants and unit factors
behind :math:`\omega=2\pi f` or :math:`\rho_a`'s conversion factors live in
:doc:`constants`, which is where to look when a page's formula needs a
concrete value rather than a symbol.

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
