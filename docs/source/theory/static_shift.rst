.. _static_shift:

Static Shift
============

Static shift is a near-surface galvanic distortion that moves apparent
resistivity curves up or down without producing the same change in phase. It
is one of the most important practical problems in MT, AMT, and CSAMT
interpretation because it can create false shallow resistivity contrasts and
can bias inversion starting models, pseudo-sections, and final sections.

In the simplest case, a station sits above a small shallow conductor or
resistor. The shallow body distorts the electric field measured at the
surface. The magnetic field and phase relationship may remain comparatively
stable, but the electric-field amplitude is scaled. Because apparent
resistivity is proportional to impedance amplitude squared, a modest
electric-field scaling can become a large apparent-resistivity shift.

This page explains the physics, diagnostics, and pyCSAMT correction workflow
for static shift. The companion tutorial :doc:`../tutorials/correct_static_shift`
shows the practical API calls.

Why Static Shift Happens
------------------------

MT, AMT, and far-field CSAMT interpretation assume that the measured
horizontal electric and magnetic fields represent regional electromagnetic
induction. Near the surface, however, small resistivity heterogeneities can
accumulate charge and distort the electric field. This is called galvanic
distortion because it is controlled by current channelling and charge
build-up in shallow conductivity contrasts.

Common causes include:

* conductive clay lenses;
* dry resistive boulders or lateritic caps;
* weathered fracture zones;
* shallow groundwater salinity changes;
* cultural conductors such as fences, pipes, cables, and grounding systems;
* local topographic or contact effects near electrodes.

The important feature is scale. Static shift is usually caused by structures
that are shallow and small compared with the induction scale of the periods
or frequencies being interpreted. They can strongly affect local electric
field amplitudes while leaving the deeper induction response broadly intact.

Mathematical Form
-----------------

For a simple scalar static shift, the observed apparent resistivity at one
station can be written as

.. math::

   \rho_{a,\mathrm{obs}}(f) = s_{\rho}\,\rho_{a,\mathrm{true}}(f),

where :math:`s_{\rho}` is a positive station-dependent scale factor. In
logarithmic units this becomes an additive offset:

.. math::

   \log_{10}\rho_{a,\mathrm{obs}}(f)
   =
   \log_{10}\rho_{a,\mathrm{true}}(f)
   + \Delta_{\rho}.

The offset :math:`\Delta_{\rho}` is approximately constant with frequency
when the distortion is truly static.

For impedance amplitude, the relation is

.. math::

   |Z_{\mathrm{obs}}| = s_Z |Z_{\mathrm{true}}|,

and because apparent resistivity is proportional to :math:`|Z|^2`,

.. math::

   s_{\rho} = s_Z^2
   \quad \hbox{and} \quad
   s_Z = \sqrt{s_{\rho}}.

pyCSAMT correction tables therefore report both ``fac_rho`` and ``fac_z``.
The impedance tensor is corrected with ``fac_z`` because the stored field
quantity is impedance, while apparent resistivity curves respond with the
square of that factor.

Tensor View
-----------

The scalar model is useful, but real galvanic distortion can mix impedance
components. A common distortion model is

.. math::

   \mathbf{Z}_{\mathrm{obs}}
   =
   \mathbf{C}\,\mathbf{Z}_{\mathrm{regional}},

where :math:`\mathbf{C}` is a real, frequency-independent distortion matrix.
Static shift is the amplitude-scaling part of this broader distortion
problem. Twist, shear, anisotropy, and component mixing can also occur.

This is why static-shift correction should be paired with tensor diagnostics.
If phase tensor skew, induction vectors, or component residuals suggest
strong 3-D mixing, a simple station-level amplitude factor may not be enough.
See :ref:`impedance_tensor` for the impedance and phase-tensor background.

Why Phase Is Less Affected
--------------------------

Apparent resistivity uses impedance amplitude:

.. math::

   \rho_a =
   \frac{1}{\mu_0 \omega}
   |Z|^2.

The impedance phase uses the ratio of imaginary and real parts:

.. math::

   \phi = \tan^{-1}
   \left(
      \frac{\mathrm{Im}(Z)}
           {\mathrm{Re}(Z)}
   \right).

If the distortion is a real scalar multiplier, both
:math:`\mathrm{Re}(Z)` and :math:`\mathrm{Im}(Z)` are scaled by the same
factor. The amplitude changes, but the phase ratio does not:

.. math::

   \tan^{-1}
   \left(
      \frac{s_Z\,\mathrm{Im}(Z)}
           {s_Z\,\mathrm{Re}(Z)}
   \right)
   =
   \tan^{-1}
   \left(
      \frac{\mathrm{Im}(Z)}
           {\mathrm{Re}(Z)}
   \right).

This is the classic diagnostic pattern: apparent resistivity curves are
offset relative to nearby stations, while phase curves remain smooth and
geologically consistent. If phase is also strongly distorted, the problem may
not be simple static shift.

Static Shift Versus Near-Surface Effect
---------------------------------------

pyCSAMT makes an important operational distinction between:

* **static shift** - a frequency-independent multiplicative shift of the
  whole apparent-resistivity curve;
* **near-surface effect** - a frequency-dependent distortion, often strongest
  at high frequencies, caused by shallow inhomogeneities or source/electrode
  effects.

Static shift can often be corrected by estimating a station-level amplitude
factor. Frequency-dependent near-surface effects should not be "fixed" with
a single static factor because the curve shape itself has changed. In that
case, better options include masking the contaminated band, improving
processing, using 2-D/3-D inversion, or explicitly modelling the near-surface
structure.

In :mod:`pycsamt.emtools.ss`, near-surface distortion classification uses
three station diagnostics:

* ``ss_delta_log10`` - the median log10 apparent-resistivity offset relative
  to an AMA spatial trend;
* ``ns_index`` - a high-frequency versus low-frequency residual spread ratio;
* ``gradient_delta`` - the difference between high-frequency and
  low-frequency log-log slopes.

The classification vocabulary is:

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * - Class
     - Meaning
     - Suggested response
   * - ``clean``
     - No strong static or near-surface flag.
     - Keep the station; inspect normally.
   * - ``static``
     - Frequency-independent offset dominates.
     - Apply static-shift correction and compare before/after.
   * - ``near_surface``
     - Frequency-dependent high-frequency distortion dominates.
     - Mask contaminated band or use modelling/inversion that can represent
       shallow structure.
   * - ``mixed``
     - Both offset and frequency-dependent distortion are present.
     - Correct cautiously; preserve diagnostics and consider excluding
       unstable bands.

How It Appears In Profiles
--------------------------

In a profile, static shift usually appears as a station-level vertical jump
in apparent resistivity. The shifted curve may preserve the same shape as
neighbouring curves but sit consistently higher or lower.

Look for:

* one or a few stations offset from a smooth spatial trend;
* apparent resistivity discontinuities without corresponding phase changes;
* TE and TM apparent resistivity shifts that are not matched by phase;
* determinant apparent resistivity offsets near one station;
* pseudo-section vertical stripes or station-centred anomalies;
* shallow inversion anomalies that disappear when the station is corrected
  or down-weighted.

Do not diagnose static shift from apparent resistivity alone. Compare phase,
tensor skew, station quality, contact resistance, field notes, and nearby
geology.

CSAMT-Specific Cautions
-----------------------

CSAMT can contain both galvanic static shift and source-related effects.
These should not be confused.

Static shift is station-local and approximately frequency independent in
log apparent resistivity. Source effects can be frequency dependent and may
vary systematically with transmitter distance, line geometry, and transition
from near field to far field. Shadow and source-overprint effects can change
curve shape, not only level.

Before applying static-shift correction to CSAMT data, check:

* whether the frequency band is far-field enough for the intended
  interpretation;
* whether the shifted station also has unusual phase behavior;
* whether neighbouring stations share the same source-distance trend;
* whether the apparent shift is actually a transmitter geometry effect;
* whether a correction improves or damages response consistency.

See :ref:`csamt_amt_mt_overview` for the broader CSAMT/AMT/MT distinction.

Correction Philosophy
---------------------

Static-shift correction should be conservative. The goal is not to make every
curve look smooth or pretty. The goal is to remove station-local amplitude
bias while preserving regional induction structure.

Good correction practice follows these principles:

* estimate factors from neighbouring stations or independent constraints;
* avoid using known bad stations as references;
* avoid correcting frequency-dependent curve-shape problems with one scalar;
* preserve before/after plots and correction factors;
* use ``inplace=False`` or copy-based workflows while exploring;
* rerun inversion diagnostics after correction.

In hydrogeological and engineering surveys, a shallow conductor may be a real
target rather than a nuisance. If a shallow body is important to the study,
document whether it is being corrected as distortion or interpreted as
geology.

AMA Correction In pyCSAMT
-------------------------

The default pyCSAMT correction path is AMA, an adaptive moving-average
spatial trend. For each station, pyCSAMT estimates a regional trend from
nearby stations and compares the station's log apparent resistivity against
that trend. The shift estimate is

.. math::

   \Delta_{\rho}
   =
   \mathrm{median}_f
   \left[
      \log_{10}\rho_{a,i}(f)
      -
      T_i(f)
   \right],

where :math:`T_i(f)` is the spatial trend estimated from neighbouring
stations at frequency :math:`f`.

The correction factor applied to apparent resistivity is

.. math::

   \mathrm{fac}_{\rho} = 10^{-\Delta_{\rho}},

and the correction factor applied to impedance is

.. math::

   \mathrm{fac}_Z = 10^{-\Delta_{\rho}/2}.

This is why the output table from :func:`pycsamt.emtools.ss.estimate_ss_ama`
contains:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Column
     - Meaning
   * - ``station``
     - Station identifier.
   * - ``delta_log10_rho``
     - Estimated log10 apparent-resistivity shift before correction.
   * - ``fac_rho``
     - Multiplicative factor for apparent resistivity.
   * - ``fac_z``
     - Multiplicative factor for impedance tensor amplitudes.
   * - ``n_used``
     - Number of samples used in the shift estimate.

Basic API example:

.. code-block:: python
   :linenos:

   from pycsamt.api import read_edis
   from pycsamt.emtools.ss import estimate_ss_ama, correct_ss_ama

   survey = read_edis("data/edis")
   sites = survey.collection

   factors = estimate_ss_ama(
       sites,
       sort_by="lon",
       half_window=3,
       weights="tri",
       pband=None,
       max_skew=6.0,
       api=True,
   )
   print(factors)

   corrected = correct_ss_ama(
       sites,
       sort_by="lon",
       half_window=3,
       weights="tri",
       inplace=False,
   )

The important parameters are:

* ``sort_by`` - station order for the spatial profile, usually ``"lon"`` for
  east-west profiles and ``"lat"`` for north-south profiles;
* ``half_window`` - number of neighbouring stations on each side;
* ``weights`` - triangular, Gaussian, or uniform spatial weights;
* ``pband`` - optional period band used to estimate the factor;
* ``max_skew`` - optional phase-tensor skew ceiling used to avoid strongly
  distorted stations in the trend.

Choosing The Estimation Band
----------------------------

The optional period band ``pband=(T_min, T_max)`` is important when only part
of the curve is trustworthy. A good band should:

* contain enough frequencies at most stations;
* avoid dead bands and noisy acquisition ranges;
* avoid strong near-field CSAMT source effects;
* avoid high-frequency near-surface curve-shape distortion;
* represent the regional trend rather than a local target.

If the chosen band is too narrow, correction factors become unstable. If it
includes contaminated frequencies, the estimated factor can absorb effects
that are not static shift.

LOESS, Reference Median, and Bilateral Options
----------------------------------------------

The pipeline catalogue includes several static-shift correction styles:

.. list-table::
   :header-rows: 1
   :widths: 18 28 54

   * - Code
     - Method
     - Typical use
   * - ``SS001``
     - AMA
     - Default profile-based correction using neighbouring stations.
   * - ``SS002``
     - LOESS
     - Smooth spatial trend when station spacing is irregular or gradual
       profile-scale drift is expected.
   * - ``SS003``
     - Reference median
     - Correction relative to a robust median or reference response.
   * - ``SS004``
     - Bilateral
     - Edge-aware spatial correction that reduces smoothing across sharp
       station-to-station changes.

The best method depends on survey geometry. AMA is simple and transparent.
LOESS can be better for uneven spacing. Reference median methods need a
credible reference. Bilateral filters can preserve sharp spatial transitions
but still require careful diagnostics.

Stratagem Workflow
------------------

The Stratagem helper wraps the same idea for AMT/CSAMT survey processing:

.. code-block:: python
   :linenos:

   from pycsamt.stratagem.process import StaticShiftCorrector

   corrector = StaticShiftCorrector(
       sort_by="lon",
       half_window=3,
       weights="tri",
       pband=None,
       max_skew=6.0,
   )

   corrected_edis = corrector.fit(edis, copy=True).out()
   factors = corrector.factors_

``copy=True`` preserves the original EDI objects. The correction is applied
to impedance amplitudes, so exported apparent resistivity and phase products
should be regenerated from the corrected impedance.

In Stratagem pipelines, static-shift correction should usually run before
frequency filtering that removes many impedance samples. The code deliberately
warns about this because AMA needs enough common frequency support to compare
neighbouring stations.

Agent And Application Workflows
-------------------------------

The same correction concept is exposed in higher-level pyCSAMT workflows:

* :class:`pycsamt.agents.static_shift.StaticShiftAgent` can run correction,
  collect summary tables, and prepare report text;
* the desktop application exposes AMA-style static-shift correction in the
  correction tools;
* the web application can call the static-shift agent as part of an assisted
  workflow;
* :mod:`pycsamt.pipeline` registers static-shift steps so they can be placed
  between QC and inversion preparation.

These interfaces should still produce the same core evidence: correction
factors, before/after curves or pseudo-sections, and enough metadata to
reproduce the correction.

Before/After Interpretation
---------------------------

After correction, inspect:

* apparent resistivity before and after correction;
* phase before and after correction;
* pseudo-sections of the correction delta;
* station-level correction factors;
* tensor skew or phase-tensor diagnostics;
* inversion residuals before and after correction.

The correction is plausible when:

* shifted apparent resistivity curves align better with neighbouring
  station trends;
* phase remains essentially unchanged;
* correction factors are not extreme;
* response fits improve without creating new residual patterns;
* the corrected section is geologically coherent and still honours known
  shallow structure.

The correction is suspicious when:

* a large factor is required for many adjacent stations;
* phase changes or phase residuals are also problematic;
* the apparent shift is frequency dependent;
* the corrected curve no longer matches nearby geology or field notes;
* inversion improves globally but worsens systematically at corrected
  stations.

Effect On Inversion
-------------------

Static shift can strongly affect inversion. If left uncorrected, a downshift
in apparent resistivity may be interpreted as a shallow conductor, while an
upshift may be interpreted as a shallow resistor. In smooth 2-D inversion,
these station-centred effects can smear laterally and vertically, producing
misleading near-surface structure.

In inversion terms, static shift changes the amplitude part of the data
vector and therefore changes the weighted residuals:

.. math::

   r_i =
   \frac{d_{\mathrm{obs},i} - d_{\mathrm{pred},i}}{\sigma_i}.

If the shifted data are given small errors, the inversion may spend model
structure trying to fit a station-local distortion. Correction, masking, or
larger error floors are all defensible responses depending on evidence.

For production inversion, keep both:

* the uncorrected data and QC record;
* the corrected data, correction table, parameters, and plots.

This makes the interpretation auditable and allows later sensitivity tests.

Common Mistakes
---------------

Avoid these mistakes:

* correcting every station simply because a method is available;
* treating frequency-dependent distortion as static shift;
* using a bad station in the neighbourhood trend;
* applying correction after aggressive frequency filtering without enough
  common bandwidth;
* accepting large correction factors without field or tensor evidence;
* comparing corrected apparent resistivity to uncorrected phase products;
* interpreting a removed shallow anomaly without checking whether it could
  be real geology;
* hiding correction factors from reports and inversion provenance.

Recommended Workflow
--------------------

A conservative workflow is:

#. Inspect raw curves and station quality.
#. Check phase, tensor skew, and obvious cultural-noise indicators.
#. Estimate static-shift factors with AMA or another spatial method.
#. Plot correction deltas and before/after curves.
#. Classify stations as clean, static, near-surface, or mixed.
#. Correct only where the evidence supports a station-level amplitude shift.
#. Prepare inversion inputs from corrected impedance data.
#. Compare inversion residuals and models before and after correction.
#. Report correction parameters and factors.

Minimal reproducible correction:

.. code-block:: python
   :linenos:

   from pycsamt.api import read_edis
   from pycsamt.emtools.ss import (
       detect_near_surface,
       correct_ss_ama,
       estimate_ss_ama,
   )

   sites = read_edis("data/edis").collection

   classes = detect_near_surface(
       sites,
       f_split=1.0,
       ns_threshold=2.0,
       ss_threshold=0.1,
       sort_by="lon",
       half_window=3,
       weights="tri",
       api=True,
   )

   factors = estimate_ss_ama(
       sites,
       sort_by="lon",
       half_window=3,
       weights="tri",
       max_skew=6.0,
       api=True,
   )

   corrected = correct_ss_ama(
       sites,
       sort_by="lon",
       half_window=3,
       weights="tri",
       inplace=False,
   )

For CLI and pipeline execution, inspect the static-shift step catalogue with:

.. code-block:: console
   :linenos:

   pycsamt pipe steps --category static_shift

Reporting Checklist
-------------------

A report should include:

* correction method, parameters, and software version;
* station ordering axis and neighbourhood window;
* period or frequency band used for estimating factors;
* maximum skew or masking criteria;
* correction factor table;
* before/after apparent resistivity and phase plots;
* stations excluded or flagged as near-surface/mixed;
* effect on inversion RMS and residuals;
* statement about interpretation uncertainty.

Next Steps
----------

For practical usage, see:

* :doc:`../tutorials/correct_static_shift`;
* :doc:`../pipeline/steps` for static-shift pipeline steps;
* :doc:`../agents/processing_agents` for ``StaticShiftAgent``;
* :ref:`impedance_tensor` for impedance and phase-tensor diagnostics;
* :ref:`inversion_concepts` for the inversion consequences of data
  weighting and residuals.

References
----------

The static-shift discussion here follows the standard galvanic-distortion
view of EM data [WardHohmann1988]_ and the practical distinction between
static and frequency-dependent near-surface effects discussed for CSAMT by
[Lei2017]_. Resistivity interpretation cautions follow [Palacky1988]_.
