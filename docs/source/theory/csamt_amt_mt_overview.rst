.. _csamt_amt_mt_overview:

CSAMT, AMT, and MT Overview
===========================

pyCSAMT works with several related electromagnetic (EM) sounding methods:
magnetotellurics (MT), audio-magnetotellurics (AMT), controlled-source
audio-frequency magnetotellurics (CSAMT), controlled-source EM (CSEM), and
time-domain EM (TDEM). They share the same physical foundation: electrical
conductivity controls how time-varying EM fields diffuse through the Earth.
They differ mainly in source type, frequency band, acquisition geometry, and
the assumptions that are safe during interpretation.

This page gives the conceptual background needed before reading the more
specific pages on :doc:`impedance_tensor`, :doc:`static_shift`, and
:doc:`inversion_concepts`.

Why Resistivity Matters
-----------------------

Most pyCSAMT workflows are ultimately trying to estimate subsurface
resistivity, commonly denoted :math:`\rho`, or its inverse conductivity,
:math:`\sigma`:

.. math::

   \rho = \frac{1}{\sigma}.

In geophysical EM methods, resistivity is not measured directly. The field
system measures electric and magnetic fields at the surface or near the
surface. Those fields are then transformed into response functions such as
impedance, apparent resistivity, phase, tipper, or transient decay curves.
The interpreter uses those responses to infer geological structure.

Resistivity is useful because rocks, fluids, alteration zones, graphite,
sulfides, clay, salinity, porosity, permeability, temperature, and fracture
connectivity can all influence bulk electrical behavior. A resistivity model
is therefore not a lithology map by itself; it is a physical property model
that must be interpreted with geology, boreholes, geochemistry, hydrology,
and survey context.

The Family Of Methods
---------------------

The methods supported by pyCSAMT can be grouped by source.

.. list-table::
   :header-rows: 1
   :widths: 18 22 24 36

   * - Method
     - Source
     - Typical domain
     - Main interpretation concern
   * - MT
     - Natural ionospheric and magnetospheric variations
     - Broad frequency band, often deep investigation
     - Assumes a plane-wave natural source at the survey scale.
   * - AMT
     - Natural audio-frequency variations
     - Higher frequencies than long-period MT, shallower targets
     - Signal strength can be weak or culturally noisy in the audio band.
   * - CSAMT
     - Controlled grounded electric dipole source
     - Audio-frequency controlled-source sounding
     - Must distinguish far-field behavior from near-field and transition-zone
       source effects.
   * - CSEM
     - Controlled electric or magnetic source
     - Frequency-domain controlled-source EM
     - Source geometry and full transmitter-receiver coupling are part of
       the forward problem.
   * - TDEM
     - Controlled transient source
     - Time-domain decay after source turn-off
     - Interpretation depends on transient diffusion and gate timing rather
       than steady harmonic impedance.

Frequency-domain methods use harmonic signals with angular frequency
:math:`\omega`:

.. math::

   \omega = 2 \pi f,

where :math:`f` is frequency in hertz. The period is:

.. math::

   T = \frac{1}{f}.

High frequencies generally sample shallower structure; low frequencies
sample deeper structure. This is not a hard boundary, because sensitivity
also depends on resistivity, dimensionality, noise, source geometry, and the
inversion model.

Maxwell Equations In The Diffusive Regime
-----------------------------------------

MT, AMT, CSAMT, and many CSEM workflows operate at frequencies where
subsurface EM behavior is usually described by a diffusive approximation.
Starting from Maxwell equations in the frequency domain:

.. math::

   \nabla \times \mathbf{E}
   =
   - i \omega \mu \mathbf{H},

.. math::

   \nabla \times \mathbf{H}
   =
   \mathbf{J}_s + \sigma \mathbf{E}
   + i \omega \epsilon \mathbf{E},

where :math:`\mathbf{E}` is electric field, :math:`\mathbf{H}` is magnetic
field, :math:`\mathbf{J}_s` is a source current density, :math:`\sigma` is
conductivity, :math:`\mu` is magnetic permeability, and :math:`\epsilon` is
dielectric permittivity.

For many geophysical EM surveys in conductive Earth materials, the
displacement current term :math:`i \omega \epsilon \mathbf{E}` is small
relative to the conduction current :math:`\sigma \mathbf{E}`. The working
equation becomes approximately:

.. math::

   \nabla \times \mathbf{H}
   \approx
   \mathbf{J}_s + \sigma \mathbf{E}.

This is why EM sounding is often described as diffusion rather than wave
propagation. EM energy diffuses downward and laterally, with a penetration
scale controlled by resistivity and frequency.

Skin Depth
----------

A useful first-order scale is the skin depth :math:`\delta`, the distance
over which a plane EM field decays by a factor of :math:`e^{-1}` in a
uniform half-space:

.. math::

   \delta
   =
   \sqrt{\frac{2 \rho}{\omega \mu}}
   =
   \sqrt{\frac{\rho}{\pi f \mu}}.

Assuming non-magnetic rocks with :math:`\mu \approx \mu_0`, a common field
approximation is:

.. math::

   \delta \approx 503 \sqrt{\frac{\rho}{f}},

where :math:`\delta` is in metres, :math:`\rho` is in ohm-m, and :math:`f`
is in hertz.

Skin depth should not be read as "the investigation depth." It is a scale
for a uniform Earth. Real sensitivity depends on the full conductivity
distribution and on the measured response. Still, it gives a useful
intuition:

* decreasing frequency increases penetration;
* increasing resistivity increases penetration;
* conductive cover can strongly attenuate high-frequency fields;
* shallow conductors may screen deeper targets.

MT And AMT
----------

MT uses naturally occurring EM fields. At suitable distances from source
regions, these fields can be approximated as plane waves incident on the
Earth. Under the plane-wave assumption, the horizontal electric and magnetic
fields at a station are related by a complex impedance tensor:

.. math::

   \begin{bmatrix}
   E_x \\
   E_y
   \end{bmatrix}
   =
   \begin{bmatrix}
   Z_{xx} & Z_{xy} \\
   Z_{yx} & Z_{yy}
   \end{bmatrix}
   \begin{bmatrix}
   H_x \\
   H_y
   \end{bmatrix}.

The tensor components vary with frequency. From them, one computes apparent
resistivity and phase. A common component-wise apparent resistivity is:

.. math::

   \rho_{a,ij}
   =
   \frac{1}{\mu_0 \omega}
   \left| Z_{ij} \right|^2,

and the impedance phase is:

.. math::

   \phi_{ij}
   =
   \tan^{-1}
   \left(
   \frac{\operatorname{Im}(Z_{ij})}
        {\operatorname{Re}(Z_{ij})}
   \right).

AMT is usually treated as the audio-frequency subset of MT. It is useful for
shallower surveys because high frequencies have shorter skin depths. AMT
surveys are often sensitive to:

* near-surface conductivity variations;
* cultural noise in the audio band;
* weak natural-field intervals;
* static shift of apparent resistivity;
* station-to-station consistency and frequency coverage.

pyCSAMT processing pages and pipeline steps often use MT/AMT language even
when the same operation is useful for CSAMT impedance-like data. The
important question is whether the response can be treated as an impedance
response under the assumptions of the workflow.

CSAMT
-----

CSAMT uses a controlled source, commonly a grounded electric dipole
transmitter, and measures fields at receivers along one or more profiles.
The controlled source improves signal repeatability and can be extremely
useful when natural AMT signal levels are weak or cultural noise is high.

The price is that the source is no longer "infinitely far away." Receiver
responses may contain source-field geometry, near-field effects, transition
zone behavior, source overprint, and shadow effects. These issues are central
to CSAMT methodology and are discussed in references such as
[Yan2004]_, [Chen2005]_, [Da2016]_, [WangLin2023]_, and [Zhang2021]_.

For a grounded electric dipole source, the source-receiver distance
:math:`r`, transmitter length, Earth resistivity, frequency, and survey
layout all affect whether a station behaves like a far-field MT-style
measurement. In the far field, CSAMT responses can often be interpreted with
MT-like apparent resistivity and phase concepts. In the near field or
transition zone, direct MT-style interpretation can be misleading.

Field Zones In CSAMT
--------------------

CSAMT interpretation often separates measurements into approximate field
zones:

.. list-table::
   :header-rows: 1
   :widths: 22 34 44

   * - Zone
     - Qualitative behavior
     - Interpretation risk
   * - Near field
     - Source geometry dominates part of the measured field.
     - MT-like apparent resistivity may reflect source coupling rather than
       subsurface layering alone.
   * - Transition zone
     - Both source geometry and diffusive Earth response matter.
     - Apparent resistivity curves may show source overprint, shadow effects,
       or frequency-dependent distortion.
   * - Far field
     - The source field approximates a plane-wave-like response at receivers.
     - MT-style impedance interpretation is more defensible, though noise,
       static shift, and dimensionality still matter.

There is no universal single number that separates these zones for all
surveys. A practical diagnostic uses the ratio between source-receiver
distance and an EM length scale such as skin depth. Let:

.. math::

   \eta = \frac{r}{\delta}.

Large :math:`\eta` values tend to be more far-field-like, while small
:math:`\eta` values tend to be more source-dominated. However, this is only
a heuristic. Transmitter configuration, anisotropy, topography, 3-D geology,
and local noise can change the practical interpretation.

In pyCSAMT, source-effect and field-zone tools should therefore be treated as
diagnostics. They help identify where a conventional MT-style workflow may
be safe, questionable, or inappropriate.

Apparent Resistivity Is Not True Resistivity
--------------------------------------------

Apparent resistivity is the resistivity of a hypothetical uniform Earth that
would produce the observed response at one frequency. It is not the actual
resistivity at a single depth. For MT-like impedance data:

.. math::

   \rho_a(f)
   =
   \frac{|Z(f)|^2}{\mu_0 \omega}.

A pseudosection plots :math:`\rho_a` or phase against station position and
period/frequency. It is useful for quality control and qualitative
geological inspection, but it is not a true cross-section. Main limitations
include:

* depth is only indirectly related to period;
* apparent resistivity blends effects from a volume of Earth;
* conductive and resistive targets have different sensitivity footprints;
* topography, source effects, and static shift can distort patterns;
* 2-D and 3-D structures can produce responses that do not map vertically
  below each station.

For this reason, pyCSAMT documentation separates quick-look products from
inversion products. Pseudosections are diagnostic and interpretive aids;
resistivity models require a forward/inverse modelling step.

Dimensionality: 1-D, 2-D, And 3-D Earth Assumptions
---------------------------------------------------

The interpretation method must match the geology closely enough to be useful.
Dimensionality describes how resistivity varies in space.

.. list-table::
   :header-rows: 1
   :widths: 16 38 46

   * - Assumption
     - Resistivity structure
     - Typical use
   * - 1-D
     - :math:`\rho = \rho(z)`
     - Layered Earth sounding, first-pass modelling, station-by-station
       inversion, synthetic tests.
   * - 2-D
     - :math:`\rho = \rho(x,z)` and approximately invariant along strike.
     - Profile interpretation when geological strike is stable and station
       spacing is along a crossing line.
   * - 3-D
     - :math:`\rho = \rho(x,y,z)`
     - Complex geology, strong lateral changes, non-profile surveys, source
       and topographic complexity.

Many CSAMT/AMT field projects start with 1-D or 2-D thinking because it is
fast and interpretable. That does not mean the Earth is actually 1-D or 2-D.
Phase tensor analysis, skew, tipper behavior, induction arrows, strike
estimates, and residual patterns help decide whether a simplified inversion
is defensible.

Source Type And Data Interpretation
-----------------------------------

A useful way to think about method choice is to ask: "What source generated
the field I am interpreting?"

Natural-source MT/AMT:

* source is not controlled by the operator;
* plane-wave assumption is central;
* impedance tensor is the primary response;
* remote reference and robust processing may be important;
* weak natural signal intervals can limit data quality.

Controlled-source CSAMT/CSEM:

* source location, orientation, waveform, and current are part of the
  experiment;
* signal strength can be high and repeatable;
* source geometry can contaminate MT-style apparent resistivity;
* near-field, transition-zone, and far-field behavior must be checked;
* full controlled-source modelling may be needed for rigorous inversion.

Time-domain TDEM:

* source is switched off or changed in time;
* measured response is a transient decay;
* time gates replace frequency samples;
* early times tend to sample shallower structure, later times deeper
  structure;
* receiver coupling, gate timing, and transmitter waveform are central.

The same field project may combine methods. For example, TDEM can constrain
shallow structure or static shift, AMT/CSAMT can provide profile-scale
coverage, and MT can extend depth of investigation.

Coordinate And Component Conventions
------------------------------------

Most MT-style processing uses horizontal components:

* :math:`E_x`, :math:`E_y`: horizontal electric fields;
* :math:`H_x`, :math:`H_y`: horizontal magnetic fields;
* :math:`H_z`: vertical magnetic field, used in tipper analysis.

Coordinate conventions must be tracked carefully. Rotating data changes the
impedance tensor:

.. math::

   \mathbf{Z}' =
   \mathbf{R}(\theta)
   \mathbf{Z}
   \mathbf{R}^{T}(\theta),

where :math:`\mathbf{R}(\theta)` is a horizontal rotation matrix. The goal is
often to align data with geological strike or profile orientation. The
details are covered in :doc:`impedance_tensor`, but the practical message is
simple: component labels only have meaning relative to a coordinate system.

Common Data Products
--------------------

pyCSAMT workflows commonly create or consume the following products.

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * - Product
     - What it shows
     - How to use it
   * - Apparent resistivity curve
     - :math:`\rho_a` versus frequency or period at one station.
     - Inspect shifts, slopes, outliers, and frequency coverage.
   * - Phase curve
     - Impedance phase versus frequency or period.
     - Identify conductive/resistive trends and dimensionality issues.
   * - Pseudosection
     - Response plotted by station and period/frequency.
     - Quick-look lateral continuity and QC; not a true depth section.
   * - Phase tensor plot
     - Distortion-resistant tensor geometry.
     - Dimensionality, strike, skew, and structural complexity.
   * - Tipper/induction vectors
     - Vertical magnetic response to horizontal magnetic fields.
     - Lateral conductivity contrasts and 3-D effects.
   * - Inversion model
     - Resistivity distribution produced by fitting data with a forward model.
     - Main quantitative interpretation product, subject to uncertainty.
   * - Residual or misfit map
     - Difference between observed and predicted data.
     - Diagnose poor fits, noisy stations, model inadequacy, or source issues.

How pyCSAMT Uses These Concepts
-------------------------------

The theory pages are not isolated from the software. They explain why the
workflow is organized the way it is.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - pyCSAMT area
     - Theoretical connection
   * - ``pycsamt.seg`` and ``pycsamt.site``
     - Store station data, fields, impedance tensors, coordinates, and survey
       context.
   * - ``pycsamt.emtools``
     - Implements processing and diagnostic operations such as frequency
       editing, tensor analysis, static shift handling, source-effect
       diagnostics, and QC plots.
   * - ``pycsamt.pipeline``
     - Chains theory-aware processing steps into reproducible workflows.
   * - ``pycsamt.models`` and ``pycsamt.inversion``
     - Prepare and interpret model-based inversion workflows such as Occam2D
       and ModEM.
   * - ``pycsamt.agents``
     - Help translate user requests into workflow plans, run diagnostics, and
       explain output products.
   * - ``pycsamt.app``
     - Provides desktop and web interfaces for the same processing concepts.

Practical Reading Of A Survey
-----------------------------

When approaching a new MT/AMT/CSAMT survey, a robust interpretation sequence
is:

1. Confirm survey type, source type, coordinate system, station spacing,
   frequency/period band, and file format.
2. Inspect station metadata and frequency coverage before interpreting any
   pseudosection.
3. Check apparent resistivity and phase curves station by station.
4. For CSAMT, inspect field-zone and source-effect diagnostics before
   treating the data as MT-like impedance data.
5. Check static shift risk, especially in rugged near-surface geology or
   resistive/conductive shallow cover.
6. Evaluate dimensionality with phase tensors, skew, strike, tipper, and
   profile consistency.
7. Choose a model class: 1-D, 2-D, 3-D, or controlled-source forward/inverse
   modelling as required.
8. Interpret inversion results with residuals, uncertainty, geology, and
   independent constraints.

This sequence is deliberately conservative. EM responses are non-unique, and
many different subsurface models can explain similar data. Good workflows
make assumptions explicit and preserve enough QC information to revisit those
assumptions later.

Common Interpretation Pitfalls
------------------------------

The most common mistakes are conceptual rather than computational:

* treating a pseudosection as a true geological section;
* ignoring source effects in CSAMT data;
* interpreting shifted apparent resistivity as a real layer without checking
  static shift;
* assuming high-frequency data always mean shallow depth in every resistivity
  environment;
* using a 2-D inversion when the data show strong 3-D behavior;
* comparing stations without checking coordinate rotation and component
  convention;
* trusting a low RMS inversion without inspecting residual distribution and
  geological plausibility.

The purpose of pyCSAMT's processing, pipeline, and agent layers is to make
these checks repeatable.

References
----------

This overview is consistent with standard EM theory and MT/CSAMT practice,
including the methodology references collected in :doc:`../references`.
Especially relevant entries include [WardHohmann1988]_, [Yan2004]_,
[Chen2005]_, [Da2016]_, [WangLin2023]_, [Zhang2021]_, [deGrootHedlin1990]_,
and [Kelbert2014]_.

Next Steps
----------

Continue with:

* :doc:`impedance_tensor` for tensor notation, apparent resistivity, phase,
  rotations, strike, and tipper;
* :doc:`static_shift` for near-surface galvanic distortion and correction
  strategies;
* :doc:`inversion_concepts` for model fitting, regularization, misfit, and
  non-uniqueness;
* :doc:`tdem_basics` for time-domain EM concepts.
