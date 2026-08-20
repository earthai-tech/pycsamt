.. _airborne_theory:

Airborne Natural-Source EM: AFMAG, ZTEM, And MobileMT
=====================================================

:doc:`csamt_amt_mt_overview` and :doc:`field_zones` already make the
key distinction that motivates this page: MT and AMT have no
controlled transmitter to be close to or far from, so the
near/transition/far dichotomy that governs CSAMT does not apply to
them -- their natural source is treated as a :term:`plane-wave field`
at essentially all usual survey scales. Every method on this page
inherits that same property, then adds a second one: the receiver is
airborne. A helicopter or fixed-wing aircraft tows one or more
magnetic-field sensors along a flight line while a fixed ground
station (or, for the historical instrument, nothing at all) provides
the reference. No transmitter is switched on anywhere in the survey.
This page explains why that combination works, what physical object
each of the three systems pyCSAMT supports actually measures, and how
those objects relate to the impedance-tensor formalism the rest of
the documentation builds around. The worked, code-driven
demonstrations live on their own pages --
:doc:`../user_guide/emtools/afmag` and
:doc:`../user_guide/emtools/ztem` -- and are cross-referenced rather
than repeated here.

Why A Natural Source Can Be Flown
---------------------------------

A grounded CSAMT transmitter radiates a field whose geometry depends
on the transmitter-receiver offset, which is exactly why
:doc:`field_zones` needs a near/transition/far classification at all.
The natural audio-frequency field that AFMAG-family methods use has
no such geometry to worry about: worldwide thunderstorm activity
(the *Schumann resonance* source band) illuminates the whole survey
area with what is, to a very good approximation, a uniform, planar,
vertically incident wave [Ward1959]_ [Legault2012]_. That is precisely
the plane-wave assumption :doc:`impedance_tensor` builds the tensor
formalism on, and it is why an AFMAG-family instrument does not need
a matched transmitter/receiver pair the way CSAMT does -- only a
receiver.

The same field diffuses downward and is attenuated on the skin-depth
scale already derived in :doc:`maxwell_forward`,

.. math::

   \delta \approx 503\sqrt{\rho/f}\ \mathrm{m}

(:eq:`eq-maxwell-skin-depth-rule`), which Legault et al. (2012, citing
Vozoff 1972) quote in exactly this form when motivating ZTEM's depth
of investigation. A lower audio frequency therefore probes deeper,
at the cost of a weaker natural signal and a coarser lateral
footprint -- the same trade-off that governs sounding-curve design in
ground MT/AMT, now applied to the choice of which of a handful of
fixed audio frequencies (typically a few tens to a few hundred Hz) an
airborne survey flies.

.. figure:: /images/theory/ztem-airbone.png
   :align: center
   :width: 70%

   Conceptual layout of a helicopter-towed natural-source airborne EM
   survey: an airborne receiver bird flies a grid of flight lines
   while natural audio-frequency energy illuminates the whole survey
   area from above. AFMAG, ZTEM, and MobileMT differ in how many field
   components the bird carries, whether a temporary ground base
   station is required, and whether that ground station measures a
   magnetic or an electric field -- the distinctions this page works
   through below.

Because there is no controlled source, an airborne AFMAG-family
survey also inherits the caution :doc:`field_zones` raises for ground
MT/AMT rather than its own: the relevant risk is not a near-field
bias but static-shift-like galvanic distortion and, uniquely to the
airborne case, motion-induced noise from the bird's own attitude
(covered below).

From A Scalar Tilt Angle To A Tensor
------------------------------------

Ward (1959) introduced AFMAG as a two-coil comparator: one coil tuned
to the local vertical, one to the horizontal, with the operator
reading a single tilt/deflection angle at each station. That angle
carries real information -- it responds to the same subsurface
conductivity contrasts a modern tensor system does -- but it is a
scalar, with no polarization-ellipse decomposition and no way to
separate the response into orthogonal components. :mod:`pycsamt.airborne.afmag`
keeps this historical reading as its own first-class shape rather
than forcing it into a tensor container it cannot honestly fill; see
:ref:`emtools_afmag` for a worked example built on a real digitised
Ward (1959) field curve.

Pedersen et al. (1994) took the conceptually next step for the VLF
band: a tensor instrument relating a fixed ground-reference horizontal
magnetic field to an airborne three-component field, giving a genuine
:term:`Interstation transfer function` rather than a bare angle,

.. math::
   :label: eq-airborne-interstation-tf

   \begin{pmatrix} H_x \\ H_y \\ H_z \end{pmatrix}_{r}
   =
   \begin{pmatrix}
   H_{xx} & H_{xy} \\
   H_{yx} & H_{yy} \\
   H_{zx} & H_{zy}
   \end{pmatrix}
   \begin{pmatrix} H_x \\ H_y \end{pmatrix}_{r_0}

(the same relation :eq:`eq-afmag-interstation-tf` states for modern
tensor AFMAG/AirMt in :ref:`emtools_afmag`). ZTEM (Z-axis Tipper
Electromagnetics; Lo and Zang 2008) specializes this idea to the
audio-frequency band and to a single output row: the airborne bird
carries only a vertical-axis coil, and the ground station only the two
horizontal coils, so the matrix above collapses to its bottom row --
the familiar tipper relation

.. math::

   H_z(f) = T_{zx}(f)\,H_x(f) + T_{zy}(f)\,H_y(f)

(:eq:`eq-ztem-tipper`). This is why pyCSAMT reuses
:attr:`~pycsamt.site.base.Site.tipper` machinery for both ground
:math:`H_z` measurements and ZTEM: mathematically, once the
horizontal ground reference is fixed, an airborne :math:`H_z` reading
*is* a tipper, no matter how it was acquired.

.. figure:: /images/theory/ztem.inJeanMLegault.png
   :align: center
   :width: 75%

   The ZTEM receiver bird in flight over the Gold Springs, Nevada
   survey, with the ground base-station coils shown alongside. The
   vertical-axis :math:`H_z` receiver coil and its GPS antenna are
   towed roughly 55-85 m below and behind the helicopter; the fixed
   ground station (bottom right) records the horizontal reference
   field :math:`H_x, H_y` throughout the flight. *Reproduced from*
   [Legault2012]_ *, Fig. 4a.*

Legault et al. (2012) flew exactly this configuration over the Gold
Springs low-sulphidation epithermal gold-silver district in
south-eastern Nevada -- 41 east-west lines at 200 m spacing plus 6
tie lines, 470 line-km, at six frequencies from 30 to 720 Hz -- and
found that resistivity highs recovered from the tipper data,
converted to a pseudo-section image, clustered around every known
gold occurrence. That real survey (in the smaller, seven-line
synthetic form committed as ``gold_springs_nv``) is the worked example
:ref:`emtools_ztem` builds on. Sattel and Witherly (2012) describe
the complementary map-view processing -- total divergence, phase
rotation, and multi-line gridding -- also demonstrated there; this
theory page only summarizes the physics: the total-divergence /
:term:`Peaker` transform converts a raw tipper crossover, which is
awkward to pick automatically from a map, into a peak centred on the
causative contact,

.. math::

   \mathrm{DT} = \frac{\partial T_{zx}}{\partial x},

a relation Lo and Zang (2008) derive and that Sattel and Witherly
(2012) note coincides, for a single flight line, with Pedersen's own
VLF-era Peaker parameter. See :term:`Total divergence` and
:ref:`emtools_ztem` for the full derivation, worked numbers, and the
map-view figures.

Why Attitude, Not Just Geology, Shows Up In The Data
----------------------------------------------------

A ground AFMAG or MT station is bolted to the earth; an airborne
receiver bird is not. As the helicopter and bird yaw, pitch, and roll
along the flight line, the sensor's own coordinate frame rotates with
respect to the (much larger) steady geomagnetic field, and that
rotation projects a spurious signal into the measured components --
signal that looks, to a naive processor, exactly like a subsurface
response. Liu et al. (2018) formalize this with a sequence of
rotation matrices carrying the field from the local geographic frame
into the sensor's own frame as a function of the platform's
instantaneous yaw, pitch, and roll angles.

.. figure:: /images/theory/afmag.liu2028_2.png
   :align: center
   :width: 55%

   Geometry of motion-induced coupling: the sensor plane's tilt angle
   :math:`\alpha(t)` with respect to the local horizontal, driven by
   the airborne platform's attitude, projects a component of the
   (otherwise DC) geomagnetic field into the measured signal.
   *Reproduced from* [Liu2018]_ *, Fig. 1.*

The key practical finding is that this motion-coupled noise
concentrates at low frequency and grows with attitude amplitude and
with the angle between the flight direction and the local geomagnetic
field -- exactly the frequency band where AFMAG-family methods also
carry their deepest, most valuable signal, since low frequency means
large skin depth in :eq:`eq-maxwell-skin-depth-rule`. There is
consequently a genuine trade-off, not just an engineering
inconvenience: gyro-stabilized bird mounts and post-flight attitude
compensation can suppress the noise, but never perfectly, so the
lowest usable frequency on a real survey is set as much by platform
stability as by target depth. :mod:`pycsamt.emtools.afmag` implements
Liu et al.'s rotation-matrix physics directly
(:func:`~pycsamt.emtools.afmag.motion_susceptibility_table`,
:func:`~pycsamt.emtools.afmag.flag_motion_susceptible_band`) so that a
survey's most exposed stations and frequency bands can be flagged
before interpretation rather than after; :ref:`emtools_afmag` and
:ref:`emtools_ztem` both show this in practice; the latter also masks
frequencies outside the vendor-declared usable band directly
(:func:`~pycsamt.emtools.ztem.mask_outside_ztem_band`).

.. figure:: /images/theory/afmag.liu2018.png
   :align: center
   :width: 65%

   A real airborne AFMAG sensor package staged for a field test, with
   the EM sensor head, inertial navigation system (INS), and a
   Mag-03 fluxgate reference magnetometer visible -- the same attitude
   and field-geometry instrumentation the motion-coupling model above
   is built from. *Reproduced from* [Liu2018]_ *, Fig. 6.*

MobileMT: A Genuinely Different Physical Object
-----------------------------------------------

Every method above is, once the ground reference is fixed, a
magnetic-field-to-magnetic-field transfer function -- structurally a
tipper or an interstation tensor, never an impedance. MobileMT
(Prikhodko et al. 2022) breaks that pattern deliberately. Instead of a
second magnetic station, MobileMT keeps a **fixed ground electric
dipole pair** (:math:`E_x, E_y`) as its reference and relates it to
the airborne bird's three magnetic coils, giving a complex
**admittance** tensor,

.. math::
   :label: eq-airborne-mobilemt-admittance

   \begin{pmatrix} H_x \\ H_y \\ H_z \end{pmatrix}
   =
   \begin{pmatrix}
   Y_{xx} & Y_{xy} \\
   Y_{yx} & Y_{yy} \\
   Y_{hzx} & Y_{hzy}
   \end{pmatrix}
   \begin{pmatrix} E_x \\ E_y \end{pmatrix}.

Compare this with :eq:`eq-airborne-interstation-tf`: the matrix shape
is identical (3 rows, 2 columns), but the input columns are now an
electric field rather than a magnetic one. That single substitution
is what makes :math:`Y` the *reciprocal* of the classical MT
impedance tensor :math:`Z` (which gives :math:`E` from :math:`H`, see
:eq:`eq-imp-matrix` in :doc:`impedance_tensor`) rather than another
tipper-shaped object. It is also why pyCSAMT deliberately refuses to
convert a MobileMT admittance transfer function into an EDI: an EDI
file has nowhere honest to put a ground electric reference alongside
an airborne magnetic response, so
:mod:`pycsamt.emtools.mobilemt` carries it instead through the
:class:`~pycsamt.airborne.AirborneEMDataset` hierarchy, never through
:class:`~pycsamt.site.base.Site`.

.. figure:: /images/theory/mobileMT.alexander.png
   :align: center
   :width: 75%

   MobileMT's two-part acquisition geometry: an airborne bird carrying
   three orthogonal magnetic-field receiver coils, referenced to a
   fixed ground electric-dipole transmitter/receiver pair
   (:math:`E_x, E_y`) rather than to a second magnetic station.
   *Reproduced from* [Prikhodko2022]_ *, Fig. 2.*

In the limit where the airborne and ground sensors are co-located --
never exactly true in practice, but the limit both Zhdanov et al.
(2024) and Sattel et al. (2019) state explicitly as the theoretical
basis for MobileMT's derived apparent-resistivity product -- the
admittance tensor reduces to the ordinary MT admittance,
:math:`Y = Z^{-1}`. Applying that identity to pyCSAMT's own
Berdichevsky-determinant convention for :math:`Z`
(:class:`pycsamt.z.resphase.ResPhase`'s ``res_det``/``phase_det``,
the same :math:`\rho_a \approx 0.2\,|Z|^2/f` shortcut already labelled
:eq:`eq-imp-rho-zonge`, generalized from a single component to the
full 2x2 determinant) and substituting :math:`\det Y = 1/\det Z`
gives, by direct algebra, a theoretical apparent conductivity and
phase for the admittance tensor:

.. math::
   :label: eq-airborne-mobilemt-sigma

   \sigma_a = 5\,f\,|\det Y|,
   \qquad
   \varphi_a = -\arg\sqrt{\det Y}.

The leading constant, :math:`5 = 1/0.2`, is the same empirical Zonge-style
factor as :eq:`eq-imp-rho-zonge`, simply inverted because the tracked
quantity is now a conductivity rather than a resistivity.
:func:`~pycsamt.emtools.mobilemt.admittance_determinant_table`
implements :eq:`eq-airborne-mobilemt-sigma` directly. It is worth
being explicit about what this quantity is *not*: it is a
theoretical, co-located-sensor-limit conductivity derived from
pyCSAMT's own determinant convention, not a reproduction of
MobileMT's proprietary processed apparent-conductivity product. When a
delivered dataset already carries that vendor-processed quantity,
:mod:`pycsamt.emtools.mobilemt` reports it alongside, unmodified, as
``apparent_conductivity_native_Sm``, precisely so the two are never
silently conflated.

A second, scale-invariant diagnostic needs no physical constant at
all and is therefore safe to compute directly from any admittance
tensor, including a non-co-located one: a Swift (1967)-style skew
ratio applied to the horizontal 2x2 admittance sub-block by direct
algebraic analogy with the impedance-tensor skew already used
elsewhere in pyCSAMT (see :doc:`dimensionality`),

.. math::
   :label: eq-airborne-mobilemt-skew

   \mathrm{skew}_Y
   =
   \frac{|Y_{xx}+Y_{yy}|}{|Y_{xy}-Y_{yx}|}.

A low value indicates a tensor close to the ideal 1-D/2-D-consistent
form; a large value signals departure from that form -- instrument
coupling error, genuine 3-D structure, or non-co-located-sensor
geometry effects that :eq:`eq-airborne-mobilemt-sigma`'s co-located
assumption does not capture.
:func:`~pycsamt.emtools.mobilemt.admittance_skew_table` implements
:eq:`eq-airborne-mobilemt-skew`.

Zhdanov et al. (2024) demonstrate the practical payoff of this
formalism on a real 3-D case: jointly inverting MobileMT admittance
data with total-magnetic-intensity (TMI) magnetic data over a
Climax-style porphyry molybdenum-copper breccia pipe target in East
Greenland, recovering a coincident resistivity-low/susceptibility-high
body consistent with the known alteration geometry. The synthetic
``flammefjeld_greenland`` sample dataset committed under
``data/mobileMT/`` is loosely inspired by that survey; the KL-22
kimberlite pipe survey of Prikhodko et al. (2022) and Sattel et al.
(2019) motivates the companion ``timiskaming_kimberlite_on`` dataset.

Comparing The Three Systems
---------------------------

.. list-table::
   :header-rows: 1
   :widths: 16 20 20 22 22

   * - System
     - Source
     - Airborne channels
     - Ground reference
     - Response object
   * - Classical AFMAG (Ward 1959)
     - Natural audio-frequency
     - 1 (tilt/deflection angle)
     - None
     - Scalar angle -- no tensor
   * - Tensor AFMAG / AirMt
     - Natural audio-frequency
     - :math:`H_x, H_y, H_z`
     - Fixed :math:`H_x, H_y` station
     - :term:`Interstation transfer function`, :eq:`eq-airborne-interstation-tf`
   * - ZTEM
     - Natural audio-frequency
     - :math:`H_z` only
     - Fixed :math:`H_x, H_y` station
     - Tipper :math:`(T_{zx}, T_{zy})`, :eq:`eq-ztem-tipper`
   * - MobileMT
     - Natural audio-frequency
     - :math:`H_x, H_y, H_z`
     - Fixed :math:`E_x, E_y` dipole pair
     - Admittance :math:`Y`, :eq:`eq-airborne-mobilemt-admittance`

Every row shares the same natural, plane-wave-illuminated source and
the same absence of a controlled transmitter; the columns that differ
are exactly what changes which pyCSAMT container carries the data
(:class:`~pycsamt.site.base.Site` and :mod:`pycsamt.emtools.ztem`/
:mod:`pycsamt.emtools.afmag` for the magnetic-only rows, versus
:class:`~pycsamt.airborne.AirborneEMDataset` and
:mod:`pycsamt.emtools.mobilemt` for the admittance row) and which
derived quantities are physically meaningful for it. A tipper or
interstation tensor has no natural resistivity/phase product because
it never involves an electric field; only the admittance tensor does,
through :eq:`eq-airborne-mobilemt-sigma`.

Practical Implications
----------------------

A few consequences follow directly from the physics above and are
worth stating explicitly before using any of the three systems'
along-line diagnostics:

.. warning::

   The total-divergence, phase-rotation, and crossover diagnostics
   summarized above (and implemented in
   :mod:`pycsamt.emtools.ztem`) all project stations onto a single
   along-line bearing. Passing a multi-line survey to any of them
   without first grouping by flight line silently differentiates or
   compares across line boundaries rather than along one coherent
   profile. :func:`~pycsamt.airborne.site.AirborneSites.select` and
   the ``LineId`` survey metadata exist precisely to make that
   grouping explicit; see :ref:`emtools_ztem` for the pattern.

.. note::

   Neither classical tilt-angle AFMAG nor a bare tipper/interstation
   tensor has a closed-form apparent-resistivity relation -- unlike
   full MT, there is no electric-field channel to form a ratio
   against. Only MobileMT's admittance tensor supports
   :eq:`eq-airborne-mobilemt-sigma`, and even then only as a
   theoretical, co-located-sensor-limit quantity, not a substitute for
   a vendor-delivered product.

.. important::

   Motion-induced noise and the finite usable frequency band are not
   optional caveats: they set the practical lower bound on
   investigation depth for every system on this page, often more
   tightly than the target's own conductivity does. Always mask or
   flag frequencies outside a system's declared usable band
   (:func:`~pycsamt.emtools.ztem.mask_outside_ztem_band`,
   :func:`~pycsamt.emtools.mobilemt.mask_outside_mobilemt_band`)
   before interpreting a total-divergence map or a resistivity
   pseudo-section.

Next Steps
----------

For worked, code-driven demonstrations built on real synthetic
surveys, continue with:

* :ref:`emtools_afmag` -- both AFMAG response families, motion-noise
  diagnostics, and the amplification parameter, on a digitised
  Ward (1959) curve and a modern tensor AirMt survey.
* :ref:`emtools_ztem` -- total divergence, phase rotation, and
  multi-line map gridding, on the Gold Springs and Forrestania
  surveys.
* :mod:`pycsamt.emtools.mobilemt` -- admittance, determinant apparent
  conductivity, and skew diagnostics on the Flammefjeld and
  Timiskaming synthetic surveys.

References
----------

The physics and figures above draw directly on [Ward1959]_,
[Liu2018]_, [Lo2008]_, [Legault2012]_, [Sattel2012]_, [Pedersen1994]_,
[Prikhodko2022]_, [Sattel2019]_, [Zhdanov2024]_, and [Swift1967]_. See
:doc:`../references` for full citations.
