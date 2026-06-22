.. _impedance_tensor:

Impedance Tensor
================

The impedance tensor is the central frequency-domain response used in
MT, AMT, and many CSAMT processing workflows. It relates horizontal
electric fields to horizontal magnetic fields at each station and each
frequency. In pyCSAMT, the tensor is also the bridge between raw field
measurements, apparent resistivity/phase products, tensor diagnostics,
static-shift correction, dimensionality analysis, and inversion input
files.

This page explains the notation, physical meaning, and practical
interpretation of the impedance tensor as used throughout pyCSAMT.

Role In MT, AMT, And CSAMT
--------------------------

For natural-source MT and AMT, the working assumption is that the incident
source field behaves approximately as a plane wave over the survey area.
Under that assumption, the horizontal electric field vector
:math:`\mathbf{E}_h` is related to the horizontal magnetic field vector
:math:`\mathbf{H}_h` by a complex frequency-dependent tensor:

.. math::

   \mathbf{E}_h(\omega)
   =
   \mathbf{Z}(\omega)
   \mathbf{H}_h(\omega).

Written component by component:

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

Each entry :math:`Z_{ij}` is complex. Its amplitude controls apparent
resistivity, and its phase describes the lag between electric and magnetic
fields. Because :math:`\mathbf{Z}` varies with frequency, it encodes how
subsurface conductivity is sampled at different EM diffusion scales.

For CSAMT, the same impedance-style quantities may be computed from
controlled-source measurements, but the interpretation must account for
source-receiver geometry. In a far-field CSAMT regime, MT-like impedance
interpretation is often useful. In near-field or transition-zone regimes,
source overprint and shadow effects may make a simple impedance
interpretation misleading. See :doc:`csamt_amt_mt_overview` for the broader
method context.

Coordinate Convention
---------------------

The tensor notation depends on a coordinate system. pyCSAMT follows the
usual two-horizontal-axis convention:

* :math:`x`: first horizontal axis, commonly north or profile-aligned;
* :math:`y`: second horizontal axis, commonly east or cross-profile;
* :math:`z`: positive vertical axis by the convention used in the input data
  or processing workflow.

The component :math:`Z_{xy}` maps :math:`H_y` into :math:`E_x`, while
:math:`Z_{yx}` maps :math:`H_x` into :math:`E_y`. Component labels are
therefore not just names; they are tied to the data rotation and survey
orientation.

In pyCSAMT code, a tensor stack is stored as an array with shape
``(n_freq, 2, 2)``:

.. list-table::
   :header-rows: 1
   :widths: 24 28 48

   * - Component
     - Array index
     - Meaning
   * - :math:`Z_{xx}`
     - ``z[:, 0, 0]``
     - :math:`H_x` contribution to :math:`E_x`.
   * - :math:`Z_{xy}`
     - ``z[:, 0, 1]``
     - :math:`H_y` contribution to :math:`E_x`.
   * - :math:`Z_{yx}`
     - ``z[:, 1, 0]``
     - :math:`H_x` contribution to :math:`E_y`.
   * - :math:`Z_{yy}`
     - ``z[:, 1, 1]``
     - :math:`H_y` contribution to :math:`E_y`.

This convention is used by :class:`pycsamt.z.z.Z`,
:func:`pycsamt.seg.ops.z_to_rho_phi`, and tensor-related processing steps.

Complex Numbers And Phase
-------------------------

Each impedance component can be written in rectangular form:

.. math::

   Z_{ij} = \operatorname{Re}(Z_{ij}) +
            i \operatorname{Im}(Z_{ij}),

or polar form:

.. math::

   Z_{ij} = |Z_{ij}| e^{i \phi_{ij}}.

The component phase is:

.. math::

   \phi_{ij}
   =
   \tan^{-1}
   \left(
   \frac{\operatorname{Im}(Z_{ij})}
        {\operatorname{Re}(Z_{ij})}
   \right).

In implementation, this is normally computed with a quadrant-aware function
such as ``atan2`` or ``np.angle``. pyCSAMT reports phase in degrees for
user-facing apparent resistivity/phase products.

Apparent Resistivity
--------------------

The apparent resistivity associated with a component :math:`Z_{ij}` is:

.. math::

   \rho_{a,ij}
   =
   \frac{|Z_{ij}|^2}{\mu_0 \omega},

where :math:`\mu_0` is magnetic permeability of free space and
:math:`\omega = 2 \pi f`.

With :math:`\mu_0 = 4\pi \times 10^{-7}` H/m, this is commonly written in
the numerical form used by pyCSAMT's internal utilities:

.. math::

   \rho_{a,ij}
   \approx
   0.2 \frac{|Z_{ij}|^2}{f},

when :math:`f` is in Hz and the impedance units follow the package
convention. The inverse relationship is:

.. math::

   |Z_{ij}|
   =
   \sqrt{5 f \rho_{a,ij}}.

This is why pyCSAMT can reconstruct complex impedance from apparent
resistivity and phase:

.. math::

   Z_{ij}
   =
   \sqrt{5 f \rho_{a,ij}}
   \left[
   \cos(\phi_{ij}) + i \sin(\phi_{ij})
   \right].

Apparent resistivity is not the true resistivity at one depth. It is the
resistivity of a uniform half-space that would produce the same component
response at that frequency. In layered or laterally varying Earth
structures, one frequency samples a broad sensitivity volume rather than a
single point.

Diagonal And Off-Diagonal Components
------------------------------------

The four tensor components carry different interpretation weight depending
on Earth dimensionality and coordinate rotation.

For a simple 1-D Earth, after using a consistent horizontal coordinate
system, the ideal tensor has zero diagonal terms and antisymmetric
off-diagonal terms:

.. math::

   \mathbf{Z}_{1D}
   =
   \begin{bmatrix}
   0 & Z_{xy} \\
   -Z_{xy} & 0
   \end{bmatrix}.

For a 2-D Earth rotated into strike coordinates, the diagonal terms are
also ideally zero, while :math:`Z_{xy}` and :math:`Z_{yx}` represent the
two principal modes:

.. math::

   \mathbf{Z}_{2D}
   =
   \begin{bmatrix}
   0 & Z_{\mathrm{TE}} \\
   Z_{\mathrm{TM}} & 0
   \end{bmatrix}.

The exact assignment of TE and TM depends on coordinate convention and
strike orientation, so pyCSAMT documentation prefers explicit component
names unless a workflow has clearly defined the strike frame.

For a 3-D Earth, all four entries may be significant:

.. math::

   Z_{xx} \ne 0,
   \quad
   Z_{yy} \ne 0,
   \quad
   Z_{xy} \ne -Z_{yx}.

Large diagonal components, inconsistent strike estimates, strong tipper
behavior, or high skew values can indicate that a simple 1-D or 2-D
interpretation is not reliable.

Tensor Rotation
---------------

Rotating the horizontal coordinate system changes the impedance tensor. In
pyCSAMT, the rotation rule used by SEG-style operations is:

.. math::

   \mathbf{Z}'
   =
   \mathbf{R}(\theta)
   \mathbf{Z}
   \mathbf{R}^{T}(\theta),

where:

.. math::

   \mathbf{R}(\theta)
   =
   \begin{bmatrix}
   \cos\theta & \sin\theta \\
   -\sin\theta & \cos\theta
   \end{bmatrix}.

The same angle can be applied to all frequencies, or a frequency-dependent
angle can be used when the processing workflow allows it. In pyCSAMT this
operation appears in functions such as
:func:`pycsamt.seg.ops.rotate_impedance` and tensor pipeline steps.

Rotation is not cosmetic. It can change:

* the size of diagonal terms;
* the relationship between :math:`Z_{xy}` and :math:`Z_{yx}`;
* whether a profile looks approximately 2-D;
* which component is used for a particular inversion mode;
* the apparent continuity of pseudosection features.

For this reason, every rotation should record the angle, convention, and
reason. A plot made after rotation must not be compared casually to a plot
made before rotation.

Strike, Principal Directions, And Ambiguity
-------------------------------------------

In 2-D interpretation, the goal is often to rotate the tensor into a
coordinate frame aligned with geological strike. In a perfect 2-D case,
there exists an angle :math:`\theta_s` such that the diagonal terms are
minimized and the off-diagonal modes are physically interpretable.

In practice, strike estimation is ambiguous because:

* MT strike often has a :math:`90^\circ` ambiguity;
* near-surface distortion can affect impedance components;
* different frequency bands may sample different structures;
* 3-D geology may not have one stable strike;
* CSAMT source geometry can bias component behavior.

pyCSAMT strike and tensor diagnostics should therefore be treated as
evidence, not as automatic truth. A robust workflow compares strike
estimates with geology, station layout, tipper vectors, phase tensors,
residuals, and inversion behavior.

Phase Tensor
------------

The phase tensor is a distortion-resistant diagnostic constructed from the
real and imaginary parts of the impedance tensor. Let:

.. math::

   \mathbf{Z} = \mathbf{X} + i \mathbf{Y},

where :math:`\mathbf{X} = \operatorname{Re}(\mathbf{Z})` and
:math:`\mathbf{Y} = \operatorname{Im}(\mathbf{Z})`. The phase tensor is:

.. math::

   \boldsymbol{\Phi}
   =
   \mathbf{X}^{-1}
   \mathbf{Y}.

It describes the phase relationship independent of galvanic amplitude
scaling under common distortion assumptions. Phase tensor ellipses are used
to inspect:

* dimensionality;
* strike direction;
* phase anisotropy;
* skew or 3-D behavior;
* station-to-station consistency;
* period-dependent structural changes.

The phase tensor does not remove every problem. It does not make poor data
good, it does not eliminate inductive 3-D effects, and it does not replace
model-based interpretation. It is a diagnostic that helps decide which
processing and inversion assumptions are defensible.

Skew And Dimensionality
-----------------------

Skew metrics summarize departures from simple dimensionality. Different
families of skew exist, and they should not be confused with the simple
linear-algebra skew :math:`Z_{xy} - Z_{yx}`.

As a practical guide:

* low skew over a continuous band supports 1-D or 2-D modelling assumptions;
* high skew may indicate 3-D structure, distortion, bad data, or source
  effects;
* isolated high-skew frequencies may reflect noise or processing artifacts;
* skew should be compared with phase tensor, tipper, and residual plots.

pyCSAMT pipeline steps can mask, select, or close gaps in low-skew frequency
bands. Those operations should be documented in the workflow report because
they affect which data enter inversion.

Tipper And Vertical Magnetic Field
----------------------------------

The impedance tensor uses horizontal electric and magnetic fields. The
vertical magnetic field is commonly described with the tipper:

.. math::

   H_z
   =
   T_x H_x + T_y H_y.

The tipper vector :math:`\mathbf{T} = [T_x, T_y]` is sensitive to lateral
conductivity contrasts. Large tipper magnitudes or coherent induction
vectors often indicate 2-D or 3-D structure. In a laterally uniform 1-D
Earth, the ideal tipper is zero.

Under horizontal rotation, pyCSAMT rotates the tipper as a vector:

.. math::

   \mathbf{T}'
   =
   \mathbf{R}(\theta)
   \mathbf{T}.

Tipper plots are especially useful when deciding whether a 2-D profile
orientation is plausible. If induction vectors suggest strong off-profile
structure, a 2-D inversion may still fit the data numerically but be
geologically misleading.

Error Propagation
-----------------

Impedance estimates should be interpreted with uncertainty. If
:math:`\Delta Z_{ij}` is an absolute component uncertainty, then the
relative amplitude error is approximately:

.. math::

   \frac{\Delta |Z_{ij}|}{|Z_{ij}|}.

Because apparent resistivity depends on :math:`|Z|^2`, the relative
resistivity error is approximately twice the relative impedance-amplitude
error:

.. math::

   \frac{\Delta \rho_a}{\rho_a}
   \approx
   2
   \frac{\Delta |Z|}{|Z|}.

Phase uncertainty depends on the relative size of the real and imaginary
parts and becomes unstable when the impedance amplitude is very small.
pyCSAMT carries absolute ``z_err`` arrays when available and propagates them
to resistivity and phase through the :class:`pycsamt.z.resphase.ResPhase`
and :class:`pycsamt.z.z.Z` containers.

Error information matters for:

* deciding whether a curve feature is real or noise;
* weighting inversion data;
* masking unstable frequencies;
* comparing observed and predicted responses;
* judging whether an apparent dimensionality feature is significant.

Static Shift And Galvanic Distortion
------------------------------------

Static shift is a near-surface galvanic distortion that multiplies apparent
resistivity without changing phase in the simplest approximation. It is not
primarily an impedance phase problem; it is an amplitude scaling problem.

If :math:`g_x` and :math:`g_y` are real electric-field gains, a simplified
distortion model can be written:

.. math::

   \mathbf{E}_{obs}
   =
   \mathbf{C}
   \mathbf{E}_{true},

where :math:`\mathbf{C}` is a real distortion matrix. Then:

.. math::

   \mathbf{Z}_{obs}
   =
   \mathbf{C}
   \mathbf{Z}_{true}.

This is why static shift can move apparent resistivity curves up or down
without producing an equivalent phase shift. The detailed interpretation is
covered in :doc:`static_shift`, but the tensor view is important: correcting
static shift changes amplitudes before apparent resistivity is interpreted
or inverted.

Determinant And Other Invariants
--------------------------------

Tensor invariants are scalar combinations of :math:`\mathbf{Z}` that can be
useful for plotting, comparison, or quality control. Common examples are:

.. math::

   \operatorname{tr}(\mathbf{Z})
   =
   Z_{xx} + Z_{yy},

.. math::

   \det(\mathbf{Z})
   =
   Z_{xx} Z_{yy} - Z_{xy} Z_{yx},

and the Frobenius norm:

.. math::

   \lVert \mathbf{Z} \rVert_F
   =
   \sqrt{
   |Z_{xx}|^2 + |Z_{xy}|^2 + |Z_{yx}|^2 + |Z_{yy}|^2
   }.

pyCSAMT exposes such quantities through the :class:`pycsamt.z.z.Z`
container. They are helpful diagnostics, but they are not a replacement for
understanding tensor orientation, uncertainty, and dimensionality.

Practical Component Choices
---------------------------

Many workflows focus on off-diagonal components because they are dominant
for ideal 1-D and 2-D MT responses. Typical choices include:

.. list-table::
   :header-rows: 1
   :widths: 24 36 40

   * - Choice
     - Use
     - Caution
   * - :math:`Z_{xy}`
     - One principal off-diagonal response.
     - Meaning depends on coordinate frame and strike convention.
   * - :math:`Z_{yx}`
     - The other principal off-diagonal response.
     - Often has opposite sign or phase convention relative to
       :math:`Z_{xy}`.
   * - Determinant response
     - Rotation-tolerant summary response for some quick-look plots.
     - Can hide mode-specific behavior and poor component quality.
   * - Average off-diagonal response
     - Simple 1-D-style sounding or quick QC.
     - Can blur TE/TM differences and 2-D structure.
   * - Phase tensor quantities
     - Dimensionality and strike diagnostics.
     - Not a replacement for amplitude-sensitive inversion.

The safest documentation and reports should always state the component or
invariant used. "Apparent resistivity" alone is incomplete; say
:math:`\rho_{a,xy}`, :math:`\rho_{a,yx}`, determinant apparent resistivity,
or the exact product used.

pyCSAMT Containers And Utilities
--------------------------------

Several pyCSAMT modules use the tensor concepts directly.

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - Object or function
     - Role
   * - :class:`pycsamt.z.z.Z`
     - High-level impedance tensor container with shape ``(n_freq, 2, 2)``,
       frequency vector, errors, rotation history, resistivity/phase
       conversion, and invariants.
   * - :class:`pycsamt.z.resphase.ResPhase`
     - Converts between complex impedance and apparent resistivity/phase,
       including error propagation.
   * - :class:`pycsamt.z.tipper.Tipper`
     - Stores complex tipper components, rotations, amplitudes, phases, and
       arrow metrics.
   * - :func:`pycsamt.seg.ops.z_to_rho_phi`
     - Computes component-wise apparent resistivity and phase from
       impedance and frequency.
   * - :func:`pycsamt.seg.ops.rho_phi_to_z`
     - Reconstructs complex impedance from apparent resistivity, phase, and
       frequency.
   * - :func:`pycsamt.seg.ops.rotate_impedance`
     - Applies :math:`\mathbf{Z}' = \mathbf{R}\mathbf{Z}\mathbf{R}^T`.
   * - :func:`pycsamt.seg.ops.rotate_tipper`
     - Rotates tipper vectors consistently with the horizontal coordinate
       frame.

Minimal Python Example
----------------------

The following example builds a simple impedance tensor, computes apparent
resistivity and phase, and rotates the tensor.

.. code-block:: python
   :linenos:

   import numpy as np

   from pycsamt.seg.ops import rotate_impedance, z_to_rho_phi

   freq = np.array([10.0, 1.0, 0.1])
   z = np.zeros((freq.size, 2, 2), dtype=complex)

   # Synthetic off-diagonal response.
   z[:, 0, 1] = 1.0 + 1.0j
   z[:, 1, 0] = -1.0 - 1.0j

   rho, phase = z_to_rho_phi(z, freq)
   z_rot = rotate_impedance(z, theta_deg=30.0)

   print(rho[:, 0, 1])
   print(phase[:, 0, 1])
   print(z_rot.shape)

The same operations can be performed through :class:`pycsamt.z.z.Z` when
you want a stateful container with errors, metadata, and derived properties.

Quality Control Questions
-------------------------

Before using impedance tensors in inversion or interpretation, ask:

1. Are frequencies positive, sorted, and consistent across stations?
2. Are tensor components in the expected coordinate convention?
3. Were rotations applied, and are the rotation angles recorded?
4. Are :math:`Z_{xy}` and :math:`Z_{yx}` internally consistent?
5. Are diagonal terms physically meaningful, or dominated by noise?
6. Are phase curves smooth enough to trust?
7. Are apparent resistivity shifts consistent with static shift?
8. Does the phase tensor suggest 1-D, 2-D, or 3-D behavior?
9. Do tipper vectors support the selected profile orientation?
10. Are errors available for inversion weighting?

Interpretation Pitfalls
-----------------------

Common tensor mistakes include:

* treating the tensor as four independent curves with no coordinate context;
* rotating data without documenting the angle and sign convention;
* assuming diagonal components are always noise;
* assuming off-diagonal antisymmetry guarantees a 1-D Earth;
* averaging :math:`Z_{xy}` and :math:`Z_{yx}` before checking mode behavior;
* interpreting apparent resistivity amplitudes without static-shift review;
* trusting a low RMS inversion while ignoring tensor residual patterns;
* using CSAMT impedance-style products without checking source effects.

Connection To Later Workflows
-----------------------------

The impedance tensor feeds directly into:

* :doc:`static_shift`, because galvanic distortion affects amplitude;
* :doc:`inversion_concepts`, because inversion fits impedance,
  resistivity/phase, or transformed tensor products;
* :doc:`../pipeline/steps`, because tensor cleanup, skew gates, rotations,
  and dimensionality diagnostics are registered pipeline steps;
* :doc:`../tutorials/correct_static_shift`, because practical correction
  begins with apparent resistivity and phase inspection;
* :doc:`../api/z`, because the public tensor containers are documented in
  the API reference.

References
----------

This page follows standard MT tensor notation and the implementation
conventions in pyCSAMT. Relevant background references are collected in
:doc:`../references`, especially [WardHohmann1988]_,
[deGrootHedlin1990]_, and [Kelbert2014]_.
