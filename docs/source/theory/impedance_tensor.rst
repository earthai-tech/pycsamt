.. _impedance_tensor:

Impedance Tensor
================

The :term:`impedance tensor` is the central frequency-domain response used in
:term:`MT`, :term:`AMT`, and many :term:`CSAMT` processing workflows. It
relates horizontal electric fields to horizontal magnetic fields at each
station and each frequency. In pyCSAMT, the tensor is also the bridge between
raw field measurements, apparent resistivity/phase products, tensor
diagnostics, static-shift correction, dimensionality analysis, and inversion
input files.

This page explains the notation, physical meaning, and practical
interpretation of the impedance tensor as used throughout pyCSAMT, continuing
from the family-of-methods picture in :doc:`csamt_amt_mt_overview` and the
numeric constants in :doc:`constants`.

Role In MT, AMT, And CSAMT
--------------------------

For natural-source MT and AMT, the working assumption is that the incident
source field behaves approximately as a :term:`plane-wave field` over the
survey area. Under that assumption, the horizontal electric field vector
:math:`\mathbf{E}_h` is related to the horizontal magnetic field vector
:math:`\mathbf{H}_h` by a complex frequency-dependent tensor:

.. math::
   :label: eq-imp-eh-zh

   \mathbf{E}_h(\omega)
   =
   \mathbf{Z}(\omega)
   \mathbf{H}_h(\omega).

Written component by component:

.. math::
   :label: eq-imp-matrix

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
:term:`source overprint` and shadow effects may make a simple impedance
interpretation misleading -- see :doc:`csamt_amt_mt_overview` for the
:func:`~pycsamt.iot.edge_csamt.classify_field_zones` diagnostic that makes
that distinction concrete.

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
   :label: eq-imp-rectangular

   Z_{ij} = \operatorname{Re}(Z_{ij}) +
            i \operatorname{Im}(Z_{ij}),

or polar form:

.. math::
   :label: eq-imp-polar

   Z_{ij} = |Z_{ij}| e^{i \phi_{ij}}.

The component phase is:

.. math::
   :label: eq-imp-phase

   \phi_{ij}
   =
   \tan^{-1}
   \left(
   \frac{\operatorname{Im}(Z_{ij})}
        {\operatorname{Re}(Z_{ij})}
   \right).

In implementation, this is normally computed with a quadrant-aware function
such as ``atan2`` or ``np.angle``. pyCSAMT reports :term:`phase` in degrees
for user-facing apparent resistivity/phase products.

Apparent Resistivity
--------------------

The apparent resistivity associated with a component :math:`Z_{ij}` is:

.. math::
   :label: eq-imp-rho-si

   \rho_{a,ij}
   =
   \frac{|Z_{ij}|^2}{\mu_0 \omega},

where :math:`\mu_0` is magnetic permeability of free space and
:math:`\omega = 2 \pi f`.

With :math:`\mu_0 = 4\pi \times 10^{-7}` H/m, this is commonly written in
the numerical form used by pyCSAMT's internal utilities:

.. math::
   :label: eq-imp-rho-zonge

   \rho_{a,ij}
   \approx
   0.2 \frac{|Z_{ij}|^2}{f},

when :math:`f` is in Hz and the impedance units follow the legacy
Zonge-style convention (:term:`apparent resistivity`'s glossary entry gives
the same shortcut). The inverse relationship is:

.. math::
   :label: eq-imp-z-from-rho

   |Z_{ij}|
   =
   \sqrt{5 f \rho_{a,ij}}.

This is why pyCSAMT can reconstruct complex impedance from apparent
resistivity and phase:

.. math::
   :label: eq-imp-z-reconstruct

   Z_{ij}
   =
   \sqrt{5 f \rho_{a,ij}}
   \left[
   \cos(\phi_{ij}) + i \sin(\phi_{ij})
   \right].

:eq:`eq-imp-rho-si` and :eq:`eq-imp-rho-zonge` are algebraically the same
physics under a unit conversion -- :doc:`constants` derives exactly this
equivalence from ``RHO_FACTOR`` and ``ZONGE_RHO_FACTOR``. What is easy to
miss is that pyCSAMT does **not** apply that conversion automatically
between its own two resistivity code paths. :func:`pycsamt.seg.ops.z_to_rho_phi`
implements :eq:`eq-imp-rho-si` directly (:math:`\rho_a = |Z|^2 \cdot
\mathrm{RHO\_FACTOR}/f`), while :class:`pycsamt.z.z.Z` -- built on
:class:`pycsamt.z.resphase.ResPhase`, the container actually populated by
the EDI, Zonge, and Jones readers -- implements :eq:`eq-imp-rho-zonge`
directly (:math:`\rho_a = 0.2\,|Z|^2/f`). Feeding the *same* raw ``z`` array
to both, without converting units first, does not raise an error; it
silently returns two resistivities that differ by a fixed factor:

.. code-block:: pycon

   >>> import numpy as np
   >>> from pycsamt.seg.ops import z_to_rho_phi
   >>> from pycsamt.z.z import Z
   >>> from pycsamt.constants import RHO_FACTOR, ZONGE_RHO_FACTOR
   >>> freq = np.array([10.0, 1.0, 0.1])
   >>> z = np.zeros((freq.size, 2, 2), dtype=complex)
   >>> z[:, 0, 1] = 1.0 + 1.0j
   >>> z[:, 1, 0] = -1.0 - 1.0j
   >>> rho_si, phase_si = z_to_rho_phi(z, freq)
   >>> rho_si[:, 0, 1]
   array([  25330.29591058,  253302.95910584, 2533029.59105844])
   >>> zc = Z(z_array=z, freq=freq)
   >>> zc.res_xy
   array([0.04, 0.4 , 4.  ])
   >>> round(float((rho_si[:, 0, 1] / zc.res_xy)[0]), 2)
   633257.4
   >>> round(RHO_FACTOR / ZONGE_RHO_FACTOR, 2)
   633257.4

The ratio is not a coincidence or a rounding artefact -- it is exactly
``RHO_FACTOR / ZONGE_RHO_FACTOR``, confirmed above to six significant
figures both ways. Unlike the skin-depth constant duplication in
:doc:`csamt_amt_mt_overview` (four numbers agreeing to a part in ten
thousand), this is a genuine convention split: :func:`~pycsamt.seg.ops.z_to_rho_phi` expects
:math:`Z` already in SI ohms, while :class:`~pycsamt.z.z.Z`/:class:`~pycsamt.z.resphase.ResPhase` expect the
legacy field-style scaling baked into the ``0.2`` factor, and pyCSAMT has no
code path that converts between the two automatically. In practice this
rarely bites, because :func:`~pycsamt.seg.ops.z_to_rho_phi`/:func:`~pycsamt.seg.ops.rho_phi_to_z` are
lightweight standalone helpers not currently called anywhere in the
EDI/Zonge/Jones ingestion pipeline -- but it means the two are only
comparable after an explicit conversion like the one worked out in
:doc:`constants`, never by passing the same array to both.

:term:`Apparent resistivity` is not the true resistivity at one depth. It is
the resistivity of a uniform :term:`half-space` that would produce the same
component response at that frequency. In layered or laterally varying Earth
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
   :label: eq-imp-z1d

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
   :label: eq-imp-z2d

   \mathbf{Z}_{2D}
   =
   \begin{bmatrix}
   0 & Z_{\mathrm{TE}} \\
   Z_{\mathrm{TM}} & 0
   \end{bmatrix}.

The exact assignment of :term:`TE mode` and :term:`TM mode` depends on
coordinate convention and strike orientation, so pyCSAMT documentation
prefers explicit component names unless a workflow has clearly defined the
strike frame.

For a 3-D Earth, all four entries may be significant:

.. math::
   :label: eq-imp-z3d

   Z_{xx} \ne 0,
   \quad
   Z_{yy} \ne 0,
   \quad
   Z_{xy} \ne -Z_{yx}.

Large diagonal components, inconsistent strike estimates, strong tipper
behavior, or high skew values can indicate that a simple 1-D or 2-D
interpretation is not reliable. The next section builds exactly these three
cases (:eq:`eq-imp-z1d`, :eq:`eq-imp-z2d`) numerically, rather than only
asserting how they behave under rotation.

Tensor Rotation
---------------

Rotating the horizontal coordinate system changes the impedance tensor. In
pyCSAMT, the rotation rule used by SEG-style operations is:

.. math::
   :label: eq-imp-rotation

   \mathbf{Z}'
   =
   \mathbf{R}(\theta)
   \mathbf{Z}
   \mathbf{R}^{T}(\theta),

where:

.. math::
   :label: eq-imp-rot-matrix

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

Rotation is not cosmetic. A genuinely 2-D tensor (:eq:`eq-imp-z2d`, unequal
:math:`Z_{\mathrm{TE}}` and :math:`Z_{\mathrm{TM}}`) picks up diagonal terms
as soon as it is rotated away from its native strike, while a perfectly 1-D
tensor (:eq:`eq-imp-z1d`) does not change at all -- rotating an isotropic
response cannot manufacture a preferred direction that was never there:

.. code-block:: pycon

   >>> import numpy as np
   >>> from pycsamt.seg.ops import rotate_impedance
   >>> z1d = np.zeros((1, 2, 2), dtype=complex)
   >>> z1d[0, 0, 1] = 1.0 + 1.0j
   >>> z1d[0, 1, 0] = -1.0 - 1.0j
   >>> z2d = np.zeros((1, 2, 2), dtype=complex)
   >>> z2d[0, 0, 1] = 1.0 + 1.0j   # Z_TE at native strike
   >>> z2d[0, 1, 0] = 3.0 + 0.5j   # Z_TM at native strike
   >>> det0 = z2d[0, 0, 0] * z2d[0, 1, 1] - z2d[0, 0, 1] * z2d[0, 1, 0]
   >>> det0
   (-2.5-3.5j)
   >>> z_rot20 = rotate_impedance(z2d, theta_deg=20.0)
   >>> z_rot20.shape
   (2, 2)
   >>> np.allclose(np.linalg.det(z_rot20), det0)
   True
   >>> np.allclose(np.trace(z_rot20), 0.0)
   True
   >>> np.allclose(z_rot20[0, 1] - z_rot20[1, 0], z2d[0, 0, 1] - z2d[0, 1, 0])
   True
   >>> angles = np.linspace(-90, 90, 181)

Note that ``z_rot20.shape`` is ``(2, 2)``, not ``(1, 2, 2)`` --
:func:`~pycsamt.seg.ops.rotate_impedance` squeezes any length-1 axis in its result, including
the frequency axis when only one frequency is given. Index the rotated
matrix directly rather than adding an extra ``[0]``.

The three ``True`` results above are not particular to this example: the
determinant, the trace, and the antisymmetric part :math:`Z_{xy}-Z_{yx}`
are each exactly invariant under :eq:`eq-imp-rotation` for *any* tensor,
for three different linear-algebra reasons -- :math:`\det(\mathbf{R}\mathbf{Z}\mathbf{R}^T)
=\det(\mathbf{Z})` because :math:`\det(\mathbf{R})=1`; the trace is
invariant under any orthogonal similarity transform; and a pure rotation
preserves a 2-D antisymmetric matrix exactly. What *does* change is visible
directly: the two diagonal terms grow individually away from strike even
though their sum never does.

.. code-block:: python
   :linenos:

   import numpy as np
   import matplotlib.pyplot as plt
   from pycsamt.seg.ops import rotate_impedance

   angles = np.linspace(-90, 90, 181)

   z1d = np.zeros((1, 2, 2), dtype=complex)
   z1d[0, 0, 1] = 1.0 + 1.0j
   z1d[0, 1, 0] = -1.0 - 1.0j
   zr1 = rotate_impedance(z1d, theta_deg=angles)
   diag1 = np.abs(zr1[:, 0, 0]) + np.abs(zr1[:, 1, 1])

   z2d = np.zeros((1, 2, 2), dtype=complex)
   z2d[0, 0, 1] = 1.0 + 1.0j
   z2d[0, 1, 0] = 3.0 + 0.5j
   zr2 = rotate_impedance(z2d, theta_deg=angles)
   diag2 = np.abs(zr2[:, 0, 0]) + np.abs(zr2[:, 1, 1])

   fig, ax = plt.subplots(1, 1, figsize=(7.0, 4.2))
   ax.plot(angles, diag2, color="#1f77b4", lw=2.2, label="2-D tensor (Zxy != Zyx)")
   ax.plot(angles, diag1, color="#d62728", lw=2.2, ls="--", label="1-D tensor (Zxy = -Zyx)")
   ax.axvline(0.0, color="0.5", lw=1.0, ls=":")
   ax.set(xlabel=r"Rotation angle $\theta$ (degrees)",
          ylabel=r"$|Z_{xx}| + |Z_{yy}|$",
          title="Diagonal growth under rotation")
   ax.legend(fontsize=9)
   ax.grid(True, alpha=0.3)
   fig.tight_layout()

.. figure:: /images/theory/impedance_rotation_diagonal.png
   :alt: Diagonal magnitude sum versus rotation angle for a 1-D and a 2-D synthetic impedance tensor
   :width: 100%

   The 2-D tensor's diagonal grows to a maximum near :math:`\pm 45^\circ`
   and returns to zero every :math:`90^\circ` -- rotating by a quarter turn
   swaps the axes and lands on another zero-diagonal frame with TE and TM
   exchanged. The 1-D tensor stays at floating-point zero for every angle:
   there is no strike to find because the response has none.

Strike, Principal Directions, And Ambiguity
-------------------------------------------

In 2-D interpretation, the goal is often to rotate the tensor into a
coordinate frame aligned with geological :term:`strike`. In a perfect 2-D
case, there exists an angle :math:`\theta_s` such that the diagonal terms
are minimized and the off-diagonal modes are physically interpretable. That
angle can be found directly by scanning :eq:`eq-imp-rotation` over a range
of trial angles and keeping the one with the smallest diagonal magnitude --
exactly the curve plotted above, read off numerically instead of by eye:

.. code-block:: pycon

   >>> z_scan = rotate_impedance(z2d, theta_deg=angles)
   >>> diag_mag = np.abs(z_scan[:, 0, 0]) + np.abs(z_scan[:, 1, 1])
   >>> round(float(angles[np.argmin(diag_mag)]), 1)
   0.0
   >>> round(float(diag_mag.min()), 6)
   0.0

The recovered strike, :math:`0.0^\circ`, matches the frame the synthetic
tensor was built in -- as expected, since :math:`\mathtt{z2d}` was
constructed already in strike coordinates. Real data will not return an
exact zero; the same grid search still applies, but the achievable minimum
reflects noise, 3-D structure, and distortion rather than floating-point
precision.

In practice, strike estimation is ambiguous because:

* MT strike often has a :math:`90^\circ` ambiguity, visible in the figure
  above as the second zero-diagonal frame near the scan's edges;
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

The :term:`phase tensor` is a distortion-resistant diagnostic constructed
from the real and imaginary parts of the impedance tensor. Let:

.. math::
   :label: eq-imp-phasetensor-def

   \mathbf{Z} = \mathbf{X} + i \mathbf{Y},

where :math:`\mathbf{X} = \operatorname{Re}(\mathbf{Z})` and
:math:`\mathbf{Y} = \operatorname{Im}(\mathbf{Z})`. The phase tensor is:

.. math::
   :label: eq-imp-phasetensor

   \boldsymbol{\Phi}
   =
   \mathbf{X}^{-1}
   \mathbf{Y}.

The defining property that makes :eq:`eq-imp-phasetensor` useful is easy to
lose sight of algebraically but obvious numerically: :math:`\mathbf{\Phi}`
is a **real**-valued matrix, even though :math:`\mathbf{Z}` is complex,
because it is built from the product of two real matrices:

.. code-block:: pycon

   >>> X, Y = z2d[0].real, z2d[0].imag
   >>> Phi = np.linalg.inv(X) @ Y
   >>> Phi
   array([[0.16666667, 0.        ],
          [0.        , 1.        ]])
   >>> Phi.dtype == np.dtype("float64")
   True

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

Skew metrics summarize departures from simple :term:`dimensionality`.
Different families of skew exist, and they should not be confused with the
simple linear-algebra skew :math:`Z_{xy} - Z_{yx}` used as a rotation
invariant above. A standard distortion-tolerant choice is Bahr skewness,
built from the symmetric and antisymmetric combinations of :math:`\mathbf{Z}`:

.. math::
   :label: eq-imp-bahr

   S_1 = Z_{xx}+Z_{yy}, \quad
   S_2 = Z_{xy}-Z_{yx}, \quad
   D_1 = Z_{xx}-Z_{yy}, \quad
   D_2 = Z_{xy}+Z_{yx},

.. math::
   :label: eq-imp-bahr-eta

   \eta_B = \sqrt{\frac{|S_1|^2 + |S_2|^2}{|D_1|^2 + |D_2|^2}}.

:func:`pycsamt.emtools.skew.bahr_skewness` implements
:eq:`eq-imp-bahr`/:eq:`eq-imp-bahr-eta` directly. Run on the same
:math:`\mathtt{z2d}` tensor before and after the :math:`20^\circ` rotation
used above, it returns the identical value both times -- Bahr skewness is a
rotational invariant by construction, just like the determinant and trace
already confirmed:

.. code-block:: pycon

   >>> from pycsamt.emtools.skew import bahr_skewness
   >>> eta0 = bahr_skewness(z2d)
   >>> eta20 = bahr_skewness(z_rot20[np.newaxis, ...])
   >>> round(float(eta0[0]), 6), round(float(eta20[0]), 6)
   (0.482573, 0.482573)

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
vertical magnetic field is commonly described with the :term:`tipper`:

.. math::
   :label: eq-imp-tipper

   H_z
   =
   T_x H_x + T_y H_y.

The tipper vector :math:`\mathbf{T} = [T_x, T_y]` is sensitive to lateral
conductivity contrasts. Large tipper magnitudes or coherent induction
vectors often indicate 2-D or 3-D structure. In a laterally uniform 1-D
Earth, the ideal tipper is zero.

Under horizontal rotation, pyCSAMT rotates the tipper as a vector:

.. math::
   :label: eq-imp-tipper-rotation

   \mathbf{T}'
   =
   \mathbf{R}(\theta)
   \mathbf{T}.

Because :eq:`eq-imp-tipper-rotation` is a plain vector rotation (not the
two-sided :eq:`eq-imp-rotation` used for the tensor), its Euclidean norm
:math:`\sqrt{|T_x|^2+|T_y|^2}` is preserved exactly -- rotation redistributes
the tipper between components without changing its overall magnitude:

.. code-block:: pycon

   >>> from pycsamt.seg.ops import rotate_tipper
   >>> t = np.array([[0.3 + 0.1j, 0.05 - 0.02j]])  # Tx, Ty at one frequency
   >>> t_rot = rotate_tipper(t, theta_deg=30.0)
   >>> t_rot.shape
   (2,)
   >>> norm_before = np.hypot(abs(t[0, 0]), abs(t[0, 1]))
   >>> norm_after = np.hypot(abs(t_rot[0]), abs(t_rot[1]))
   >>> round(float(norm_before), 6) == round(float(norm_after), 6)
   True

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
   :label: eq-imp-z-rel-err

   \frac{\Delta |Z_{ij}|}{|Z_{ij}|}.

Because apparent resistivity depends on :math:`|Z|^2`, the relative
resistivity error is approximately twice the relative impedance-amplitude
error:

.. math::
   :label: eq-imp-rho-rel-err

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

:term:`Static shift` is a near-surface :term:`galvanic distortion` that
multiplies apparent resistivity without changing phase in the simplest
approximation. It is not primarily an impedance phase problem; it is an
amplitude scaling problem.

If :math:`g_x` and :math:`g_y` are real electric-field gains, a simplified
distortion model can be written:

.. math::
   :label: eq-imp-distortion-e

   \mathbf{E}_{obs}
   =
   \mathbf{C}
   \mathbf{E}_{true},

where :math:`\mathbf{C}` is a real :term:`distortion matrix`. Then:

.. math::
   :label: eq-imp-distortion-z

   \mathbf{Z}_{obs}
   =
   \mathbf{C}
   \mathbf{Z}_{true}.

Because :math:`\mathbf{C}` is real, left-multiplying by it rescales each
row of :math:`\mathbf{Z}` without touching phase -- exactly the "amplitude
only" behavior the name implies:

.. code-block:: pycon

   >>> C = np.diag([1.4, 0.6])  # real, positive static-shift gains
   >>> z_true = z1d[0]
   >>> z_obs = C @ z_true
   >>> np.round(np.angle(z_obs, deg=True), 6)
   array([[   0.,   45.],
          [-135.,    0.]])
   >>> np.round(np.angle(z_true, deg=True), 6)
   array([[   0.,   45.],
          [-135.,    0.]])
   >>> np.abs(z_obs) / np.abs(z_true)
   array([[nan, 1.4],
          [0.6, nan]])

Phase is identical before and after distortion (the ``nan`` entries above
are the zero diagonal, ``0/0``); only the off-diagonal amplitudes are
scaled, by :math:`g_x=1.4` and :math:`g_y=0.6` respectively. This is why
static shift can move apparent resistivity curves up or down without
producing an equivalent phase shift. The detailed interpretation is covered
in :doc:`static_shift`, but the tensor view is important: correcting static
shift changes amplitudes before apparent resistivity is interpreted or
inverted.

Determinant And Other Invariants
--------------------------------

Tensor invariants are scalar combinations of :math:`\mathbf{Z}` that can be
useful for plotting, comparison, or quality control. Common examples are:

.. math::
   :label: eq-imp-trace

   \operatorname{tr}(\mathbf{Z})
   =
   Z_{xx} + Z_{yy},

.. math::
   :label: eq-imp-det

   \det(\mathbf{Z})
   =
   Z_{xx} Z_{yy} - Z_{xy} Z_{yx},

and the Frobenius norm:

.. math::
   :label: eq-imp-frobenius

   \lVert \mathbf{Z} \rVert_F
   =
   \sqrt{
   |Z_{xx}|^2 + |Z_{xy}|^2 + |Z_{yx}|^2 + |Z_{yy}|^2
   }.

:eq:`eq-imp-trace` and :eq:`eq-imp-det` are exactly the two invariants
already confirmed under rotation in the worked example above -- "invariant"
here is not a qualitative claim, it is the ``np.allclose(..., True)`` result
from that block. pyCSAMT exposes such quantities through the
:class:`pycsamt.z.z.Z` container (``.trace``, ``.det``, ``.norm``). They are
helpful diagnostics, but they are not a replacement for understanding tensor
orientation, uncertainty, and dimensionality.

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

The safest documentation and reports should always state the component,
invariant, **and** resistivity convention used. "Apparent resistivity"
alone is incomplete; say :math:`\rho_{a,xy}`, :math:`\rho_{a,yx}`,
determinant apparent resistivity, or the exact product used -- and, per the
worked example above, whether it came from :func:`~pycsamt.seg.ops.z_to_rho_phi`'s SI
convention or :class:`~pycsamt.z.z.Z`'s legacy Zonge-style convention, since the two are
not interchangeable on the same raw array.

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
       conversion (legacy Zonge-style :eq:`eq-imp-rho-zonge` convention --
       this is what the EDI, Zonge, and Jones readers populate), and
       invariants (``.trace``, ``.det``, ``.norm``, ``.skew``).
   * - :class:`pycsamt.z.resphase.ResPhase`
     - Converts between complex impedance and apparent resistivity/phase
       using :eq:`eq-imp-rho-zonge`, including error propagation.
   * - :class:`pycsamt.z.tipper.Tipper`
     - Stores complex tipper components, rotations, amplitudes, phases, and
       arrow metrics.
   * - :func:`pycsamt.seg.ops.z_to_rho_phi`
     - Computes component-wise apparent resistivity and phase from
       impedance and frequency using the SI convention
       (:eq:`eq-imp-rho-si`) -- not currently wired into the EDI/Zonge/Jones
       ingestion pipeline, see the worked comparison above.
   * - :func:`pycsamt.seg.ops.rho_phi_to_z`
     - Reconstructs complex impedance from apparent resistivity, phase, and
       frequency (same SI convention as :func:`~pycsamt.seg.ops.z_to_rho_phi`).
   * - :func:`pycsamt.seg.ops.rotate_impedance`
     - Applies :math:`\mathbf{Z}' = \mathbf{R}\mathbf{Z}\mathbf{R}^T`;
       squeezes length-1 axes in its output (see the rotation example
       above).
   * - :func:`pycsamt.seg.ops.rotate_tipper`
     - Rotates tipper vectors consistently with the horizontal coordinate
       frame; same squeeze behavior for a single frequency.
   * - :func:`pycsamt.emtools.skew.bahr_skewness`
     - Computes the rotation-invariant Bahr skew (:eq:`eq-imp-bahr-eta`)
       used for dimensionality screening.

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
11. Which resistivity convention (SI or legacy Zonge-style) produced the
    numbers being compared?

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
* using CSAMT impedance-style products without checking source effects;
* mixing :func:`~pycsamt.seg.ops.z_to_rho_phi`'s SI convention with :class:`~pycsamt.z.z.Z`'s legacy
  Zonge-style convention on the same raw impedance values.

Connection To Later Workflows
-----------------------------

The impedance tensor feeds directly into:

* :doc:`static_shift`, because galvanic distortion affects amplitude;
* :doc:`inversion_concepts`, because inversion fits impedance,
  resistivity/phase, or transformed tensor products;
* :doc:`../user_guide/pipeline/steps`, because tensor cleanup, skew gates, rotations,
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
