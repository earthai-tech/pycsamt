.. _dimensionality:

Dimensionality, Distortion, And The Phase Tensor
================================================

:doc:`impedance_tensor` introduces the phase tensor as a real, distortion-
resistant object built from :math:`\mathbf{Z}`, and :doc:`csamt_amt_mt_overview`
and :doc:`inversion_concepts` both lean on :term:`Dimensionality` as a
modelling choice without saying exactly how it is measured. This page closes
that gap: it derives the phase tensor's orientation and skew angles, uses
them to classify 1-D/2-D/3-D behavior station by station, and covers the
:term:`Groom-Bailey decomposition`, the classical model for separating local
galvanic :term:`distortion matrix` effects from the :term:`regional tensor`
that :doc:`static_shift` only touches briefly through its own tensor view.

Every number below comes from the same 28-station bundled AMT survey used in
:doc:`static_shift`, ``data/AMT/WILLY_DATA/L18PLT``
(:func:`pycsamt.api.read_edis`), processed with
:mod:`pycsamt.emtools.tensor`, :mod:`pycsamt.emtools.dimensionality`, and
:mod:`pycsamt.emtools.gb`.

Alpha, Beta, And Ellipticity
----------------------------

Recall from :doc:`impedance_tensor` that the phase tensor is
:math:`\boldsymbol{\Phi}=\mathbf{X}^{-1}\mathbf{Y}` (:eq:`eq-imp-phasetensor`),
a real 2x2 matrix even though :math:`\mathbf{Z}=\mathbf{X}+i\mathbf{Y}` is
complex. Write its entries as :math:`\boldsymbol{\Phi}=\begin{bmatrix}a & b
\\ c & d\end{bmatrix}`. Caldwell, Bibby & Brown (2004) [Caldwell2004]_ define
two angles from these entries. The orientation angle:

.. math::
   :label: eq-dim-alpha

   \alpha = \tfrac{1}{2}\operatorname{atan2}(b + c,\, a - d),

and the skew angle:

.. math::
   :label: eq-dim-beta

   \beta = \tfrac{1}{2}\operatorname{atan2}(b - c,\, a + d).

The two look almost interchangeable on the page, but they behave completely
differently under coordinate rotation. :math:`\alpha` rotates together with
the survey's coordinate frame -- it is a *direction*, like :term:`strike` --
while :math:`\beta` (the phase-tensor :term:`Skew`) is built from the
combination of entries that stays fixed no matter which way the horizontal
axes are turned, which is exactly what makes it useful as a distortion and
dimensionality diagnostic instead of just another orientation angle.

pyCSAMT computes both in :func:`pycsamt.emtools.tensor.build_phase_tensor_table`.
While researching this page, rotating an arbitrary synthetic phase tensor
through the implementation caught a real bug: a previous version of the
underlying ``_angles_deg`` helper had swapped the two formulas above (and
additionally sign-flipped the skew one), so the dataframe's ``"alpha"``
column was actually rotation-invariant and its ``"beta"``/``"skew"`` column
tracked coordinate rotation -- exactly backwards from what a skew diagnostic
needs. The check that exposed it is simple enough to reproduce here: rotate
:math:`\boldsymbol{\Phi}=\begin{bmatrix}0.6 & 0.35\\0.10 & 0.45\end{bmatrix}`
by :math:`0^\circ`, :math:`20^\circ`, and :math:`45^\circ` and recompute
:eq:`eq-dim-alpha`/:eq:`eq-dim-beta` at each angle:

.. code-block:: pycon

   >>> import numpy as np
   >>> from pycsamt.emtools.tensor import _angles_deg
   >>> def rotate_phi(Phi, deg):
   ...     th = np.radians(deg)
   ...     R = np.array([[np.cos(th), np.sin(th)], [-np.sin(th), np.cos(th)]])
   ...     return R @ Phi @ R.T
   >>> Phi = np.array([[0.6, 0.35], [0.10, 0.45]])
   >>> for ang in (0.0, 20.0, 45.0):
   ...     Pr = rotate_phi(Phi, ang)
   ...     a, b = _angles_deg(
   ...         np.array([Pr[0, 0]]), np.array([Pr[0, 1]]),
   ...         np.array([Pr[1, 0]]), np.array([Pr[1, 1]]),
   ...     )
   ...     print(ang, round(float(a[0]), 3), round(float(b[0]), 3))
   0.0 35.783 6.696
   20.0 15.783 6.696
   45.0 -9.217 6.696

:math:`\alpha` tracks the applied rotation exactly (it drops by precisely the
rotation angle each time: :math:`35.783-20=15.783`,
:math:`15.783-25=-9.217`), while :math:`\beta=6.696^\circ` never moves. That
asymmetry -- one angle rotation-variant, the other not -- is the whole test;
a function that gets the two formulas backwards passes every ordinary unit
test that only checks numeric ranges, and only fails a check that rotates
the input and compares. The fix (already released) restores
:eq:`eq-dim-alpha`/:eq:`eq-dim-beta` to the standard assignment, and the one
downstream pyCSAMT page that quoted the old, swapped skew values --
:doc:`../tutorials/map_porphyry_mineralization_from_noisy_amt` -- has been
regenerated with the corrected numbers. :func:`~pycsamt.emtools.strike.estimate_strike_consensus`
and everything built on it were unaffected, because strike estimation uses
the tensor's SVD-based orientation (below), not ``alpha``/``beta``.

The tensor's shape, independent of either angle, comes from its singular
values. If :math:`\boldsymbol{\Phi}=\mathbf{U}\mathbf{S}\mathbf{V}^{T}` with
:math:`\mathbf{S}=\operatorname{diag}(\phi_{\max}, \phi_{\min})`, the
:term:`Ellipticity` is:

.. math::
   :label: eq-dim-ellipt

   \lambda = \frac{\phi_{\max} - \phi_{\min}}{\phi_{\max} + \phi_{\min}}.

A circular phase tensor (:math:`\phi_{\max}=\phi_{\min}`, :math:`\lambda=0`)
is the 1-D signature: the phase relationship between :math:`\mathbf{E}` and
:math:`\mathbf{H}` is the same in every direction. Growing :math:`\lambda`
means the phase response depends on direction, which is what 2-D and 3-D
structure produce. :func:`pycsamt.emtools.tensor.build_phase_tensor_table`
returns ``alpha``, ``beta`` (aliased as ``skew``), ``theta`` (the SVD major-axis
orientation used by strike estimation), and ``ellipt`` together, one row per
station-period:

.. code-block:: pycon

   >>> from pycsamt.api import read_edis
   >>> from pycsamt.emtools.tensor import build_phase_tensor_table
   >>> survey = read_edis(
   ...     "data/AMT/WILLY_DATA/L18PLT", recursive=False, strict=False,
   ...     progress=False,
   ... )
   >>> sites = survey.collection
   >>> len(sites)
   28
   >>> pt = build_phase_tensor_table(sites)
   >>> len(pt)
   1484
   >>> pt.head(3)[["freq", "alpha", "beta", "theta", "ellipt"]].round(3)
         freq   alpha   beta    theta  ellipt
   0  10400.0 -56.701  2.612  120.688   0.195
   1   8707.0 -54.693  1.964  123.342   0.210
   2   7289.0 -51.452  1.804  126.744   0.214
   >>> round(pt["beta"].abs().median(), 2), round(pt["beta"].abs().quantile(0.9), 2)
   (14.92, 62.28)
   >>> round(pt["ellipt"].abs().median(), 3)
   0.674

A median :math:`|\beta|` of :math:`14.9^\circ` is well above the few-degree
level usually taken as "clean" 2-D, and the 90th percentile of
:math:`62.3^\circ` shows a long high-skew tail -- consistent with the
station-period counts below, where 3-D is the majority label across this
line, not the exception.

Classifying Dimensionality
--------------------------

:func:`pycsamt.emtools.dimensionality.classify_dimensionality` turns
:math:`|\beta|` and :math:`|\lambda|` into a three-way label using two
thresholds, ``skew_th`` (default :math:`3.0^\circ`) and ``ellipt_th``
(default :math:`0.2`):

.. list-table::
   :header-rows: 1
   :widths: 12 30 58

   * - Label
     - Condition
     - Meaning
   * - 0 (1-D)
     - :math:`|\beta|\le\text{skew\_th}` and :math:`|\lambda|\le\text{ellipt\_th}`
     - Low skew, low ellipticity -- isotropic phase response.
   * - 1 (2-D)
     - :math:`|\beta|\le\text{skew\_th}` and :math:`|\lambda|>\text{ellipt\_th}`
     - Low skew, but the phase tensor is elongated -- direction-dependent
       phase, no strong 3-D skew.
   * - 2 (3-D)
     - :math:`|\beta|>\text{skew\_th}`
     - High skew regardless of ellipticity -- the diagnostic that most
       directly signals 3-D or unresolved distortion.

Running the classifier on the same 28 stations:

.. code-block:: pycon

   >>> from pycsamt.emtools.dimensionality import classify_dimensionality
   >>> dimdf = classify_dimensionality(sites)
   >>> len(dimdf)
   1484
   >>> dimdf["dim"].value_counts().to_dict()
   {2: 1271, 1: 155, 0: 58}

Out of 1484 station-period rows, 1271 (85.6%) are labelled 3-D, 155 (10.4%)
2-D, and only 58 (3.9%) 1-D. This is a real, if somewhat sobering, result for
a line often processed with a 2-D workflow: the *majority* of periods do not
pass the default 2-D skew test. It does not mean a 2-D inversion of this line
is meaningless -- 2-D inversion routinely tolerates some 3-D contamination,
and the low-skew 58+155 rows show that genuinely clean bands do exist -- but
it does mean the default thresholds should be read as a screening tool, not
a verdict, and that residual 3-D behavior should be expected in the
inversion misfit rather than treated as a surprise.

The figure below combines the phase-tensor ellipse pseudo-section (top,
:func:`pycsamt.emtools.tensor.plot_phase_tensor_psection`, fill colour =
:math:`\beta`) with the resulting dimensionality pseudo-section (bottom,
:func:`pycsamt.emtools.tensor.plot_dimensionality_psection`) for the same 28
stations and the same period range:

.. code-block:: python
   :linenos:

   import matplotlib.pyplot as plt
   from pycsamt.emtools.tensor import (
       plot_phase_tensor_psection,
       plot_dimensionality_psection,
   )

   fig, axes = plt.subplots(2, 1, figsize=(10.5, 8.5))
   plot_phase_tensor_psection(
       sites, ax=axes[0], title="Phase-tensor ellipses (fill = skew beta)"
   )
   plot_dimensionality_psection(sites, ax=axes[1])
   axes[1].set_title("Rule-based dimensionality (0=1D, 1=2D, 2=3D)")
   fig.tight_layout()

.. figure:: /images/theory/dim_phase_tensor_psection.png
   :alt: Phase-tensor ellipse pseudo-section and rule-based dimensionality pseudo-section for the L18PLT AMT line
   :width: 100%

   Top: ellipse shape encodes :math:`\phi_{\max}/\phi_{\min}` (ellipticity)
   and orientation encodes :math:`\theta`; fill colour is skew :math:`\beta`.
   Bottom: the same data reduced to the three-way dimensionality label. The
   yellow (3-D) majority is concentrated station-wide, with pockets of
   teal (2-D) and a handful of purple (1-D) cells mostly at the shortest
   periods -- exactly the pattern a near-surface, laterally heterogeneous
   overburden would produce on top of more consistent deeper structure.

Groom-Bailey Galvanic Distortion
--------------------------------

:doc:`impedance_tensor` already writes the static-shift distortion model as
:math:`\mathbf{Z}_{obs}=\mathbf{C}\mathbf{Z}_{true}` (:eq:`eq-imp-distortion-z`)
for a diagonal gain-only :math:`\mathbf{C}`. Groom & Bailey (1989)
[GroomBailey1989]_ generalize this to a full real 2x2
:term:`distortion matrix` :math:`\mathbf{D}` acting on a :term:`regional
tensor` :math:`\mathbf{Z}_{2D}` that is assumed anti-diagonal after rotation
to strike:

.. math::
   :label: eq-dim-gb-model

   \mathbf{Z}_{obs}(f) \approx \mathbf{D}\,\mathbf{Z}_{2D}(f),
   \qquad
   \mathbf{Z}_{2D}(f) =
   \begin{bmatrix} 0 & Z_{\mathrm{TE}}(f) \\ Z_{\mathrm{TM}}(f) & 0 \end{bmatrix}.

Because :math:`\mathbf{D}` is real and frequency-independent while
:math:`\mathbf{Z}_{2D}` carries all the frequency dependence, fitting
:eq:`eq-dim-gb-model` across a band of frequencies over-determines
:math:`\mathbf{D}` from a single station -- that is the leverage that makes
the decomposition possible at all.
:func:`pycsamt.emtools.gb._fit_gb_distortion` solves it iteratively: given a
trial :math:`\mathbf{D}`, solve for the anti-diagonal
:math:`\mathbf{Z}_{2D}` that best explains the observed :math:`\mathbf{Z}`,
then refit :math:`\mathbf{D}` given that :math:`\mathbf{Z}_{2D}`, alternating
until the residual stabilizes.

The fitted :math:`\mathbf{D}` is decomposed into four named parameters. Gain
is its scale:

.. math::
   :label: eq-dim-gb-gain

   \mathrm{gain} = \sqrt{|\det \mathbf{D}|},

twist is the antisymmetric (rotational) part of the gain-normalized
:math:`\mathbf{D}_n = \mathbf{D}/\mathrm{gain}`:

.. math::
   :label: eq-dim-gb-twist

   \mathrm{twist} =
   \operatorname{atan2}\!\big(D_{n,12}-D_{n,21},\, D_{n,11}+D_{n,22}\big),

and shear and anisotropy come from rotating :math:`\mathbf{D}_n` back by
:math:`-\mathrm{twist}` into :math:`\mathbf{M}=\mathbf{R}(-\mathrm{twist})
\mathbf{D}_n` and reading its own symmetric/antisymmetric split:

.. math::
   :label: eq-dim-gb-shear

   \mathrm{shear} = \frac{M_{12}+M_{21}}{M_{11}+M_{22}},
   \qquad
   \mathrm{anisotropy} = \frac{M_{11}-M_{22}}{M_{11}+M_{22}}.

There is a real, easy-to-miss consequence of :eq:`eq-dim-gb-gain` combined
with how :func:`~pycsamt.emtools.gb._fit_gb_distortion` normalizes
:math:`\mathbf{D}` at every iteration (to unit determinant, so the fit
doesn't drift to an arbitrary scale): ``gain`` in
:func:`pycsamt.emtools.gb.groom_bailey_table` is not a free parameter at all,
it is pinned to :math:`1.0` by construction, for every station, every time:

.. code-block:: pycon

   >>> from pycsamt.emtools.gb import groom_bailey_table
   >>> gb = groom_bailey_table(sites)
   >>> len(gb), (gb["status"] == "ok").sum()
   (28, 28)
   >>> round(float(gb["gain"].min()), 6), round(float(gb["gain"].max()), 6)
   (1.0, 1.0)

This is not a bug -- it is the classical Groom-Bailey scale ambiguity: an
overall multiplicative gain is indistinguishable from :term:`static shift`
in this decomposition, so the fit deliberately normalizes it away and
reports only the shape parameters (twist, shear, anisotropy) that static
shift alone cannot produce. Any real amplitude scaling still present in the
data has to be handled separately, by the static-shift tools in
:doc:`static_shift`, not read off this table's ``gain`` column.

The other three parameters do vary station to station and carry real
information:

.. code-block:: pycon

   >>> gb = gb.sort_values("station")
   >>> round(float(gb["twist_deg"].abs().median()), 2)
   11.72
   >>> round(float(gb["shear"].abs().median()), 3)
   0.236
   >>> row0 = gb[gb["station"] == sites[0].station].iloc[0]
   >>> [round(float(row0[k]), 4) for k in
   ...  ("twist_deg", "shear", "anisotropy", "rms_fit")]
   [9.6171, -0.5186, 0.0756, 0.1607]
   >>> round(float(gb["diagonal_ratio_before"].median()), 3), round(float(gb["diagonal_ratio_after"].median()), 3)
   (0.423, 0.367)

Station ``18-001A`` has a modest :math:`9.6^\circ` twist and a shear of
:math:`-0.52`; across the whole line, median :math:`|{\rm twist}|` is
:math:`11.7^\circ` and median :math:`|{\rm shear}|` is :math:`0.24`, both
comfortably nonzero -- this line is not distortion-free. The median
diagonal/off-diagonal amplitude ratio drops from :math:`0.423` before
correction to :math:`0.367` after, a real but partial improvement, not a
clean collapse to zero.

"Partial" undersells how uneven that improvement is station by station.
Applying the fitted correction with
:func:`pycsamt.emtools.gb.groom_bailey_decomposition` and comparing
before/after diagonal ratios one station at a time shows the correction
working best exactly where distortion is mild, and making the diagonal
ratio *worse* at a cluster of high-twist, high-shear stations:

.. code-block:: python
   :linenos:

   import matplotlib.pyplot as plt
   import numpy as np

   x = np.arange(len(gb))
   st = gb["station"].to_list()
   fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))

   axes[0].bar(x, gb["twist_deg"], color="#1f77b4")
   axes[0].set_title("Twist (degrees)")

   axes[1].bar(x, gb["shear"], color="#d62728")
   axes[1].set_title("Shear (dimensionless)")

   axes[2].plot(x, gb["diagonal_ratio_before"], "o-", color="0.5", label="before GB")
   axes[2].plot(x, gb["diagonal_ratio_after"], "o-", color="#2ca02c", label="after GB")
   axes[2].set_title("Diagonal / off-diagonal ratio")
   axes[2].legend(fontsize=8)

   for ax in axes:
       ax.set_xticks(x[::3])
       ax.set_xticklabels([st[i] for i in x[::3]], rotation=45, ha="right", fontsize=7)
       ax.grid(True, alpha=0.3)
   fig.tight_layout()

.. figure:: /images/theory/dim_groom_bailey.png
   :alt: Groom-Bailey twist, shear, and before/after diagonal-ratio bar and line plots for the L18PLT AMT line
   :width: 100%

   Twist and shear both grow sharply for the same run of stations
   (``18-020A`` through ``18-024U``, :math:`|\mathrm{twist}|` up to
   :math:`58^\circ`, :math:`|\mathrm{shear}|` up to :math:`0.92`). At those
   same stations the green "after" curve sits
   *above* the grey "before" curve -- correcting for a fitted distortion
   matrix made the residual diagonal terms larger, not smaller. 13 of the 28
   stations end up worse by this measure, and it is not simply the stations
   with the worst fit residual (``rms_fit``): the correlation between
   ``rms_fit`` and the before/after change is only 0.17.

The honest reading is that :eq:`eq-dim-gb-model`'s core assumption -- a
single frequency-independent real :math:`\mathbf{D}` explains the data on
top of a purely anti-diagonal regional tensor -- holds unevenly across this
line. Where it holds, the correction genuinely suppresses diagonal energy.
Where the underlying structure is not well described by *any* real
distortion matrix over an anti-diagonal regional response (large, unstable
twist/shear are the warning sign, not the fit residual alone), forcing the
correction can inject as much diagonal energy as it removes. This is a
standard, documented limitation of classical galvanic-distortion
decomposition, not a pyCSAMT-specific defect -- but it means station-by-
station before/after comparison, not just the survey median, should guide
whether to trust and apply a Groom-Bailey correction.

.. code-block:: pycon

   >>> from pycsamt.emtools.gb import groom_bailey_decomposition
   >>> res = groom_bailey_decomposition(sites, apply=True)
   >>> res.n_station, res.applied
   (28, True)

:class:`~pycsamt.emtools.gb.GroomBaileyResult` returns the corrected sites
alongside the fitted parameter table, ready for the same downstream steps
(static-shift review, strike estimation, dimensionality reclassification,
inversion) used elsewhere in this documentation.

Connection To Strike And Later Workflows
----------------------------------------

The SVD orientation angle :math:`\theta` returned alongside ``alpha``/
``beta`` in :func:`~pycsamt.emtools.tensor.build_phase_tensor_table` is one
of the inputs :func:`pycsamt.emtools.strike.estimate_strike_consensus`
combines with other strike estimators; this page does not repeat that
material, which belongs to strike estimation specifically rather than
dimensionality classification. What this page adds is upstream of strike: a
reason to trust (or distrust) a 2-D/strike-based workflow in the first
place, and a distortion correction that can be applied before strike
rotation and inversion.

Dimensionality and Groom-Bailey diagnostics feed directly into:

* :doc:`inversion_concepts`, because 1-D/2-D/3-D dimensionality is a
  modelling choice made partly from the counts and figure above;
* :doc:`static_shift`, because Groom-Bailey's gain-normalization means
  amplitude (static-shift) correction has to be handled separately, not
  read off the same table;
* :doc:`csamt_amt_mt_overview`, which introduces phase tensor and skew as
  quick-look diagnostics without deriving them -- this page is where that
  derivation lives.

Practical Guidance And Pitfalls
-------------------------------

* Do not read ``alpha`` as a skew or ``beta`` as an orientation -- confirm
  which is which for any pyCSAMT version by rotating a synthetic tensor, the
  same check that caught the swap described above.
* Treat ``classify_dimensionality``'s default thresholds
  (:math:`\text{skew\_th}=3^\circ`, :math:`\text{ellipt\_th}=0.2`) as a
  screening rule, not a verdict -- a majority-3-D result over a real survey
  is common and does not automatically disqualify a 2-D inversion.
* Do not average phase-tensor quantities across periods before checking
  whether the dimensionality label is stable with period; shallow and deep
  structure can classify differently at the same station.
* Groom-Bailey's ``gain`` column is always 1.0 by construction; look to
  static-shift tools, not this table, for amplitude correction.
* Always inspect the before/after diagonal ratio per station, not only the
  survey median, before trusting a Groom-Bailey correction -- it can make
  individual stations worse even when the median improves.
* Large, unstable twist or shear is itself diagnostic: it suggests the
  anti-diagonal regional-tensor assumption is breaking down, which is
  useful information even when the fit residual (``rms_fit``) looks
  unremarkable.

References
----------

This page follows Caldwell, Bibby & Brown's phase-tensor formulation
[Caldwell2004]_ and Groom & Bailey's galvanic-distortion decomposition
[GroomBailey1989]_. See also :doc:`../references`.
