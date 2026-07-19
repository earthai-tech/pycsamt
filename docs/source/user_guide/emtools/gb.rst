.. _emtools_gb:

Groom-Bailey Galvanic Distortion
================================

``pycsamt.emtools.gb`` estimates and optionally removes
frequency-independent :term:`galvanic distortion` from MT/AMT/CSAMT
:term:`impedance tensor` data. It is designed as an auditable
preprocessing step before 2-D interpretation and inversion.

Galvanic distortion is a near-surface effect. Small local conductivity
heterogeneities deflect and scale the electric field before the receiver
records it, so the measured tensor can contain diagonal leakage and mode
mixing that are not part of the deeper regional induction response. The
:term:`Groom-Bailey decomposition` is one common way to describe that
effect when the regional response is close enough to 2-D.

The fitted model is:

.. math::

   Z_{obs}(f) \approx D Z_{2D}(f)

where:

- ``Z_obs`` is the observed :term:`impedance tensor`;
- ``D`` is a real, frequency-independent 2 x 2
  :term:`distortion matrix`;
- ``Z_2D`` is the best anti-diagonal :term:`regional tensor` at each
  frequency.

Full callable signatures live in the :doc:`API reference <../../api/emtools>`.
This page explains how to fit the table, read the distortion parameters,
apply the correction, and record the result in a pre-2D workflow.

When To Use Groom-Bailey
------------------------

Use this workflow when the data appear close enough to 2-D for
:term:`galvanic distortion` correction to be meaningful, but the
:term:`impedance tensor` has diagonal leakage or station-dependent
distortion that should be documented before inversion.

Good use cases include:

- preparing a 2-D inversion input after dimensionality and strike checks;
- testing whether diagonal tensor leakage is reduced after correction;
- documenting :term:`twist`, :term:`shear`, and anisotropy-style
  distortion parameters;
- comparing corrected and uncorrected impedance curves at the same
  station.

Poor use cases include:

- strongly 3-D data with no stable strike or 2-D period band;
- too few valid frequencies in the selected band;
- using the fitted gain as a unique static-shift solution;
- applying correction without saving the fit diagnostics.

Core Assumptions
----------------

The implementation fits a real distortion matrix that is constant over
the selected period band. That is the galvanic assumption: the
distortion is local and frequency-independent, while the regional tensor
varies with frequency.

This assumption is powerful but narrow. It says the shallow distortion
changes the measured electric-field axes, while the deeper
:term:`regional tensor` still carries the frequency-dependent induction
physics. If the data are strongly 3-D across the whole band, the fitted
matrix may become only a mathematical approximation, not a meaningful
correction.

The regional tensor is forced to be anti-diagonal:

.. math::

   Z_{2D}(f) =
   \begin{bmatrix}
   0 & u(f) \\
   v(f) & 0
   \end{bmatrix}

The observed tensor is then approximated by multiplying this 2-D tensor
by ``D``. Because the model :math:`Z_{obs} = D\,Z_{2D}` is bilinear —
linear in ``D`` for fixed :math:`Z_{2D}`, and linear in :math:`u, v` for
fixed ``D`` — each half of the fit has a closed-form least-squares
solution, and the iteration just alternates between them. With ``D``
held fixed, the best anti-diagonal tensor at each frequency is the
projection

.. math::

   u(f) = \frac{D_{xx} Z_{xy}(f) + D_{yx} Z_{yy}(f)}{D_{xx}^2 + D_{yx}^2},
   \qquad
   v(f) = \frac{D_{xy} Z_{xx}(f) + D_{yy} Z_{yx}(f)}{D_{xy}^2 + D_{yy}^2},

and with :math:`u, v` held fixed across all frequencies, each row of
``D`` is refit by ordinary least squares against the corresponding
tensor components. The loop repeats until the relative RMS residual
(defined below) stops improving by more than ``tol``, or ``max_iter``
is reached.

After each row-solve, the fitted matrix is rescaled so that
:math:`|\det D| = 1` — with :math:`u, v` rescaled inversely so the
product :math:`D\,Z_{2D}` is unchanged — before it is summarized as
gain, twist, shear, and anisotropy-style parameters below.

Fit A Distortion Table
----------------------

Start by estimating parameters without changing the data.

.. code-block:: python
   :linenos:

   from pycsamt.emtools.gb import groom_bailey_table

   survey = "data/AMT/WILLY_DATA/L18PLT"

   table = groom_bailey_table(
       survey,
       band=(1e-3, 10.0),
       rotate_deg=None,
       min_freq=4,
       max_iter=30,
       tol=1e-6,
       robust=True,
   )

   print(
       table[
           [
               "station",
               "status",
               "n_freq",
               "twist_deg",
               "shear",
               "anisotropy",
               "rms_fit",
               "diagonal_ratio_before",
               "diagonal_ratio_after",
           ]
       ].head()
   )

.. code-block:: text

      station status  ...  diagonal_ratio_before  diagonal_ratio_after
   0  18-001A     ok  ...               0.614153              0.273213
   1  18-002U     ok  ...               0.434659              0.261954
   2  18-003A     ok  ...               0.373588              0.379570
   3  18-004A     ok  ...               0.454007              0.323360
   4  18-005U     ok  ...               0.496585              0.305364

   [5 rows x 9 columns]

The ``band`` argument is in period seconds, not hertz. Choose a band
that is justified by dimensionality, strike stability, and data quality.

Table Columns
-------------

Successful rows have ``status == "ok"`` and include:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Column
     - Meaning
   * - ``station``
     - Station name.
   * - ``n_freq``
     - Number of valid frequencies used in the fit.
   * - ``period_min_s`` / ``period_max_s``
     - Period range actually used after band selection.
   * - ``rotate_deg``
     - Rotation angle applied before fitting, or ``NaN``.
   * - ``distortion_xx`` ... ``distortion_yy``
     - Entries of the fitted real 2 x 2 :term:`distortion matrix`.
   * - ``gain``
     - :math:`\sqrt{|\det D|}` of the matrix handed to the twist/shear
       decomposition — see the note below.
   * - ``twist_deg``
     - :term:`Twist` angle inferred from the normalized matrix.
   * - ``shear``
     - Dimensionless :term:`shear`-style parameter, clipped to
       ``[-0.99, 0.99]``.
   * - ``shear_angle_deg``
     - ``atan(shear)`` in degrees.
   * - ``anisotropy``
     - Dimensionless :term:`anisotropy`-style parameter, clipped to
       ``[-0.99, 0.99]``.
   * - ``rms_fit``
     - Relative fit residual.
   * - ``diagonal_ratio_before``
     - Median diagonal/off-diagonal tensor ratio before correction.
   * - ``diagonal_ratio_after``
     - Median diagonal/off-diagonal tensor ratio after applying the
       fitted inverse matrix to the fitted band.
   * - ``robust``
     - Whether robust residual weighting was used.
   * - ``method``
     - Current method label, ``gb_real_distortion_2d``.

Rows with too few valid frequencies have
``status == "insufficient_frequencies"`` and include the available
``n_freq``. Increase the band, lower ``min_freq`` only with care, or
exclude that station from correction.

In the terminology used by this table, :term:`twist` is rotational
mixing, :term:`shear` is non-orthogonal mixing, and the reported
:term:`anisotropy` parameter is the directional scaling left after the
twist part has been removed. These are distortion parameters, not a
complete geological interpretation by themselves.

The twist/shear/anisotropy decomposition takes whatever matrix
``distortion_xx`` ... ``distortion_yy`` reports — call it :math:`D` —
and unwinds it into a rotation, a shear, and an anisotropy, the same
way a 2x2 real matrix decomposes in general:

.. math::

   \mathrm{gain} = \sqrt{|\det D|}, \qquad D_n = D / \mathrm{gain},

.. math::

   \mathrm{twist\_deg} = \operatorname{atan2}\!\bigl(
   D_{n,xy} - D_{n,yx},\ D_{n,xx} + D_{n,yy}\bigr) \times
   \frac{180}{\pi},

.. math::

   M = R(-\mathrm{twist\_deg})\, D_n, \qquad
   \mathrm{shear} = \frac{M_{xy}+M_{yx}}{M_{xx}+M_{yy}}, \qquad
   \mathrm{anisotropy} = \frac{M_{xx}-M_{yy}}{M_{xx}+M_{yy}},

with ``shear`` and ``anisotropy`` clipped to :math:`[-0.99, 0.99]`
against the case :math:`M_{xx}+M_{yy}\approx 0`. The residual and
diagonal-ratio columns are:

.. math::

   \mathrm{rms\_fit} =
   \frac{\sqrt{\left\langle |Z_{obs}-D\,Z_{2D}|^2\right\rangle}}
        {\sqrt{\left\langle |Z_{obs}|^2\right\rangle}},
   \qquad
   \mathrm{diagonal\_ratio} = \mathrm{median}\!\left(
   \frac{\sqrt{|Z_{xx}|^2+|Z_{yy}|^2}}{\sqrt{|Z_{xy}|^2+|Z_{yx}|^2}}
   \right),

where :math:`\langle \cdot \rangle` averages over every frequency and
tensor component in the fitted band. ``diagonal_ratio`` is computed the
same way before correction (on the raw tensor) and after (on the fitted
inverse matrix applied to the fitted band) — it is *not* the same
quantity as the ``diagonal`` confidence score in :ref:`emtools_qc`,
which reports a leakage *fraction* rather than a raw ratio.

The twist formula is structurally the same
:math:`\arctan2`-of-off-diagonal-over-diagonal shape as the phase-tensor
skew :math:`\beta` in :ref:`emtools_skew` — both come from decomposing a
real 2x2 matrix into a rotation plus a symmetric remainder — but the two
quantities are not interchangeable: one describes distortion in ``D``,
the other describes the regional tensor itself.

One honest caveat about ``gain``: the iterative fit already rescales its
working matrix to :math:`|\det D| = 1` after every row-solve (see
*Core Assumptions* above), and that rescaled matrix is what gets stored
in ``distortion_xx`` ... ``distortion_yy``. Recomputing
:math:`\sqrt{|\det D|}` from that already-normalized matrix returns
``1.0`` for every successful fit — which is exactly what you will see if
you print the column. Read ``gain`` as confirmation that the stored
matrix is in unit-determinant form, not as a per-station distortion
amplitude; the classic gain/static-shift ambiguity this page already
warns about is not something this column resolves.

Reading The Parameters
----------------------

The most useful diagnostic columns are usually:

- ``rms_fit``: lower values indicate that the fitted model describes the
  selected band better.
- ``diagonal_ratio_before`` and ``diagonal_ratio_after``: correction is
  behaving sensibly when the after value is lower.
- ``twist_deg``: large :term:`twist` can imply strong
  :term:`galvanic distortion` or a poor 2-D assumption.
- ``shear`` and ``anisotropy``: large absolute values deserve station
  inspection.
- ``n_freq``: low values make the fit less stable.

``gain`` is not one of them — as explained above, it is always ``1.0``
by construction, not a per-station static-shift estimate.

Rank Stations For Review
------------------------

Use the table to find stations with poor fits or strong residual
diagonal leakage.

.. code-block:: python
   :linenos:

   from pycsamt.emtools.gb import groom_bailey_table

   table = groom_bailey_table(
       "data/AMT/WILLY_DATA/L18PLT",
       band=(1e-3, 10.0),
       robust=True,
   )

   ok = table.loc[table["status"] == "ok"].copy()
   ok["diag_reduction"] = (
       ok["diagonal_ratio_before"] - ok["diagonal_ratio_after"]
   )

   ranked = ok.sort_values(
       ["rms_fit", "diagonal_ratio_after"],
       ascending=[False, False],
   )

   print(
       ranked[
           [
               "station",
               "n_freq",
               "rms_fit",
               "diag_reduction",
               "twist_deg",
               "shear",
               "anisotropy",
           ]
       ].head(10)
   )

.. code-block:: text

       station  n_freq   rms_fit  diag_reduction  twist_deg     shear  anisotropy
   22  18-022U      39  0.552447        0.322411  22.141866 -0.117145   -0.008538
   21  18-021U      39  0.537890       -0.419458 -63.216543  0.990000    0.404853
   20  18-021B      39  0.430286       -0.604539 -37.911526  0.866936    0.211520
   17  18-018A      39  0.427220       -0.564491  56.099182  0.990000   -0.234017
   24  18-023A      39  0.412376        0.298014  19.376516 -0.123178    0.012153
   27  18-025A      39  0.390159       -0.384625 -62.091367  0.504772    0.012596
   18  18-019U      39  0.363720       -0.016977  17.861709 -0.015653    0.008229
   8   18-009A      39  0.335742       -0.014381  17.904448  0.118619   -0.006990
   19  18-020A      39  0.326943       -0.391908 -59.746000  0.870488    0.525820
   12  18-013U      39  0.322409        0.166075   8.996415 -0.388653    0.079927

Stations with high ``rms_fit`` or little diagonal reduction should be
reviewed before applying correction automatically.

Use A Strike Rotation
---------------------

If you have selected a strike angle, pass it as ``rotate_deg`` before
fitting.

.. code-block:: python
   :linenos:

   from pycsamt.emtools.gb import groom_bailey_table

   strike_deg = 35.0

   table = groom_bailey_table(
       "data/AMT/WILLY_DATA/L18PLT",
       band=(1e-3, 10.0),
       rotate_deg=strike_deg,
       robust=True,
   )

The rotation is applied to the tensor before fitting the distortion
matrix. Use a strike that has been justified by the strike and
dimensionality workflows, not one chosen to improve the GB fit alone.

Apply A Precomputed Table
-------------------------

Use ``apply_groom_bailey`` when you have already inspected and accepted
a table.

.. code-block:: python
   :linenos:

   from pycsamt.emtools.gb import apply_groom_bailey, groom_bailey_table

   survey = "data/AMT/WILLY_DATA/L18PLT"

   table = groom_bailey_table(
       survey,
       band=(1e-3, 10.0),
       robust=True,
   )

   accepted = table.loc[
       (table["status"] == "ok")
       & (table["rms_fit"] < 0.25)
       & (table["diagonal_ratio_after"] < table["diagonal_ratio_before"])
   ].copy()

   corrected = apply_groom_bailey(
       survey,
       table=accepted,
       inplace=False,
   )

Only stations present in the accepted table are corrected. Stations
missing from the table or with invalid matrices are left unchanged.

Estimate And Apply In One Step
------------------------------

Use ``groom_bailey_decomposition`` when you want a result container with
the fitted table and optionally corrected sites.

.. code-block:: python
   :linenos:

   from pycsamt.emtools.gb import groom_bailey_decomposition

   result = groom_bailey_decomposition(
       "data/AMT/WILLY_DATA/L18PLT",
       apply=True,
       band=(1e-3, 10.0),
       rotate_deg=None,
       robust=True,
       inplace=False,
   )

   print(result.summary())
   corrected_sites = result.sites
   gb_table = result.table

.. code-block:: text

   GroomBaileyResult(stations=28, applied=True, median_rms=0.2797)

The result container records:

- ``sites``: corrected sites when ``apply=True``, otherwise loaded sites.
- ``table``: fitted parameter table.
- ``applied``: whether correction was applied.
- ``method``: method label.
- ``n_station``: number of fitted table rows.

Compare Robust And Non-Robust Fits
----------------------------------

With ``robust=True``, every iteration after the first re-weights
frequencies by a Huber-style rule built from the per-frequency residual
norm :math:`r_i = \sqrt{\langle |Z_{obs,i}-D\,Z_{2D,i}|^2\rangle}`
(averaged over the four tensor components at frequency :math:`i`):

.. math::

   c = 1.4826\,\mathrm{med}\bigl(|r_i - \mathrm{med}(r)|\bigr)
       \times 1.345,
   \qquad
   w_i = \mathrm{clip}\!\left(\min\!\left(1,\ \frac{c}{r_i}\right),\
   0.05,\ 1\right),

so a frequency at or below the robust scale :math:`c` keeps full
weight, and one far above it is downweighted roughly in proportion to
how far it overshoots — never dropped to zero, since that could make an
already-unstable fit rank-deficient. The constant ``1.4826`` converts a
median absolute deviation to a normal-equivalent standard deviation, and
``1.345`` is the standard Huber tuning constant for 95% efficiency under
Gaussian residuals — this is the same M-estimator recipe used for
robust regression generally, not something specific to galvanic
distortion. Compare both modes when outliers are suspected.

.. code-block:: python
   :linenos:

   from pycsamt.emtools.gb import groom_bailey_table

   survey = "data/AMT/WILLY_DATA/L18PLT"

   robust = groom_bailey_table(
       survey,
       band=(1e-3, 10.0),
       robust=True,
       api=False,
   )
   plain = groom_bailey_table(
       survey,
       band=(1e-3, 10.0),
       robust=False,
       api=False,
   )

   compare = robust.merge(
       plain,
       on="station",
       suffixes=("_robust", "_plain"),
   )

   print(
       compare[
           [
               "station",
               "rms_fit_robust",
               "rms_fit_plain",
               "twist_deg_robust",
               "twist_deg_plain",
           ]
       ].head()
   )

.. code-block:: text

      station  rms_fit_robust  rms_fit_plain  twist_deg_robust  twist_deg_plain
   0  18-001A        0.138687       0.138478         10.728457        10.602753
   1  18-002U        0.212166       0.211317          8.600826         9.526479
   2  18-003A        0.278134       0.277152          4.295350         4.190446
   3  18-004A        0.276268       0.274535         17.547552        16.011004
   4  18-005U        0.249210       0.248697          4.692963         3.924469

If robust and non-robust parameters differ strongly, inspect the station
for outlier frequencies, poor dimensionality, or unstable strike.

Synthetic Sanity Check
----------------------

For development and training, it is useful to test the decomposition on
a known distorted 2-D tensor. This example constructs a small synthetic
site-like object with a known distortion matrix and checks whether the
correction reduces diagonal leakage.

.. code-block:: python
   :linenos:

   import numpy as np

   from pycsamt.emtools.gb import groom_bailey_table

   class ZBlock:
       def __init__(self, z, freq):
           self.z = z
           self.freq = freq
           self.z_err = None

   class Site:
       station = "SYN001"

       def __init__(self, z, freq):
           self.Z = ZBlock(z, freq)

   freq = np.logspace(0, 3, 12)
   regional = np.zeros((freq.size, 2, 2), dtype=complex)
   regional[:, 0, 1] = 1.0 + 0.2j
   regional[:, 1, 0] = -0.8 + 0.1j

   D = np.array([[1.0, 0.25], [-0.15, 1.1]])
   observed = D[None, :, :] @ regional
   site = Site(observed, freq)

   table = groom_bailey_table([site], robust=False)

   print(table[["station", "rms_fit", "diagonal_ratio_before", "diagonal_ratio_after"]])

.. code-block:: text

     station       rms_fit  diagonal_ratio_before  diagonal_ratio_after
   0  SYN001  3.678917e-16               0.187225          4.388355e-17

This pattern is useful when you need to verify behavior after changing
preprocessing code. Real surveys should still be assessed with their
own dimensionality and strike diagnostics.

Integrate With Pre-2D Assessment
--------------------------------

The dimensionality guide includes ``pre2d_inversion_assessment``. After
running Groom-Bailey, record whether it was attempted and applied.

.. code-block:: python
   :linenos:

   from pycsamt.emtools.dimensionality import pre2d_inversion_assessment
   from pycsamt.emtools.gb import groom_bailey_decomposition

   survey = "data/AMT/WILLY_DATA/L18PLT"
   band = (1e-3, 10.0)

   gb = groom_bailey_decomposition(
       survey,
       apply=True,
       band=band,
       robust=True,
   )

   assessment = pre2d_inversion_assessment(
       gb.sites,
       band=band,
       rotation_applied=False,
       groom_bailey_attempted=True,
       groom_bailey_applied=gb.applied,
       groom_bailey_reason="Applied pycsamt.emtools.gb real 2-D distortion fit.",
   )

   print(assessment[["station", "frac_3d", "groom_bailey_applied", "recommendation"]].head())

.. code-block:: text

      station   frac_3d  groom_bailey_applied               recommendation
   0  18-001A  0.974359                  True  review_3d_effects_before_2d
   1  18-002U  0.974359                  True  review_3d_effects_before_2d
   2  18-003A  0.974359                  True  review_3d_effects_before_2d
   3  18-004A  0.974359                  True  review_3d_effects_before_2d
   4  18-005U  0.923077                  True  review_3d_effects_before_2d

This makes the correction auditable in reports and manuscripts.

Reading The Results
-------------------

Use this interpretation order:

1. Confirm that dimensionality and strike are acceptable in the selected
   period band.
2. Fit ``groom_bailey_table`` without applying correction.
3. Inspect ``status``, ``n_freq``, ``rms_fit``, and diagonal ratios.
4. Compare robust and non-robust fits if outliers are likely.
5. Apply correction only to stations with acceptable fits.
6. Save the table and pre-2D assessment with the inversion inputs.

Common Failure Modes
--------------------

Insufficient frequencies
   The selected period band has fewer than ``min_freq`` valid tensor
   rows. Widen the band or skip the station.

High residual fit
   The station may not be well described by a frequency-independent
   real distortion matrix times a 2-D regional tensor.

Diagonal ratio does not improve
   Correction may not be useful for that station. Review strike,
   dimensionality, and the period band.

Very large twist, shear, or anisotropy
   Large parameters may indicate strong galvanic distortion, but they
   can also indicate a poor model assumption.

Treating gain as static shift
   The scalar gain ambiguity is not uniquely solved here. Use static
   shift workflows and independent constraints when gain matters.

Applying correction globally
   Do not apply every fitted row blindly. Filter by status and quality
   diagnostics first.

Saving A Reproducible Bundle
----------------------------

Save the fitted table, accepted subset, and pre-2D assessment.

.. code-block:: python
   :linenos:

   from pathlib import Path

   from pycsamt.emtools.dimensionality import pre2d_inversion_assessment
   from pycsamt.emtools.gb import apply_groom_bailey, groom_bailey_table

   survey = "data/AMT/WILLY_DATA/L18PLT"
   band = (1e-3, 10.0)
   out = Path("outputs/gb_l18plt")
   out.mkdir(parents=True, exist_ok=True)

   table = groom_bailey_table(survey, band=band, robust=True)
   accepted = table.loc[
       (table["status"] == "ok")
       & (table["rms_fit"] < 0.25)
       & (table["diagonal_ratio_after"] < table["diagonal_ratio_before"])
   ].copy()

   corrected = apply_groom_bailey(survey, table=accepted, inplace=False)
   assessment = pre2d_inversion_assessment(
       corrected,
       band=band,
       groom_bailey_attempted=True,
       groom_bailey_applied=True,
       groom_bailey_reason="Applied accepted Groom-Bailey station fits.",
   )

   table.to_csv(out / "groom_bailey_table.csv", index=False)
   accepted.to_csv(out / "groom_bailey_accepted.csv", index=False)
   assessment.to_csv(out / "pre2d_assessment_after_gb.csv", index=False)

Worked Workflow
---------------

There is currently no dedicated Sphinx-Gallery ``plot_gb.py`` example in
``docs/examples/emtools``. Until one is added, use the examples in this
page as the worked workflow:

- estimate the table without applying correction;
- filter stations by fit quality;
- apply correction to accepted stations;
- record the result in the pre-2D assessment.
