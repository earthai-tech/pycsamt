.. _ai_processing_classify:

Learned Dimensionality Classification
========================================

:func:`~pycsamt.emtools.dimensionality.classify_dimensionality` labels
each (station, frequency) observation 1-D, 2-D, or 3-D from two
thresholds on phase-tensor :term:`skew` and ellipticity -- transparent,
but a hard boundary at a fixed skew value treats a sample just below the
cutoff identically to one far below it, and ignores the other three
phase-tensor-derived quantities that
:func:`~pycsamt.emtools.dimensionality.phase_features_table` already
computes.
:class:`~pycsamt.ai.processing.classify.DimensionalityClassifier` is a
:term:`multi-layer perceptron` (MLP) trained on all five features at
once --
:math:`[\,|\beta|,\ \text{ellipticity},\ \log_{10}\rho_{\det},\
\phi_{\det},\ |T|\,]` -- with two output heads sharing one backbone: a
3-class softmax for dimensionality, and a linear regression head for
:term:`geoelectric strike` on the 2-D subset.

Self-training from the rule-based labels
-------------------------------------------

:meth:`~pycsamt.ai.processing.classify.DimensionalityClassifier.from_features_table`
builds training labels *from* the same rule-based thresholds it is meant
to refine when no ``label_col`` is supplied -- the network's job is not
to invent a new definition of dimensionality, but to interpolate the
rule-based boundary more smoothly across all five features jointly,
rather than two features and a hard cutoff.

.. code-block:: pycon

   >>> from pycsamt.emtools.dimensionality import phase_features_table
   >>> feats = phase_features_table(sites)
   >>> len(feats), feats["station"].nunique()
   (1484, 28)

.. warning::

   The strike head is a *separate* branch that only receives gradient
   updates through its own loss term, added only when real strike
   targets are supplied via ``strike=`` (in :meth:`fit`) or
   ``strike_col=`` (in :meth:`from_features_table`). Call
   ``from_features_table(feats)`` with neither, as the self-training
   description above might suggest, and the strike head trains on
   nothing at all -- :meth:`predict_strike` still returns numbers, but
   they carry no information and should not be interpreted physically.
   Always supply real strike targets when the predicted strike matters,
   as done next.

Supervising the strike head with a classical estimate
---------------------------------------------------------

:func:`~pycsamt.emtools.strike.estimate_strike_phase_tensor` computes one
classical, per-station consensus strike angle from the same phase-tensor
table. Broadcasting it back onto every frequency row of its station
turns it into a real ``strike_col`` teacher signal -- the strike-head
counterpart of the rule-based dimensionality labels:

.. code-block:: pycon

   >>> from pycsamt.emtools.strike import estimate_strike_phase_tensor
   >>> from pycsamt.ai.processing import DimensionalityClassifier
   >>> consensus = estimate_strike_phase_tensor(sites)
   >>> consensus[["station", "ang", "iqr"]].head(3)
     station        ang        iqr
   0 18-001A -50.071764  28.974476
   1 18-002U -43.757578  16.329767
   2 18-003A -38.949283  18.993311
   >>> feats = feats.merge(
   ...     consensus[["station", "ang"]], on="station", how="left",
   ... ).rename(columns={"ang": "strike_target"})
   >>> clf = DimensionalityClassifier.from_features_table(
   ...     feats, strike_col="strike_target", epochs=120, seed=0,
   ...     verbose=False,
   ... )
   DimensionalityClassifier(n_classes=3, torch, fitted)
   >>> table = clf.predict_table(sites)
   >>> table["dim_label"].value_counts()
   3D    1368
   2D     116
   Name: count, dtype: int64
   >>> n2d = int((table["dim"] == 1).sum())
   >>> table.loc[table["dim"] == 1, "strike"].mean(), \
   ...     table.loc[table["dim"] == 1, "strike"].std()
   (-42.9, 3.3)
   >>> consensus["ang"].mean()
   -39.4

Averaged over the 116 samples the network itself classified as 2-D, the
predicted strike (:math:`-42.9^\circ \pm 3.3^\circ`) lands within four
degrees of the classical per-station consensus averaged over all 28
stations (:math:`-39.4^\circ`) -- direct evidence that the strike head
learned a physically meaningful regression rather than collapsing to an
arbitrary constant, which is precisely the failure mode the warning
above describes for the unsupervised case.

.. figure:: /images/user_guide/ai_processing/dimensionality_summary.png
   :align: center
   :width: 95%

   :func:`~pycsamt.ai.processing.plot.plot_dimensionality_summary` --
   the station x frequency class map, class proportions, and the 2-D
   strike rose (folded to the standard 180-degree line-ambiguity
   convention, with the circular mean overlaid as a dashed line).

L18PLT classifies as predominantly 3-D (92%), with a 2-D minority (8%)
concentrated at short period along most of the line and a handful of
long-period 2-D windows near the profile's last few stations. This
matches the picture built independently in :doc:`qc`: the stations whose
median QC score collapsed under the hard skew rule are exactly the ones
a genuinely 3-D subsurface response would be expected to produce. 1-D
observations are rare in the rule-based training labels themselves --
well under 5% of the L18PLT dataset -- and this run predicted none at
all, which is not a bug: with training labels this imbalanced, a
self-trained network is not guaranteed to recover the minority class
every time. If 1-D detection matters for your survey, check the raw
rule-based label balance first and consider class weighting.

.. note::

   The exact counts above will drift a little from run to run --
   :class:`DimensionalityClassifier`'s PyTorch weight initialization is
   not seeded (the same caveat noted for
   :class:`~pycsamt.ai.processing.anomaly.AnomalyDetector` in
   :doc:`anomaly`), so ``fit(seed=0)`` reproduces the data split
   exactly but not the network's starting point. The qualitative
   picture -- overwhelmingly 3-D, a real 2-D minority near short
   period, an unreliable 1-D minority, and a predicted strike within a
   few degrees of the classical consensus -- is stable across runs
   even when the exact tallies are not.

.. figure:: /images/user_guide/ai_processing/dimensionality_training.png
   :align: center
   :width: 70%

   :func:`~pycsamt.ai.processing.plot.plot_training_history` -- both
   curves report classification cross-entropy only, on the same scale,
   so they can be compared directly even though the strike MSE term
   also contributes to the gradient during training.

Parameters
-----------

``hidden`` (default ``(128, 64)``) sets the shared-backbone width;
``dropout`` (default ``0.2``) regularizes it. The strike loss weight is
fixed at 0.1 relative to the classification cross-entropy, low enough
that strike supervision refines the shared features without dominating
the primary classification objective. Falls back to a random-forest
classifier (scikit-learn) when neither PyTorch nor TensorFlow is
installed -- :meth:`predict_strike` then always returns ``NaN``, since
the fallback has no regression head at all.

Full signatures:
:class:`~pycsamt.ai.processing.classify.DimensionalityClassifier`,
:func:`~pycsamt.emtools.dimensionality.phase_features_table`,
:func:`~pycsamt.emtools.strike.estimate_strike_phase_tensor`,
:func:`~pycsamt.ai.processing.plot.plot_dimensionality_map`,
:func:`~pycsamt.ai.processing.plot.plot_predicted_strike_rose`, and
:func:`~pycsamt.ai.processing.plot.plot_dimensionality_summary`.

.. note::

   :func:`~pycsamt.ai.processing.plot.plot_predicted_strike_rose` draws
   a rose from an *already-computed* array of predicted angles (as
   returned by :meth:`predict_table`) and is intentionally distinct from
   :func:`~pycsamt.emtools.strike.plot_strike_rose`, which recomputes
   strike from a site collection using the classical estimators
   described in :doc:`../emtools/strike`.

.. code-dropdown:: ../../../scripts/generate_user_guide_ai_processing_classify_figures.py
   :language: python
   :pyobject: run_classify
   :linenos:
   :title: View the complete DimensionalityClassifier example source

.. container:: pyc-download-cta

   :download:`Point load_l18() at your own EDI folder to adapt this script <../../../scripts/generate_user_guide_ai_processing_classify_figures.py>`
