.. _ai_processing_uncertainty:

Learned Uncertainty Calibration
===================================

:doc:`qc` answers "is this (station, frequency) cell trustworthy?"
with a score in ``[0, 1]``. It does not answer "how large should the
error bar on this cell be?" -- yet that second question is exactly
what an inversion error floor needs, and today pyCSAMT, like most MT
software, uses whatever ``z_err`` field the processing software wrote,
or a flat percentage-of-|Z| floor. Both are frequently miscalibrated:
field error propagation is often too optimistic in quiet bands and too
pessimistic near noise spikes.
:class:`~pycsamt.ai.processing.uncertainty.UncertaintyCalibrator`
re-estimates that error directly, regressing a smoothed, internally
consistent fractional error from the same kind of per-cell diagnostic
features :class:`~pycsamt.ai.processing.qc.EMQCScorer` already
extracts.

There is no pre-existing "correct error bar" to train against, the way
:doc:`qc`'s rule-based labels or :doc:`classify`'s skew/ellipticity
thresholds exist. ``UncertaintyCalibrator`` trains in *recalibration*
mode: the regression target is the field-processing ``z_err`` itself,
expressed as a fractional error rather than an absolute value,

.. math::
   :label: eq-ai-uncertainty-target

   e_\text{frac} = \frac{|z_\text{err}|}{|Z|}

so it stays comparable across stations and frequencies with very
different :math:`|Z|`. On L18PLT the raw ``z_err`` column itself spans
three decades (1.3 to 1625) purely from :math:`|Z|` scale, while
:math:`e_\text{frac}` is a tight, physically meaningful range:

.. code-block:: pycon

   >>> from pycsamt.ai.processing import build_uncertainty_features_table
   >>> from pycsamt.emtools import ensure_sites
   >>> sites = ensure_sites("data/AMT/WILLY_DATA/L18PLT", recursive=True)
   >>> feats = build_uncertainty_features_table(sites)
   >>> len(feats), feats["station"].nunique()
   (1484, 28)
   >>> feats["z_err_frac"].min(), feats["z_err_frac"].max(), \
   ...     feats["z_err_frac"].median()
   (0.0152, 0.6702, 0.0547)

The network only ever sees four diagnostic features -- Swift skew,
off-diagonal amplitude asymmetry, and both off-diagonal phases -- so
its output is a smoothed re-estimate driven by signal quality, not a
copy of ``z_err_frac``.

.. note::

   :class:`~pycsamt.ai.processing.qc.EMQCScorer`'s fifth feature,
   ``snr``, is deliberately **not** one of those four, even though
   ``build_uncertainty_features_table`` still reports it as a
   diagnostic column. ``snr`` is built from the exact same pooled
   amplitude/error terms as ``z_err_frac`` -- verified on this data,
   ``corr(1/snr, z_err_frac) = 1 - 2e-16``, identical to
   floating-point precision -- so feeding it in would let the network
   "recalibrate" by trivially inverting one input column rather than
   genuinely re-estimating error from signal *shape*. This is a
   *target-leakage* pitfall -- a feature that deterministically encodes
   the answer -- a different failure mode from the *degenerate*-feature
   pitfall caught in
   :class:`~pycsamt.ai.processing.distortion.DistortionTypeClassifier`
   (Groom-Bailey's ``gain`` column carried no information at all,
   always exactly 1.0), but both were caught the same way: by checking
   a candidate feature's actual behaviour against real data -- its
   variance, or here its correlation with the target -- before
   trusting it.

Fitting and a held-out quantitative check
----------------------------------------------

A field ``z_err`` is itself only an estimate, so validating the
calibrator means checking how well it predicts *held-out* cells'
``z_err_frac`` -- the same masked-reconstruction logic as
:doc:`imputer`'s validation, but simpler here since there is only one
target unit to track rather than several feature components:

.. code-block:: pycon

   >>> import numpy as np
   >>> from pycsamt.ai.processing import UncertaintyCalibrator
   >>> rng = np.random.default_rng(0)
   >>> idx = rng.permutation(len(feats))
   >>> n_val = int(0.2 * len(feats))
   >>> val_df = feats.iloc[idx[:n_val]].reset_index(drop=True)
   >>> train_df = feats.iloc[idx[n_val:]].reset_index(drop=True)
   >>> cal = UncertaintyCalibrator(hidden=(32, 16))
   >>> cal.fit(train_df, epochs=150, seed=0, verbose=False)
   UncertaintyCalibrator(n_features=4, torch, fitted)
   >>> y_true = val_df["z_err_frac"].to_numpy()
   >>> y_pred = cal.transform(val_df)
   >>> rmse = float(np.sqrt(np.mean((y_true - y_pred) ** 2)))
   >>> rmse
   0.0298

Against a "predict the median" baseline (RMSE 0.0369) that is a real,
if modest, ~19% error reduction -- :math:`R^2 = 0.317` over the 296
held-out cells, meaning the four diagnostic features explain roughly a
third of the field error's cell-to-cell variance. The rest is
plausibly genuine measurement-specific noise no station-level feature
can predict, not a sign the calibrator has failed.

.. warning::

   As with every network-based estimator in this package, ``seed``
   fixes the held-out split but not PyTorch's own weight-
   initialization RNG -- rerunning the block above reproduces the same
   *protocol*, not bit-identical numbers. See the equivalent caveat in
   :doc:`denoise`.

.. figure:: /images/user_guide/ai_processing/uncertainty_summary.png
   :align: center
   :width: 95%

   :func:`~pycsamt.ai.processing.plot.plot_uncertainty_summary` -- (a)
   the station x frequency calibrated-error map (log colour scale);
   (b) median fractional error per station, original versus
   calibrated; (c) the held-out validation scatter; (d) the training
   curve.

Panel (b) shows the calibration's main visible effect: the
station-to-station spread shrinks -- the standard deviation of the
28 per-station medians drops from 0.0117 to 0.0062, roughly halved.
The stations whose field ``z_err`` was most extreme move the most:
``18-018A`` (the lowest original median, 0.0388) rises to 0.0516, and
``18-022U`` (the highest, 0.0861) falls to 0.0713, both pulled toward
a common, narrower band -- exactly what "recalibration" should do when
the field estimates are noisier than the signal-quality features
underneath them.

Working directly with site collections
---------------------------------------

:meth:`~pycsamt.ai.processing.uncertainty.UncertaintyCalibrator.predict_table`
scores an entire site collection at once, and
:meth:`~pycsamt.ai.processing.uncertainty.UncertaintyCalibrator.apply`
closes the sites-in / sites-out loop the way
:meth:`~pycsamt.ai.processing.qc.EMQCScorer.apply` does -- but instead
of masking or replacing data outright, it *rescales* ``Z.z_err``: for
every (station, frequency) cell it computes the ratio between the
calibrated fractional error and that station's own original pooled
fractional error, then multiplies *all four* tensor components'
existing ``z_err`` by that one ratio, preserving whatever
per-component structure the field error already carries (for example
a larger diagonal-term error) while shifting the overall level to the
calibrated estimate:

.. code-block:: pycon

   >>> cal = UncertaintyCalibrator(hidden=(32, 16))
   >>> cal.fit(feats, epochs=150, seed=0, verbose=False)
   UncertaintyCalibrator(n_features=4, torch, fitted)
   >>> table = cal.predict_table(sites)
   >>> table[["station", "freq", "z_err_frac",
   ...        "z_err_frac_calibrated"]].head(3)
     station     freq  z_err_frac  z_err_frac_calibrated
   0 18-001A  10400.0    0.029769                0.044529
   1 18-001A   8707.0    0.033041                0.044105
   2 18-001A   7289.0    0.032470                0.043674
   >>> calibrated_sites = cal.apply(sites, inplace=False)
   >>> from pycsamt.emtools._core import _get_z_block, _iter_items
   >>> ed_in = next(_iter_items(sites))
   >>> ed_out = next(_iter_items(calibrated_sites))
   >>> _, z_in, _, ze_in = _get_z_block(ed_in, with_errors=True)
   >>> _, z_out, _, ze_out = _get_z_block(ed_out, with_errors=True)
   >>> bool(np.allclose(z_in, z_out, equal_nan=True))   # Z.z: untouched
   True
   >>> ze_in[:3, 0, 1].round(2)
   array([59.08, 57.87, 51.19])
   >>> ze_out[:3, 0, 1].round(2)
   array([88.37, 77.25, 68.85])

``Z.z`` is byte-identical before and after -- this method only ever
touches ``Z.z_err`` -- and the first station's ``Zxy`` error grows by
roughly 35-50%, consistent with its calibrated median fractional error
(0.0461, panel (b) above) sitting above its own original pooled value
(0.0405).

.. note::

   ``floor_frac`` (default ``None``) enforces a minimum fractional
   error after calibration, e.g. ``cal.apply(sites,
   floor_frac=0.02)`` for a 2% floor -- applied before the rescaling
   ratio, so it interacts with the *calibrated* estimate rather than
   the original one.

Parameters and limitations
-----------------------------

``hidden`` (default ``(64, 32)``) sets the regressor backbone width,
narrowed to ``(32, 16)`` above for a 28-station survey the same way
:doc:`denoise` narrows its own network for L18PLT. ``log_target``
(default ``True``) trains in log10 space -- fractional error is
right-skewed -- and :meth:`transform` always de-logs back to the
original linear scale regardless.

.. warning::

   ``fit(X, y=None)`` self-trains only when ``X`` is a
   ``build_uncertainty_features_table`` DataFrame carrying its own
   ``z_err_frac`` column. There is no rule-based fallback the way the
   classifiers in this package synthesise class labels from
   thresholds -- a plain feature array requires an explicit ``y=``.
   This also means an *empirical* target (a genuine bootstrap or
   repeat-measurement spread computed independently, for example from
   :mod:`pycsamt.seg.spectra`) plugs into the exact same
   :meth:`fit(X, y=...) <pycsamt.ai.processing.uncertainty.UncertaintyCalibrator.fit>`
   call without any API change -- only the automatic ``y=None``
   default is recalibration-mode.

.. note::

   When neither PyTorch nor TensorFlow is installed,
   :class:`~pycsamt.ai.processing.uncertainty.UncertaintyCalibrator`
   falls back to a random-forest regressor (scikit-learn).

Full signatures:
:class:`~pycsamt.ai.processing.uncertainty.UncertaintyCalibrator`,
:func:`~pycsamt.ai.processing.uncertainty.build_uncertainty_features_table`,
:func:`~pycsamt.ai.processing.plot.plot_uncertainty_map`,
:func:`~pycsamt.ai.processing.plot.plot_uncertainty_validation`, and
:func:`~pycsamt.ai.processing.plot.plot_uncertainty_summary`.

.. code-dropdown:: ../../../scripts/generate_user_guide_ai_processing_uncertainty_figures.py
   :language: python
   :pyobject: run_uncertainty
   :linenos:
   :title: View the complete UncertaintyCalibrator example source

.. container:: pyc-download-cta

   :download:`Point load_l18() at your own EDI folder to adapt this script <../../../scripts/generate_user_guide_ai_processing_uncertainty_figures.py>`
