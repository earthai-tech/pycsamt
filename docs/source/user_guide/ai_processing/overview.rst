.. _user_guide_ai_processing_overview:

AI Processing Overview
=========================

Every tool in this section shares the ``fit`` / ``transform`` estimator
pattern, and -- for five of the seven -- can be called directly on a
site collection and chained the same way as a rule-based correction
function, returning a corrected, masked, filtered, filled, or
recalibrated site collection rather than a bare array. ``transform``
is the estimator-interface method name used throughout -- the trained
model applied to new data, in the same sense as a scikit-learn
transformer.

.. admonition:: When to reach for the learned version
   :class: important

   A hard threshold (SNR below 3, skew above 0.3) is easy to justify and
   audit, but it treats every station the same way regardless of how its
   features co-vary. A learned model can flag observations that fail no
   single rule yet are still unusual in combination -- at the cost of
   needing training data, and of predictions that are only as good as the
   examples they were fit on. Use the rule-based tools in
   :doc:`../emtools/qc` and :doc:`../emtools/dimensionality` when you need
   a fully transparent, reproducible criterion; reach for the tools in
   this section when the rules alone leave too many borderline cases, or
   when you want a self-consistent screen across many correlated features
   at once.

:doc:`qc`, :doc:`denoise`, :doc:`anomaly`, and :doc:`classify` follow
the same real dataset -- pyCSAMT's bundled AMT line
``data/AMT/WILLY_DATA/L18PLT`` (28 stations, the line already used
throughout :doc:`../emtools/qc`) -- so those four tools can be read as
one continuous walk through a single survey: score data quality,
suppress noise, flag whole-spectrum outliers, and classify geoelectric
dimensionality with a geoelectric strike estimate for the 2-D subset.
:doc:`imputer` switches to a different real survey,
``data/MT/broken-hill``, for a concrete reason: L18PLT has no
genuinely missing cells to fill, and Broken Hill does.
:doc:`uncertainty` and :doc:`distortion` return to L18PLT -- their
inputs (the field ``z_err`` every station carries, and a real,
physically meaningful spread of Groom-Bailey twist and shear) are
already there. :doc:`ts_denoise` works on raw field time series rather
than impedance spectra, so none of L18PLT's EDI files apply; it uses a
synthetic record with a known clean reference for a genuine before/after
number, then applies the same pipeline to a real burst in pyCSAMT's
bundled long-period recording ``data/MT/TS/kap103as.ts``.

Shared building blocks
-------------------------

* a scikit-learn-style ``fit(X, **kwargs)`` / :meth:`transform(X)
  <pycsamt.ai.processing.qc.EMQCScorer.transform>` interface (see
  :class:`~pycsamt.ai._base.BaseEMProcessor`);
* a PyTorch or TensorFlow backend selected automatically through
  :mod:`pycsamt.backends`, with a scikit-learn or SciPy fallback (Isolation
  Forest, PCA, random forest, or Gaussian smoothing) when neither deep
  learning framework is installed, so the same script runs -- more slowly,
  and with a printed notice -- on a minimal install;
* a ``history_`` property on every network-based estimator (
  :class:`~pycsamt.ai.processing.denoise.EMDenoiser`,
  :class:`~pycsamt.ai.processing.anomaly.AnomalyDetector`,
  :class:`~pycsamt.ai.processing.classify.DimensionalityClassifier`,
  :class:`~pycsamt.ai.processing.imputer.EMImputer`,
  :class:`~pycsamt.ai.processing.uncertainty.UncertaintyCalibrator`),
  plotted with the shared
  :func:`~pycsamt.ai.processing.plot.plot_training_history` function
  introduced in :doc:`denoise` -- :class:`EMQCScorer
  <pycsamt.ai.processing.qc.EMQCScorer>`'s Isolation Forest has no
  training epochs, so it has no ``history_``;
* an ``apply(sites, ...)`` method, introduced in :doc:`denoise`, that
  takes a site collection and hands one back -- the same convention
  rule-based functions such as
  :func:`~pycsamt.emtools.remove_noise.notch_powerline` already use.
  :doc:`ts_denoise` is the one exception: its
  :meth:`~pycsamt.ai.processing.tsdenoise.TimeSeriesDenoiser.apply`
  takes and returns a :class:`~pycsamt.ts.TSData` raw time-series
  record instead, since it works before a site collection's impedance
  tensor even exists. Elsewhere, ``apply``'s shape follows what each
  tool actually does to the data:
  :meth:`EMDenoiser.apply
  <pycsamt.ai.processing.denoise.EMDenoiser.apply>` corrects ``Z.z``
  in place (``inplace=True/False``), :meth:`EMQCScorer.apply
  <pycsamt.ai.processing.qc.EMQCScorer.apply>` masks or interpolates
  over flagged rows (same ``inplace`` option), and
  :meth:`AnomalyDetector.apply
  <pycsamt.ai.processing.anomaly.AnomalyDetector.apply>` *drops*
  flagged stations via :meth:`Sites.select()
  <pycsamt.site.base.Sites.select>` and so has no ``inplace`` option
  at all -- :class:`~pycsamt.site.base.Sites` has no in-place removal
  to mirror, and :meth:`EMImputer.apply
  <pycsamt.ai.processing.imputer.EMImputer.apply>` writes its
  reconstruction back only into cells that were genuinely missing on
  the input, leaving every observed measurement untouched.
  :meth:`UncertaintyCalibrator.apply
  <pycsamt.ai.processing.uncertainty.UncertaintyCalibrator.apply>`
  never touches ``Z.z`` at all -- it *rescales* ``Z.z_err`` toward the
  calibrated level, preserving whatever per-component error structure
  the field estimate already had. :doc:`classify` and :doc:`distortion`
  have no ``apply``: neither produces corrected impedance data, only a
  classification table --
  :meth:`~pycsamt.ai.processing.classify.DimensionalityClassifier.predict_table`
  or
  :meth:`~pycsamt.ai.processing.distortion.DistortionTypeClassifier.predict_table`
  -- which the *user* then routes to the right correction tool, shown
  explicitly in :doc:`distortion`;
* :meth:`save` / :meth:`load` checkpointing to a single ``.npz`` file
  that bundles hyperparameters and weights together.

.. warning::

   None of these tools substitutes for looking at the data. Treat a
   learned score, denoised spectrum, or predicted class the same way you
   would treat any other automated flag: as a prioritization aid for
   review, not as a silent accept/reject gate ahead of inversion.
