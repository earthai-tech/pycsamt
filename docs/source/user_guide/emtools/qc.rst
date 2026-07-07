.. _emtools_qc:

pycsamt.emtools.qc — Quality-Control Confidence Scoring
=============================================================

:mod:`pycsamt.emtools.qc` turns "does this transfer function look
trustworthy?" into numbers: per-station and per-frequency composite
confidence scores built from data coverage, tensor uncertainty,
off-diagonal consistency, diagonal leakage, phase smoothness, and
spatial coherence with neighboring stations — plus a wide family of
plots built on top of those scores, from a Kouadio et al. (2024)-style
along-line profile to an uncertainty-propagated apparent-resistivity
fan chart.

Functions
---------

- :func:`~pycsamt.emtools.qc.build_qc_table`
- :func:`~pycsamt.emtools.qc.qc_flags`
- :func:`~pycsamt.emtools.qc.confidence_ratio`
- :func:`~pycsamt.emtools.qc.station_confidence_table`
- :func:`~pycsamt.emtools.qc.frequency_confidence_table`
- :func:`~pycsamt.emtools.qc.plot_confidence_profile`
- :func:`~pycsamt.emtools.qc.plot_frequency_confidence_psection`
- :func:`~pycsamt.emtools.qc.plot_station_confidence_spectrum`
- :func:`~pycsamt.emtools.qc.plot_station_confidence_dashboard`
- :func:`~pycsamt.emtools.qc.plot_confidence_band_summary`
- :func:`~pycsamt.emtools.qc.plot_qc_quicklook`
- :func:`~pycsamt.emtools.qc.plot_consistency_fan`
- :func:`~pycsamt.emtools.qc.plot_xyyx_crossover_map`
- :func:`~pycsamt.emtools.qc.overlay_noise_cone`
- :func:`~pycsamt.emtools.qc.overlay_spectral_holes`

Full signatures and parameter documentation live in the
:doc:`API reference <../api/emtools>`.

Confidence-Ratio Formulation
----------------------------

The composite confidence ratio (CR) is a bounded weighted mean of finite
diagnostic scores:

.. math::

   \mathrm{CR}_{i,f} =
   \frac{\sum_k w_k s_{k,i,f}\,\mathbf{1}_{s_{k,i,f}\ \mathrm{finite}}}
        {\sum_k w_k\,\mathbf{1}_{s_{k,i,f}\ \mathrm{finite}}},
   \qquad 0 \leq s_k \leq 1.

The default score vector is ``coverage``, ``uncertainty``, ``offdiag``,
``diagonal``, ``phase``, and ``spatial`` with weights
``0.35, 0.20, 0.15, 0.10, 0.10, 0.10``.  The component named
``uncertainty`` is derived from the median transfer-function relative
error, so it is the direct bridge between CR and impedance/coherence-style
uncertainty: higher relative error lowers CR.  When cross-spectra are
available, :mod:`pycsamt.emtools.spectra` provides squared coherence
(:math:`\gamma^2`) and the coherence SNR proxy
:math:`\gamma^2 / (1-\gamma^2)`, which can be used alongside or upstream
of the CR table.

The default manuscript classes are:

* ``CR >= 0.95`` - safe / retained;
* ``0.85 <= CR < 0.95`` - recoverable or marginal, review before inversion;
* ``CR < 0.85`` - reject or down-weight.

The reported ``confidence_err`` is the spread of the available component
scores.  If only one component is available, it falls back to the binomial
standard error :math:`\sqrt{\mathrm{CR}(1-\mathrm{CR})/n}`.

Propagation to Inversion
------------------------

For MARE2DEM exports created from EDI data, CR-derived uncertainty
propagation can be enabled with
``pycsamt.models.mare2dem.edi.make_mt_data_from_edi(...,
confidence_weighting=True)`` or through ``Mare2DEMAgent`` input
``{"confidence_weighting": True}``.  For each station-frequency datum,
pyCSAMT inflates the relative impedance error before writing inversion
standard errors:

.. math::

   \epsilon_{Z,\mathrm{eff}} =
   \epsilon_Z
   \left[\frac{1}{\max(\mathrm{CR}, \mathrm{CR}_{\min})}\right]^p .

The defaults are :math:`\mathrm{CR}_{\min}=0.05` and :math:`p=1`.  The
usual impedance-to-observable propagation is then applied:

.. math::

   \sigma_{\rho_a,\mathrm{eff}} =
   2\rho_a\,\epsilon_{Z,\mathrm{eff}},
   \qquad
   \sigma_{\phi,\mathrm{eff}} =
   \frac{180}{\pi}\epsilon_{Z,\mathrm{eff}} .

The inversion weight is therefore reduced as confidence decreases
(:math:`w \propto 1/\sigma_{\mathrm{eff}}^2`).  Existing error floors
remain active after CR inflation, so CR never lowers a datum's stated
uncertainty.

Worked Example
--------------

Built on **L18PLT** (``data/AMT/WILLY_DATA/``), whose real ``z_err``
tensors make the uncertainty-aware scores and the Monte Carlo fan
chart genuine rather than illustrative. It moves from simple summary
tables (where a surprising finding shows up immediately: every station
gets flagged ``high_skew`` under the default threshold, since this
line's real Bibby skew runs 30-50 degrees — a structural signal, not a
data defect) through the ``presence``-vs-``composite`` confidence
contrast, the along-line profile and period-vs-station pseudo-section,
a worst-station/best-station dashboard comparison that separates *data
quality* from the *structural complexity* already explored in the
``anisotropy``/``impedance`` examples, a survey-wide period-band
summary, a coverage/SNR quicklook, and finishes with genuine
Monte Carlo error propagation and crossover/hole-detection diagnostics
on the same strongly anisotropic station featured throughout the
gallery.

.. include:: ../examples/emtools/plot_qc.rst
   :start-after: .. _sphx_glr_examples_emtools_plot_qc.py:
   :end-before: .. _sphx_glr_download_examples_emtools_plot_qc.py:
