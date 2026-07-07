.. _emtools_remove_noise:

pycsamt.emtools.remove_noise — Noise Removal and Spatial Filtering
=====================================================================

:mod:`pycsamt.emtools.remove_noise` is the largest ``emtools`` module:
power-line notching, log-frequency and rho/phase trend smoothing,
outlier and spatial denoising (Hampel, spatial median, low-rank
RPCA-style), off-diagonal consistency enforcement, frequency masking,
group-trend shrinkage, static-shift/EMAP spatial filters
(AMA/FLMA/TMA), a confidence-gated EMAP filter that reuses
:mod:`pycsamt.emtools.qc`'s composite confidence scores, a full
pipeline, and a family of dedicated QC plots.

Functions
---------

- :func:`~pycsamt.emtools.remove_noise.snr_table`
- :func:`~pycsamt.emtools.remove_noise.emi_mitigation_report`
- :func:`~pycsamt.emtools.remove_noise.notch_powerline`
- :func:`~pycsamt.emtools.remove_noise.smooth_logfreq`
- :func:`~pycsamt.emtools.remove_noise.smooth_rho_phase`
- :func:`~pycsamt.emtools.remove_noise.shrink_to_group_trend`
- :func:`~pycsamt.emtools.remove_noise.remove_noise_pipeline`
- :func:`~pycsamt.emtools.remove_noise.hampel_filter_freq`
- :func:`~pycsamt.emtools.remove_noise.spatial_median_filter`
- :func:`~pycsamt.emtools.remove_noise.rpca_offdiag_denoise`
- :func:`~pycsamt.emtools.remove_noise.enforce_offdiag_consistency`
- :func:`~pycsamt.emtools.remove_noise.mask_incoherent_freqs`
- :func:`~pycsamt.emtools.remove_noise.drop_freqs_manual`
- :func:`~pycsamt.emtools.remove_noise.correct_static_shift`
- :func:`~pycsamt.emtools.remove_noise.apply_emap_filter`
- :func:`~pycsamt.emtools.remove_noise.fixed_length_moving_average`
- :func:`~pycsamt.emtools.remove_noise.trimmed_moving_average`
- :func:`~pycsamt.emtools.remove_noise.confidence_gated_emap_filter`
- :func:`~pycsamt.emtools.remove_noise.emap_filter_report`
- :func:`~pycsamt.emtools.remove_noise.plot_emap_filter_profile`
- :func:`~pycsamt.emtools.remove_noise.plot_emap_filter_psection`
- :func:`~pycsamt.emtools.remove_noise.nr_qc_delta_offdiag_psection`
- :func:`~pycsamt.emtools.remove_noise.nr_qc_snr_gain_profile`
- :func:`~pycsamt.emtools.remove_noise.nr_qc_harmonic_waterfall`
- :func:`~pycsamt.emtools.remove_noise.nr_qc_station_offdiag_curves`
- :class:`~pycsamt.emtools.remove_noise.EMAPFilterResult`

Full signatures and parameter documentation live in the
:doc:`API reference <../api/emtools>`.

Remote Reference and EMI Mitigation
-----------------------------------

The ``emtools`` noise-removal layer operates on estimated transfer
functions stored in EDI/Sites objects.  It does not perform time-series
remote-reference transfer-function estimation internally.  If remote
reference was performed by field software or an external processor, the
resulting remote-referenced EDIs can still be loaded and documented.  If
remote time series were not available, use
:func:`~pycsamt.emtools.remove_noise.emi_mitigation_report` to record that
choice and the post-estimation mitigation applied.

For power-line and EMI contamination, pyCSAMT provides:

* :func:`~pycsamt.emtools.remove_noise.notch_powerline` for 50/60 Hz
  harmonics using masking or interpolation;
* :func:`~pycsamt.emtools.remove_noise.mask_incoherent_freqs` and
  :mod:`pycsamt.emtools.spectra` coherence tables for low-coherence
  samples when cross-spectra are available;
* :func:`~pycsamt.emtools.remove_noise.hampel_filter_freq`,
  :func:`~pycsamt.emtools.remove_noise.spatial_median_filter`, and
  :func:`~pycsamt.emtools.remove_noise.rpca_offdiag_denoise` for isolated,
  spatially coherent, or low-rank off-diagonal noise;
* :func:`~pycsamt.emtools.remove_noise.nr_qc_harmonic_waterfall` and
  :func:`~pycsamt.emtools.remove_noise.nr_qc_snr_gain_profile` for
  before/after QC figures.

For the manuscript response, this means remote reference should be stated
as "not attempted" when no independent remote time series were acquired;
the reproducible mitigation is then the notch/mask/interpolation and
robust denoising chain documented by the report table and QC figures.

Worked Example
--------------

Built on **L18PLT** (``data/AMT/WILLY_DATA/``) throughout, with a
small, honestly constructed dense-frequency synthetic survey wherever
a function specifically needs mains-harmonic contamination that this
real, sparsely log-sampled 53-frequency CSAMT line does not have. It
moves from the SNR diagnostic table through power-line notching,
log-frequency vs. rho/phase trend smoothing, three different
denoising philosophies (per-frequency Hampel, cross-station spatial
median, survey-wide low-rank RPCA), consistency enforcement and
frequency masking, group-trend shrinkage, static-shift and EMAP
spatial filters with their dedicated before/after plots, a
confidence-gated EMAP filter that reuses the ``qc`` example's
composite scores directly (and resurfaces the same lowest-confidence
station, ``18-022U``, through a completely different mechanism), and
finishes with the full pipeline and its four QC diagnostic plots.
Several of this module's conservative defaults turn out to have little
or nothing to do on this particular real, comparatively clean survey —
each such case is verified and reported honestly rather than hidden,
alongside a deliberately loosened setting that confirms the mechanism
still works.

.. include:: ../examples/emtools/plot_remove_noise.rst
   :start-after: .. _sphx_glr_examples_emtools_plot_remove_noise.py:
   :end-before: .. _sphx_glr_download_examples_emtools_plot_remove_noise.py:
