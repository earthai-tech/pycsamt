.. _ai_processing_imputer:

Learned Gap Filling
======================

:doc:`denoise` assumes every (station, frequency) cell has *some*
value and improves its quality. Real surveys also have cells with *no*
value at all -- a frequency dropped mid-acquisition, a component that
never converged, a band skipped for poor coherence -- and nothing else
in pyCSAMT fills those gaps with anything better than carrying ``NaN``
straight through to inversion.
:class:`~pycsamt.ai.processing.imputer.EMImputer` targets that gap
directly: reconstructing the *missing* cells from the *observed* ones
using spectral smoothness -- the same station, nearby frequencies --
a matrix-completion problem, not a denoising one.

The recipe is the masked-reconstruction family, the same idea behind
BERT-style masked-language modelling or image inpainting:

.. math::
   :label: eq-ai-imputer-objective

   \hat{\mathbf{x}} = f_\theta(\mathbf{x} \odot \mathbf{m},\ \mathbf{m}),
   \qquad
   \mathcal{L} = \frac{\big\| (\hat{\mathbf{x}} - \mathbf{x}) \odot
   \mathbf{m}_\text{syn} \big\|_2^2}{\sum \mathbf{m}_\text{syn}}

where :math:`\mathbf{m} \in \{0, 1\}` is the mask fed to the network
(1 = "this cell is observed") and :math:`\mathbf{m}_\text{syn}` marks
the subset of *genuinely observed* cells synthetically hidden on that
training step -- exactly where the loss is evaluated. On every step,
:meth:`~pycsamt.ai.processing.imputer.EMImputer.fit` hides a random
``mask_frac`` fraction of the truly-observed cells (any real ``NaN``
gaps already in the training data are always hidden and never used as
a target, since their true value is unknown), mirroring the way
:class:`~pycsamt.ai.processing.denoise.EMDenoiser` adds synthetic
noise. Training data therefore only needs to be *mostly* complete,
not to come with real missing-data examples of known ground truth.

The network is the same 1-D convolutional encoder-decoder family as
:class:`~pycsamt.ai.processing.denoise.EMDenoiser`, with one
structural difference: the encoder takes ``2 * n_components`` input
channels -- the zeroed-where-hidden values, concatenated with the
observed-cell mask along the channel axis -- so it can distinguish
"zero because hidden" from "zero because that is genuinely the
value", while the decoder still only reconstructs ``n_components``
channels; the mask itself is never reconstructed.

A survey with real gaps
---------------------------

L18PLT, used throughout the rest of this section, has no genuinely
missing cells to fill -- every earlier page here either works with
complete data or synthesises noise on top of it. This page instead
uses ``data/MT/broken-hill``: 21 real "ultra-wide-band" magnetotelluric
soundings recorded by the University of Sydney around the Broken Hill
Pb-Zn-Ag deposit, New South Wales (AlQahtani et al., 2026, open-access
CC-BY). Several of its stations really do have missing frequency rows
-- dropped estimates, most likely from reduced coherence at the
long-period end of a station's recording, where a wideband instrument
has proportionally less time to average down noise.

.. dropdown:: Data and model citation
   :animate: fade-in
   :color: secondary

   AlQahtani, Y., Ozaydin, S., Chatzaras, V., Rey, P. F., & Passos, T.
   (2026). Why does the Broken Hill deposit sit in resistive crust?
   Magnetotelluric evidence for metamorphic decoupling of a
   world-class mineral system. *Journal of Geophysical Research: Solid
   Earth*, **131**, e2026JB035666.
   `<https://doi.org/10.1029/2026JB035666>`__

   AlQahtani, Y., Ozaydin, S., Chatzaras, V., Rey, P. F., & Passos, T.
   (2026). Open data for the article "Why does the Broken Hill deposit
   sit in resistive crust?". *Zenodo*.
   `<https://doi.org/10.5281/zenodo.21272106>`__

   See :doc:`/tutorials/model_broken_hill_mt_3d` for the full survey
   carried through to a published 3-D inversion model, and that
   tutorial's own citation dropdown for the complete provenance.

.. code-block:: pycon

   >>> from pycsamt.ai.processing import EMImputer, prepare_z_features
   >>> from pycsamt.emtools import ensure_sites
   >>> sites = ensure_sites("data/MT/broken-hill/edis", recursive=True)
   >>> X = prepare_z_features(sites, n_components=4)
   >>> X.shape
   (21, 4, 96)
   >>> import numpy as np
   >>> missing = np.isnan(X).any(axis=1)   # (n_sites, n_freqs)
   >>> int(missing.sum())
   43

``prepare_z_features`` resamples every station onto one common
frequency grid (``freq_ref``, the first site's own grid by default) --
43 (station, frequency) rows come out with at least one ``NaN``
component after that resampling, spread across 5 of the 21 stations:

.. code-block:: pycon

   >>> per_site = missing.sum(axis=1)
   >>> from pycsamt.emtools._core import _iter_items, _name
   >>> labels = [_name(ed, i) for i, ed in enumerate(_iter_items(sites))]
   >>> {labels[i]: int(per_site[i]) for i in range(len(labels))
   ...  if per_site[i] > 0}
   {'BH_11_imp': 18, 'BH_19_imp': 10, 'BH_1_imp': 3, 'BH_4_imp': 9, 'BH_5_imp': 3}

.. figure:: /images/user_guide/ai_processing/imputer_summary.png
   :align: center
   :width: 95%

   :func:`~pycsamt.ai.processing.plot.plot_imputer_summary` -- (a) the
   station x frequency missing-data map (note the gaps cluster almost
   entirely at long period, consistent with the coherence explanation
   above); (b) two of the five gapped stations, filled curve in blue
   with reconstructed cells marked in red; (c) the training curve.

Quantitative validation on held-out cells
----------------------------------------------

A genuinely missing cell has no known value to check a reconstruction
against, so the standard masked-reconstruction validation protocol
does the opposite of what a real gap is: it hides a known fraction of
the *already observed* cells, reconstructs them, and compares the
reconstruction to the real, held-back values.

.. code-block:: pycon

   >>> rng = np.random.default_rng(0)
   >>> observed = np.isfinite(X)
   >>> holdout = (rng.random(X.shape) < 0.12) & observed
   >>> X_holdout = X.copy()
   >>> X_holdout[holdout] = np.nan
   >>> int(holdout.sum())
   963
   >>> imp = EMImputer(channels=(32, 64, 32))
   >>> imp.fit(X_holdout, mask_frac=0.15, epochs=150, seed=0, verbose=False)
   EMImputer(n_freqs=96, n_components=4, fitted)
   >>> X_recon = imp.transform(X_holdout)
   >>> y_true, y_pred = X[holdout], X_recon[holdout]

.. warning::

   As with every network-based estimator in this package, ``seed``
   fixes the held-out split and the synthetic training mask but not
   PyTorch's own weight-initialization RNG -- rerunning the block
   above reproduces the same *protocol*, not bit-identical numbers.
   See the equivalent caveat in :doc:`denoise`.

A single pooled RMSE across all four feature channels would repeat
the mistake :doc:`denoise` warns about: log-amplitude and
degree-valued phase channels sit on completely different scales, and
whichever one has the larger numeric range would dominate. Broken
down by channel against a channel-mean "predict nothing" baseline,
every channel improves substantially:

.. list-table::
   :header-rows: 1
   :widths: 30 15 20 20

   * - Channel
     - n cells
     - RMSE, EMImputer
     - RMSE, baseline
   * - :math:`\log_{10}|Z_{xy}|`
     - 234
     - 0.179
     - 0.978
   * - :math:`\phi_{xy}` (deg)
     - 215
     - 8.098
     - 16.480
   * - :math:`\log_{10}|Z_{yx}|`
     - 279
     - 0.182
     - 1.026
   * - :math:`\phi_{yx}` (deg)
     - 235
     - 6.081
     - 17.983

.. figure:: /images/user_guide/ai_processing/imputer_validation.png
   :align: center
   :width: 90%

   :func:`~pycsamt.ai.processing.plot.plot_imputer_validation`, one
   panel per feature component -- true value on the x-axis,
   reconstructed value on the y-axis, dashed line at :math:`y = x`.
   Both log-amplitude channels track the diagonal closely
   (:math:`R^2` 0.966 and 0.968); the phase channels are noisier but
   still clearly informative (:math:`R^2` 0.755 and 0.885), consistent
   with phase generally being the harder quantity to reconstruct from
   neighbouring frequencies alone.

.. figure:: /images/user_guide/ai_processing/imputer_reconstruction.png
   :align: center
   :width: 95%

   :func:`~pycsamt.ai.processing.plot.plot_imputer_reconstruction` at
   full resolution for the four most-gapped stations. The reconstructed
   points (red) sit smoothly on the continuation of each curve's real,
   observed trend rather than jumping off it -- direct visual evidence
   that the network is extrapolating the station's own spectral shape,
   not producing arbitrary values.

Sites in, filled sites out
-------------------------------

:meth:`~pycsamt.ai.processing.imputer.EMImputer.apply` closes the same
sites-in / sites-out loop as :doc:`denoise`'s
:meth:`~pycsamt.ai.processing.denoise.EMDenoiser.apply`: pass a site
collection, get one back, with an ``inplace`` flag controlling whether
the input is mutated. The version fitted here trains on the *full*
array -- the real gaps stay ``NaN`` and excluded from the loss, none
of the synthetic held-out cells from the validation above are
involved:

.. code-block:: pycon

   >>> imp = EMImputer(channels=(32, 64, 32))
   >>> imp.fit(X, mask_frac=0.15, epochs=150, seed=0, verbose=False)
   EMImputer(n_freqs=96, n_components=4, fitted)
   >>> filled_sites = imp.apply(sites, inplace=False)

Checking station ``BH_11_imp`` directly confirms both halves of that
contract: every cell that was already observed comes back unchanged,
and every genuinely missing row now has a reconstructed value.

.. code-block:: pycon

   >>> from pycsamt.emtools._core import _get_z_block
   >>> gi = labels.index("BH_11_imp")
   >>> ed_in = list(_iter_items(sites))[gi]
   >>> ed_out = list(_iter_items(filled_sites))[gi]
   >>> _, z_in, _ = _get_z_block(ed_in, with_errors=False)[:3]
   >>> _, z_out, _ = _get_z_block(ed_out, with_errors=False)[:3]
   >>> obs = ~np.isnan(z_in).any(axis=(1, 2))
   >>> bool(np.allclose(z_in[obs], z_out[obs], equal_nan=True))
   True
   >>> was_missing = np.isnan(z_in).any(axis=(1, 2))
   >>> int((was_missing & ~np.isnan(z_out[:, 0, 1])
   ...      & ~np.isnan(z_out[:, 1, 0])).sum())
   20

Every observed cell in ``BH_11_imp`` is untouched, and all 20 of its
missing rows now have :math:`Z_{xy}` and :math:`Z_{yx}` filled.

.. warning::

   With the default ``n_components=4``, only :math:`Z_{xy}` and
   :math:`Z_{yx}` are modelled -- a real missing measurement usually
   means the *whole* frequency row (all four components) is absent,
   and the diagonal :math:`Z_{xx}, Z_{yy}` at those rows is left
   exactly as it was, ``NaN`` if it was ``NaN``, the same
   "not modelled -> untouched" contract
   :meth:`~pycsamt.ai.processing.denoise.EMDenoiser.apply` uses. Fit
   and apply with ``n_components=8`` to fill every component.

Parameters and limitations
-----------------------------

``mask_frac`` (default ``0.15``) controls how much of the observed
data is hidden on every training step -- too low and the network sees
little missing-data signal to learn from; too high and too little
context remains around a hidden cell for it to reconstruct anything.
``channels`` sets the encoder/decoder width at each stage, narrowed to
``(32, 64, 32)`` above for a 21-station survey the same way
:doc:`denoise` narrows it for L18PLT. ``n_components`` (``4`` or
``8``) trades model complexity against diagonal-component coverage,
as shown above.

.. note::

   When neither PyTorch nor TensorFlow is installed,
   :class:`~pycsamt.ai.processing.imputer.EMImputer` falls back to
   per-sample, per-channel linear interpolation along the frequency
   axis -- no cross-station information, unlike the network path, and
   a channel with no observed value anywhere in the training batch
   falls back further still, to that channel's training mean.

Full signatures:
:class:`~pycsamt.ai.processing.imputer.EMImputer`,
:func:`~pycsamt.ai.processing.denoise.prepare_z_features`,
:func:`~pycsamt.ai.processing.plot.plot_imputer_gaps`,
:func:`~pycsamt.ai.processing.plot.plot_imputer_validation`,
:func:`~pycsamt.ai.processing.plot.plot_imputer_reconstruction`, and
:func:`~pycsamt.ai.processing.plot.plot_imputer_summary`.

.. code-dropdown:: ../../../scripts/generate_user_guide_ai_processing_imputer_figures.py
   :language: python
   :pyobject: run_imputer
   :linenos:
   :title: View the complete EMImputer example source

.. container:: pyc-download-cta

   .. container:: pyc-download-cta-text

      **Adapt this to your own survey** -- download the script below
      and point ``load_broken_hill()`` at your own EDI folder.

   :download:`Download script <../../../scripts/generate_user_guide_ai_processing_imputer_figures.py>`
