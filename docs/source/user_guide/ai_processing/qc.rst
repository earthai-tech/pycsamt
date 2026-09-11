.. _ai_processing_qc:

Machine-Learning Quality Scoring
==================================

:doc:`../emtools/qc` builds its confidence ratio from named, individually
interpretable component scores -- coverage, uncertainty, tensor shape,
phase smoothness. That is a strength for auditability, but it also means
every component is scored independently: nothing in the rule-based table
notices that a station's SNR, skew, and phase asymmetry are jointly
unusual even when no single one crosses its own threshold.
:class:`~pycsamt.ai.processing.qc.EMQCScorer` adds that joint view. It
keeps two hard, auditable rules from the rule-based world -- a minimum
:term:`SNR` and a maximum Swift :term:`skew` -- and combines them with an
:term:`isolation forest` fitted on five per-(station, frequency) features:

* **SNR** -- :math:`|\bar Z| / \sigma(Z)`, from the impedance error tensor
  when available, or from local spectral smoothness (a
  median-filter residual) otherwise;
* **Swift skew** --
  :math:`|\beta_\text{Swift}| = |(Z_{xx}+Z_{yy})/(Z_{xy}-Z_{yx})|`, the
  same regional-distortion diagnostic used throughout
  :doc:`../emtools/tensor`;
* **off-diagonal asymmetry** -- :math:`\log_{10}(|Z_{xy}|/|Z_{yx}|)`;
* **phase of** :math:`Z_{xy}` **and** :math:`Z_{yx}`, in degrees.

Every observation that fails the hard SNR or skew rule is scored 0
outright. Surviving observations get a data-driven score in ``[0, 1]``
from the Isolation Forest's decision function, min-max normalized across
the batch and then combined with the rule outcome by geometric mean, so a
rule pass cannot be overridden by the anomaly model but a rule pass alone
does not guarantee a high score either.

Fitting and scoring
--------------------

Load the bundled L18PLT line the same way :doc:`../emtools/qc` does, then
fit and score in two calls:

.. code-block:: pycon

   >>> from pathlib import Path
   >>> from pycsamt.emtools import ensure_sites
   >>> from pycsamt.ai.processing import EMQCScorer
   >>> sites = ensure_sites(
   ...     Path("data/AMT/WILLY_DATA/L18PLT"),
   ...     recursive=True, on_dup="replace", strict=True,
   ... )
   >>> scorer = EMQCScorer(random_state=0)
   >>> scorer.fit(sites)
   EMQCScorer(ml+rules)
   >>> table = scorer.score_table(sites)
   >>> table[["station", "freq", "score", "flag"]].head()
     station      freq     score  flag
   0 18-001A  10400.0  0.940062     1
   1 18-001A   8707.0  0.954844     1
   2 18-001A   7289.0  0.949452     1
   3 18-001A   6102.0  0.963360     1
   4 18-001A   5108.0  0.982671     1

``score_table`` mirrors
:func:`~pycsamt.emtools.qc.station_confidence_table`'s row shape --
``station``, ``freq``, one column per feature, ``score``, and a binary
``flag`` -- so it can be filtered, merged, or exported the same way. On
this line, 556 of the 1484 station-frequency cells (37%) fall below the
default ``score_threshold=0.5``:

.. code-block:: pycon

   >>> (table["flag"] == 0).sum(), len(table)
   (556, 1484)

.. warning::

   ``random_state=None`` by default, matching
   :class:`sklearn.ensemble.IsolationForest`'s own default: every
   unseeded call to :meth:`fit` trains a slightly different forest, and
   individual scores can shift by a few percent between runs. Pass an
   explicit ``random_state`` (as above) for a reproducible score table;
   the hard SNR/skew rules and their resulting ``flag`` values do not
   depend on the random state at all.

Reading the score chart
------------------------

.. figure:: /images/user_guide/ai_processing/qc_summary.png
   :align: center
   :width: 95%

   :func:`~pycsamt.ai.processing.plot.plot_qc_summary` on the L18PLT
   score table: (a) per-station bars are the *median* score over all 53
   frequencies at that station, with individual per-frequency scores
   shown as translucent dots; (b) the pooled score distribution; (c) the
   same distribution as a violin.

Seven stations -- ``18-004A``, ``18-018A``, ``18-019U``, ``18-022U``,
``18-022V``, ``18-023A``, ``18-024U`` -- sit at exactly zero in panel (a)
even though their scattered per-frequency dots range up to 0.9 or higher.
A median collapses to zero whenever more than half of a station's
frequencies are hard-flagged, and for these seven stations that is
exactly what happens: more than half of their 53 frequencies carry a
Swift skew above the ``skew_threshold=0.3`` cutoff. That is a real
electromagnetic signal, not scoring noise -- these are the stations
:doc:`classify` later classifies as predominantly 3-D -- and it is a
direct consequence of summarizing a per-frequency quantity with a single
per-station number: always check the scatter, not just the bar, before
excluding a station outright.

.. figure:: /images/user_guide/ai_processing/qc_heatmap.png
   :align: center
   :width: 95%

   :func:`~pycsamt.ai.processing.plot.plot_qc_heatmap` -- the same
   scores resolved by station *and* period, with the ``score=0.5``
   contour traced in white.

The heat-map makes the period dependence explicit: red regions cluster
at the long-period (low-frequency) end for most stations, the classic
signature of dead-band and cultural-noise contamination in AMT data, and
the run of solid red across stations 18-018U through 18-023A at short
period corresponds to the same skew-driven rejection seen in the summary
bars.

.. figure:: /images/user_guide/ai_processing/qc_feature_heatmap.png
   :align: center
   :width: 85%

   :func:`~pycsamt.ai.processing.plot.plot_qc_feature_heatmap` -- one
   panel per input feature, letting you check *why* a region scored low
   rather than only *that* it did. The Swift-skew panel (second from
   top) shows the same stations lighting up pale yellow (skew well above
   the default 0.3 rule) that collapsed to zero in the summary bar
   chart.

Sites in, masked sites out
----------------------------

Everything above operates on a ``score_table()`` DataFrame. For
everyday use,
:meth:`~pycsamt.ai.processing.qc.EMQCScorer.apply` follows the same
sites-in / sites-out convention as :doc:`denoise`'s
:meth:`~pycsamt.ai.processing.denoise.EMDenoiser.apply` and the
rule-based :func:`~pycsamt.emtools.remove_noise.notch_powerline`: pass
a site collection, get one back with every hard-flagged
station-frequency row masked (``mode="mask"``, the default) or
replaced by a linear interpolation from its nearest good neighbours
(``mode="interp"``) -- the same two modes ``notch_powerline`` itself
offers, reusing its own interpolation helper.

.. code-block:: pycon

   >>> masked_sites = scorer.apply(sites, mode="mask")
   >>> len(list(masked_sites)) == len(list(sites))
   True

Masking changes ``Z.z``, not the station count: all 28 stations
survive, but every flagged row across them now reads ``NaN`` in all
four impedance components, ready to be dropped or down-weighted by
whatever inversion input builder consumes ``masked_sites`` next.

.. code-block:: pycon

   >>> n_nan_before = 0  # sites: 0 masked rows anywhere
   >>> n_nan_after = 556  # masked_sites: one NaN row per flagged cell
   >>> n_nan_before, n_nan_after
   (0, 556)

.. warning::

   A whole frequency row is masked together -- all four impedance
   components -- even though only the off-diagonal components feed the
   score. This mirrors ``notch_powerline``'s own convention: a
   hard-flagged frequency (failing the SNR or skew rule) is treated as
   unreliable across the board, not selectively for the components
   that happened to trigger the rule.

Parameters
-----------

``contamination`` (default ``0.05``) is the Isolation Forest's expected
outlier fraction and controls how aggressively it separates the "normal"
region of feature space; raise it if you expect a noisier survey,
lower it for a clean, well-instrumented one. ``snr_threshold`` (default
``3.0``) and ``skew_threshold`` (default ``0.3``) are the two hard rules
described above -- tightening either one increases the number of
stations forced to a median of zero the way panel (a) shows.
``score_threshold`` (default ``0.5``) only affects the binary ``flag``
column and the plotted "Good" / "Review" split; it does not change the
underlying score. Set ``use_ml=False`` to disable the Isolation Forest
entirely and fall back to the two hard rules alone (scores are then
strictly 0 or 1) -- useful as a sanity baseline before trusting the
combined score, or when scikit-learn is unavailable.

Full parameter and return-value documentation is in the API reference:
:class:`~pycsamt.ai.processing.qc.EMQCScorer` and
:func:`~pycsamt.ai.processing.plot.plot_qc_scores`,
:func:`~pycsamt.ai.processing.plot.plot_qc_heatmap`,
:func:`~pycsamt.ai.processing.plot.plot_qc_feature_heatmap`,
:func:`~pycsamt.ai.processing.plot.plot_qc_score_distribution`,
:func:`~pycsamt.ai.processing.plot.plot_qc_score_spread`, and
:func:`~pycsamt.ai.processing.plot.plot_qc_summary`. Every plotting
function accepts a plain ``ndarray``, a ``{profile: ndarray}`` dict, or
the ``score_table()`` DataFrame directly, and every station-axis
function accepts a
:class:`~pycsamt.ai.plot.StationTickConfig` for fine control over tick
spacing on long profiles.

.. code-dropdown:: ../../../scripts/generate_user_guide_ai_processing_qc_figures.py
   :language: python
   :pyobject: run_qc
   :linenos:
   :title: View the complete EMQCScorer example source

.. container:: pyc-download-cta

   .. container:: pyc-download-cta-text

      **Adapt this to your own survey** -- download the script below
      and point ``load_l18()`` at your own EDI folder.

   :download:`Download script <../../../scripts/generate_user_guide_ai_processing_qc_figures.py>`
