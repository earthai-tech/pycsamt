.. _user_guide_stratagem_quality_control:

Quality Control
===============

Every station in :doc:`coordinates`'s injected K2 batch now has a real
position, but position alone says nothing about whether its impedance
tensor is trustworthy. :class:`~pycsamt.stratagem.qc.QualityController`
and :class:`~pycsamt.stratagem.qc.FrequencyFilter` are the two answers
to that: one reports, station by station, how good the data actually
is; the other acts on that judgement, removing frequency bins that
should not be trusted. Neither reimplements its statistics from
scratch -- both delegate to :mod:`pycsamt.emtools.qc`,
:mod:`pycsamt.emtools.frequency`, and
:mod:`pycsamt.emtools.remove_noise`, the same functions any other
:mod:`pycsamt.emtools` workflow uses. What they add is Stratagem-specific
wiring: aligning a raw hardware reader's per-frequency SNR mask onto
each station by number (the same :meth:`~pycsamt.stratagem.io.StratagemRawReader.match_to_edis`
lookup from :doc:`loading`, not by position), and persisting the result
for inspection rather than only mutating data in place.

Per-Station QC Report
---------------------

:meth:`~pycsamt.stratagem.qc.QualityController.fit` takes the injected
EDI list from :doc:`coordinates` and, optionally, a fitted
:class:`~pycsamt.stratagem.io.StratagemRawReader` to enrich the report
with hardware coverage:

.. code-block:: pycon

   >>> from pycsamt.stratagem import EDIBatch, CoordinateInjector, StratagemRawReader
   >>> from pycsamt.stratagem.qc import QualityController
   >>> import pandas as pd
   >>> from pathlib import Path
   >>> from tempfile import TemporaryDirectory

   >>> aligned = pd.read_csv("data/stratagem/K2/k2-gps-aligned.csv")
   >>> survey_coords = aligned[aligned["use_for_survey"]]
   >>> with TemporaryDirectory() as tmp:
   ...     coord_csv = Path(tmp) / "coords.csv"
   ...     survey_coords.to_csv(coord_csv, index=False)
   ...     batch = EDIBatch("data/stratagem/K2/k2-edi").fit()
   ...     edi_objects = [e for i, e in enumerate(batch.edi_objects_) if i != 0]
   ...     injector = CoordinateInjector(epsg=32649, order="forward").fit(
   ...         edi_objects, coord_csv,
   ...         easting_col="easting", northing_col="northing",
   ...         elev_col="elev", station_col="edi_file",
   ...     )

   >>> rdr = StratagemRawReader("data/stratagem/K2/k2-HX", component="X").fit()
   >>> qc = QualityController().fit(injector.edi_objects_, raw_reader=rdr)
   >>> cols = ["station", "n_freq", "snr_med", "skew_med", "hw_coverage"]
   >>> qc.report_[cols].round(2).head(3)
      station  n_freq  snr_med  skew_med  hw_coverage
   0  Z2HX002      39   687.98     17.48         0.83
   1  Z2HX003      39  1491.65     20.06         0.83
   2  Z2HX004      39   932.54     28.87         0.83

``hw_coverage`` is the same per-station fraction from
:meth:`~pycsamt.stratagem.io.StratagemRawReader.station_frame` in
:doc:`loading` (0.83 there too, for the same stations), joined in by
station number rather than assumed to line up positionally. ``snr_med``
is a different, complementary signal -- an impedance-domain SNR
computed from each frequency row's ``Z``/``Z.z_err`` ratio, not the
hardware stack count. A station can have excellent stack coverage and
still produce a noisy impedance estimate, or the reverse, which is
exactly why the report keeps both rather than merging them into one
score.

QC Flag Thresholds
------------------

:meth:`~pycsamt.stratagem.qc.QualityController.summary` condenses the
report into flag counts against three thresholds --
``min_frac_ok``, ``min_snr_med``, ``max_skew_med``:

.. code-block:: pycon

   >>> print(qc.summary())
   QualityController: 86 stations
     flagged  : 86 (100%)
     frac_ok  : 1.00 mean
     snr_med  : 2767.1 median
     skew_med : 38.8° median
     flag breakdown:
       high_skew: 86

Every station passes ``frac_ok``/``min_snr_med`` comfortably -- K2's
impedance coverage is complete and its median SNR is in the thousands,
nowhere near the default ``min_snr_med=2.0`` floor -- and yet every
single station is flagged, all for the same reason:
:term:`Phase tensor` skew. The default ``max_skew_med=6.0`` degrees is
nowhere close to what this line actually produces:

.. code-block:: pycon

   >>> qc.report_["skew_med"].describe().round(2)
   count    86.00
   mean     37.19
   std      14.16
   min       8.84
   25%      26.14
   50%      38.82
   75%      45.89
   max      67.34
   Name: skew_med, dtype: float64

Not one station comes close to 6°; the *minimum* across the whole line
is 8.84°. A threshold that flags literally every station is not telling
you the survey is unusable -- it is telling you the default does not
fit this survey. A skew this large and this consistent across 86
stations reads as genuine 2-D/3-D geoelectric structure along the line
(plausible for the hillside terrain the field report describes), not
measurement noise, and treating "flagged" as "bad data" here would be a
misreading of what the flag means. Picking a threshold that actually
discriminates within this dataset -- rather than adopting the module
default unexamined -- changes the picture substantially:

.. code-block:: pycon

   >>> for thresh in (6.0, 40.0):
   ...     qc_t = QualityController(max_skew_med=thresh).fit(injector.edi_objects_)
   ...     print(thresh, len(qc_t.flagged_stations()))
   6.0 86
   40.0 40

Raising ``max_skew_med`` to 40° -- roughly this line's own median --
flags 40 of 86 stations instead of all of them: the ones whose skew
sits meaningfully above the rest of the line, which is a more useful
signal than a threshold every station fails identically.
:meth:`~pycsamt.stratagem.qc.QualityController.flagged_stations` returns
the station names behind either count, ready to feed into a manual
review or an exclusion list.

Frequency Filtering Stages
--------------------------

:class:`~pycsamt.stratagem.qc.FrequencyFilter` acts on the same EDI
list, removing frequency bins rather than just reporting on them. It
runs three passes in a fixed order -- hardware mask, then band
selection, then a statistical incoherent-frequency mask -- and keeps a
running count of how many (station, frequency) cells each stage
actually touched:

.. code-block:: pycon

   >>> from pycsamt.stratagem.qc import FrequencyFilter

   >>> filt = FrequencyFilter(fmin=10.0, fmax=10000.0)
   >>> filt.fit(injector.edi_objects_, raw_reader=rdr)
   FrequencyFilter(fmin=10.0, fmax=10000.0, snr_thresh=2.5, n_masked_hw_=460, n_dropped_band_=2832)
   >>> filt.n_masked_hw_, filt.n_dropped_band_, filt.n_masked_stat_
   (460, 2832, 0)

``fmin``/``fmax`` here match the band :class:`~pycsamt.stratagem.survey.StratagemSurvey`'s
own K2 example uses in :doc:`pipeline`. 460 cells fail hardware SNR,
2832 more fall outside the 10-10 000 Hz band -- and the third,
statistical stage finds nothing left to remove. That is not the stage
doing nothing: the first two passes already removed the frequencies
most likely to fail a coherence check, so by the time
``mask_incoherent_freqs`` runs, what remains genuinely clears its
``snr_thresh``/``min_frac`` bar. On a delivery without a raw reader to
supply the hardware mask, this third stage is where most of that
cleanup work would fall instead.

Underlying Station Matching
---------------------------

Both :meth:`~pycsamt.stratagem.qc.QualityController.fit`'s hardware
enrichment and :meth:`~pycsamt.stratagem.qc.FrequencyFilter.fit`'s
hardware-mask stage call
:meth:`~pycsamt.stratagem.io.StratagemRawReader.match_to_edis`
internally, for the same reason :doc:`loading` introduced it: raw and
EDI sequences do not necessarily start at the same station. K2's own
raw reader still has all 87 stations, but ``injector.edi_objects_`` here
has only the 86 that survived dropping the calibration shot in
:doc:`coordinates` -- so EDI list index 0 is station 2, not station 1:

.. code-block:: pycon

   >>> mapping = rdr.match_to_edis(injector.edi_objects_)
   >>> mapping[0], rdr.stations_[mapping[0]], injector.edi_objects_[0].station
   (1, 'X2HX.002', 'Z2HX002')

EDI list index 0 (``Z2HX002``) correctly resolves to raw index 1
(``X2HX.002``) -- the one-station offset from dropping station 1 is
absorbed automatically, exactly as it needs to be for ``hw_coverage``
and the hardware-mask stage to line up the right raw station with the
right EDI, not merely the one sitting at the same list position.
:doc:`processing` picks up from here, applying static-shift correction
and noise removal to the same frequency-filtered data.
