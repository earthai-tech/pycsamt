.. _user_guide_stratagem_pipeline:

The StratagemSurvey Pipeline
============================

Every page from :doc:`concepts` through :doc:`export_rename` worked one
stage of the K2 survey at a time -- loading, coordinate injection, QC,
static-shift correction, frequency filtering, noise removal, export --
each rebuilt from scratch to keep that page self-contained.
:class:`~pycsamt.stratagem.survey.StratagemSurvey` is the same eight
stages composed into a single object with a fluent, chainable API:
:meth:`~pycsamt.stratagem.survey.StratagemSurvey.fit` loads and injects
coordinates, and every step after it -- :meth:`~pycsamt.stratagem.survey.StratagemSurvey.run_qc`,
:meth:`~pycsamt.stratagem.survey.StratagemSurvey.remove_static_shift`,
:meth:`~pycsamt.stratagem.survey.StratagemSurvey.drop_frequencies`,
:meth:`~pycsamt.stratagem.survey.StratagemSurvey.remove_noises`,
:meth:`~pycsamt.stratagem.survey.StratagemSurvey.export`,
:meth:`~pycsamt.stratagem.survey.StratagemSurvey.rename` -- returns
``self``, so the whole survey reads as one expression. Nothing about
composing them changes what each stage does; the point of this page is
showing that the numbers already verified piece by piece across this
section come out identically when the pipeline runs them end to end.

The fit() Step
--------------

:meth:`~pycsamt.stratagem.survey.StratagemSurvey.fit` is the one step
every other method depends on: it loads the EDI batch (and, if
``raw_dir`` is given, the raw hardware reader), drops any
``drop_stations`` before coordinate injection, and runs
:class:`~pycsamt.stratagem.gis_correct.CoordinateInjector`. K2's
station 1 -- the calibration/test shot from :doc:`concepts` -- is
dropped by its 0-based ``EDIBatch`` index, exactly the way
:doc:`coordinates` did it by hand:

.. code-block:: pycon

   >>> import pandas as pd
   >>> from pathlib import Path
   >>> from pycsamt.stratagem import StratagemSurvey

   >>> aligned = pd.read_csv("data/stratagem/K2/k2-gps-aligned.csv")
   >>> survey_coords = aligned[aligned["use_for_survey"]]
   >>> survey_coords.to_csv("k2_coords_for_injection.csv", index=False)

   >>> sv = StratagemSurvey(
   ...     edi_dir="data/stratagem/K2/k2-edi",
   ...     coord_file="k2_coords_for_injection.csv",
   ...     raw_dir="data/stratagem/K2/k2-HX",
   ...     epsg=32649,
   ...     order="forward",
   ...     drop_stations=[0],
   ...     easting_col="easting",
   ...     northing_col="northing",
   ...     elev_col="elev",
   ...     station_col="edi_file",
   ... ).fit()

   >>> sv.batch_.n_stations_, sv.raw_reader_.n_stations_, sv.n_stations_
   (87, 87, 86)
   >>> sv.edi_objects_[0].station
   'Z2HX002'

``batch_``/``raw_reader_`` still report all 87 stations -- they are the
raw loaders, untouched by ``drop_stations`` -- while ``n_stations_``
and ``edi_objects_`` reflect the 86 that actually went into coordinate
injection, starting at station 2 exactly as :doc:`coordinates` and
:doc:`loading` found by hand. ``coord_file`` still needs one row per
*surviving* station, which is why ``survey_coords`` is filtered by
``use_for_survey`` before being written, the same filtering
:doc:`coordinates` did explicitly.

Chained Processing Steps
------------------------

Every subsequent call mutates ``edi_objects_`` in place and returns
``self``, so the full pipeline -- QC, static-shift correction,
frequency filtering, noise removal -- reads as one continued
expression:

.. code-block:: pycon

   >>> sv = (
   ...     sv
   ...     .run_qc()
   ...     .remove_static_shift()
   ...     .drop_frequencies(fmin=10.0, fmax=10000.0)
   ...     .remove_noises()
   ... )
   >>> print(sv.summary())
   StratagemSurvey
     edi_dir    : data/stratagem/K2/k2-edi
     coord_file : k2_coords_for_injection.csv
     n_stations : 86
     epsg / zone: 32649 / 49N
     raw_reader : yes (87 stations)
     coord order: forward
     QC flags   : 86 / 86 stations
     SS fac_z   : median=0.952
     freq filter: hw=460, band=2832, incoh=0
     noise rm   : 86 stations

Every number here is one already seen in isolation: 86 of 86 stations
flagged at the default QC skew threshold (:doc:`quality_control`), a
median static-shift factor of 0.952 (:doc:`processing`), 460 hardware-
masked and 2832 band-dropped frequency cells with nothing left for the
incoherence stage to catch (:doc:`quality_control` again). Running the
whole chain through one object rather than by hand changes nothing
about the computation -- it is the same
:class:`~pycsamt.stratagem.gis_correct.CoordinateInjector`,
:class:`~pycsamt.stratagem.qc.QualityController`,
:class:`~pycsamt.stratagem.process.StaticShiftCorrector`,
:class:`~pycsamt.stratagem.qc.FrequencyFilter`, and
:class:`~pycsamt.stratagem.process.NoiseRemover` underneath, only
wired together and stored on ``self`` instead of tracked by hand
across separate variables.

Export and Rename
-----------------

:meth:`~pycsamt.stratagem.survey.StratagemSurvey.export` and
:meth:`~pycsamt.stratagem.survey.StratagemSurvey.rename` close the
chain out the same way :doc:`export_rename` did directly with
:class:`~pycsamt.stratagem.rename.EDIWriter`/:class:`~pycsamt.stratagem.rename.EDIRenamer`:

.. code-block:: pycon

   >>> sv = sv.export("k2_corrected").rename(basename="T2.", dst_path="k2_renamed")
   >>> exported = sorted(Path("k2_corrected").glob("*.edi"))
   >>> renamed = sorted(Path("k2_renamed").glob("*.edi"))
   >>> len(exported), exported[0].name, exported[-1].name
   (86, 'Z2HX002.edi', 'Z2HX087.edi')
   >>> len(renamed), renamed[0].name, renamed[-1].name
   (86, 'T2.000.edi', 'T2.085.edi')

``export`` keeps each station's original filename (the same
``EDIWriter``-default behaviour from :doc:`export_rename`), and
``rename`` imposes the sequential ``T2.0NN.edi`` convention on top of
that -- ``export``'s output is untouched by it, since ``rename``
defaults to reading from the most recent export directory rather than
overwriting it. The same reminder from :doc:`export_rename` applies
here too: if the destination were :attr:`~pycsamt.stratagem.survey.StratagemSurvey.sites_`'s
generic ``Sites.write()`` instead of ``rename()``, the ``T2.0NN``
identity would not carry over to it.

Alternative EDI Sources
-----------------------

``edi_dir`` does not require a directory path. Anything
:func:`~pycsamt.emtools.ensure_sites` already accepts -- a
:class:`~pycsamt.site.base.Sites`, an ``EDICollection``, or a plain
list of :class:`~pycsamt.seg.edi.EDIFile` -- works too, which is useful
when a batch is already loaded and processed by other
:mod:`pycsamt.emtools` code before it reaches Stratagem-specific
handling:

.. code-block:: pycon

   >>> from pycsamt.stratagem import EDIBatch
   >>> from pycsamt.emtools import ensure_sites

   >>> batch = EDIBatch("data/stratagem/K2/k2-edi").fit()
   >>> sites = ensure_sites(batch.edi_objects_)
   >>> survey_coords.to_csv("k2_coords_for_injection.csv", index=False)

   >>> sv2 = StratagemSurvey(
   ...     edi_dir=sites,
   ...     coord_file="k2_coords_for_injection.csv",
   ...     epsg=32649,
   ...     order="forward",
   ...     drop_stations=[0],
   ...     easting_col="easting",
   ...     northing_col="northing",
   ...     elev_col="elev",
   ...     station_col="edi_file",
   ... ).fit()
   >>> sv2.batch_ is None, sv2.n_stations_, sv2.edi_objects_[0].station
   (True, 86, 'Z2HX002')

``batch_`` stays ``None`` in this form -- there was no directory listing
to report, only whatever order ``sites`` already had -- but every other
attribute behaves identically, including ``drop_stations`` counting
against ``sites``' own order the same way it counts against
``EDIBatch``'s natural-sorted order.

Validating a Full Run
---------------------

Composing every stage into one call does not remove the reasons
:doc:`processing` and :doc:`quality_control` gave for checking
intermediate results before trusting a correction -- those results are
still sitting on ``sv``, not hidden behind the fluent interface:

.. code-block:: pycon

   >>> len(sv._ss_corrector_.factors_)
   74
   >>> round(sv._ss_corrector_.factors_["fac_z"].max(), 2)
   93.03
   >>> len(sv.qc_.flagged_stations())
   86

Station 61's ``fac_z=93`` outlier from :doc:`processing` is exactly as
reachable here, through ``sv._ss_corrector_.factors_``, as it was from
a standalone :class:`~pycsamt.stratagem.process.StaticShiftCorrector`.
Nothing about running the whole survey as one expression is a reason to
skip reading ``factors_``, ``qc_.flagged_stations()``, or
:attr:`~pycsamt.stratagem.survey.StratagemSurvey.coordinate_frame`
before treating the exported EDIs as final -- the one-call form saves
bookkeeping, not judgement. With that, the K2 line is fully
reconciled, corrected, and exported: from raw hardware files and a
misleading GPS table in :doc:`concepts` to 86 renamed, denoised EDIs
ready for :doc:`../emtools/index` or inversion.
