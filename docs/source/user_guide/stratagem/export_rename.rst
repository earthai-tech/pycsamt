.. _user_guide_stratagem_export_rename:

Export and Renaming
===================

:doc:`processing` leaves a batch of in-memory
:class:`~pycsamt.seg.edi.EDIFile` objects, statically shift-corrected
and denoised, that still only exist in Python. Getting them onto disk
in a form ready for :doc:`../emtools/index` or inversion is the last
step, and Stratagem gives two purpose-built ways to do it --
:class:`~pycsamt.stratagem.rename.EDIWriter` and
:class:`~pycsamt.stratagem.rename.EDIRenamer` -- alongside the generic
:meth:`~pycsamt.site.base.Sites.write` already introduced in
:doc:`loading`. All three ultimately call the same
:meth:`~pycsamt.seg.edi.EDIFile.write`; what differs is how each one
decides a station's output filename and whether it touches ``>HEAD``
metadata on the way out.

Writing EDI Files
-----------------

:class:`~pycsamt.stratagem.rename.EDIWriter` is the plain case: write
each object using the name it already carries, changing nothing else.
Continuing from :doc:`processing`'s corrected, denoised ``nr.edi_objects_``
-- rebuilt here the same way that page built it, K2's calibration shot
dropped, coordinates injected, static shift corrected, then denoised:

.. code-block:: pycon

   >>> import pandas as pd
   >>> from pathlib import Path
   >>> from tempfile import TemporaryDirectory
   >>> from pycsamt.stratagem import EDIBatch, CoordinateInjector
   >>> from pycsamt.stratagem.process import StaticShiftCorrector, NoiseRemover
   >>> from pycsamt.stratagem.rename import EDIWriter, EDIRenamer

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
   >>> sc = StaticShiftCorrector(sort_by="lon").fit(injector.edi_objects_)
   >>> nr = NoiseRemover(mains_hz=50.0).fit(sc.edi_objects_)

   >>> with TemporaryDirectory() as out1:
   ...     wr = EDIWriter().fit(nr.edi_objects_, out1)
   ...     print(wr.n_written_, wr.written_[0].name, wr.written_[-1].name)
   86 Z2HX002.edi Z2HX087.edi

Every station keeps its original ``Z2HX0NN.edi`` name and ``DATAID``,
because :meth:`~pycsamt.stratagem.rename.EDIWriter.fit` falls back to
``edi.path.name`` whenever the object still carries the path it was
loaded from -- which every K2 station here does, since none of
:doc:`coordinates`, :doc:`quality_control`, or :doc:`processing`
touched the filename, only the data. ``dataid_prefix`` is what actually
changes identity, standardising every station to a sequential label
regardless of what it was called before:

.. code-block:: pycon

   >>> with TemporaryDirectory() as out2:
   ...     wr2 = EDIWriter(dataid_prefix="S", zero_pad=3).fit(nr.edi_objects_, out2)
   ...     print(wr2.n_written_, wr2.written_[0].name, wr2.written_[-1].name)
   86 S000.edi S085.edi
   >>> nr.edi_objects_[0].station, nr.edi_objects_[-1].station
   ('S000', 'S085')

The prefix does not just change the *filename* -- ``nr.edi_objects_``
itself now reports ``S000``/``S085`` as each station's identity, because
:meth:`~pycsamt.stratagem.rename.EDIWriter.fit` sets
``edi.station = f"{prefix}{i:0{zero_pad}d}"`` before writing, and that
setter is the same one covered next.

Standardised Renaming
---------------------

:class:`~pycsamt.stratagem.rename.EDIRenamer` goes further: rather than
writing under whatever name a station already has, it imposes
``{basename}{index:0Nd}{trailer}.edi`` on every station, from a
directory, a list of paths, or -- as here -- a list of already-loaded
:class:`~pycsamt.seg.edi.EDIFile` objects. Renaming twice into the same
directory, the second time without ``overwrite=True``, shows both the
rename itself and what happens when the destination already exists:

.. code-block:: pycon

   >>> out3 = TemporaryDirectory()
   >>> rn = EDIRenamer(basename="T2.", zero_pad=3).fit(nr.edi_objects_, out3.name)
   >>> rn.n_renamed_, rn.dst_paths()[0].name, rn.dst_paths()[-1].name
   (86, 'T2.000.edi', 'T2.085.edi')
   >>> src, dst = rn.renamed_pairs_[0]
   >>> src.name, "->", dst.name
   ('Z2HX002.edi', '->', 'T2.000.edi')

   >>> rn2 = EDIRenamer(basename="T2.", zero_pad=3).fit(nr.edi_objects_, out3.name)
   >>> rn2.n_renamed_, len(rn2.skipped_)
   (0, 86)
   >>> out3.cleanup()

``renamed_pairs_[0]``'s source name is still ``Z2HX002.edi``, not the
``S000.edi`` the previous section's ``dataid_prefix`` gave the station's
*identity* -- ``EDIRenamer`` resolves each source entry from
``edi.path``, which only ever reflects where the object was originally
loaded from, not its current ``station``/``DATAID`` value. Renaming
changes what a station is called; it does not change where
``pycsamt`` remembers it came from, so ``renamed_pairs_`` stays a
genuine provenance record back to the real WinGLink export even after
several relabelling passes. The second call recognises every destination file
already exists and records each one as skipped rather than repeating
the work or raising: ``n_renamed_`` drops to ``0`` and ``skipped_``
picks up all 86. That is the same idempotent-by-default behaviour
:meth:`~pycsamt.stratagem.gis_correct.CoordinateInjector.export` and
:meth:`~pycsamt.stratagem.qc.FrequencyFilter.out` use, and it is why
re-running a pipeline script against an existing output directory is
safe rather than silently duplicating or clobbering work -- pass
``overwrite=True`` explicitly when clobbering is actually intended.

``update_dataid=True`` (the default) means the new name from the first
call above was not cosmetic: ``EDIFile.station``'s setter writes the
same string into ``>HEAD.DATAID`` *and* every linked section's
``SECTID`` (``MTSECT``, ``SPECTRA``, ``TIMESERIES``) that the object
has, so the file's own internal identity matches its filename on disk,
not just the name on the outside:

.. code-block:: pycon

   >>> nr.edi_objects_[0].station
   'T2.000'

Sites.write() and Identity Persistence
--------------------------------------

:doc:`loading` already showed :func:`~pycsamt.emtools.ensure_sites`
wrapping an ``EDIBatch``'s objects into a
:class:`~pycsamt.site.base.Sites` container; that same container's
:meth:`~pycsamt.site.base.Sites.write` is a perfectly usable third way
to reach disk, and :attr:`~pycsamt.stratagem.survey.StratagemSurvey.sites_`
(covered in :doc:`pipeline`) is built on exactly this path. Wrapping
the *already-renamed* ``nr.edi_objects_`` from the previous section is
worth trying directly rather than assuming it picks up where
``EDIRenamer`` left off:

.. code-block:: pycon

   >>> from pycsamt.emtools import ensure_sites

   >>> sites = ensure_sites(nr.edi_objects_)
   >>> with TemporaryDirectory() as out4:
   ...     paths = sites.write(out4, exist_ok=True)
   ...     names = sorted(p.name for p in paths)
   ...     print(len(paths), names[:2], names[-2:])
   86 ['Z2HX002.edi', 'Z2HX003.edi'] ['Z2HX086.edi', 'Z2HX087.edi']

The ``T2.0NN`` identity is gone -- ``Sites.write`` wrote the *original*
K2 filenames, even though every station's ``DATAID`` reads ``T2.0NN``
at this point (confirmed above). :class:`~pycsamt.site.base.Site`
resolves its own station identity independently, preferring an
in-memory marker (documented as being set by
:meth:`~pycsamt.site.base.Site.rename`) over ``DATAID``, and falling
back to the object's original file stem when that marker is unset --
and ``EDIRenamer``/``EDIWriter`` only ever write ``DATAID`` and the
linked ``SECTID`` fields, never that marker. So the two naming systems
are simply independent: a rename made through
:mod:`pycsamt.stratagem.rename` is invisible to
:class:`~pycsamt.site.base.Sites`, and ``Sites.write`` falls back to
the original on-disk name as if the rename had never happened. This is
worth verifying directly before trusting it either way -- setting the
marker by hand on a couple of test objects did make a small,
hand-built :class:`~pycsamt.site.base.Sites` list respect it, but doing
the same across the full 86-station batch through ``ensure_sites`` did
not reliably carry through, so it is not a dependable workaround. The
safe rule is not to mix the two: if the destination is
``Sites.write``/``StratagemSurvey.sites_.write``, rename through
``Site``/``Sites`` itself rather than through ``EDIRenamer`` beforehand;
if the destination is a plain directory of files, ``EDIWriter``/``EDIRenamer``
already do exactly what they say and there is no need to involve
``Sites`` at all. :doc:`pipeline` closes this section out, composing
every stage from :doc:`concepts` through here into the one fluent
:class:`~pycsamt.stratagem.survey.StratagemSurvey` call.
