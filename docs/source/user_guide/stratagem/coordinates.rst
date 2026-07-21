.. _user_guide_stratagem_coordinates:

Coordinate Injection
====================

:doc:`loading` ended with every K2 station still sitting at
``0.0, 0.0, 0.0`` -- the WinGLink placeholder from :doc:`concepts`. The
real positions live in a separate GPS deliverable, almost always in a
projected :term:`CRS` rather than WGS84 directly: metre-scale easting
and northing on some local zone, not longitude and latitude.
:class:`~pycsamt.stratagem.gis_correct.CoordinateInjector` is what turns
that table into real ``LAT``/``LONG``/``ELEV`` values in each EDI's
``>HEAD`` section, converting through
:func:`~pycsamt.gis.utils.project_point_utm2ll`. Getting there for a
real delivery takes two separate problems solved in order: figuring out
which columns of the table actually hold the coordinates, and making
sure the table has exactly one row per station before injection can run
at all. K2's own GPS table fails both, in ways worth walking through
directly rather than assuming a clean input.

Coordinate Column Detection
---------------------------

K2's raw coordinate export, ``k2-gps-ref.csv``, is a
:term:`Gauss-Kruger` table with the columns literally named ``lon`` and
``lat`` -- and neither one is a longitude or a latitude:

.. code-block:: pycon

   >>> import pandas as pd

   >>> gps_ref = pd.read_csv("data/stratagem/K2/k2-gps-ref.csv", encoding="utf-8-sig")
   >>> gps_ref.columns = [c.strip() for c in gps_ref.columns]
   >>> gps_ref[["lon", "lat", "elev", "step"]].head(3)
              lon          lat      elev     step
   0  2850828.525  362607.5696  257.9511   0.0000
   1  2850821.854  362626.9799  257.4241  19.7991
   2  2850815.405  362646.2191  257.4001  40.3187

Those are metre-scale Gauss-Kruger values, not decimal degrees, and the
``lon`` column -- 2.85 million -- is actually the northing;
``lat`` -- 362 thousand -- is the easting. This is a real, common
habit in Chinese field deliveries, not a K2-specific mistake, which is
why :class:`~pycsamt.stratagem.gis_correct.CoordinateInjector` never
trusts a column name: its private ``_detect_coord_cols`` helper compares
the *median magnitude* of the numeric candidate columns instead, on the
reasoning that northing (distance from the equator) is reliably larger
than easting for any standard UTM/Gauss-Kruger zone in the northern
hemisphere:

.. code-block:: pycon

   >>> from pycsamt.stratagem.gis_correct import _detect_coord_cols

   >>> e_col, n_col = _detect_coord_cols(gps_ref, None, None, exclude={"elev", "step"})
   >>> e_col, n_col
   ('lat', 'lon')

The detector correctly picks ``lat`` as the easting column and ``lon``
as the northing column -- the exact reverse of what the names claim --
because it never looks at the names at all. That heuristic only has a
sound basis when exactly two numeric candidates remain once
non-coordinate columns are excluded. The reconciled coordinate table
built later in this page, ``k2-gps-aligned.csv``, has five numeric
columns rather than two, and detection refuses to guess among them:

.. code-block:: pycon

   >>> aligned_raw = pd.read_csv("data/stratagem/K2/k2-gps-aligned.csv")
   >>> aligned_raw.select_dtypes(include="number").columns.tolist()
   ['station_num', 'chainage_design_m', 'northing', 'easting', 'elev']
   >>> _detect_coord_cols(aligned_raw, None, None, exclude={"elev"})
   Traceback (most recent call last):
       ...
   pycsamt.exceptions.ValidationError: Cannot auto-detect easting/northing unambiguously: 4 numeric candidate columns remain after exclusion: ['station_num', 'chainage_design_m', 'northing', 'easting'].  Picking the min/max-median pair among more than two candidates is a guess (e.g. a station index or chainage column could be mistaken for a coordinate). Pass explicit easting_col and northing_col to CoordinateInjector.fit(), or add the extra columns to the exclusion set.

``station_num`` and ``chainage_design_m`` are both numeric and both
plausible-looking by magnitude, so guessing between them and the real
``easting``/``northing`` pair would be exactly the silent
misidentification the module docstring warns about. Passing
``easting_col="easting", northing_col="northing"`` explicitly -- shown
later in this page -- sidesteps the ambiguity entirely rather than
trying to make the heuristic smarter.

Station-Count Reconciliation
----------------------------

:class:`~pycsamt.stratagem.gis_correct.CoordinateInjector` requires
*exactly* one coordinate row per EDI station -- and ``k2-gps-ref.csv``
has 80 rows against 87 K2 stations. Cross-referencing the field report
against the raw table (done once by
``data/stratagem/K2/build_gps_aligned.py``, reproducible with
``python data/stratagem/K2/build_gps_aligned.py``) turns up four
distinct real-world reasons for that gap, not one:

* Station 1 is a calibration/"parallel test" shot at chainage 1 m, off
  the real profile -- a file exists, but there is no corresponding
  survey position to inject.
* Stations 43, 69, 75, and 87 are "check point" repeats: QC
  measurements taken again at the same physical location as the
  station immediately before them, so they share a chainage with an
  existing GPS row rather than needing a new one.
* One physical point at chainage 560 m was GPS-shot but the instrument
  never recorded a file there (too much road traffic, per the field
  report) -- an orphan GPS row with nothing to attach to.
* Three further chainages (~1080 m, ~1260 m, ~1300 m) were simply never
  GPS-shot.

The raw table hides a second problem underneath the first: its
``station`` column is shifted one row relative to the trustworthy
``step`` (chainage) column. Row 0's station label encodes chainage
19.7991 m, which is row 1's ``step`` value, not row 0's:

.. code-block:: pycon

   >>> gps_ref[["station", "step"]].head(3)
          station     step
   0  K0+019.7991   0.0000
   1  K0+040.3187  19.7991
   2  K0+060.6005  40.3187

and the table's last row is a bad shot outright -- its ``step`` value
claims chainage 1640 m (the far end of the line), but its coordinates
sit almost exactly on top of row 0's, at the profile's *start*:

.. code-block:: pycon

   >>> gps_ref[["station", "lon", "lat", "step"]].iloc[[0, -1]]
           station          lon          lat       step
   0   K0+019.7991  2850828.525  362607.5696     0.0000
   79  K0+000.7991  2850835.864  362589.1739  1640.4414

``build_gps_aligned.py`` treats ``step`` as the only trustworthy
chainage (never ``station``), drops that last row as a mislabeled
re-shot rather than trusting its claimed position, and reconstructs the
true endpoint (chainage 1640 m, stations 86 and 87) by a linear fit
through the last five good points instead. The three ungapped chainages
above are filled by linear interpolation between their chainage-adjacent
neighbours. Every row in the resulting 87-row table records how its
coordinate was actually obtained:

.. code-block:: pycon

   >>> aligned = pd.read_csv("data/stratagem/K2/k2-gps-aligned.csv")
   >>> aligned["source"].value_counts()
   source
   gps_exact           78
   checkpoint_dup       4
   interpolated         3
   test_shot_approx     1
   extrapolated         1
   Name: count, dtype: int64
   >>> aligned[aligned["station_num"].isin([86, 87])][["station_num", "source", "northing", "easting"]]
       station_num          source      northing        easting
   85           86    extrapolated  2.850260e+06  364146.054375
   86           87  checkpoint_dup  2.850260e+06  364146.054375

Stations 86 and 87 land on identical coordinates for a documented
reason, not a bug in the reconciliation: 87 is the check-point repeat
of 86, both shot at the same chainage, so 87 legitimately inherits 86's
extrapolated position rather than getting one of its own.

Running the Coordinate Injection
--------------------------------

With one coordinate row per station secured, injection itself is
direct. Station 1's calibration shot has no survey position -- its
``use_for_survey`` flag is ``False`` -- so it is filtered out of both
the coordinate table and the EDI batch before
:meth:`~pycsamt.stratagem.gis_correct.CoordinateInjector.fit` runs, the
same 0-based ``drop_stations`` idea :class:`~pycsamt.stratagem.survey.StratagemSurvey`
automates in :doc:`pipeline`:

.. code-block:: pycon

   >>> from pycsamt.stratagem import EDIBatch, CoordinateInjector

   >>> survey_coords = aligned[aligned["use_for_survey"]]
   >>> survey_coords.to_csv("k2_coords_for_injection.csv", index=False)
   >>> len(survey_coords)
   86

   >>> batch = EDIBatch("data/stratagem/K2/k2-edi").fit()
   >>> edi_objects = [e for i, e in enumerate(batch.edi_objects_) if i != 0]
   >>> len(edi_objects)
   86

   >>> injector = CoordinateInjector(epsg=32649, order="forward")
   >>> injector.fit(
   ...     edi_objects, "k2_coords_for_injection.csv",
   ...     easting_col="easting", northing_col="northing",
   ...     elev_col="elev", station_col="edi_file",
   ... )
   CoordinateInjector(epsg=32649, utm_zone='49N', order='forward', n_stations_=86, reversed_=False)

   >>> injector.coordinate_frame().head(3)
          station   latitude   longitude  elevation
   0  Z2HX002.edi  25.769107  109.629891   257.9511
   1  Z2HX003.edi  25.769048  109.630085   257.4241
   2  Z2HX004.edi  25.768992  109.630278   257.4001

``epsg=32649`` (UTM Zone 49N WGS84) overrides the class's own default
of ``15921`` (Beijing 1954 / Gauss-Kruger Zone 49) -- neither default is
guaranteed right for a survey outside the specific belt each one covers,
and 32649 is what actually reproduces plausible ground truth here: the
injected line sits around 25.77°N, 109.63°E with elevations in the
200-350 m range, consistent with the hillside terrain the field report
describes for this part of Guangxi. That agreement is corroborating
evidence for the EPSG choice, not proof of it -- an EDI's ``>HEAD``
section has no field for "which EPSG was this projected from", so a
wrong-but-plausible zone would inject silently. Confirm the EPSG against
your own project records before treating injected coordinates as final,
particularly outside a survey area you already know.

Station Ordering
----------------

Stratagem numbers its raw files from the first measurement taken;
nothing guarantees the GPS table was recorded walking the same
direction. :class:`~pycsamt.stratagem.gis_correct.StationLocator`,
called internally by ``CoordinateInjector.fit``, is what maps EDI index
*i* to the correct GPS row -- ``order="forward"`` for ``i -> i``,
``order="reversed"`` for ``i -> n-1-i``, or an explicit ``order="mapping"``
list for anything else. K2's own coordinate table already runs in
acquisition order, which is why the injection above used
``order="forward"`` and reported ``reversed_=False``.

``order="auto"`` exists for the case where you do not already know the
relationship, but it is worth being precise about what it currently
does: it accepts a candidate direction and returns a valid
permutation of GPS rows, but its resolution today always chooses
forward, regardless of the GPS table's actual geographic direction --
the median-northing comparison the class docstring describes is
computed but not yet used to decide between the two:

.. code-block:: pycon

   >>> from pycsamt.stratagem.gis_correct import StationLocator
   >>> import numpy as np

   >>> lats_north_to_south = np.linspace(26.0, 25.5, 5)  # decreasing latitude
   >>> loc = StationLocator(order="auto").fit([None] * 5, lats_north_to_south, np.zeros(5))
   >>> loc.index_map_, loc.reversed_
   ([0, 1, 2, 3, 4], False)

A GPS table running clearly north-to-south while station numbers
increase south-to-north is exactly the case ``order="reversed"`` exists
for, but ``order="auto"`` returns the forward map regardless. Until that
heuristic is filled in, treat ``"auto"`` as a synonym for ``"forward"``
and pass ``order="reversed"`` (or an explicit ``mapping``) whenever you
have independent evidence -- a field sketch, chainage direction, a known
survey convention -- that acquisition ran opposite the table's row
order. :doc:`quality_control` picks up from here, now that every station
has a real position to report QC results against.
