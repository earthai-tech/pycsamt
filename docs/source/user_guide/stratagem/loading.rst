.. _user_guide_stratagem_loading:

Loading Hardware and EDI Data
=============================

:doc:`concepts` introduced the shape of both file worlds a Stratagem
survey lives in -- the raw 19-column tables and the WinGLink EDI export
-- and used :class:`~pycsamt.stratagem.io.StratagemRawReader` and
:class:`~pycsamt.stratagem.io.EDIBatch` just enough to make that shape
concrete. This page is the fuller tour of both classes: every
diagnostic the raw reader exposes beyond a single coverage number, what
loading actually gives you once an EDI batch is in memory, and the one
place the two worlds have to be reconciled by station *number* rather
than by position -- which is also where a WinGLink export that quietly
skipped a station turns into a silently mismatched dataset if you are
not careful.

Raw Reader Diagnostics
----------------------

:meth:`~pycsamt.stratagem.io.StratagemRawReader.station_frame` answers
"which stations are the problem"; its counterpart,
:meth:`~pycsamt.stratagem.io.StratagemRawReader.freq_frame`, answers the
transposed question, "which frequencies are the problem", one row per
frequency bin instead of one row per station:

.. code-block:: pycon

   >>> from pycsamt.stratagem import StratagemRawReader

   >>> rdr = StratagemRawReader("data/stratagem/K2/k2-HX", component="X").fit()
   >>> rdr.freq_frame().head(8)
      freq_hz  stations_ok  frac_ok  med_stacks
   0     6.25            0      0.0         NaN
   1     7.50            0      0.0         NaN
   2     8.75            0      0.0         NaN
   3    10.00            0      0.0         NaN
   4    11.30           87      1.0        24.0
   5    12.50            0      0.0         NaN
   6    13.80            0      0.0         NaN
   7    15.00           87      1.0        24.0

This is the low-frequency striping from the coverage figure in
:doc:`concepts`, quantified rather than just seen: 6.25, 7.5, 8.75, and
10 Hz are dead across every one of the 87 stations (``stations_ok=0``,
hence a ``NaN`` median with no valid stack to take it over), then 11.3
Hz works everywhere, then two more dead bins, then 15 Hz works
everywhere. That is a property of the *band*, not of any one station --
which is exactly why the earlier figure showed it as clean vertical
stripes rather than scattered noise.

Both ``station_frame()`` and ``freq_frame()`` collapse the full
(station, frequency) grid along one axis.
:meth:`~pycsamt.stratagem.io.StratagemRawReader.stack_audit` keeps both
axes, returning the raw stack-count matrix itself with station names as
the row index and frequency values as columns, for the times a summary
statistic is not enough and you need the actual numbers behind it:

.. code-block:: pycon

   >>> audit = rdr.stack_audit()
   >>> audit.shape
   (87, 292)
   >>> audit.loc["X2HX.005"].iloc[:6]
   6.25      0
   7.50      0
   8.75      0
   10.00     0
   11.30    24
   12.50     0
   Name: X2HX.005, dtype: int32

Station 5 reproduces the same dead/live pattern from ``freq_frame()``
at the level of one station's actual stack counts: zero stacks at the
dead bins, 24 at the live one -- ``station_frame()``'s ``med_stacks``
column for this station is exactly the median across a row like this
one, taken only over the non-zero entries.

Component Selection
-------------------

Every example so far has read ``component="X"``. Passing
``component="ALL"`` reads X, Y, and Z together and keeps each
component's mask separately in
:attr:`~pycsamt.stratagem.io.StratagemRawReader.component_masks_`,
which is convenient for a delivery where every component is trustworthy
-- but K2's is not, and ``"ALL"`` does not protect you from that:

.. code-block:: pycon

   >>> rdr_all = StratagemRawReader("data/stratagem/K2/k2-HX", component="ALL").fit()
   >>> rdr_all.n_stations_, rdr_all.n_freqs_
   (87, 292)
   >>> for comp, (mask, stacks) in rdr_all.component_masks_.items():
   ...     print(comp, mask.shape, round(mask.mean(), 3))
   X (87, 292) 0.794
   Y (87, 4093) 0.01
   Z (87, 116) 0.424

``rdr_all.n_freqs_`` and the top-level ``snr_mask_`` still come from X,
the "primary" component read first, so nothing above looks obviously
wrong at a glance. But ``component_masks_["Y"]`` carries the same 4093-bin
garbage grid from :doc:`concepts` -- ``"ALL"`` reads every component
that exists on disk, it does not check that what it read was the
19-column table it expected. And Z, which *is* genuine ASCII data, has
its own legitimate 116-bin grid rather than X's 292 -- the vertical
magnetic-field channel was simply recorded over fewer frequencies on
this instrument, a real hardware asymmetry rather than a parsing
problem. Neither fact is visible from ``n_stations_``/``n_freqs_``
alone, which is the point of keeping this reader's raw per-component
diagnostics around rather than trusting only the top-level summary.

EDIBatch Sequence Protocol
--------------------------

:class:`~pycsamt.stratagem.io.EDIBatch` wraps its natural-sorted file
list in the plain Python sequence protocol, so it behaves like a list of
:class:`~pycsamt.seg.edi.EDIFile` wherever a list would do -- indexing,
negative indexing, and iteration all work directly against the
natural-sorted order established when :meth:`~pycsamt.stratagem.io.EDIBatch.fit`
ran:

.. code-block:: pycon

   >>> from pycsamt.stratagem import EDIBatch

   >>> batch = EDIBatch("data/stratagem/K2/k2-edi").fit()
   >>> len(batch)
   87
   >>> batch[0].station, batch[-1].station
   ('Z2HX001', 'Z2HX087')
   >>> [edi.station for edi in batch][:3]
   ['Z2HX001', 'Z2HX002', 'Z2HX003']

Nothing here re-reads the directory or re-sorts anything -- ``batch[-1]``
is ordinary Python list indexing against ``edi_objects_``, which is why
it lands on station 87 rather than whatever file a filesystem happened
to list last.

Station-Number Matching
-----------------------

K2's raw and EDI directories happen to cover the exact same 87 stations,
so pairing ``rdr.stations_[i]`` with ``batch[i]`` by shared index ``i``
would work by coincidence. It is still the wrong tool to reach for,
because nothing guarantees a WinGLink export covers the same stations
as the raw delivery it came from -- WinGLink can skip a station the
hardware recorded (a failed conversion, an operator excluding a bad
run), and once that happens every EDI *after* the gap silently pairs
with the wrong raw station under plain index alignment. The following
reconstructs exactly that situation from real K2 files -- raw data for
stations 1, 2, 4, 5 (station 3 missing) against an EDI export for
stations 1, 2, 3, 5 (station 4 missing), a deliberately built gap on
each side using genuine file content, not a hypothetical:

.. code-block:: pycon

   >>> import shutil
   >>> from pathlib import Path
   >>> from tempfile import TemporaryDirectory

   >>> with TemporaryDirectory() as tmp:
   ...     raw_dir, edi_dir = Path(tmp) / "raw", Path(tmp) / "edi"
   ...     raw_dir.mkdir(); edi_dir.mkdir()
   ...     for n in [1, 2, 4, 5]:
   ...         _ = shutil.copy(f"data/stratagem/K2/k2-HX/X2HX.{n:03d}", raw_dir / f"X2HX.{n:03d}")
   ...     for n in [1, 2, 3, 5]:
   ...         _ = shutil.copy(f"data/stratagem/K2/k2-edi/Z2HX{n:03d}.edi", edi_dir / f"Z2HX{n:03d}.edi")
   ...     gap_rdr = StratagemRawReader(raw_dir, component="X").fit()
   ...     gap_batch = EDIBatch(edi_dir).fit()
   ...     naive = [
   ...         (gap_batch[i].station, gap_rdr.stations_[i])
   ...         for i in range(len(gap_batch))
   ...     ]
   ...     mapping = gap_rdr.match_to_edis(gap_batch.edi_objects_)
   ...     print(naive)
   ...     print(mapping)
   [('Z2HX001', 'X2HX.001'), ('Z2HX002', 'X2HX.002'), ('Z2HX003', 'X2HX.004'), ('Z2HX005', 'X2HX.005')]
   {0: 0, 1: 1, 3: 3}

Plain index alignment pairs EDI station 3 with raw station 4 -- silently
wrong, and every station after the gap inherits the same one-station
offset. :meth:`~pycsamt.stratagem.io.StratagemRawReader.match_to_edis`
looks up each EDI's station number in the raw file's own numbering
instead of trusting position, so it recovers ``{0: 0, 1: 1, 3: 3}``:
EDI batch index 2 (station 3) is simply absent, because no raw file for
station 3 exists to map it to, rather than being guessed into the
nearest available one. This is exactly the lookup
:class:`~pycsamt.stratagem.qc.QualityController` and
:class:`~pycsamt.stratagem.qc.FrequencyFilter` use internally whenever a
raw reader is supplied alongside an EDI batch, covered in
:doc:`quality_control`.

Sites Interoperability
----------------------

``EDIBatch.edi_objects_`` is a plain list of
:class:`~pycsamt.seg.edi.EDIFile`, which means it is already valid input
anywhere :mod:`pycsamt.emtools` expects EDI-like data.
:func:`~pycsamt.emtools.ensure_sites` is the single normalisation point
the rest of that stack uses, converting a path, a list of ``EDIFile``,
or an existing :class:`~pycsamt.site.base.Sites` into one canonical
``Sites`` container:

.. code-block:: pycon

   >>> from pycsamt.emtools import ensure_sites

   >>> sites = ensure_sites(batch.edi_objects_)
   >>> len(sites)
   87
   >>> sites["Z2HX001"].name
   'Z2HX001'
   >>> sites[0].summary()
   {'name': 'Z2HX001', 'nfreq': 34, 'lat': 0.0, 'lon': 0.0, 'elev': 0.0, 'components': ['Zxx', 'Zxy', 'Zyx', 'Zyy'], 'tipper': False}

The ``lat``/``lon``/``elev`` of ``0.0`` here are the same WinGLink
placeholders from :doc:`concepts` -- ``Sites`` does not invent
coordinates any more than ``EDIBatch`` does, it just gives the batch a
name-indexable, station-aware container that every other
:mod:`pycsamt.emtools` tool already understands. This is the same door
:class:`~pycsamt.stratagem.survey.StratagemSurvey` walks through in the
other direction: its ``edi_dir`` parameter accepts a directory path
(loaded internally exactly the way ``EDIBatch`` does above) or an
already-built ``Sites``/list of ``EDIFile``, run through this same
``ensure_sites`` call, and its
:attr:`~pycsamt.stratagem.survey.StratagemSurvey.sites_` property hands
the pipeline's current state back out the same way at any point --
before or after coordinate injection, static-shift correction, or
noise removal. :doc:`pipeline` picks that up in full; the next page,
:doc:`coordinates`, is what turns the ``0.0, 0.0, 0.0`` placeholders
above into real station positions.
