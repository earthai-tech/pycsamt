.. _user_guide_data:

Loading electromagnetic data
============================

Loading is the boundary between field files and the scientific workflow.
Before correcting, plotting, or inverting a survey, verify that pyCSAMT found
the intended stations, preserved their order, and read a usable frequency grid
for each station. A wrong directory or an unnoticed duplicate can otherwise
change every downstream table while still producing plausible-looking output.

The public entry points are :func:`pycsamt.api.read_edi` for one EDI file and
:func:`pycsamt.api.read_edis` for one or more files. ``read_edis`` returns an
:class:`~pycsamt.api.APISurvey`: a stable survey view with station metadata, a
compact inventory, and access to the underlying EDI collection when an
advanced processing function requires it.

Reading a survey directory
--------------------------

The examples use the bundled WILLY L18PLT line, so their output can be checked
without downloading data. The directory contains 28 EDI files from one AMT
survey line. From the repository root, load it with:

.. code-block:: pycon

   >>> from pycsamt.api import read_edis
   >>> survey = read_edis(
   ...     "data/AMT/WILLY_DATA/L18PLT",
   ...     recursive=False,
   ...     strict=True,
   ...     on_dup="replace",
   ...     progress=False,
   ... )
   >>> type(survey).__name__, survey.n_sites
   ('APISurvey', 28)
   >>> survey.stations[:3]
   ['23-18-001A', '23-18-002U', '23-18-003A']

``recursive=False`` restricts discovery to that directory. Set it to ``True``
only when the supplied directory contains survey-line subdirectories that
should all be included. Keeping the search boundary explicit prevents a
parent field folder from silently combining unrelated lines.

``strict`` controls failure when no usable EDI data can be resolved. Use
``strict=False`` during initial inspection of a mixed field drop, where
skipping unrelated or malformed files is useful. Use ``strict=True`` in a
reproducible workflow so an empty or invalid input fails before processing.
It does not mean that every optional EDI section must be present.

.. important::

   The duplicate policies supported by ``read_edis`` are ``"replace"`` and
   ``"keep"``.
   ``replace`` retains the later occurrence of a repeated station identity;
   ``keep`` retains the earlier occurrence. Neither policy establishes which
   field file is scientifically authoritative. Inspect duplicate files and
   construct an explicit input list when that choice matters.

Inspecting what was loaded
--------------------------

An :class:`~pycsamt.api.APISurvey` is deliberately small. Inspect its station
inventory before passing it downstream:

.. code-block:: pycon

   >>> summary = survey.summary()
   >>> print(summary)
   APIFrame: edi_survey_summary
   kind: edi.summary
   shape: 28 rows x 6 columns
   columns: station, path, n_freq, tipper, spectra, ts
   numeric: 1 columns
   missing: 0.0%
   source: data/AMT/WILLY_DATA/L18PLT

The inventory reports one row per station. ``n_freq`` is the number of
frequencies in the impedance section; ``tipper``, ``spectra``, and ``ts``
indicate whether those optional sections are available. A station can load
successfully without them, so check these flags before calling a method that
depends on a particular response.

For station :math:`i`, the frequency grid is

.. math::
   :label: eq-data-loading-frequency-grid

   \mathbf{f}_i = (f_{i,0}, f_{i,1}, \ldots, f_{i,m_i-1}),

where :math:`f_{i,j}` is frequency in hertz and :math:`m_i` is the value
reported by ``n_freq``. Stations need not have identical grids, but differences
affect frequency matching, station-to-station comparisons, and inversion
preparation. Unexpected variation in ``n_freq`` is therefore a useful early
warning, not merely a formatting detail.

Convert the inventory to pandas only when a library or notebook operation
requires a plain dataframe:

.. code-block:: pycon

   >>> inventory = summary.to_pandas(copy=True)
   >>> inventory.shape
   (28, 6)
   >>> inventory.loc[0, ["station", "n_freq"]].to_dict()
   {'station': '23-18-001A', 'n_freq': 53}

The copied dataframe can be edited without changing the survey view.

Reading one file or an explicit selection
-----------------------------------------

Use :func:`pycsamt.api.read_edi` when a calculation genuinely concerns one
station and the raw :class:`~pycsamt.seg.edi.EDIFile` is the desired object:

.. code-block:: pycon

   >>> from pycsamt.api import read_edi
   >>> edi = read_edi("data/AMT/WILLY_DATA/L18PLT/18-001A.edi")
   >>> edi.station
   '23-18-001A'
   >>> len(edi.Z.freq)
   53

For reusable code, ``read_edis`` can also load one file and preserve the same
survey interface used for a full line. To select several known stations, pass
their paths explicitly:

.. code-block:: pycon

   >>> from pathlib import Path
   >>> line = Path("data/AMT/WILLY_DATA/L18PLT")
   >>> selected = read_edis(
   ...     [line / "18-001A.edi", line / "18-002U.edi"],
   ...     recursive=False,
   ...     strict=True,
   ...     progress=False,
   ... )
   >>> selected.n_sites
   2

An explicit list is preferable when scratch exports, corrected copies, or
neighboring lines share a directory. It records exactly which observations
entered the analysis.

Choosing the right representation
---------------------------------

Keep the public survey view for inspection, reporting, and most user-facing
workflows. Its ``collection`` attribute exposes the parsed
:class:`~pycsamt.seg.collection.EDICollection` for lower-level functions:

.. code-block:: pycon

   >>> type(survey.collection).__name__
   'EDICollection'
   >>> survey.to_collection() is survey.collection
   True

Pass ``survey.collection`` when a documented function explicitly requests an
EDI collection or a site container. For the processing-oriented ``Sites``
representation, use the public normalization function described next.

Normalizing inputs
--------------------

:func:`pycsamt.emtools.ensure_sites` is the public normalization boundary used
by the electromagnetic processing tools. It accepts paths, an ``EDIFile`` or
``EDICollection``, a ``Site`` or ``Sites`` instance, or an iterable of EDI-like
objects, and always returns :class:`pycsamt.site.base.Sites`. This is useful
when reusable processing code should accept several input forms without
implementing its own coercion rules.

.. code-block:: pycon

   >>> from itertools import islice
   >>> from pycsamt.emtools import ensure_sites
   >>> sites = ensure_sites(
   ...     survey.collection,
   ...     recursive=False,
   ...     order_by="input",
   ...     strict=True,
   ... )
   >>> type(sites).__name__, len(sites)
   ('Sites', 28)
   >>> [site.name for site in islice(sites, 3)]
   ['18-001A', '18-002U', '18-003A']
   >>> first = sites[0].summary()
   >>> first["name"], first["nfreq"], first["tipper"]
   ('18-001A', 53, False)
   >>> first["components"]
   ['Zxx', 'Zxy', 'Zyx', 'Zyy']

The shorter names in ``Sites`` come from the site-normalization layer, whereas
the public survey view reports the EDI header identities. Keep the chosen
representation stable within a table or workflow rather than joining outputs
on an assumed naming convention.

``order_by`` controls the returned station sequence. ``"input"`` preserves
discovery order; ``"station"``, ``"latitude"``, and ``"longitude"`` sort by
the corresponding property; ``"chainage"`` follows distance along the
profile. ``"auto"`` uses coordinate-derived chainage only when the geometry
looks like a reliable single line. Passing ``None`` uses the package-wide
ordering configuration. Choose an explicit value when reproducible row order
is important.

Unlike ``read_edis``, ``ensure_sites`` supports ``"replace"``,
``"keep_first"``, ``"keep_last"``, and ``"raise"`` for ``on_dup``. Use
``"raise"`` in controlled processing when a repeated normalized site name
should stop the workflow. The policy is applied at the ``Sites`` boundary and
cannot recover duplicates that an earlier loader has already collapsed.

``strict=True`` raises when no input can be normalized; it does not reject a
station merely because an optional component is absent. A positive ``verbose``
value reports coercion and duplicate warnings. Import the function from
``pycsamt.emtools`` as shown above—the private implementation module is not a
user-facing import path.

Preparing data for downstream science
-------------------------------------

Before quality control or inversion, confirm at least four properties:

* the station count agrees with the field manifest;
* station identifiers are unique and ordered as expected along the line;
* coordinates and elevation use the expected reference and units;
* frequency coverage and required impedance or tipper components are present.

Station order matters because downstream rows should retain a stable mapping
to physical locations. If the loaded sequence is

.. math::
   :label: eq-data-loading-station-sequence

   \mathcal{S} = (S_0, S_1, \ldots, S_{n-1}),

then row :math:`i` in a station-level QC or correction table should continue
to describe :math:`S_i`. Sort or subset deliberately and retain the resulting
station identifiers with every exported table.

.. warning::

   Successful parsing does not prove scientific validity. EDI files can have
   swapped coordinates, inconsistent units, sparse frequency coverage, or
   physically implausible transfer functions while remaining syntactically
   valid. Loading should always be followed by inventory and quality-control
   checks.

For an end-to-end inspection with the same dataset, figures, dataframe output,
and command-line equivalents, continue with
:doc:`../tutorials/read_edi_survey`. Quality-control workflows are introduced
in :doc:`emtools/index`.

Common loading problems
-----------------------

No stations are returned
   Confirm that the path exists and contains ``.edi`` files. Enable
   ``recursive=True`` only if files are nested below the supplied directory.
   During diagnosis, use ``strict=False`` and a nonzero ``verbose`` value;
   restore ``strict=True`` for the reproducible pipeline.

Fewer stations are returned than expected
   Compare ``survey.stations`` and ``survey.paths`` with the field manifest.
   Repeated station identities may have been collapsed by the selected
   duplicate policy, or individual files may have failed parsing.

Frequency counts differ between stations
   Inspect the EDI impedance sections and decide whether the difference is an
   expected acquisition condition or missing data. Do not assume that arrays
   can be compared element by element until their frequency grids have been
   aligned.

Tipper, spectra, or time-series flags are false
   These sections are optional. Their absence is acceptable for methods that
   use only impedance, but any workflow requiring them must exclude or repair
   the affected stations explicitly.

See also
--------

* :doc:`../getting_started/quickstart` — the shortest complete first workflow.
* :doc:`../tutorials/read_edi_survey` — detailed EDI survey inspection.
* :doc:`../api/api` — public loading functions and survey-view API.
* :doc:`workflow_overview` — choose the route that matches your input and goal.
