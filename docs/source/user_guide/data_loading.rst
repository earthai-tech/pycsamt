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
advanced processing function requires it. Both are deliberately EDI-only --
narrow, predictable, and unchanged below. Ground stations delivered as
EMTF-XML instead of EDI, and airborne surveys with no EDI-shaped impedance at
all, load through two further, equally stable boundaries introduced later on
this page: :func:`pycsamt.emtools.ensure_sites` (:ref:`data-loading-xml`) and
:func:`pycsamt.airborne.site.ensure_asites` (:ref:`data-loading-airborne`).
Nothing about ``read_edi``/``read_edis``/``APISurvey`` changes because of
either.

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
------------------

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

.. _data-loading-xml:

Loading EMTF-XML the same way
-----------------------------

``ensure_sites`` is not an EDI-only boundary. A directory, path, or object
built from the modern EMTF-XML transfer-function format normalizes into the
exact same :class:`~pycsamt.site.base.Sites` container, with no separate
function to learn and no change to what a returned :class:`Site` can do.
The examples below use the real Gabbs Valley station ``gv100`` (public-domain
USGS data, `doi:10.5066/P9GZ9Z56 <https://doi.org/10.5066/P9GZ9Z56>`__),
committed both as EDI (``data/gv_data/gv_final_edi/``) and as the three-
station EMTF-XML subset pyCSAMT converted it to (``data/gv_data/xml/``):

.. code-block:: pycon

   >>> xml_sites = ensure_sites("data/gv_data/xml")
   >>> type(xml_sites).__name__, len(xml_sites)
   ('Sites', 3)
   >>> sorted(s.name for s in xml_sites)
   ['gv100', 'gv101', 'gv102']
   >>> gv100 = xml_sites[0]
   >>> gv100.backend
   'xml'
   >>> xml_summary = gv100.summary()
   >>> xml_summary["name"], xml_summary["nfreq"], xml_summary["tipper"]
   ('gv100', 48, True)
   >>> xml_summary["components"]
   ['Zxx', 'Zxy', 'Zyx', 'Zyy']

:attr:`~pycsamt.site.base.Site.backend` is the only visible difference --
``gv100``'s impedance (``Zxx``...``Zyy``) and tipper are both filled exactly
as they would be from an EDI file, because every numeric accessor on
:class:`~pycsamt.site.base.Site` (``z``, ``tipper``, ``rho``, ``phase``,
:meth:`~pycsamt.site.base.Site.to_dataframe`) works identically regardless of
which format is native. EMTF-XML also carries richer typed metadata than an
EDI header can (site layout, provenance, processing history, per-component
variance) through :attr:`~pycsamt.site.base.Site.site_meta`,
:attr:`~pycsamt.site.base.Site.provenance`, and related properties -- the
full lazy-materialization and metadata story, including how a
:class:`Site` converts between the two representations on first access, is
covered in :doc:`site/containers`; see :doc:`metadata/index` for what each
metadata object holds.

Writing back out is symmetric, not a one-way import. A :class:`Sites`
collection converts to either format regardless of which one is native,
because pyCSAMT round-trips EDI and EMTF-XML rather than treating one as
authoritative:

.. code-block:: pycon

   >>> import tempfile
   >>> from pathlib import Path
   >>> with tempfile.TemporaryDirectory() as tmp:
   ...     xml_paths = xml_sites.write_xml(tmp)
   ...     edi_paths = xml_sites.write(tmp)
   ...     xml_names = sorted(Path(p).name for p in xml_paths)
   ...     edi_names = sorted(Path(p).name for p in edi_paths)
   >>> xml_names
   ['gv100.xml', 'gv101.xml', 'gv102.xml']
   >>> edi_names
   ['gv100.edi', 'gv101.edi', 'gv102.edi']

``xml_sites`` above is XML-native, yet :meth:`~pycsamt.site.base.Sites.write`
still produces genuine, re-readable EDI files, not placeholders -- the same
:func:`~pycsamt.emtf.converters.edi.edi_to_emtf`/``emtf_to_edi`` conversion
:doc:`site/containers` documents in full runs automatically underneath both
calls. Nothing downstream needs to know or care which file format a survey
originally arrived in.

.. _data-loading-airborne:

Loading airborne surveys
------------------------

Airborne EM (ZTEM, AFMAG, MobileMT) does not fit the ground-station model at
all -- a genuine airborne measurement has no fixed electric-field channel, so
there is no impedance to bridge into an EDI-shaped :class:`Site`, and forcing
one would mean either bending ``Site`` until it stops meaning "ground EDI
station" or fabricating data that was never actually measured. pyCSAMT keeps
these as a parallel, EMTF-XML-only boundary instead:
:func:`pycsamt.airborne.site.ensure_asites`, returning
:class:`~pycsamt.airborne.site.AirborneSites` rather than ``Sites``:

.. code-block:: pycon

   >>> from pycsamt.airborne.site import ensure_asites
   >>> asites = ensure_asites("data/ZTEM/forrestania_wa")
   >>> type(asites).__name__, len(asites), asites.technologies
   ('AirborneSites', 15, ('ztem',))
   >>> station = asites[0]
   >>> station.name, station.technology
   ('FO_001', 'ztem')
   >>> station.z is None
   True
   >>> station.tipper.shape
   (6, 1, 2)

``station.z`` is ``None`` -- not a loading failure, but an honest
reflection of what a ZTEM system actually records: a tipper relating the
airborne vertical field to a fixed-ground horizontal reference, never an
electric field. Reading synthetic sample surveys for all three airborne
technologies, the technology-per-station auto-detection ``asites.technologies``
performs above, and the structural quality-control layer this boundary feeds
into are covered in full in :doc:`airborne/index`, starting with
:doc:`airborne/site`.

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
   Confirm that the path exists and contains files of the expected kind for
   the boundary in use: ``.edi`` for ``read_edis``, ``.edi``/``.xml`` for
   ``ensure_sites``, ``.xml`` for ``ensure_asites``. Enable
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

An ``AirborneSite``'s ``z`` is always ``None``
   This is expected, not a loading failure -- see :ref:`data-loading-airborne`.
   ZTEM, AFMAG, and MobileMT have no electric-field channel to build an
   impedance from; check ``has_component("tipper")``/``"admittance"`` instead
   of assuming ``z`` will ever populate.

See also
--------

* :doc:`../getting_started/quickstart` — the shortest complete first workflow.
* :doc:`../tutorials/read_edi_survey` — detailed EDI survey inspection.
* :doc:`site/containers` — full ``Site``/``Sites`` EDI/EMTF-XML dual-backend
  reference, including lazy materialization and the rich metadata properties.
* :doc:`emtf/xml` — the EMTF-XML format itself, independent of ``Site``.
* :doc:`airborne/index` — the airborne data model, starting with
  :doc:`airborne/site` for reading real ZTEM/AFMAG/MobileMT sample surveys.
* :doc:`../api/api` — public loading functions and survey-view API.
* :doc:`../api/airborne` — the complete airborne callable reference.
* :doc:`workflow_overview` — choose the route that matches your input and goal.
