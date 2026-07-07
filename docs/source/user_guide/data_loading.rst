.. _user_guide_data:

Data Loading
============

Data loading is the first boundary in a pyCSAMT workflow. The goal is to
turn files, directories, parsed EDI objects, or existing site containers into
one predictable object that the science tools can use.

For processing workflows, the canonical loader is
:func:`pycsamt.emtools._core.ensure_sites`. It returns a
:class:`pycsamt.site.base.Sites` container. Most functions in
``pycsamt.emtools`` call this same helper internally, so using it in notebooks,
scripts, agents, and tests gives you the same behavior as the package itself.

Use :func:`pycsamt.api.read_edis` when you want a public API survey view for
notebooks, reports, and API-style tables. The survey view still wraps the same
EDI-level data, but ``ensure_sites`` is the direct science/workflow input.

Loading EDI Directories
-----------------------

The common case is a directory containing one or more ``.edi`` files. The
example below uses one of the repository's example survey lines.

.. code-block:: python
   :linenos:

   from pycsamt.emtools import ensure_sites

   edi_dir = "data/AMT/WILLY_DATA/L18PLT"
   sites = ensure_sites(
       edi_dir,
       recursive=True,
       strict=False,
       on_dup="replace",
       verbose=0,
   )

   print(f"Loaded {len(sites)} stations")
   print([site.name for site in list(sites)[:3]])

Line 1 imports the loader used by the processing modules. Line 3 points at a
directory rather than at a single file. Lines 4-10 normalize the input into a
``Sites`` object. ``recursive=True`` searches subdirectories, ``strict=False``
allows malformed or unrelated files to be skipped, and ``on_dup="replace"``
keeps the last station encountered when duplicated station names are found.

The returned ``sites`` object is iterable. Each item is a
:class:`pycsamt.site.base.Site` wrapper around one parsed EDI object.

Loading Explicit Files
----------------------

When the desired stations are known, pass a list of files. This avoids
accidentally loading scratch files or neighboring survey lines.

.. code-block:: python
   :linenos:

   from pathlib import Path

   from pycsamt.emtools import ensure_sites

   line = Path("data/AMT/WILLY_DATA/L18PLT")
   selected = [
       line / "18-001A.edi",
       line / "18-002U.edi",
       line / "18-003A.edi",
   ]

   sites = ensure_sites(selected, recursive=False, strict=True)

Line 6 creates an explicit list of path-like objects. Line 13 disables
directory recursion because the input is already a file list. ``strict=True``
is useful here: if none of the requested inputs can be resolved into valid
stations, loading fails immediately instead of returning an empty container.

Normalizing Existing Objects
----------------------------

``ensure_sites`` is intentionally permissive at the boundary. It accepts:

- a directory, a single EDI path, a glob-like path, or a list of paths;
- an ``EDIFile`` or an ``EDICollection``;
- a single ``Site``;
- an existing ``Sites`` object;
- an iterable containing EDI-like objects.

This lets you normalize once and then write the rest of the workflow against
``Sites``.

.. code-block:: python
   :linenos:

   from pycsamt.api import read_edis
   from pycsamt.emtools import ensure_sites

   survey = read_edis("data/AMT/WILLY_DATA/L18PLT", recursive=True)
   sites = ensure_sites(survey.collection, recursive=False)

   same_sites = ensure_sites(sites)
   assert same_sites is sites

Lines 4-5 show the bridge from the public API reader to the science container.
``read_edis`` returns an ``APISurvey`` view, while ``survey.collection`` is the
underlying EDI collection. Line 7 is a no-op normalization: if the input is
already a ``Sites`` instance, the same object is returned.

Inspecting Loaded Stations
--------------------------

Before applying QC, static-shift correction, dimensionality analysis, or
inversion preparation, inspect what actually loaded.

.. code-block:: python
   :linenos:

   first = sites[0]
   same = sites[first.name]
   maybe = sites.get("missing_station")

   print(first.summary())
   print(first.has_component("Zxy"))
   print(first.has_component("tipper"))

   for site in sites:
       row = site.summary()
       print(row["name"], row["nfreq"], row["lat"], row["lon"])

Line 1 retrieves by zero-based position. Line 2 retrieves by station name;
name lookup is case-insensitive and raises ``KeyError`` if the station is not
present. Line 3 uses safe lookup and returns ``None`` when the station is
missing.

``Site.summary()`` returns a dictionary with the station name, number of
frequencies, coordinates, present impedance components, and whether tipper data
are available. ``has_component("Zxy")`` and ``has_component("tipper")`` are
quick checks before running tools that require those arrays.

Selecting Stations
------------------

Use ``Sites.select`` to keep a named subset or stations matching a predicate.
The method returns a new ``Sites`` container and does not mutate the original.

.. code-block:: python
   :linenos:

   named = sites.select(names=["18-001A", "18-002U"])

   with_tipper = sites.select(
       predicate=lambda site: site.has_component("tipper")
   )

   dense = sites.select(
       predicate=lambda site: site.summary()["nfreq"] >= 20
   )

Line 1 keeps specific station names. Lines 3-5 keep only stations with tipper
data. Lines 7-9 keep stations with at least 20 frequency samples. Predicates
receive a ``Site`` object, so they can use station names, coordinates,
frequency counts, and component checks.

Duplicate Stations
------------------

Field directories sometimes contain repeated station names: an original file,
a corrected export, and a manually edited copy may all describe the same
station. Choose the policy explicitly.

.. code-block:: python
   :linenos:

   latest = ensure_sites("data/edis", on_dup="replace")
   first = ensure_sites("data/edis", on_dup="keep_first")
   last = ensure_sites("data/edis", on_dup="keep_last")

   checked = ensure_sites(
       "data/edis",
       on_dup="raise",
       strict=True,
   )

``replace`` is the default and keeps the latest duplicate found by the lower
level collection loader. ``keep_first`` keeps the first station name
encountered. ``keep_last`` keeps the last station name encountered. ``raise``
is the safest production setting when duplicates should be treated as a data
management error.

For path inputs, pyCSAMT maps these policies onto the EDI collection parser
and performs any extra duplicate check after parsing. For already parsed
objects, the same policy is applied while building the ``Sites`` container.

Strict Versus Exploratory Loading
---------------------------------

Use ``strict=False`` while exploring mixed folders. Use ``strict=True`` when a
pipeline should fail early if no valid stations can be loaded.

.. code-block:: python
   :linenos:

   exploratory = ensure_sites(
       "field_drop/raw_export",
       recursive=True,
       strict=False,
       verbose=1,
   )

   production = ensure_sites(
       "field_drop/approved_edi",
       recursive=True,
       strict=True,
       on_dup="raise",
       verbose=1,
   )

In the exploratory call, unrelated files or broken EDI files can be skipped so
you can inspect what is salvageable. In the production call, an empty result
raises ``ValueError`` and duplicated station names raise before downstream
processing can silently use the wrong station.

Using Loaded Data In Processing
-------------------------------

Once data are normalized to ``Sites``, pass the same object to QC,
static-shift, plotting, frequency tools, and inversion preparation.

.. code-block:: python
   :linenos:

   from pycsamt.emtools.qc import build_qc_table, qc_flags
   from pycsamt.emtools.ss import estimate_ss_ama

   qc_table = build_qc_table(sites)
   flags = qc_flags(sites, min_frac_ok=0.6, min_snr_med=2.0)
   ss_table = estimate_ss_ama(
       sites,
       half_window=3,
       sort_by="lon",
       verbose=0,
   )

   print(qc_table.head())
   print(flags.head())
   print(ss_table)

Lines 4-11 reuse the same loaded container. This is the expected pattern:
normalize at the boundary, then pass ``sites`` through the workflow. You do
not need to reload the EDI directory for each processing step.

Loading Through The Agent Layer
-------------------------------

The loader agent is useful when a workflow should return an
``AgentResult`` with status, summary text, and a quality table.

.. code-block:: python
   :linenos:

   from pycsamt.agents.loader import MTLoaderAgent

   agent = MTLoaderAgent(recursive=True, on_dup="replace")
   result = agent.execute({
       "path": "data/AMT/WILLY_DATA/L18PLT",
   })

   if result.status != "success":
       raise RuntimeError(result.summary)

   print(result["n_stations"])
   print(result["quality_table"])

Lines 3-6 delegate loading to :class:`pycsamt.agents.loader.MTLoaderAgent`.
Internally, the agent calls ``ensure_sites`` and then performs a lightweight
quality scan. Use this route for conversational workflows, dashboards, and
pipeline orchestration. Use ``ensure_sites`` directly for lower-level scripts
and notebooks.

Writing And Reloading EDI Files
-------------------------------

``Sites.write`` writes one EDI file per station. This is useful after
renaming, topography alignment, frequency slicing, denoising, or correction
steps that return a new ``Sites`` object.

.. code-block:: python
   :linenos:

   from pathlib import Path

   outdir = Path("results/loaded_line")
   written = sites.write(
       outdir,
       template="{station}.edi",
       exist_ok=True,
   )

   reloaded = ensure_sites(outdir, recursive=False, strict=True)
   print([path.name for path in written])
   print(len(reloaded))

Lines 4-8 create or reuse the output directory and write station files with
names derived from normalized station IDs. Line 10 reloads the written files
as a quick integrity check before handing them to another tool.

Troubleshooting
---------------

``ensure_sites: got None``
    The loader received ``None``. Check that your script is passing a real
    path, EDI object, ``Site``, ``Sites``, or iterable of those objects.

``ensure_sites(strict=True): no sites were resolved``
    No valid EDI-like items were found. Confirm the path, enable
    ``recursive=True`` for nested folders, and inspect whether files contain
    impedance data.

``Duplicate site ... encountered with on_dup='raise'``
    Two or more inputs have the same station identity. Remove the duplicate
    files, rename stations intentionally, or use ``keep_first``/``keep_last``
    during exploratory work.

An empty or tiny station count
    Make sure you are pointing to the survey-line directory, not a parent
    folder with unrelated files. Use an explicit file list when only selected
    stations should be loaded.

Missing ``Z`` or tipper components
    Loading can succeed even when optional components are absent. Check
    ``site.summary()`` and ``site.has_component(...)`` before running a method
    that requires specific transfer-function components.

See Also
--------

:doc:`../getting_started/quickstart`
    End-to-end first workflow using real v2 entry points.
:doc:`../tutorials/read_edi_survey`
    Public ``read_edis`` tutorial and API survey view examples.
:doc:`../api/api`
    Public API helpers, including ``read_edis``.
:mod:`pycsamt.emtools._core`
    Canonical ``ensure_sites`` loader used by processing tools.
:mod:`pycsamt.site.base`
    ``Site`` and ``Sites`` container implementation.
