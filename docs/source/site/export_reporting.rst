.. _site_export_reporting:

Export And Reporting
====================

The export and reporting tools are the last mile of a site workflow. Export
functions write cleaned EDI-like objects to disk or package them into an
archive. Report classes summarize one station or a whole survey for terminal
inspection, notebooks, QA logs, and downstream tables.

Use this page when you need to:

* write one corrected station to a known file path;
* batch-write a selected survey with deterministic file names;
* create a manifest CSV for a delivery directory or zip archive;
* package cleaned sites into a zip file;
* generate single-station and survey-level report tables;
* understand how terminal reports, dictionaries, and DataFrames differ.

If the goal is not only to write existing sites, but also to rotate,
frequency-filter, recompute resistivity/phase, and regenerate EDI files from
another software source, use :doc:`recompute`. The export helpers documented
here are simpler writer utilities for objects that have already been prepared.

Tool Map
--------

.. list-table::
   :header-rows: 1
   :widths: 24 32 44

   * - Tool
     - Module
     - Main purpose
   * - :func:`pycsamt.site.export.write_site`
     - :mod:`pycsamt.site.export`
     - Write one EDI-like object to a target path.
   * - :func:`pycsamt.site.export.write_sites`
     - :mod:`pycsamt.site.export`
     - Batch-write many sites to an output directory with templated file
       names and an optional manifest CSV.
   * - :func:`pycsamt.site.export.pack_zip`
     - :mod:`pycsamt.site.export`
     - Stage sites in a temporary directory, compress them, and optionally
       write a manifest.
   * - :class:`pycsamt.site.report.SiteReport`
     - :mod:`pycsamt.site.report`
     - Compute and display statistics for one
       :class:`pycsamt.site.base.Site`.
   * - :class:`pycsamt.site.report.SitesReport`
     - :mod:`pycsamt.site.report`
     - Compute and display survey-level and per-station statistics for a
       collection.

Export Workflow
---------------

Export normally happens after loading, selecting, editing, and quality checks.

.. code-block:: python
   :linenos:

   from pycsamt.seg.collection import EDICollection
   from pycsamt.site.export import write_sites
   from pycsamt.site.selection import drop_empty, keep_finite_z

   sites = EDICollection.from_sources("data/edi")

   sites = drop_empty(sites)
   sites = keep_finite_z(sites)

   paths = write_sites(
       sites,
       outdir="deliveries/clean_edi",
       template="{index:03d}_{station}.edi",
       manifest_csv="deliveries/clean_edi_manifest.csv",
   )

   print(paths)

The export functions accept:

* a :class:`pycsamt.site.base.Sites` object;
* an ``EDICollection``;
* a list or other iterable of EDI-like objects;
* a single EDI-like object.

An object is EDI-like for export when it supports at least one common writer
method: ``write(new_edifn=...)``, ``write(path)``, ``to_file(path)``, or
``save(path)``.

Writing One Site
----------------

Use :func:`pycsamt.site.export.write_site` when the output path is already
known and you are only writing one object.

.. code-block:: python
   :linenos:

   from pycsamt.seg.edi import EDIFile
   from pycsamt.site.export import write_site

   edi = EDIFile("data/edi/S01.edi")

   out = write_site(
       edi,
       "deliveries/single/S01_corrected.edi",
   )

   print(out)

The parent directory is created automatically. If the target file already
exists, overwrite behavior depends on the underlying EDI writer. For explicit
collision control in batch exports, use :func:`pycsamt.site.export.write_sites`
with ``exist_ok``.

Batch Writing Sites
-------------------

Use :func:`pycsamt.site.export.write_sites` for a selected survey or station
group.

.. code-block:: python
   :linenos:

   from pycsamt.site.export import write_sites

   written = write_sites(
       sites,
       outdir="deliveries/line_01",
       template="{index:03d}_{station}",
       exist_ok=False,
       manifest_csv="deliveries/line_01_manifest.csv",
   )

   for path in written:
       print(path.name)

If the rendered file name does not end with ``.edi``, the extension is added
automatically. In the example above, a station named ``S01`` at index ``0`` is
written as ``000_S01.edi``.

Filename Templates
------------------

Filename templates are rendered from a context collected from each EDI header
and its input order.

.. list-table::
   :header-rows: 1
   :widths: 24 28 48

   * - Template key
     - Source
     - Notes
   * - ``{station}``
     - Station name helper.
     - Uses the same station-name resolution as the site utilities.
   * - ``{index}``
     - Input iteration order.
     - Zero-based integer. Supports Python format specifiers such as
       ``{index:03d}``.
   * - ``{lat}``
     - EDI ``HEAD`` coordinate.
     - Missing values become ``NaN``.
   * - ``{lon}``
     - EDI ``HEAD`` coordinate.
     - Also accepts ``long`` and ``longitude`` internally.
   * - ``{elev}``
     - EDI ``HEAD`` coordinate.
     - Missing values become ``NaN``.
   * - ``{chainage}``
     - EDI object attribute.
     - Reads ``edi.chainage`` when present.

Examples:

.. code-block:: python
   :linenos:

   write_sites(
       sites,
       "out/by_station",
       template="{station}.edi",
   )

   write_sites(
       sites,
       "out/by_order",
       template="{index:03d}_{station}",
   )

   write_sites(
       sites,
       "out/by_line_position",
       template="{index:03d}_{station}_{chainage:.0f}m",
   )

Unknown template keys are rendered as empty strings. Keep templates simple and
avoid path separators inside station names when building files for external
delivery.

Collision Policy
----------------

By default, :func:`pycsamt.site.export.write_sites` refuses to overwrite an
existing file.

.. code-block:: python
   :linenos:

   from pycsamt.site.export import write_sites

   # Raises FileExistsError if any rendered file already exists.
   write_sites(
       sites,
       "deliveries/line_01",
       template="{station}.edi",
       exist_ok=False,
   )

Pass ``exist_ok=True`` when an export directory is intentionally being
refreshed.

.. code-block:: python
   :linenos:

   write_sites(
       sites,
       "deliveries/line_01",
       template="{station}.edi",
       exist_ok=True,
   )

Manifest CSV
------------

Both :func:`pycsamt.site.export.write_sites` and
:func:`pycsamt.site.export.pack_zip` can write a manifest CSV. The manifest is
a compact audit table for exported stations.

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Column
     - Meaning
   * - ``index``
     - Zero-based input order used during export.
   * - ``station``
     - Resolved station name.
   * - ``lat``
     - Latitude from the EDI header, or ``NaN``.
   * - ``lon``
     - Longitude from the EDI header, or ``NaN``.
   * - ``elev``
     - Elevation from the EDI header, or ``NaN``.
   * - ``chainage``
     - ``edi.chainage`` when present, otherwise ``NaN``.
   * - ``filename``
     - File name written in the output directory or archive.
   * - ``path``
     - Full exported file path for directory exports, or the zip archive path
       for archive exports.

.. code-block:: python
   :linenos:

   import pandas as pd

   from pycsamt.site.export import write_sites

   write_sites(
       sites,
       "deliveries/clean_edi",
       template="{index:03d}_{station}",
       manifest_csv="deliveries/clean_edi_manifest.csv",
   )

   manifest = pd.read_csv("deliveries/clean_edi_manifest.csv")
   print(manifest[["index", "station", "filename"]])

The manifest is only written when at least one row exists.

Zip Packaging
-------------

Use :func:`pycsamt.site.export.pack_zip` when you need a single archive for a
delivery, web upload, or reproducible artifact.

.. code-block:: python
   :linenos:

   from pycsamt.site.export import pack_zip

   archive = pack_zip(
       sites,
       out_zip="deliveries/line_01_clean.zip",
       template="{station}.edi",
       manifest_csv="deliveries/line_01_clean_manifest.csv",
   )

   print(archive)

The function writes each site to a temporary directory, then stores those
files inside the archive with ``zipfile.ZIP_DEFLATED`` compression. The
original EDI sources are not deleted or modified by packaging.

You can inspect the archive names with the standard library:

.. code-block:: python
   :linenos:

   from zipfile import ZipFile

   with ZipFile("deliveries/line_01_clean.zip", "r") as zf:
       print(sorted(zf.namelist()))

Reporting Workflow
------------------

Reports are read-only summaries. They do not modify site data.

.. code-block:: python
   :linenos:

   from pycsamt.site.base import Sites
   from pycsamt.site.report import SitesReport

   sites = Sites.from_any("data/edi")

   report = SitesReport(sites)
   report.report(top=10)

   table = report.to_dataframe(api=False)
   print(table.head())

If the optional :mod:`rich` package is installed, terminal reports use rich
panels and tables. Without ``rich``, pyCSAMT falls back to plain-text output.

Single-Site Reports
-------------------

:class:`pycsamt.site.report.SiteReport` computes statistics for one
:class:`pycsamt.site.base.Site`-like object.

.. code-block:: python
   :linenos:

   from pycsamt.site.base import Site
   from pycsamt.site.report import SiteReport

   site = Site(edi)
   report = SiteReport(site)

   print(report.summary())
   report.report()

   info = report.to_dict()
   z_frame = report.to_dataframe("z", api=False)

The dictionary returned by :meth:`pycsamt.site.report.SiteReport.to_dict`
contains:

* ``name``, ``lat``, ``lon``, and ``elev``;
* ``nfreq``, ``freq_min``, ``freq_max``, ``per_min``, and ``per_max``;
* component availability for ``Zxx``, ``Zxy``, ``Zyx``, and ``Zyy``;
* ``has_tipper``;
* mean and standard deviation for apparent resistivity and phase summaries of
  the ``xy`` and ``yx`` components;
* per-component finite-data quality fractions.

The DataFrame returned by
:meth:`pycsamt.site.report.SiteReport.to_dataframe` delegates to the
underlying site ``to_dataframe`` method. Use the ``kind`` argument to choose
the station table representation supported by :class:`pycsamt.site.base.Site`.

Survey Reports
--------------

:class:`pycsamt.site.report.SitesReport` computes one station record per site
and a survey-level summary.

.. code-block:: python
   :linenos:

   from pycsamt.site.report import SitesReport

   report = SitesReport(sites)

   print(report.summary())
   report.report(top=20)

   records = report.to_dict()
   frame = report.to_dataframe(api=False)

   print(len(records))
   print(frame.columns)

The survey DataFrame contains one row per station with columns:

.. code-block:: text
   :linenos:

   station
   lat
   lon
   elev
   nfreq
   freq_min
   freq_max
   has_Zxx
   has_Zxy
   has_Zyx
   has_Zyy
   has_tipper
   rho_xy
   rho_xy_std
   rho_yx
   rho_yx_std
   phi_xy
   phi_xy_std
   phi_yx
   phi_yx_std

Use ``top`` to limit only the terminal display. It does not limit
``to_dict`` or ``to_dataframe``.

API View Wrapping
-----------------

Report DataFrame methods accept ``api``:

``api=False``
   Return a plain pandas DataFrame.

``api=True``
   Return a pyCSAMT API view wrapper around the DataFrame.

``api=None``
   Defer to the global API-view configuration.

.. code-block:: python
   :linenos:

   from pycsamt.site.report import SiteReport, SitesReport

   one_station = SiteReport(site).to_dataframe("resphase", api=False)
   survey_view = SitesReport(sites).to_dataframe(api=True)

   print(type(one_station))
   print(survey_view.kind)

Use plain DataFrames for local analysis and API views when returning data from
public-facing pyCSAMT APIs.

End-To-End Delivery Example
---------------------------

The following pattern creates a cleaned station export, a manifest, an archive,
and a survey report table.

.. code-block:: python
   :linenos:

   from pathlib import Path

   from pycsamt.seg.collection import EDICollection
   from pycsamt.site.edit import select_freq_all
   from pycsamt.site.export import pack_zip, write_sites
   from pycsamt.site.report import SitesReport
   from pycsamt.site.selection import by_freq, drop_empty, keep_finite_z

   root = Path("deliveries/line_01")
   root.mkdir(parents=True, exist_ok=True)

   sites = EDICollection.from_sources("data/edi")
   sites = drop_empty(sites)
   sites = keep_finite_z(sites)
   sites = by_freq(sites, fmin=0.1, fmax=100.0)
   sites = select_freq_all(
       sites,
       fmin=0.1,
       fmax=100.0,
       inplace=False,
   )

   write_sites(
       sites,
       root / "edi",
       template="{index:03d}_{station}",
       manifest_csv=root / "manifest.csv",
   )

   pack_zip(
       sites,
       root / "edi.zip",
       template="{index:03d}_{station}.edi",
       manifest_csv=root / "zip_manifest.csv",
   )

   report = SitesReport(sites)
   report.report(top=15)
   summary = report.to_dataframe(api=False)
   summary.to_csv(root / "site_report.csv", index=False)

Common Mistakes
---------------

Assuming :func:`pycsamt.site.export.write_site` enforces overwrite policy
   The single-site writer delegates to the backend writer. Use
   :func:`pycsamt.site.export.write_sites` with ``exist_ok=False`` when you
   need explicit collision protection.

Forgetting that ``{index}`` is zero-based
   The index in filenames and manifests is the input iteration position,
   starting from ``0``.

Using :func:`pycsamt.site.export.pack_zip` when you also need individual files
   The zip function stages files in a temporary directory. Use
   :func:`pycsamt.site.export.write_sites` first when the directory tree itself
   is part of the delivery.

Expecting reports to clean data
   Reports summarize the current objects. Run selection and editing before
   generating reports.

Confusing ``top`` with filtering
   ``SitesReport.report(top=10)`` limits terminal display only. It does not
   remove stations from the report object.

Next Pages
----------

* :doc:`selection` for choosing which stations to export or report;
* :doc:`editing` for changing station data before writing;
* :doc:`recompute` for regenerating EDI files from raw or foreign EDI inputs;
* :doc:`computed_diagnostics` for analysis tables that complement site
  reports;
* :doc:`containers` for the ``Site`` and ``Sites`` wrappers used by reporting.
