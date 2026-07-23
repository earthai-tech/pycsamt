.. _getting-started-first-survey:

First Survey
============

This page walks through the first practical pyCSAMT workflow: load a field
survey, inspect the station inventory, check parser errors, build a quality
table, save a first quality-control figure, and run a small reproducible
pipeline.

The examples use EDI files because EDI is the normal exchange format for
MT/AMT/CSAMT transfer functions.  If your first dataset is a Zonge AVG,
J-format, TEM, or intermediate inversion file, start with
:doc:`data_formats` and then return here once you have a collection of EDI
sites or a pyCSAMT site collection.

What This First Workflow Should Answer
--------------------------------------

A first survey pass is not an inversion and it is not a final interpretation.
It is the short diagnostic pass that tells you whether the dataset is ready
for deeper processing.

By the end of this page you should know:

* how many stations were loaded;
* which stations have coordinates and tipper data;
* whether any EDI files failed to parse;
* whether the station/frequency coverage is usable;
* where pyCSAMT writes first-pass outputs;
* which documentation page to open next.

Recommended Project Layout
--------------------------

Use a small project directory before you begin.  Keeping raw input, generated
tables, figures, and pipeline outputs separate avoids confusion later when
you begin filtering, static-shift correction, or inversion preparation.

.. code-block:: text
   :linenos:

   my_survey/
   |-- data/
   |   `-- edis/
   |       |-- S001.edi
   |       |-- S002.edi
   |       `-- ...
   |-- results/
   |   `-- first_survey/
   `-- notebooks/

The examples below assume the EDI directory is ``data/edis`` and the output
directory is ``results/first_survey``.  Replace these paths with your own
survey paths.

Configure The Session
---------------------

pyCSAMT has separate configuration layers for API views, plotting style, and
pipeline output.  For a first run, keep the configuration explicit at the top
of the script or notebook.

.. code-block:: python
   :linenos:

   from pycsamt.api import configure_api_view, configure_pipe, use_style

   configure_api_view(backend="pycsamt")
   configure_pipe(
       output_root="results/first_survey",
       show_progress=True,
       plot_dpi=200,
       plot_fmt="png",
   )
   use_style("publication")

The configuration call is intentionally lightweight:

``configure_api_view(backend="pycsamt")``
    Keeps returned tables as pyCSAMT ``APIFrame`` objects.  These objects still
    expose the underlying pandas dataframe through ``.df``.

``configure_pipe(output_root=...)``
    Sets the default root directory for pipeline outputs when a pipeline run
    does not receive an explicit ``outdir``.

``use_style("publication")``
    Applies the package plotting defaults for readable saved figures.

For a deeper explanation of these layers, see :doc:`configuration`.

Read The Survey
---------------

The simplest survey reader is :func:`pycsamt.api.read_edis`.  It accepts a
single file, a directory, or a sequence of files/directories and returns an
``APISurvey`` object.

.. code-block:: python
   :linenos:

   from pycsamt.api import read_edis

   survey = read_edis(
       "data/edis",
       recursive=True,
       strict=False,
       on_dup="replace",
       progress="auto",
   )

   print(survey)

The important options are:

``recursive=True``
    Search subdirectories below ``data/edis``.  This is useful when field
    exports are grouped by line, date, or crew.

``strict=False``
    Load every readable file and collect parser failures for review.  Use
    ``strict=True`` when you want the read operation to fail immediately on
    the first malformed file.

``on_dup="replace"``
    If duplicate station names are found, keep the last station encountered.
    Use ``on_dup="keep"`` only when you know that the first occurrence is the
    one you want.

``progress="auto"``
    Shows progress when the current environment supports it.

Inspect The Loaded Survey
-------------------------

Start with the survey object itself.  It exposes a compact public interface
for common first-pass questions.

.. code-block:: python
   :linenos:

   print("Number of stations:", survey.n_sites)
   print("First stations:", survey.stations[:10])
   print("First paths:", survey.paths[:3])

   if survey.n_sites == 0:
       raise RuntimeError("No EDI stations were loaded.")

``survey.stations`` is the station-name inventory.  ``survey.paths`` is the
source-file inventory.  If the station count is lower than expected, check
for duplicate names, failed parses, or an input path that does not include
all exported EDI files.

Check Parser Errors
-------------------

When ``strict=False``, pyCSAMT keeps readable stations and stores parser
errors for the files that could not be loaded.  Always check these errors
before trusting a first survey summary.

.. code-block:: python
   :linenos:

   errors = survey.errors()

   if errors:
       print(f"{len(errors)} file(s) could not be read.")
       for path, exc in errors[:5]:
           print(f"- {path}: {type(exc).__name__}: {exc}")
   else:
       print("All discovered EDI files were read successfully.")

If only a few files fail, you can continue with QC while you repair or remove
the bad files.  If many files fail, stop and inspect the export format before
running processing steps.

Build A Station Summary
-----------------------

``survey.summary()`` returns an ``APIFrame``.  It prints cleanly in terminals,
but it also behaves like a dataframe wrapper.

.. code-block:: python
   :linenos:

   summary = survey.summary()
   print(summary)

   station_table = summary.df
   print(station_table.head())
   print(station_table.columns)

For first-pass review, the most useful columns are usually:

``station``
    Station name resolved from the EDI metadata or file identity.

``n_freq``
    Number of frequency samples available for the station.

``period_min`` and ``period_max``
    Period range represented by the station.

``has_tipper``
    Whether vertical-field transfer-function information is available.

``lat`` and ``lon``
    Geographic coordinates when present in the EDI metadata.

You can request a smaller station table when preparing a report or debugging
specific metadata.

.. code-block:: python
   :linenos:

   compact = survey.summary(
       fields=("station", "n_freq", "period_min", "period_max", "lat", "lon")
   )
   print(compact.df)

Inspect One Station
-------------------

Before plotting the whole survey, inspect one station object.  This confirms
that station naming, frequency arrays, impedance tensors, and optional tipper
metadata are reachable.

.. code-block:: python
   :linenos:

   station_name = survey.stations[0]
   site = survey.get_site(station_name)

   print("Station:", station_name)
   print(site)

   # Many EDI-like objects expose these attributes.  Use getattr so the
   # first workflow remains tolerant of partial metadata.
   print("Path:", getattr(site, "path", None))
   print("Frequencies:", getattr(site, "freq", None))
   print("Has tipper:", getattr(site, "tipper", None) is not None)

If ``survey.get_site(name)`` returns ``None``, the station name you supplied
does not match a loaded station, filename stem, or source path.

Use The Public Table Helper
---------------------------

The survey object is convenient when you want to keep a loaded collection in
memory.  For a quick script, the public table helper can read a directory and
return a station summary in one call.

.. code-block:: python
   :linenos:

   from pycsamt.api import sites_summary

   sites = sites_summary(
       "data/edis",
       recursive=True,
       strict=False,
       on_dup="replace",
   )

   print(sites)
   print(sites.df.head())

Use ``read_edis`` when the next step will reuse the same loaded survey.  Use
``sites_summary`` when you only need a station inventory table.

Build A First Quality Table
---------------------------

The station inventory tells you what was loaded.  The quality table starts to
answer whether the transfer-function data are suitable for processing.

.. code-block:: python
   :linenos:

   from pycsamt.api import quality_dataframe

   quality = quality_dataframe(survey.to_collection())
   print(quality)
   print(quality.df.head())

Depending on the metadata available in the loaded EDI files, the quality
table may include coverage, missing-value, uncertainty, and station-level
diagnostic columns.  Treat it as a screening table: it helps you identify
stations that deserve closer plotting, not as a replacement for geophysical
judgment.

Save A First QC Figure
----------------------

A confidence profile is a compact first figure because it reduces the station
quality assessment to one value per station along the profile.  It is a good
way to spot stations that require cleaning, masking, or removal before
inversion.

.. code-block:: python
   :linenos:

   from pathlib import Path

   import matplotlib.pyplot as plt

   from pycsamt.emtools.qc import plot_confidence_profile

   outdir = Path("results/first_survey")
   outdir.mkdir(parents=True, exist_ok=True)

   ax = plot_confidence_profile(
       survey.to_collection(),
       method="composite",
       ci_hi=0.95,
       ci_lo=0.50,
       station_label_step=None,
   )
   fig = ax.get_figure()
   fig.savefig(outdir / "station_confidence_profile.png", dpi=200,
               bbox_inches="tight")
   plt.close(fig)

Use ``method="presence"`` for a very simple completeness check.  Use
``method="composite"`` when you want a more informative first-pass score that
combines several transfer-function diagnostics.

For a station-frequency view, save a pseudo-section as well.

.. code-block:: python
   :linenos:

   from pycsamt.emtools.qc import plot_frequency_confidence_psection

   ax = plot_frequency_confidence_psection(
       survey.to_collection(),
       method="composite",
   )
   fig = ax.get_figure()
   fig.savefig(outdir / "frequency_confidence_psection.png", dpi=200,
               bbox_inches="tight")
   plt.close(fig)

The station profile shows where weak stations are located along the line.
The frequency pseudo-section shows whether the weak data are isolated at
particular periods or distributed through an entire station.

Run A First Pipeline
--------------------

Once the survey can be read and inspected, run the built-in ``basic_qc``
preset.  This is the recommended first pipeline because it is small enough
for fast feedback and structured enough to produce reproducible outputs.

.. code-block:: python
   :linenos:

   from pycsamt.pipeline import Pipeline

   pipe = Pipeline.from_preset("basic_qc")
   print(pipe)

   result = pipe.run(
       survey.to_collection(),
       outdir="results/first_survey/basic_qc",
       save_plots=True,
       save_edis=True,
       save_report=True,
   )

   print(result.summary())
   print("Pipeline ok:", result.ok)
   print("Step errors:", result.n_errors)
   print("Output directory:", result.outdir)

The pipeline writes a reproducible run directory.  The exact files depend on
the enabled output settings and the steps in the preset, but the first things
to inspect are:

``pipeline.yaml``
    The pipeline configuration that produced the run.

``reports/``
    Text or HTML summaries of the run, including step status.

``plots/``
    QC plots produced by the pipeline steps.

``processed/``
    Processed EDI files when ``save_edis=True``.

If ``result.ok`` is ``False``, inspect the reported step errors before using
the processed files.  With the default warning policy, a failed step may leave
the previous station collection in place so that the rest of the run can
continue.

Command-Line Equivalent
-----------------------

The same first pass can be started from the command line.  These commands are
useful for quick validation, automation, and sharing a workflow with another
developer.

.. code-block:: console
   :linenos:

   pycsamt edi info data/edis
   pycsamt edi stations data/edis --format text
   pycsamt edi validate data/edis

   pycsamt pipe run data/edis \
       --preset basic_qc \
       --out results/first_survey/basic_qc \
       --on-error warn \
       --dpi 200 \
       --plot-fmt png

Use ``pycsamt pipe run ... --dry-run`` before a long run to confirm that the
survey path, preset, step count, and output directory resolve as expected.

.. code-block:: console
   :linenos:

   pycsamt pipe run data/edis \
       --preset basic_qc \
       --out results/first_survey/basic_qc \
       --dry-run

First-Survey Checklist
----------------------

Before moving to filtering, static-shift correction, dimensionality analysis,
or inversion preparation, confirm the following:

* ``survey.n_sites`` matches the expected number of field stations.
* ``survey.errors()`` is empty or every failed file is understood.
* Station names are unique and ordered as expected.
* Coordinates are present when profile distance or maps matter.
* Frequency/period ranges are comparable across neighbouring stations.
* Missing tipper data are expected, or downstream steps do not require tipper.
* The confidence profile does not show unexplained clusters of weak stations.
* The pipeline report has no unexpected failed steps.
* The output directory contains the run configuration and generated reports.

Troubleshooting
---------------

No stations are loaded
    Check that the input path exists and contains files with an EDI extension.
    If the files are nested by survey line or date, keep ``recursive=True``.

The station count is too low
    Look for duplicate station names and parser errors.  Re-run with
    ``on_dup="keep"`` if the first duplicate should win, or fix station names
    in the source files before processing.

Many files appear in ``survey.errors()``
    Confirm the files are EDI files rather than SEG blocks, Zonge AVG files,
    or another export format.  See :doc:`data_formats` for format-specific
    entry points.

Coordinates are missing
    The EDI metadata may not include latitude and longitude.  You can still
    inspect frequencies and transfer functions, but maps and distance-based
    plots require station-location metadata.

The plot is empty
    Confirm that the survey contains stations and that the station objects
    expose frequency and impedance arrays.  If the dataset is partial, try a
    simpler table first with ``survey.summary().df``.

The pipeline writes somewhere unexpected
    Pass ``outdir=...`` to ``Pipeline.run`` or ``--out ...`` to the CLI.  The
    pipeline otherwise falls back to its own configured output directory and
    then to the global pipeline configuration.

Where To Go Next
----------------

After a successful first survey pass, continue with the page that matches
your task:

* :doc:`data_formats` for converting or understanding non-EDI inputs.
* :doc:`../tutorials/read_edi_survey` for a deeper EDI-reading tutorial.
* :doc:`../tutorials/inspect_and_qc_survey` for richer station diagnostics.
* :doc:`../user_guide/pipeline/presets` for choosing a processing preset.
* :doc:`../user_guide/pipeline/outputs` for understanding generated pipeline
  files.
* :doc:`../tutorials/prepare_occam2d_inversion` when the data are ready for
  Occam2D preparation.
* :doc:`../user_guide/agents/overview` if you want the AI-assisted agents to
  help review a survey, plan processing, or summarize results.

In Short
--------

A good first pyCSAMT survey pass is deliberately simple:

.. code-block:: python
   :linenos:

   from pycsamt.api import read_edis, quality_dataframe
   from pycsamt.pipeline import Pipeline

   survey = read_edis("data/edis", strict=False)
   print(survey)
   print(survey.summary().df.head())
   print(quality_dataframe(survey.to_collection()).df.head())

   pipe = Pipeline.from_preset("basic_qc")
   result = pipe.run(survey.to_collection(), outdir="results/first_survey")
   print(result.summary())

If this short script reads the survey, prints sensible station tables, and
produces a successful ``basic_qc`` run, the project is ready for deeper
processing.
