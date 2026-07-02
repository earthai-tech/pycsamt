.. _getting-started:

Getting Started
===============

This section is the shortest path from a clean Python environment to a first
pyCSAMT survey workflow.  It is written for users who want to install the
package, understand which data format they have, configure a working session,
and inspect a first field survey before moving into processing, inversion, or
AI-assisted workflows.

Start here when you are new to pyCSAMT v2, when you are setting up a new
machine, or when you need a reproducible first-pass workflow for a survey
directory.

Fast Path
---------

For a standard Python environment and a folder of EDI files, the complete
first path is:

.. code-block:: console
   :linenos:

   python -m pip install pycsamt
   pycsamt --help
   pycsamt edi info data/edis
   pycsamt edi validate data/edis

Then in Python:

.. code-block:: python
   :linenos:

   from pycsamt.api import read_edis
   from pycsamt.pipeline import Pipeline

   survey = read_edis("data/edis", strict=False)
   print(survey)
   print(survey.summary().df.head())

   pipe = Pipeline.from_preset("basic_qc")
   result = pipe.run(survey.to_collection(), outdir="results/first_survey")
   print(result.summary())

If that works, pyCSAMT is installed, the EDI reader can see your data, the
public API is available, and the processing pipeline can produce a first
quality-control run.

For a condensed end-to-end demonstration — processing, Occam2D inversion,
geological interpretation, and AI-based inversion in a few code blocks —
see :doc:`quickstart`.

Recommended Reading Order
-------------------------

Read the pages in this order for the smoothest first experience:

.. list-table::
   :header-rows: 1
   :widths: 22 38 40

   * - Step
     - Page
     - Use it to
   * - 1
     - :doc:`installation`
     - Choose an environment, install core pyCSAMT or optional extras, and
       verify the CLI, API, AI backends, GUI, and documentation tools.
   * - 2
     - :doc:`data_formats`
     - Identify whether your inputs are EDI, Zonge AVG/AMTAVG, J-format,
       spectra, TDEM/TEM, site tables, or inversion/model files.
   * - 3
     - :doc:`configuration`
     - Configure API table behavior, plotting style, pipeline output, AI
       backend selection, and agent/LLM settings.
   * - 4
     - :doc:`first_survey`
     - Load a first EDI survey, inspect station metadata, check parser
       errors, build QC tables, save figures, and run ``basic_qc``.

The order matters because pyCSAMT workflows become much easier when you know
which optional dependencies are installed, which field format you are using,
and where outputs will be written before you begin processing data.

What Each Page Covers
---------------------

:doc:`installation`
    Installation recipes for PyPI, source checkouts, development, docs,
    geospatial workflows, desktop/web apps, AI backends, and LLM provider
    SDKs.  It also includes verification commands and troubleshooting notes.

:doc:`data_formats`
    A practical guide to supported input families.  It explains the standard
    EDI path, legacy/field formats, conversion routes into EDI, site metadata,
    inversion files, duplicate-station policies, and unit conventions.

:doc:`configuration`
    The core runtime settings users usually touch first: APIFrame versus
    pandas output, pipeline output directories, plotting defaults, AI backend
    selection, and agent configuration.

:doc:`first_survey`
    The first end-to-end survey walkthrough.  It demonstrates
    :func:`pycsamt.api.read_edis`, station summaries, parser-error review,
    quality tables, confidence plots, the ``basic_qc`` pipeline preset, and
    CLI equivalents.

:doc:`quickstart`
    A condensed tour of the core v2 workflow: EDI processing, 2-D inversion
    with Occam2D, geological interpretation and export, and AI-based 1-D
    inversion — each in a short, copy-paste-ready code block.

Choose Your Entry Point
-----------------------

Different users arrive with different needs.  Use this table to jump directly
to the right page.

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - If you want to...
     - Start with...
   * - install pyCSAMT on a new machine
     - :doc:`installation`
   * - decide between ``venv``, conda, core install, and full extras
     - :doc:`installation`
   * - enable the desktop GUI or Dash-based app components
     - :doc:`installation`
   * - enable PyTorch, TensorFlow, or LLM-backed agents
     - :doc:`installation` and :doc:`../agents/llm_configuration`
   * - understand whether your field files are supported
     - :doc:`data_formats`
   * - convert AVG, J-format, spectra, or TEM/TDEM data into an EDI workflow
     - :doc:`data_formats`
   * - configure where pipeline results and figures are written
     - :doc:`configuration`
   * - load an EDI directory and inspect stations
     - :doc:`first_survey`
   * - run the first quality-control pipeline
     - :doc:`first_survey`

Before You Continue
-------------------

Before moving into tutorials, pipeline customization, or inversion
preparation, confirm these basics:

* ``python -c "import pycsamt; print(pycsamt.__version__)"`` prints a version.
* ``pycsamt --help`` shows the command groups.
* You know whether your input data are EDI, AVG, J-format, spectra, TDEM/TEM,
  site metadata, or model/inversion files.
* Your output directory is explicit, for example ``results/first_survey``.
* ``read_edis(..., strict=False)`` loads the expected number of EDI stations.
* Parser errors are understood before processing results are trusted.
* The first ``basic_qc`` pipeline run writes a report and plots.

Where To Go Next
----------------

Once the getting-started workflow is working, continue with:

* :doc:`../tutorials/read_edi_survey` for a deeper EDI loading tutorial;
* :doc:`../tutorials/inspect_and_qc_survey` for richer survey diagnostics;
* :doc:`../pipeline/index` for reproducible processing pipelines;
* :doc:`../pipeline/presets` for built-in workflow presets;
* :doc:`../tutorials/prepare_occam2d_inversion` for a first inversion
  preparation path;
* :doc:`../agents/overview` for AI-assisted survey review, processing plans,
  and result summaries;
* :doc:`../development/index` if you are contributing to pyCSAMT itself.

.. toctree::
   :maxdepth: 1
   :caption: Getting Started Pages

   installation
   data_formats
   configuration
   first_survey
   quickstart
