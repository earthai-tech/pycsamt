.. _tutorials:

Tutorials
=========

Tutorials are practical, end-to-end recipes for common pyCSAMT survey tasks.
They are written in a workflow style: each page states the input assumptions,
shows runnable Python snippets, includes CLI equivalents where available,
describes expected outputs, and points to the related API and user-guide
sections.

The tutorial sequence follows the normal path from raw EDI files to reviewed,
processed data that can be used for modelling or inversion.

How to Use These Tutorials
--------------------------

If you are new to pyCSAMT v2, read the tutorials in order. The first two pages
teach the basic survey object and QC workflow. The later pages build on that
foundation.

If you already know the data-loading API, use the table below to jump to the
task you need.

.. list-table::
   :header-rows: 1
   :widths: 26 44 30

   * - Tutorial
     - Use It When
     - Main Outputs
   * - :doc:`read_edi_survey`
     - You need to load one EDI file, a line directory, or a recursive survey
       tree and inspect what pyCSAMT found.
     - ``APISurvey``, station inventory, parser diagnostics.
   * - :doc:`inspect_and_qc_survey`
     - You want to review station coverage, usable frequencies, skew proxies,
       signal quality, and stations that need manual attention.
     - QC tables, confidence tables, review CSV files, diagnostic plots.
   * - :doc:`correct_static_shift`
     - You have inspected the survey and want a conservative static-shift
       correction workflow before inversion preparation.
     - Static-shift factors, corrected site collection, corrected EDIs.
   * - :doc:`prepare_occam2d_inversion`
     - You have cleaned and reviewed a line survey and want to prepare a 2-D
       Occam2D inversion workspace.
     - Occam2D data, model, startup, and run-directory files.
   * - :doc:`run_pipeline_from_config`
     - You want a repeatable processing workflow stored in YAML, JSON, or
       Python config files and runnable from Python or the CLI.
     - Processed EDIs, plots, ``pipeline.yaml``, ``summary.txt``,
       ``report.html``.

Recommended Learning Path
-------------------------

1. Start with :doc:`read_edi_survey`.

   This page introduces :func:`pycsamt.api.read_edis`, the ``APISurvey``
   object, station summaries, duplicate policies, parser errors, and the
   basic EDI CLI commands.

2. Continue with :doc:`inspect_and_qc_survey`.

   Use this after loading the survey. It helps you decide whether the data are
   ready for correction, whether some stations should be reviewed, and which
   frequency bands look reliable.

3. Move to :doc:`correct_static_shift` when the QC pass is understood.

   Static-shift correction should not be the first operation applied blindly.
   The tutorial shows how to estimate, inspect, apply, and export correction
   factors conservatively.

4. Use :doc:`prepare_occam2d_inversion` for a cleaned profile.

   This tutorial connects the processing workflow to inversion preparation. It
   focuses on the files and checks needed before handing data to Occam2D.

5. Use :doc:`run_pipeline_from_config` when the workflow should be repeated.

   Once the steps are stable, move them into a config file. This gives you a
   reproducible processing chain that can be reviewed, rerun, and archived with
   the results.

Before You Start
----------------

The tutorials assume that pyCSAMT is installed and that you have a folder of
EDI files. If you are still setting up the environment or learning the data
formats, read these pages first:

- :doc:`../getting_started/installation`
- :doc:`../getting_started/data_formats`
- :doc:`../getting_started/first_survey`
- :doc:`../getting_started/configuration`

Most examples use paths such as ``data/edis`` and ``results/first_qc``. Replace
them with your local project paths.

Workflow Conventions
--------------------

The tutorial pages use the same conventions:

``survey``
    Public survey object returned by :func:`pycsamt.api.read_edis`.

``sites``
    Lower-level site or EDI collection passed into QC, editing, pipeline, and
    export helpers. In many examples it is obtained from ``survey.collection``.

``result``
    Returned object from an operation, such as a pipeline run or inversion
    preparation step.

``outdir``
    Root output directory for generated files.

``strict=False``
    Used during first inspection so the reader can continue past recoverable
    EDI issues.

``strict=True`` or ``--on-error raise``
    Used for production validation when the workflow should fail immediately
    on malformed input or processing errors.

Choosing Python or the CLI
--------------------------

Use Python when you need to keep objects in memory, combine tables, make custom
plots, or integrate pyCSAMT with notebooks and analysis scripts.

Use the CLI when you need quick inspection, reproducible batch processing, or a
shell-friendly workflow:

.. code-block:: bash
   :linenos:

   pycsamt edi info data/edis
   pycsamt edi validate data/edis
   pycsamt pipe init --preset basic_qc --name first_qc
   pycsamt pipe run --config first_qc.yaml --survey data/edis --out results/first_qc

The tutorials usually show both forms when a CLI command exists.

Tutorial Index
--------------

.. toctree::
   :maxdepth: 1

   read_edi_survey
   inspect_and_qc_survey
   correct_static_shift
   prepare_occam2d_inversion
   run_pipeline_from_config

Related Sections
----------------

- :doc:`../pipeline/index`
- :doc:`../site/index`
- :doc:`../models/index`
- :doc:`../forward/index`
- :doc:`../cli/index`
- :doc:`../api/api`
