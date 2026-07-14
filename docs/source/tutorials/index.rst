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
   * - :doc:`compare_survey_lines_for_qc`
     - You need to decide whether two line surveys can start from the same
       first-pass QC and processing configuration.
     - Line-comparison tables, frequency overlap, confidence distributions,
       processing decision table.
   * - :doc:`correct_static_shift`
     - You have inspected the survey and want a conservative static-shift
       correction workflow before inversion preparation.
     - Static-shift factors, corrected site collection, corrected EDIs.
   * - :doc:`condition_mt_line_with_tipper_and_rotation`
     - You have MT data with tipper and full impedance tensors and need a
       transparent conditioning workflow before inversion.
     - Raw tensor/tipper plots, bad-frequency screen, filters, static shift,
       strike rose, phase tensor grid, rotated tensors.
   * - :doc:`prepare_occam2d_inversion`
     - You have cleaned and reviewed a line survey and want to prepare a 2-D
       Occam2D inversion workspace.
     - Occam2D data, model, startup, and run-directory files.
   * - :doc:`ai_inversion_from_corrected_edis`
     - You have corrected EDIs and want to use AI inversion instead of, or
       alongside, classical inversion.
     - AI-mode decision table, training coverage check, 1-D/2-D/3-D code
       paths, prediction and validation figures.
   * - :doc:`essential_3d_ai_inversion`
     - You want a focused L18PLT 3-D AI inversion workflow using real station
       geometry, real EDI topography, and agent-returned depth layers.
     - L18 station-topography profile and embedded-topography 3-D AI inversion
       block.
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

3. Use :doc:`compare_survey_lines_for_qc` when a project has several lines.

   This page helps you decide whether one first-pass config can be reused
   across related lines, or whether each line needs separate QC parameters.

4. Move to :doc:`correct_static_shift` when the QC pass is understood.

   Static-shift correction should not be the first operation applied blindly.
   The tutorial shows how to estimate, inspect, apply, and export correction
   factors conservatively.

5. Use :doc:`condition_mt_line_with_tipper_and_rotation` for advanced MT data.

   KP-style MT data with tipper need a fuller pre-inversion review: raw tensor
   curves, tipper response, weak-frequency handling, static-shift review,
   phase tensors, strike, and rotation.

6. Use :doc:`prepare_occam2d_inversion` for a cleaned profile.

   This tutorial connects the processing workflow to inversion preparation. It
   focuses on the files and checks needed before handing data to Occam2D.

7. Use :doc:`ai_inversion_from_corrected_edis` when choosing AI inversion.

   This page starts from corrected EDIs and explains when to use 1-D AI, when
   to bypass it for 2-D profile inversion, and when 3-D graph AI is justified.

8. Use :doc:`essential_3d_ai_inversion` for a focused real-topography 3-D AI run.

   This page treats L18PLT as corrected EDI, runs 3-D AI, and drapes the result
   on package-extracted station topography.

9. Use :doc:`run_pipeline_from_config` when the workflow should be repeated.

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

.. grid:: 1 1 2 3
   :gutter: 3
   :class-container: cta-tiles

   .. grid-item-card:: Read an EDI Survey
      :link: read_edi_survey
      :link-type: doc
      :img-top: ../_static/icons/user-guide-data-loading.svg
      :class-card: pycsamt-card sd-text-center

      Load one file, a line directory, or a survey tree and inspect
      what pyCSAMT found.

   .. grid-item-card:: Inspect & QC a Survey
      :link: inspect_and_qc_survey
      :link-type: doc
      :img-top: ../_static/icons/diagnostic.svg
      :class-card: pycsamt-card sd-text-center

      Review coverage, usable frequencies, signal quality, and stations
      that need attention.

   .. grid-item-card:: Compare Survey Lines
      :link: compare_survey_lines_for_qc
      :link-type: doc
      :img-top: ../_static/icons/external-validation.svg
      :class-card: pycsamt-card sd-text-center

      Decide whether two lines can share one first-pass QC and
      processing configuration.

   .. grid-item-card:: Correct Static Shift
      :link: correct_static_shift
      :link-type: doc
      :img-top: ../_static/icons/static-shift.svg
      :class-card: pycsamt-card sd-text-center

      Estimate, inspect, apply, and export a conservative static-shift
      correction.

   .. grid-item-card:: Condition an MT Line
      :link: condition_mt_line_with_tipper_and_rotation
      :link-type: doc
      :img-top: ../_static/icons/transfer-tipper.svg
      :class-card: pycsamt-card sd-text-center

      Tipper review, bad-frequency screening, phase tensors, strike,
      and rotation before inversion.

   .. grid-item-card:: Prepare Occam2D Inversion
      :link: prepare_occam2d_inversion
      :link-type: doc
      :img-top: ../_static/icons/user-guide-inversion.svg
      :class-card: pycsamt-card sd-text-center

      Build the data, model, startup, and run-directory files for a
      2-D Occam2D workspace.

   .. grid-item-card:: AI Inversion from EDIs
      :link: ai_inversion_from_corrected_edis
      :link-type: doc
      :img-top: ../_static/icons/user-guide-ai.svg
      :class-card: pycsamt-card sd-text-center

      Choose between 1-D, 2-D, and 3-D AI code paths starting from
      corrected EDIs.

   .. grid-item-card:: Essential 3-D AI Inversion
      :link: essential_3d_ai_inversion
      :link-type: doc
      :img-top: ../_static/icons/model.svg
      :class-card: pycsamt-card sd-text-center

      A focused L18PLT 3-D AI run with real station geometry and
      embedded topography.

   .. grid-item-card:: Run a Pipeline from Config
      :link: run_pipeline_from_config
      :link-type: doc
      :img-top: ../_static/icons/workflow.svg
      :class-card: pycsamt-card sd-text-center

      Store the workflow in YAML, JSON, or Python and rerun it from
      Python or the CLI.

.. toctree::
   :maxdepth: 1
   :hidden:

   read_edi_survey
   inspect_and_qc_survey
   compare_survey_lines_for_qc
   correct_static_shift
   condition_mt_line_with_tipper_and_rotation
   prepare_occam2d_inversion
   ai_inversion_from_corrected_edis
   essential_3d_ai_inversion
   run_pipeline_from_config

Related Sections
----------------

- :doc:`../user_guide/pipeline/index`
- :doc:`../user_guide/site/index`
- :doc:`../user_guide/models/index`
- :doc:`../user_guide/forward/index`
- :doc:`../cli/index`
- :doc:`../api/api`
