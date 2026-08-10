.. _user-guide-workflow-overview:
.. _user_guide_prerequisites:

Choose a workflow
=================

The user guide explains how to move electromagnetic observations from field
files to reviewed scientific products. It assumes that pyCSAMT is already
installed and importable. If that is not yet true, follow the guided
:doc:`../getting_started/installation` page or consult the complete
:doc:`../installation` reference.

Most workflows share the same backbone:

.. code-block:: text

   identify the input format
       -> load and inspect the survey
       -> check data quality
       -> correct or transform the responses
       -> model or invert
       -> interpret and report

The exact route depends on the instrument output and the scientific question.
Choose the nearest starting point below rather than reading every section in
order.

Starting from EDI files
-----------------------

Begin with :doc:`data_loading`. It explains the public ``read_edis`` survey
view, the ``ensure_sites`` processing boundary, station ordering, duplicate
handling, and frequency-grid checks. Continue to :doc:`emtools/index` for
quality control, noise removal, frequency editing, dimensionality analysis,
and static-shift correction.

Map the inspected stations with :doc:`map/index`. If elevation must influence
a section or model, prepare it with :doc:`topo/index` before interpreting the
result. A syntactically valid EDI file is not necessarily scientifically
usable, so do not skip the inventory and quality-control stages.

Starting from instrument or intermediate formats
------------------------------------------------

Use :doc:`transformers` when the input is AVG, J-format, spectra, time-series,
or transient electromagnetic data that must be converted or normalized before
the standard EDI workflow.

Geometrics/EMI Stratagem surveys have additional hardware-specific steps.
:doc:`stratagem/index` covers raw import, coordinate injection, WinGLink
interchange, and validation before the data enter the general processing
tools. For connected acquisition systems, :doc:`iot/index` covers device
configuration, telemetry, timing, edge quality control, and provenance.

Preparing an inversion
----------------------

After reviewing the observations, use :doc:`inversion/index` to choose an
inversion path. That guide separates three related layers:

* :doc:`models/index` configures and runs classical external solvers;
* :doc:`ai_inversion/index` covers learned and physics-informed inversion;
* :doc:`forward/index` generates responses from known earth models for
  sensitivity studies, validation, and synthetic training data.

An inversion result is not the end of the workflow. Continue to
:doc:`interpretation/index` to evaluate uncertainty, geological context, and
the claims that the data can support. Use the mapping and topography guides to
prepare spatially defensible figures rather than treating model colours as a
standalone interpretation.

Working with stations directly
------------------------------

:doc:`site/index` documents station-level selection, coordinate and metadata
editing, profile construction, normalization, and export. Use it when the
scientific task concerns which stations enter the workflow or how station
metadata are represented. Use :doc:`emtools/index` when the task concerns the
electromagnetic responses themselves.

Automating a reviewed workflow
------------------------------

First establish the workflow with explicit API calls and inspect its
intermediate products. Then use :doc:`pipeline/index` to record the sequence
as a named, versioned recipe. :doc:`agents/index` adds higher-level automation
and orchestration; it does not remove the need to validate data, configuration,
and scientific assumptions.

.. warning::

   Automation can reproduce an invalid workflow as reliably as a valid one.
   Record station identities, input paths, coordinate reference information,
   processing parameters, software version, and random seeds where applicable.

Finding the right kind of documentation
---------------------------------------

Use the user guide for concepts, decisions, and reusable workflows. The other
documentation sections serve different purposes:

* :doc:`../getting_started/index` installs pyCSAMT and completes a first run;
* :doc:`../tutorials/index` provides longer end-to-end worked examples;
* :doc:`../theory/index` explains the physical and mathematical background;
* :doc:`../api/index` gives signatures and object-level API details.

When you know the scientific goal but not the function name, stay in the user
guide. When you already know the object and need its exact parameters or return
type, switch to the API reference.
