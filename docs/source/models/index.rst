.. _models:

Model Integrations
==================

The :mod:`pycsamt.models` package contains pyCSAMT's interfaces to external
electromagnetic modelling and inversion ecosystems. These modules are not
only data containers. They prepare engine-specific input files, validate
working directories, run external binaries when requested, read native
outputs, and convert results back into pyCSAMT objects for plotting,
interpretation, and reporting.

This section documents the model-integration layer separately from the
generated API reference because users often need conceptual guidance before
they need class signatures. The generated API remains available at
:doc:`../api/models`.

What This Section Covers
------------------------

Model integrations sit at the boundary between pyCSAMT and external modelling
codes. They are intentionally practical: each page explains what the backend
is for, which native files matter, how configuration is organized, how runs
are prepared, and which diagnostics should be reviewed before interpretation.

Use this section when you need to:

* choose between Occam2D, ModEM, and MARE2DEM;
* prepare native input files in a reproducible run directory;
* connect an external executable to a pyCSAMT workflow;
* load completed inversion products back into Python;
* understand how backend-specific files relate to the common inversion API;
* document model settings for a project report or review.

Where Models Fit
----------------

pyCSAMT has three related layers:

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * - Layer
     - Purpose
     - Typical modules
   * - Inversion workflow
     - Backend-neutral configuration and result handling.
     - :mod:`pycsamt.inversion`
   * - Model integrations
     - Engine-specific files, runners, logs, meshes, and plotting.
     - :mod:`pycsamt.models`
   * - Pipeline and agents
     - Repeatable orchestration, QC, reports, and assisted workflows.
     - :mod:`pycsamt.pipeline`, :mod:`pycsamt.agents`

For example, an :class:`pycsamt.inversion.config.InversionConfig` can select
``backend="occam2d"`` or ``backend="modem"``. The inversion workflow then
delegates engine-specific preparation and loading to the corresponding model
subpackage.

Backend Summary
---------------

.. list-table::
   :header-rows: 1
   :widths: 18 24 32 26

   * - Package
     - External family
     - Best suited for
     - Start here
   * - :mod:`pycsamt.models.occam2d`
     - Occam2DMT-style smooth 2-D inversion
     - MT, AMT, and CSAMT profile inversions where a smooth 2-D section is
       desired and the survey geometry is line based.
     - :doc:`occam2d`
   * - :mod:`pycsamt.models.modem`
     - ModEM modular EM inversion
     - 2-D and 3-D MT or AMT workflows, covariance-controlled 3-D models,
       ModEM file management, and native ModEM project review.
     - :doc:`modem`
   * - :mod:`pycsamt.models.mare2dem`
     - MARE2DEM finite-element modelling and inversion
     - 2.5-D MT and CSEM workflows, finite-element meshes, topography-aware
       model studies, and MARE2DEM-native project folders.
     - :doc:`mare2dem`
   * - :mod:`pycsamt.models.config_io`
     - Shared config readers and writers
     - Reproducible model templates and backend options
     - :doc:`configuration_and_io`

Quick Selection Guide
---------------------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Situation
     - Recommended starting page
   * - You have a 2-D survey line and want a smooth resistivity section.
     - Start with :doc:`occam2d`.
   * - You have a station grid or off-profile structure is important.
     - Start with :doc:`modem` and review the 3-D input-set workflow.
   * - You need finite-element modelling with 2.5-D geometry or CSEM-style
       control over source/receiver files.
     - Start with :doc:`mare2dem`.
   * - You are unsure which backend fits the survey.
     - Read :doc:`choosing_backend` before configuring files.
   * - You already have native files from another team or an older project.
     - Read the engine page, then :doc:`configuration_and_io` for run-folder
       and template policy.

Design Goals
------------

The model layer follows a few practical rules:

* external binaries are not hidden from the user;
* generated files are written into explicit working directories;
* native file formats remain inspectable;
* runners are separate from builders and result loaders;
* validation should happen before expensive external execution;
* pyCSAMT result objects should preserve provenance and native files.

This separation matters because field projects often need to rerun only one
stage: rebuild inputs, launch a binary on another machine, reload finished
outputs, or regenerate figures without rerunning inversion.

Recommended Reading Order
-------------------------

For a new project:

#. Read :doc:`choosing_backend` to choose an engine.
#. Read the engine page: :doc:`occam2d`, :doc:`modem`, or :doc:`mare2dem`.
#. Read :doc:`configuration_and_io` for reproducible templates, run folders,
   and file naming.
#. Connect the model page to :doc:`../pipeline/index` when the run should be
   repeated or automated.
#. Use :doc:`../api/models` when you need class-level details.

For an existing native run:

#. Identify the engine and native file set.
#. Use the matching engine page to load or validate files.
#. Review result diagnostics before reusing the model.
#. Archive the configuration, command, logs, and final interpretation notes.

Relationship To Theory
----------------------

Model pages assume the scientific background described in
:doc:`../theory/index`. In particular:

* :doc:`../theory/impedance_tensor` explains the impedance components used by
  ModEM and Occam-style MT/AMT workflows.
* :doc:`../theory/inversion_concepts` explains RMS, regularization, model
  roughness, and uncertainty floors.
* :doc:`../theory/static_shift` explains a common distortion that can strongly
  affect resistivity models.
* :doc:`../theory/csamt_amt_mt_overview` explains method differences that
  should influence backend choice.

Contents
--------

.. toctree::
   :maxdepth: 1

   choosing_backend
   configuration_and_io
   occam2d
   modem
   mare2dem
