.. _changelog:

Changelog
=========

A terse, version-by-version record of every notable change. For the
narrative behind each release — highlights, migration guidance, and known
limitations — see the :ref:`Release Notes <release-notes>`.

This project follows `Keep a Changelog <https://keepachangelog.com>`__
conventions and `semantic versioning <https://semver.org>`__. Every entry
carries a colour-coded badge so the log can be scanned at a glance:

.. rst-class:: changelog-legend

* |Feature| |New| — a new capability, module, or application
* |Enhancement| |Perf| — improved behaviour or performance
* |Fix| — a bug fix
* |API Change| |Deprecated| — a change to the public API surface
* |Breaking| |Security| — read carefully before upgrading
* |Docs| |Build| |Tests| — documentation and project tooling

----

.. _changelog-2-0-0rc1:

2.0.0rc1 |Feature| |API Change| |Breaking| |Docs|
-------------------------------------------------

*Release candidate — the first pre-stable tag of the v2 line. The public
API is* **not** *backwards-compatible with v1 and may still change before
2.0.0 final; see* :ref:`migration` *for a v1 → v2 name map.*

Added
~~~~~

* |Feature| **Multi-backend AI engine** — every neural model (FCN,
  ResNet, CNN1D, UNet2D, DRCNN) dispatches transparently to PyTorch
  **or** TensorFlow through a unified abstraction layer
  (:mod:`pycsamt.backends`).
* |Feature| **Named-step processing pipeline** — :mod:`pycsamt.pipeline`
  chains catalogued steps (QC, notch, band selection, static shift,
  rotation, …) into reproducible runs with reports, presets, and YAML
  round-tripping.
* |Feature| **EM analytics toolbox** — :mod:`pycsamt.emtools`: QC tables,
  dimensionality, skew, strike, phase-tensor diagnostics, corrections,
  and one-call plotting for every view.
* |Feature| **Modeling stack** — synthetic forward modelling
  (:mod:`pycsamt.forward`), external engines Occam2D / ModEM / MARE2DEM
  under :mod:`pycsamt.models`, and deep-learning inversion in
  :mod:`pycsamt.ai` (1-D/2-D/3-D nets, ensembles, uncertainty).
* |Feature| **Geological interpretation** — :mod:`pycsamt.interp`
  supersedes ``geodrill``: method-agnostic model calibration, a 25-rock
  EM database, and Oasis Montaj / LAS 2.0 / VTK export.
* |Feature| **AI agents** — specialised loader / QC / processing /
  inversion agents plus a workflow orchestrator that turns a plain-text
  request into a previewable agent chain (:mod:`pycsamt.agents`).
* |Feature| **IoT field telemetry** — :mod:`pycsamt.iot`: edge QC on the
  recorder, telemetry schemas and transports, power budgeting, GPS-sync
  audit, and a live field-session dashboard.
* |New| **Four applications, one engine** — the desktop GUI, the web app,
  the Agent Master chat surface, and the MapView workbench
  (see :doc:`/applications/index`).
* |New| **Offline agent example gallery** — seven runnable, zero-cost
  examples covering routing, context parsing, the model zoo, forward
  modelling, plan validation, coordination, and the ``AgentMaster`` front
  door (see :doc:`/examples/index`).

Changed
~~~~~~~

* |API Change| **Configuration is explicit** — runtime behaviour (output
  directories, plot styles, axis conventions, agent budgets) is set
  through the documented ``configure_*`` / ``reset_*`` families instead
  of scattered keyword arguments; see the
  :doc:`API configuration guide </api_guide/index>`.
* |API Change| **EDI is the source of truth** — the impedance
  :class:`~pycsamt.z.z.Z` object and
  :class:`~pycsamt.io.EDICollection` now drive every processing and
  inversion pipeline.

Fixed
~~~~~

* |Fix| **Encoding-safe figure saving** — saving a figure or writing the
  plot config no longer raises ``UnicodeEncodeError`` on a legacy Windows
  (cp1252) console; the file is already on disk, so the save is no longer
  reported as failed (:mod:`pycsamt.api.plot`).
* |Fix| **Robust agent previews** — ``AgentCoordinator`` dry-run previews
  route their output through an encoding-safe printer, so a preview never
  crashes on non-ASCII glyphs (:mod:`pycsamt.agents`).

Removed
~~~~~~~

* |Breaking| ``pycsamt.geodrill`` — removed; use :mod:`pycsamt.interp`.
* |Breaking| ``pycsamt.ff`` — removed; use :mod:`pycsamt.emtools`.
* |Breaking| ``pycsamt.viewer`` — removed; use module-level ``.plot()``
  methods on result objects.
* |Breaking| ``pycsamt.modeling`` — removed; use :mod:`pycsamt.models`.
* |Breaking| **All v1 CLI entry points** — removed and replaced by the v2
  command set.

Docs & tooling
~~~~~~~~~~~~~~~

* |Docs| **Documentation rebuilt** — executed example galleries for every
  tool family (:doc:`/examples/index`), task-based user guides,
  application manuals, and a full autodoc reference.
* |Docs| **Reference hub** — a new :doc:`Resources <resources>` landing
  page groups the :doc:`Glossary <glossary>` of MT/AMT/CSAMT terms and the
  :doc:`bibliography <references>`.
* |Build| **PyPI-installable** — pure ``pyproject.toml``; ``setup.py``
  removed.

.. _migration:

Migration from v1
-----------------

The v1 module tree is replaced by the subsystem packages above. v1 scripts
do **not** run unchanged against v2 — start from the
:doc:`user guide </user_guide/index>` equivalents of your workflow, and use
this map for the most common names:

.. list-table::
   :header-rows: 1
   :widths: 45 55
   :class: migration-table

   * - v1
     - v2 equivalent
   * - ``pycsamt.geodrill.geocore.GeoStratigraphy``
     - :class:`pycsamt.interp.ModelCalibrator`
   * - ``pycsamt.geodrill.geocore.Geodrill``
     - :class:`pycsamt.interp.ResistivityModel` + :class:`~pycsamt.interp.ModelCalibrator`
   * - ``pycsamt.geodrill.geodatabase.GeoDataBase``
     - :class:`pycsamt.interp.RockDatabase`
   * - ``pycsamt.modeling.occam2d``
     - :mod:`pycsamt.models.occam2d`
   * - ``pycsamt.viewer.plot``
     - Module-level ``plot_*()`` methods on result objects
   * - ``pycsamt.ff.processing``
     - :mod:`pycsamt.emtools`

----

.. _changelog-1-2-1:

1.2.1 |Docs|
------------

Last stable v1 release. See ``README_v1.md`` for the full v1 changelog.
