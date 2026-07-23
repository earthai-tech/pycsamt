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

.. _changelog-2-0-0:

2.0.0 |Feature| |API Change| |Breaking|
---------------------------------------

*First stable release of the v2 line. It consolidates* :ref:`2.0.0rc1
<changelog-2-0-0rc1>` *and* :ref:`2.0.0rc2 <changelog-2-0-0rc2>` *— see those
entries for the full v2 feature set — plus the changes below. The public API
is now stable; v1 users should read* :ref:`migration` *for the v1 → v2 name
map.*

Fixed
~~~~~

* |Fix| |Tests| **Python 3.9 Qt-toolbar segfault** — stub matplotlib's Qt
  navigation toolbar under the offscreen platform so the CI *interfaces* shard
  no longer crashes on the Python 3.9 / matplotlib 3.9.x combination.

Changed
~~~~~~~

* |Docs| **PyPI summary** rescoped to the v2 feature set — "Scientific Python
  for electromagnetic geophysics — processing, inversion, AI agents, and apps".

----

.. _changelog-2-0-0rc2:

2.0.0rc2 |Fix| |Build| |Docs|
-----------------------------

*Second pre-release of the v2 line, published to PyPI and TestPyPI for
community testing ahead of the 2.0.0 final. Bug fixes, packaging, and docs
only — no public API change since* :ref:`2.0.0rc1 <changelog-2-0-0rc1>`.

Fixed
~~~~~

* |Fix| |Tests| **Python 3.9 interpreter segfaults** — bounded runaway
  matplotlib figure accumulation across the test session (a per-test
  ``plt.close("all")``) and stubbed matplotlib's Qt navigation toolbar under
  the offscreen platform. Both crashed the process on the Python 3.9 /
  matplotlib 3.9.x combination the CI *interfaces* shard runs.
* |Fix| **SEG parsing and plotting** (:mod:`pycsamt.seg`) — consistent mixin
  MRO discovery, ``>=``-prefixed section-name normalisation, DMS hemisphere
  parsing, and trend carry-through on topography resample.

Changed
~~~~~~~

* |Build| **Lean distributions** — the source distribution no longer bundles
  the documentation tree (which pushed it past PyPI's 100 MB per-file limit)
  or the test suites; the sdist drops ~153 MB → ~4 MB and the wheel
  5.8 MB → ~4.4 MB, with every runtime resource retained.
* |Docs| **Hosted-applications status** — the Applications cards now reflect
  the in-progress hosted rollout and link each app's installation guide
  rather than promising a live instance.

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
* |Fix| **Qt-free headless imports** — matplotlib Qt backends are imported
  lazily, so CI runners and documentation builds never pull Qt in.
* |Fix| **Docs sidebar handler crash on full rebuilds** — pages whose HTML
  context carries ``meta=None`` no longer break the ``html-page-context``
  hook (``'NoneType' object has no attribute 'pop'``).

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
* |Docs| **Documentation moved to** `pycsamt.org <https://pycsamt.org>`__ —
  hosted on Netlify with cached incremental builds (only changed gallery
  examples re-execute). Canonical URLs, the version switcher, README
  badges, and every in-app *Documentation* link now point at the new
  domain; the Read the Docs configuration is retired and its site remains
  as the v1 legacy archive.
* |Docs| **scikit-learn-style landing page** — full-bleed hero carousel
  with a rotating survey-method keyword and clickable workflow strip; the
  six capability cards show *real package output* (noise removal, forward
  responses, a ModEM section beside a learned inversion, pseudosection and
  stratigraphic fence, pipeline timings with an agent dry-run chain, and a
  QC coverage audit), each regenerable via ``scripts/home_card_*.py``; the
  "Code in action" panel is a themed editor window in light and dark mode.
* |Docs| **Navigation polish** — the right "On this page" sidebar is back
  on ordinary pages (only the home page, the API reference tables, and the
  gallery index pages stay full-width); Map Tools and Tutorials indexes
  gained icon card grids; Site Tools uses a two-column grid; the desktop and
  web application entries point directly to their deep guides; header icon
  links now include the issue tracker and Stack Overflow.
* |Build| **Faster documentation builds** — third-party module
  highlighting disabled in ``viewcode``, the Sphinx environment and the
  executed gallery persist between Netlify deploys, and the docs
  environment installs CPU-only torch so the AI examples run without
  CUDA wheels.
* |Build| **Test-coverage lift** — offline agent contract batteries,
  TDEM/AVG parser suites, and map topography/export tests; coverage now
  omits externally sourced shims and static configuration modules
  (``compat/``, ``_typing.py``, ``config.py``, ``projection.py``).
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
