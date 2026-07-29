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

.. When cutting a release, add a "*Released YYYY-MM-DD.*" line directly
   under that version's title (before the summary paragraph). It is the
   single source of truth the "New" badge on the site navbar reads at
   build time (see conf.py: _write_whats_new_json) -- an entry with no
   such line is treated as unreleased and never shows a badge. See the
   Pre-release checklist in development/documentation_build.rst.

----

.. _changelog-2-0-1:

2.0.1 |Feature| |Fix| |API Change| |Docs| |Tests|
-------------------------------------------------

*Maintenance release for deterministic, spatially correct station ordering
and validated, auditable station-metadata editing, plus the first
physics-grounded building blocks of AI inversion. See*
:ref:`release_v2_0_1` *for upgrade guidance and migration examples.*

Added
~~~~~

* |Feature| **Staged inversion loss package** -- added
  :mod:`pycsamt.ai.losses` (:mod:`~pycsamt.ai.losses.model`,
  :mod:`~pycsamt.ai.losses.spatial`, :mod:`~pycsamt.ai.losses.response`,
  :mod:`~pycsamt.ai.losses.boundary`, :mod:`~pycsamt.ai.losses.uncertainty`)
  implementing the staged objective ``L = w_m*L_model + lambda_x*L_grad_x +
  lambda_z*L_grad_z + lambda_tv*L_TV + lambda_d*L_response``.
* |Feature| **Scientific validation package** -- added
  :mod:`pycsamt.ai.validation` (:mod:`~pycsamt.ai.validation.recovery`,
  :mod:`~pycsamt.ai.validation.residuals`,
  :mod:`~pycsamt.ai.validation.calibration`,
  :mod:`~pycsamt.ai.validation.ood`) for recovery metrics, response
  residuals, predictive calibration, and out-of-distribution screening.
* |Feature| **2-D Maxwell training-data pipeline** -- added
  :func:`pycsamt.ai.training.dataset2d.generate_2d_maxwell_dataset` and
  :class:`~pycsamt.ai.training.dataset2d.Maxwell2DDatasetConfig`, generating
  spatially correlated 2-D geological realizations and solving them with
  :class:`~pycsamt.forward.maxwell.mt2d.MT2DAdapter` into a versioned,
  realization-split, resumable-cache dataset.
* |API Change| **Inv2DAgent physics mode** -- added ``physics="mt2d"`` to
  :class:`pycsamt.agents.Inv2DAgent`, training on the new dataset generator
  with the staged spatial-regularization loss and reporting a held-out
  recovery check. ``physics="mt1d"`` (tiled 1-D forward models) remains the
  unchanged default. :meth:`pycsamt.ai.inversion.inv2d.EMInverter2D.fit`
  gained ``lambda_x``, ``lambda_z``, and ``lambda_tv`` staged-loss weights
  (PyTorch backend only).
* |API Change| **Global station-ordering policy** -- added
  :func:`pycsamt.api.configure_ordering`,
  :func:`pycsamt.api.reset_ordering`, ``PYCSAMT_ORDERING``, and
  :class:`pycsamt.api.SiteOrderingConfig`. Configure ``mode="auto"`` once and
  APIs using :func:`pycsamt.emtools.ensure_sites` inherit it.
* |API Change| **Sites ordering API** -- added
  :meth:`pycsamt.site.Sites.ordered` with ``auto``, ``chainage``, ``input``,
  natural ``station``, ``latitude``, and ``longitude`` strategies, plus
  :attr:`pycsamt.site.Sites.ordering` diagnostics.
* |API Change| **Transactional metadata API** -- added
  :class:`pycsamt.site.metadata.SiteMetadataEditor`,
  :class:`pycsamt.site.metadata.MetadataChange`,
  :func:`pycsamt.site.metadata.rename_sites`,
  :func:`pycsamt.site.metadata.update_metadata`, and
  :func:`pycsamt.site.metadata.update_metadata_all`. The editor supports
  copy-on-write staging, dry-run plans, batch validation, atomic in-place
  commits, audit records, and metadata-aware export.
* |API Change| **Container metadata conveniences** -- added
  :meth:`pycsamt.site.Site.update_metadata`,
  :meth:`pycsamt.site.Sites.update_metadata`, and
  :meth:`pycsamt.site.Sites.rename`. Renaming changes station identity without
  reordering the collection or changing station coordinates.
* |Tests| **Spatial regression coverage** -- added synthetic edge cases and
  real-data tests for the bundled L18PLT and L22PLT survey lines, including a
  combined multiple-line collection.
* |Tests| **Metadata regression coverage** -- added mapping, sequence,
  callable, DataFrame, and CSV sources; nested ``HEAD``/``INFO`` changes;
  coordinate and duplicate-name validation; error policies; atomic commits;
  planning, auditing, export/reload; and real L18PLT and L22PLT workflows.

Fixed
~~~~~

* |Fix| **Oblique survey lines** -- order is derived from both geographic
  coordinates and projected chainage instead of longitude, latitude, or
  lexical station names alone.
* |Fix| **Unsafe geometry guesses** -- ``auto`` validates coordinate coverage,
  linearity, and cross-track spread and preserves input order when the spatial
  evidence is insufficient.
* |Fix| **Multiple-line interleaving** -- separated parallel profiles are
  ordered independently.
* |Fix| **Pseudosection station order** -- dataframe pivots no longer replace
  canonical profile order with lexical column order.
* |Fix| **Processing consistency** -- static-shift and near-surface methods,
  field-zone and CS/AMT sections, strike profiles, and pipeline presets now
  inherit the shared ordering strategy by default.
* |Fix| **2-D mesh/receiver air-layer mismatch** -- the 2-D training-data
  mesh builder no longer relies on ``build_solver_mesh``'s air layers, which
  were incompatible with :class:`~pycsamt.forward.maxwell.mt2d.MT2DAdapter`
  (receivers must sit exactly at ``z=0``) and produced spurious 393x errors
  despite an apparently converged solve.
* |Fix| **Depth-zero float precision** -- ``Maxwell2DDatasetConfig``'s
  surface-depth check now uses a numerical tolerance instead of strict
  equality, so floating-point-derived grid spacings no longer raise
  spuriously.
* |Fix| **BatchNorm crash on a trailing batch of size 1** --
  :meth:`pycsamt.ai.inversion.inv2d.EMInverter2D.fit`'s PyTorch training
  loop now drops a trailing size-1 batch instead of crashing with "Expected
  more than 1 value per channel when training".
* |Fix| **TM-mode interface averaging** --
  :func:`pycsamt.forward.em2d._assemble_tm` now uses a thickness-weighted
  harmonic mean, not an arithmetic one, for the resistivity coefficient at
  a cell interface, matching the parallel-current-path physics there.
  ``Maxwell2DDatasetConfig`` generates and trains on both the TE-mode
  (``zxy``) and TM-mode (``zyx``) response by default again.
* |Fix| **Field-survey impedance units** --
  :func:`pycsamt.ai.domain_gap.survey_fit.survey_data_from_sites` now
  converts :attr:`~pycsamt.site.base.Site.z` from the EDI-native
  ``[mV/km]/[nT]`` field convention to the SI convention
  :class:`~pycsamt.ai.data.contracts.SurveyData` declares, instead of
  copying it unconverted.

Changed
~~~~~~~

* |API Change| :mod:`pycsamt.ai.domain_gap.willy_fit` is renamed to
  :mod:`pycsamt.ai.domain_gap.survey_fit` -- the field-survey-to-domain-gap
  bridge works for any AMT, CSAMT, or MT survey, not only the bundled WILLY
  line.
* |API Change| Processing ``sort_by=None`` and
  :func:`pycsamt.emtools.ensure_sites` ``order_by=None`` now mean "use the
  global policy". Explicit arguments still override it; select ``input`` to
  retain loader order.
* |Docs| Added configuration, migration, override, threshold, reset, and
  compatibility guidance for the ordering policy. Reworked the complete Site
  Tools guide with executable ``pycon`` transcripts, captured outputs, labeled
  equations, generated figures, and dedicated metadata guidance including the
  L22PLT rename-versus-order comparison.
* |API Change| :meth:`pycsamt.site.Sites.map` remains a callable mapper and
  does not accept a dictionary as a rename table. Use
  :meth:`pycsamt.site.Sites.rename` for explicit station-name mappings.

Docs & tooling
~~~~~~~~~~~~~~

* |Docs| **Code of Conduct and Netlify attribution** -- added a project
  Code of Conduct and a "powered by Netlify" link in the footer of every
  documentation page, meeting the requirements of Netlify's Open Source
  Plan.
* |Build| **Release-gated docs deploys** -- ``pycsamt.org`` now rebuilds
  only on a ``vX.Y.Z`` tag push (a GitHub Actions workflow calling a
  Netlify build hook), rather than on every commit to ``master``. The
  previous per-commit pattern exhausted the hosting team's free-plan
  credits and suspended the site.
* |Docs| **AI-inversion user guide, written in full** --
  :doc:`/user_guide/ai_inversion/roadmap`,
  :doc:`/user_guide/ai_inversion/data_contracts`,
  :doc:`/user_guide/ai_inversion/experiments`,
  :doc:`/user_guide/ai_inversion/forward_physics`,
  :doc:`/user_guide/ai_inversion/dataset2d`,
  :doc:`/user_guide/ai_inversion/domain_gap`,
  :doc:`/user_guide/ai_inversion/losses`, and
  :doc:`/user_guide/ai_inversion/scientific_validation` are now complete
  guides covering :mod:`pycsamt.ai.data`, :mod:`pycsamt.ai.domain_gap`,
  :mod:`pycsamt.ai.experiments`, :mod:`pycsamt.ai.geology`,
  :mod:`pycsamt.ai.losses`, :mod:`pycsamt.ai.training.dataset2d`,
  :mod:`pycsamt.ai.validation`, and :mod:`pycsamt.forward.maxwell`
  end to end, with real, externally captured code output, generated
  figures, and labeled equations throughout.
  :doc:`/user_guide/ai_inversion/geology_priors` remains a stub.
* |Tests| **AI-inversion regression coverage** -- added tests for
  :mod:`pycsamt.ai.losses`, :mod:`pycsamt.ai.validation`,
  :mod:`pycsamt.ai.training.dataset2d` (including analytic half-space
  regression checks for both the TE- and TM-mode response),
  :func:`pycsamt.forward.em2d._assemble_tm`'s interface-averaging
  coefficients, the staged-loss ``EMInverter2D`` fit path, and the
  ``Inv2DAgent`` ``physics="mt2d"`` end-to-end path.

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
