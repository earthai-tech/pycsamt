.. _coverage_todo:

Coverage improvement TODO
=========================

Purpose
-------

This TODO is a working guide for improving pyCSAMT's coverage in a way that
future developers can maintain.  The goal is not to chase a percentage by
touching easy lines.  The goal is to make the existing deterministic tests
visible in CI, add explicit tests for public behavior, and protect scientific
contracts that matter to users.

Treat coverage as a diagnostic.  A useful test should normally assert one of
these things:

* a scientific invariant, such as finite positive static-shift factors;
* a public API contract, such as stable columns, result fields, or exit codes;
* a regression path, especially for parsing, validation, and numerical edge
  cases;
* an offline workflow behavior that CI can run without credentials, external
  services, or developer-machine paths.

When starting a new coverage pass, first establish the current project
baseline from the branch being worked on, then improve it with focused,
reviewable slices.  Do not hard-code historical percentages into test goals:
the baseline will change as CI starts running more of the existing suite.

Definition of done
------------------

* Tests live beside their package in ``pycsamt/<package>/tests/test_*.py``.
* The default suite is deterministic, offline, and excludes ``live`` tests.
* Every new test checks an observable result; import-only and call-only tests
  are insufficient unless they protect a public API contract.
* Numerical tests use justified tolerances and small synthetic fixtures.
* Filesystem tests use ``tmp_path`` and do not write into the repository.
* Plot tests use a non-interactive backend and inspect returned figures,
  axes, artists, or saved output rather than image pixels where possible.
* Optional-dependency behavior is tested both when unavailable (mocked) and,
  in a suitable CI job, when installed.
* Coverage is reported with branches and missing lines, then uploaded from
  one canonical Python version.

Phase 0 -- establish a trustworthy baseline
--------------------------------------------

Priority: immediate.

.. checklist items intentionally use normal bullets so they render in all
   supported Sphinx themes.

* [ ] Save the current XML/JSON coverage artifact as a baseline.
* [x] Run ``python -m pytest --collect-only -m "not live"`` and resolve
  collection errors, duplicate module names, and unintended skips.
  Missing ``__init__.py`` files in seven ``tests`` packages made duplicate
  basenames (``test_view.py``, ``test_forward.py``, ``test_backends.py``)
  abort collection for the whole suite; all tests dirs are packages now.
* [ ] Inventory tests by package and classify them as unit, integration,
  optional-dependency, GUI, large-data, or live.
* [ ] Confirm that tests using WILLY data resolve repository-relative paths
  and do not depend on a developer machine.
* [x] Add markers for slow, integration, GUI, and optional-dependency tests.
* [x] Remove the overlap between ``.coveragerc`` and ``pyproject.toml`` or
  document which file is authoritative.  ``.coveragerc`` is authoritative
  (it wins the config-discovery order); the duplicate ``[tool.coverage.*]``
  tables were removed from ``pyproject.toml``.  Test code and ``conftest``
  files no longer count toward coverage.

Phase 1 -- expose existing tests in CI
-------------------------------------

Priority: highest return for the least new code.

* [x] Replace the two-path pytest command with package groups or the full
  ``pycsamt`` test tree, always using ``-m "not live"``.
* [x] Split CI into fast unit shards so failures remain attributable.  A
  practical initial grouping is:

  * science: ``emtools``, ``site``, ``seg``, ``z``, ``forward``;
  * workflows: ``agents``, ``assistant``, ``pipeline``, ``ai``;
  * formats/models: ``api``, ``models``, ``inversion``, ``zonge``, ``jones``;
  * interfaces: ``cli``, ``app``, ``map``, ``iot``, ``tdem``;
  * support: ``core``, ``compat``, ``gis``, ``interp``, ``topo``, ``utils``.

* [x] Combine shard coverage data before generating ``coverage.xml``; do not
  upload several incomplete reports as though each represented the project.
* [x] Upload shard coverage artifacts with ``if: always()`` and run the
  combine job even when a shard is red: a failing shard used to silently
  drop every one of its packages to 0% in the combined report.
* [x] Guard the suite with ``pytest-timeout`` (``--timeout=300``) and force
  the headless platform (``QT_QPA_PLATFORM=offscreen``, ``MPLBACKEND=Agg``)
  in CI so a hung GUI test cannot stall a shard.  Two such hangs existed:
  the ``MainWindow`` quit-confirmation ``QMessageBox`` and the global
  ``QApplication.setStyleSheet`` restyle; both are neutralized by autouse
  fixtures in ``pycsamt/app/tests/conftest.py``.
* [ ] Keep GUI, torch, GIS, and other heavy optional stacks in dedicated jobs
  if they make the fast suite fragile.
* [ ] Publish test counts, skipped-test reasons, and coverage artifacts on
  every pull request.

Milestone: existing offline tests are collected in CI and the reported score
reflects what is actually exercised.

Phase 2 -- protect the core scientific workflows
------------------------------------------------

Priority: correctness before broad line coverage.

``pycsamt/emtools/tests``
~~~~~~~~~~~~~~~~~~~~~~~~~

* [ ] Add ``test_core_loader.py`` for ``ensure_sites`` with a ``Sites``
  instance, a single EDI file, a directory, recursive discovery, invalid
  input, and EDI files without valid impedance.
* [ ] Add explicit QC contract tests for ``build_qc_table`` and ``qc_flags``:
  empty inputs, missing components, NaNs, low coverage, low SNR, stable column
  names, and deterministic station ordering.  Initial synthetic coverage now
  checks stable columns, invalid-Z skipping, low coverage, and low SNR.
* [ ] Add ``test_static_shift_ama.py`` for ``estimate_ss_ama``,
  ``correct_ss_ama``, and ``apply_ss_factors``.  Assert finite positive
  factors, unchanged phase, expected resistivity scaling, station alignment,
  rejection of NaN/zero/negative factors, and the valid empty-table result on
  strongly 3-D data.
* [ ] Test ``plot_ss_summary`` and ``plot_ss_1d_curves`` with synthetic data,
  existing axes, empty results, and non-interactive saving.  Initial synthetic
  tests now cover both plotting functions under the non-interactive backend.

``pycsamt/site/tests``, ``pycsamt/seg/tests``, and ``pycsamt/z/tests``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* [ ] Test impedance shape/frequency validation, component availability,
  derived resistivity and phase, error propagation, and frequency ordering.
* [ ] Test ``Sites`` selection, station-name stability, iteration, copying,
  serialization, and mixed valid/invalid EDI collections.
* [ ] Add EDI parse/write/parse round trips using minimal fixtures and assert
  numerical values plus metadata, not only successful parsing.

``pycsamt/forward/tests`` and ``pycsamt/ai/tests``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* [ ] Test MT1D limiting cases, layer validation, output shapes, finite
  responses, reproducibility, and known synthetic models.
* [ ] Test ``EMInverter1D`` preprocessing, fit/predict contracts, invalid
  dimensions, reproducible seeds, serialization, and backend-unavailable
  messages with tiny synthetic datasets.

Milestone: ``emtools``, ``site``, ``seg``, ``z``, and ``forward`` each reach a
meaningful branch baseline, with all safety invariants covered.

Phase 3 -- test agents and end-to-end orchestration
--------------------------------------------------

Priority: public workflow reliability.

``pycsamt/agents/tests``
~~~~~~~~~~~~~~~~~~~~~~~~

* [ ] Parameterize the shared keyword registry tests so
  ``ContextInputAgent`` and ``WorkflowOrchestratorAgent`` are checked against
  ``pycsamt.agents._workflows`` as the single source of truth.  Initial
  contract tests now exercise registry classification through both agents.
* [ ] Test ``MTLoaderAgent``, ``DataQCAgent``, ``StaticShiftAgent``,
  ``PhaseAnalysisAgent``, ``AIInversionAgent``, and ``CodeGenerationAgent``
  for successful results and actionable validation failures.
* [ ] Assert the ``AgentResult`` contract: ``status``, ``summary``, ``data``,
  warnings/errors, and output paths.  Initial dry-run and coordinator failure
  tests now assert status, summary, data, errors, hints, and cost.
* [ ] Add orchestration tests for load -> QC, load -> QC -> static shift, and
  load -> QC -> phase analysis using small fixtures and ``tmp_path``.
* [ ] Verify failure propagation and cleanup when a middle workflow step
  declines or fails.  Required-step abort and optional-step continuation are
  now covered with fake agents.

``pycsamt/assistant/**/tests`` and ``pycsamt/pipeline/tests``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* [ ] Test project-registry survey-line resolution instead of hard-coded
  paths, including unknown and ambiguous lines.
* [ ] Test ``validate_generated_code`` against real and invented imports.
* [ ] Test offline RAG ingest/retrieve/context assembly with a tiny temporary
  index; keep remote embeddings and LLM calls under ``live``.
* [ ] Test session/project memory isolation and workflow trace persistence.
* [ ] Test pipeline resume, partial output, overwrite policy, and deterministic
  manifests.

Milestone: every documented high-level workflow has one successful integration
test and at least one failure-path test.

Phase 4 -- expand format, model, CLI, and UI coverage
----------------------------------------------------

Priority order follows risk and statement count.

* [ ] ``models/tests``: configuration validation, mesh/data/result
  round trips, malformed files, unit conversion, and backend failure modes.
* [ ] ``inversion/tests``: retain the strong existing API coverage and add
  branch tests for invalid observations, bounds, convergence status, and
  backend selection.
* [ ] ``api/tests``: request validation, serialization, status/error mapping,
  and temporary-file handling without starting external services.
* [ ] ``cli/tests``: use the CLI runner for help, invalid options, successful
  small workflows, exit codes, and machine-readable output.  Initial root
  interface contracts now cover registered command groups, help/version,
  unknown commands, JSON config output, and forward-model error output.
* [ ] ``app/tests``: isolate controllers/models from Qt or Dash rendering;
  exercise callbacks and state transitions with mocked workers.  Put true GUI
  smoke tests in a dedicated job.
* [ ] ``iot/tests`` and ``tdem/tests``: packet/schema validation, numerical
  edge cases, provenance, transforms, and offline store/forward behavior.
* [ ] ``zonge/tests`` and ``jones/tests``: malformed records, optional blocks,
  missing values, and parse/write/parse round trips.

Milestone: no user-facing package remains at 0% merely because its tests were
excluded from CI; critical interfaces have explicit contract tests.

Phase 5 -- cover support packages selectively
---------------------------------------------

* [ ] Add focused tests for ``core``, ``compat``, ``gis``, ``interp``,
  ``topo``, ``transformers``, and ``utils`` based on public behavior and known
  regressions.  Initial ``utils`` contracts now cover unit conversion,
  dataframe cleaning, imputation, constant-axis pruning, NaN filling, and
  missing-value alignment.
* [ ] Avoid spending effort on trivial constants, typing-only modules,
  generated UI assets, and unreachable platform guards solely to raise the
  percentage.
* [x] Exclude third-party-derived and deprecated modules from the coverage
  denominator (see ``.coveragerc``): ``utils/funcutils.py`` and
  ``utils/exmath.py`` are BSD-3 collections imported from watex, and
  ``utils/em.py`` is v1-legacy code slated for deprecation.  These are not
  tested directly; behavior needed by v2 code should migrate into focused
  modules with their own tests.
* [ ] Add ``pragma: no cover`` only for a documented, genuinely unreachable
  branch; do not use it to hide untested business logic.
* [ ] Review exclusions that currently ignore broad ``raise ...Error`` lines:
  error paths should remain visible when they represent user-facing behavior.

Phase 6 -- stabilize the new coverage baseline
----------------------------------------------

* [x] Run the new focused slices together:
  ``emtools/tests/test_core_qc_contracts.py``,
  ``agents/tests/test_step3_orchestration_contracts.py``,
  ``cli/tests/test_step4_interface_contracts.py``, and
  ``utils/tests/test_step5_utils_contracts.py``.
* [ ] Run one CI-style shard locally at a time and fix collection/runtime
  failures before adding a hard coverage gate.
* [ ] Generate ``coverage.xml`` and ``coverage.json`` from the combined
  offline suite, then save the result as the first trustworthy baseline.
* [ ] Compare the new report against the original 10% baseline and list the
  packages that remain at or near 0% because they still lack executable
  offline tests.
* [ ] Add a small non-blocking coverage summary to pull requests before
  enforcing a minimum percentage.
* [ ] Only after the baseline is stable, set the first gate near the observed
  value rather than jumping directly to the long-term 25% milestone.

Coverage gates
--------------

Raise gates only after the baseline is stable.  Suggested project milestones
are 25%, 40%, 55%, and then 65%, with small ratcheting increases that never
allow a pull request to reduce the established baseline.  In addition to the
global score, track the core packages separately so a well-tested small module
cannot hide an untested scientific workflow.

For each pull request:

* changed production lines should normally be covered;
* new public APIs require success, boundary, and failure tests;
* bug fixes require a regression test that fails without the fix;
* coverage drops need an explanation, even before coverage is a hard gate.

Recommended local commands
--------------------------

.. code-block:: bash

   # Fast, offline collection check
   python -m pytest pycsamt --collect-only -m "not live"

   # Full offline suite with branch coverage
   python -m pytest pycsamt -m "not live" \
       --cov=pycsamt --cov-branch \
       --cov-report=term-missing --cov-report=xml --cov-report=json

   # A package while developing focused tests
   python -m pytest pycsamt/emtools/tests -m "not live" \
       --cov=pycsamt.emtools --cov-branch --cov-report=term-missing
