Continuous Integration
======================

pyCSAMT uses GitHub Actions to keep the package installable, tested,
documentable, and releasable. The current automation exercises focused
map and inversion checks, builds the Sphinx documentation, validates
Python package artifacts, and publishes release artifacts to PyPI through
trusted publishing.

This page is for contributors who need to understand what runs on pull
requests, how to reproduce failures locally, and how to extend CI
coverage without making the workflow slow or fragile.

Workflow Map
------------

The repository currently has two GitHub Actions workflow files.

.. list-table::
   :header-rows: 1
   :widths: 24 24 52

   * - Workflow
     - File
     - Purpose
   * - CI
     - ``.github/workflows/ci.yml``
     - Runs linting, selected tests, coverage upload, documentation
       build, and package validation.
   * - Publish to PyPI
     - ``.github/workflows/python-publish.yml``
     - Builds release distributions, stores them as an artifact, and
       publishes them to PyPI from the protected ``pypi`` environment.

The CI workflow is the everyday quality gate. The publish workflow is a
release workflow and should be treated as part of the packaging and
distribution process.

Main CI Workflow
----------------

The main workflow is named ``CI`` and runs for:

* pushes to ``main``, ``master``, ``develop``, and ``v2``;
* pull requests;
* manual ``workflow_dispatch`` runs from GitHub.

The workflow has read-only repository contents permission:

.. code-block:: yaml

   permissions:
     contents: read

That keeps the default CI jobs conservative. They can read the source
tree, install dependencies, run checks, and upload coverage through the
configured action, but they do not receive broad write permissions.

Test Job
--------

The ``test`` job is the main compatibility check. It runs on
``ubuntu-latest`` against this Python matrix:

* Python 3.9
* Python 3.10
* Python 3.11
* Python 3.12

The matrix uses ``fail-fast: false``. A failure on one Python version
does not cancel the others, which is useful when a dependency, typing
change, or language-version behavior only affects one runtime.

Each matrix entry performs the same core steps:

1. Check out the repository with full history available to the job.
2. Set up the selected Python version and enable pip caching.
3. Upgrade ``pip``.
4. Install pyCSAMT in editable mode with the ``full`` extra.
5. Run the configured Ruff lint command.
6. Run the selected pytest suite, excluding tests marked ``live``.
7. Upload coverage to Codecov from the Python 3.12 job only.

The install command is:

.. code-block:: bash

   python -m pip install -e ".[full]"

The ``full`` extra expands to the project extras used by development,
documentation, app, agent, geospatial, and optional torch-backed
features. This makes CI heavier than a minimal install, but it is useful
because pyCSAMT spans scientific APIs, CLIs, documentation builds, and
apps.

Lint Scope
----------

The current lint command is deliberately scoped:

.. code-block:: bash

   python -m ruff check pycsamt/map pycsamt/inversion/mt

Ruff itself is configured in ``pyproject.toml``. Important settings
include:

* target Python version: ``py39``;
* line length: ``62``;
* enabled rule families: ``E``, ``F``, ``W``, ``I``, ``UP``, and ``B``;
* ignored rule: ``B008``;
* ``__init__.py`` ignores unused-import warnings used for public API
  re-export patterns.

When extending lint coverage, prefer adding one stable package at a
time. That makes failures easier to understand and avoids mixing a real
regression with a large historical cleanup.

Test Scope
----------

The current pytest command is:

.. code-block:: bash

   python -m pytest pycsamt/map/tests pycsamt/inversion/tests/test_inversion_api.py -m "not live" --cov=pycsamt --cov-report=term-missing --cov-report=xml

This checks the map package tests and the inversion API tests while
excluding tests marked ``live``. The ``live`` marker is defined in
``pyproject.toml`` for tests that require an LLM API key or another
external live service.

Default CI tests should be:

* deterministic;
* reasonably fast;
* independent of private credentials;
* independent of network availability unless the workflow explicitly
  exists to test an online integration;
* safe to run on pull requests from contributors.

If a test needs an API key, a live remote service, or large local data,
mark it as ``live`` or move it into a dedicated workflow. Do not let
default pull-request CI depend on secrets or changing external state.

Coverage Upload
---------------

Every Python matrix entry writes coverage output, including
``coverage.xml``. Only the Python 3.12 job uploads to Codecov:

.. code-block:: yaml

   if: matrix.python-version == '3.12'

The Codecov upload uses:

* ``coverage.xml`` as the explicit coverage file;
* the ``unit`` flag;
* a per-version name based on the matrix Python version;
* ``fail_ci_if_error: false``.

Because ``fail_ci_if_error`` is disabled, Codecov upload problems should
not block normal development. A coverage upload failure is still worth
investigating, especially on release branches, but it is not treated the
same as a lint, test, docs, or packaging failure.

Docs Job
--------

The ``docs`` job verifies that the Sphinx documentation builds. It runs
on ``ubuntu-latest`` with Python 3.12.

The job installs the project with:

.. code-block:: bash

   python -m pip install -e ".[full]"

Then it builds HTML documentation:

.. code-block:: bash

   python -m sphinx -b html docs/source docs/_build/html

Two environment variables are set for the docs build:

.. code-block:: yaml

   MPLCONFIGDIR: ${{ runner.temp }}/matplotlib
   PYCSAMT_DOCS_BUILD: "1"

``MPLCONFIGDIR`` keeps Matplotlib cache and configuration files inside
the GitHub Actions runner's temporary directory. ``PYCSAMT_DOCS_BUILD``
allows package code to detect documentation imports and avoid
interactive behavior or expensive side effects where necessary.

Package Job
-----------

The ``package`` job verifies that the project can produce valid Python
distribution artifacts. It runs on Python 3.12 and executes:

.. code-block:: bash

   python -m pip install --upgrade pip build twine
   python -m build
   python -m twine check dist/*

``python -m build`` creates source and wheel distributions in
``dist/``. ``twine check`` validates package metadata and long
description rendering before any release workflow attempts to publish.

Publish Workflow
----------------

The publish workflow is named ``Publish to PyPI``. It runs when:

* a GitHub release is published;
* a maintainer starts it manually with ``workflow_dispatch``.

The ``build`` job creates and validates distributions, then uploads the
``dist/`` directory as a GitHub Actions artifact named
``python-package-distributions``.

The ``publish`` job depends on ``build``. It downloads that exact
artifact and publishes it to PyPI using:

.. code-block:: yaml

   pypa/gh-action-pypi-publish@release/v1

The job is attached to the ``pypi`` environment and declares:

.. code-block:: yaml

   permissions:
     id-token: write

That is the permission used for PyPI trusted publishing through OpenID
Connect. The workflow does not need a long-lived PyPI API token checked
into repository secrets.

Local Reproduction
------------------

The fastest way to debug CI is to run the same command that failed.
Inside the project root, use the repository's development environment if
it is available:

.. code-block:: powershell

   conda activate pycsamt-v2
   python -m pip install -e ".[full]"

Run the CI lint command:

.. code-block:: powershell

   python -m ruff check pycsamt/map pycsamt/inversion/mt

Run the CI test command:

.. code-block:: powershell

   python -m pytest pycsamt/map/tests pycsamt/inversion/tests/test_inversion_api.py -m "not live" --cov=pycsamt --cov-report=term-missing --cov-report=xml

Build the docs locally:

.. code-block:: powershell

   $env:MPLCONFIGDIR = "$env:TEMP\pycsamt-matplotlib"
   $env:PYCSAMT_DOCS_BUILD = "1"
   python -m sphinx -b html docs/source docs/_build/html

Check packaging:

.. code-block:: powershell

   python -m pip install --upgrade build twine
   python -m build
   python -m twine check dist/*

For a quick check of this page, use Docutils:

.. code-block:: powershell

   python -m docutils --halt=warning docs/source/development/ci.rst NUL

Failure Triage
--------------

Install failures usually mean a dependency constraint changed upstream,
an optional dependency in ``full`` lacks a wheel for the selected Python
version, or project metadata in ``pyproject.toml`` is inconsistent. Read
the first resolver error before changing constraints.

Ruff failures are usually direct. Run the same ``ruff check`` command
locally and either fix the code or update the lint scope intentionally.
Prefer narrow per-file ignores over broad rule removals.

Test failures should be handled from the smallest failing test outward.
Run the specific failing test with ``-vv`` and then widen back to the CI
command once the local failure is understood.

Docs failures may be syntax errors, import-time errors, broken
references, missing optional dependencies, or examples that assume files
not present on the runner. If the failure happens while importing a
module, decide whether the import should be lighter during docs builds.

Package failures are often metadata or README-rendering issues. Run
``python -m build`` and ``python -m twine check dist/*`` locally before
changing the publish workflow.

Codecov upload failures are currently non-blocking. They should be
investigated when coverage disappears from pull requests, but they do
not necessarily indicate that tests failed.

Publish failures should be treated carefully. Check whether the
distribution artifact was built correctly, whether the ``pypi``
environment requires approval, and whether PyPI trusted publishing is
configured for the repository, workflow file, and environment.

Adding CI Coverage
------------------

When adding more code to CI, prefer small, explainable expansions.

Good candidates for default CI are:

* pure unit tests for scientific transformations;
* parser tests using small fixture files;
* CLI tests that use temporary directories and local sample data;
* docs pages that import stable public APIs;
* packaging checks that do not publish anything.

Poor candidates for default CI are:

* tests requiring private credentials;
* tests requiring a user-specific data directory;
* long-running inversion jobs over large survey datasets;
* tests depending on a live LLM, remote API, or mutable online service;
* tests that write outside the workspace.

For expensive or live checks, use a dedicated marker, a scheduled
workflow, or a manually triggered workflow. Keep the normal pull-request
path deterministic and fast enough that contributors trust it.

Release Checklist
-----------------

Before publishing a release, confirm that:

* the main CI workflow is passing on the release branch;
* docs build successfully;
* package validation succeeds with ``twine check``;
* version metadata in ``pyproject.toml`` is correct;
* release notes or changelog entries are ready;
* the GitHub release is intentional and approved;
* PyPI trusted publishing is configured for the repository and
  ``pypi`` environment.

The publish workflow builds the release artifacts itself. Do not assume
that a locally built ``dist/`` directory is what will be published.

Design Principles
-----------------

pyCSAMT CI should be practical. The best checks are the ones that catch
real regressions in workflows users rely on: loading data, transforming
survey products, running inversion APIs, generating figures, building
documentation, and producing installable packages.

As the project grows, CI should grow in layers:

* keep a fast pull-request path for everyday development;
* move slow scientific workflows into explicit jobs;
* mark live or credentialed tests clearly;
* preserve packaging and docs checks on important branches;
* publish only from controlled release events.

