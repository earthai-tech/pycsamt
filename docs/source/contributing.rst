.. _contributing:

Contributing
============

pyCSAMT is developed in the open at
`github.com/earthai-tech/pycsamt <https://github.com/earthai-tech/pycsamt>`_,
and contributions of every kind move it forward: bug reports with a
reproducible trace, documentation fixes, new gallery examples, field
wisdom about a data format, or code. You do not need to be a
geophysicist *and* a software engineer — either half is valuable.

.. grid:: 1 1 2 3
   :gutter: 3
   :class-container: cta-tiles

   .. grid-item-card:: Report an issue
      :link: https://github.com/earthai-tech/pycsamt/issues
      :img-top: _static/icons/stackoverflow.svg
      :class-card: pycsamt-card sd-text-center

      Found a bug, a confusing result, or a gap in the docs? The issue
      tracker is the front door — see the checklist below.

   .. grid-item-card:: Improve the documentation
      :link: development/documentation_build
      :link-type: doc
      :img-top: _static/icons/doc-build.svg
      :class-card: pycsamt-card sd-text-center

      Docs and gallery examples are code too: build them locally,
      fix a page, or add an executed example.

   .. grid-item-card:: Contribute code
      :link: development/index
      :link-type: doc
      :img-top: _static/icons/developer-notes-icon-tools.svg
      :class-card: pycsamt-card sd-text-center

      Set up a development install and follow the workflow below;
      the development section covers style, CI, and API policy.

Development Setup
-----------------

Clone the repository and install in editable mode with the development
extras:

.. code-block:: bash

   git clone https://github.com/earthai-tech/pycsamt.git
   cd pycsamt
   pip install -e ".[full,docs,dev]"

Verify the install by running a fast subset of the test suite:

.. code-block:: bash

   python -m pytest pycsamt/emtools/tests -q

Development Workflow
--------------------

#. **Branch** — create a feature branch from the active development
   branch rather than committing to it directly.

#. **Change with tests** — every subpackage keeps its tests next to the
   code (``pycsamt/<subpackage>/tests/``). A bug fix comes with a test
   that fails before the fix; a feature comes with tests that document
   its contract.

#. **Lint and test locally** — the same checks CI runs:

   .. code-block:: bash

      python -m ruff check pycsamt/
      python -m pytest pycsamt/ -v --tb=short

#. **Build the docs when they are affected** — see
   :ref:`contributing-docs` below.

#. **Open a pull request** — small and focused beats large and mixed.
   Describe what changed and why, mention the affected subpackages, and
   include before/after screenshots for documentation or application UI
   changes.

Code Style
----------

* Linting is enforced with `ruff <https://docs.astral.sh/ruff/>`_; the
  configuration lives in ``pyproject.toml``, so ``python -m ruff check``
  picks it up automatically.
* Public functions, classes, and modules follow the numpydoc
  conventions documented in :doc:`development/docstring_style`.
* What counts as public API — and what a change to it requires — is
  defined in :doc:`development/api_policy`.
* New classes inherit from the package base objects
  (``PyCSAMTObject``, ``CoreObject``, ``MTBase``) rather than plain
  ``object`` — :doc:`development/extending` explains what each base
  provides.
* Match the style of the file you are editing: naming, comment density,
  and import grouping are kept locally consistent.

Running The Tests
-----------------

The full suite:

.. code-block:: bash

   python -m pytest pycsamt/ -v --tb=short

One subpackage or one file while iterating:

.. code-block:: bash

   python -m pytest pycsamt/iot/tests -q
   python -m pytest pycsamt/emtools/tests/test_emtools_tensor.py -q

Tests should not write inside the repository — use ``tmp_path`` or the
system temp directory for generated files.

.. _contributing-docs:

Documentation Contributions
---------------------------

The docs build with Sphinx from ``docs/``:

.. code-block:: bash

   cd docs
   pip install -r requirements-docs.txt
   python -m sphinx -M html source build
   # open build/html/index.html

Two frequent contributions have a well-worn footprint:

* **Gallery examples** — every page under :doc:`Examples </examples/index>`
  is an executed script. Add a ``plot_*.py`` to an existing section
  under ``docs/examples/<section>/``, or create a new section with its
  own ``README.txt`` and register it in ``subsection_order`` in
  ``docs/source/conf.py``. Scripts run at build time, so they must be
  self-contained and use bundled data or synthetics.
* **Guide pages** — task-based narrative lives under
  ``docs/source/user_guide/``; the :doc:`development/documentation_build`
  page covers the build environment and layout in detail.

Continuous Integration
----------------------

Every push and pull request runs the GitHub Actions pipeline: a
Python 3.9–3.12 test matrix on Linux with a ruff lint step, the pytest
suite, and coverage upload to Codecov. Releases are published to PyPI
automatically from version tags. :doc:`development/ci` describes the
pipeline and how to reproduce it locally.

Commit Style
------------

* Descriptive imperative subject line (≤ 72 characters), with the body
  explaining *why* when it is not obvious.
* Group related changes into one commit rather than many fragments.
* Co-author line:
  ``Co-Authored-By: earthai-tech <earthai-tech@users.noreply.github.com>``.

Reporting Issues
----------------

Use the `GitHub issue tracker
<https://github.com/earthai-tech/pycsamt/issues>`_. A report that can be
acted on quickly includes:

* the pyCSAMT version
  (``python -c "import pycsamt; print(pycsamt.__version__)"``),
  Python version, and operating system;
* the smallest command or snippet that reproduces the problem —
  ideally against the bundled sample data (``data/``);
* the full traceback or the wrong output, and what you expected
  instead;
* for data-format issues: a minimal anonymised sample file when
  sharing is possible.

Feature requests are welcome through the same tracker — describe the
workflow you are trying to achieve rather than only the API you expect.

License And Citation
--------------------

pyCSAMT is distributed under the LGPL-3.0 license; contributions are
accepted under the same terms. If pyCSAMT supports published work,
please cite the project — see :doc:`references` for the citable papers.
