.. _contributing:

Contributing
============

pyCSAMT is developed in the open at
`github.com/earthai-tech/pycsamt <https://github.com/earthai-tech/pycsamt>`_,
and contributions of every kind move it forward: bug reports with a
reproducible trace, documentation fixes, new gallery examples, field wisdom
about a data format, or code. You do not need to be a geophysicist *and* a
software engineer — either half is valuable.

This page is the front door: how to set up, how a change travels from your
clone to a merged pull request, and how to report a problem. It deliberately
does not restate the rules it links to — the rest of this section is the
reference, and each page below goes deeper than a summary here could:

.. list-table::
   :header-rows: 1
   :widths: 38 62

   * - When you need
     - Read
   * - Docstring conventions for anything public
     - :doc:`development/docstring_style`
   * - What counts as public API, and what changing it costs
     - :doc:`development/api_policy`
   * - Which base class a new object should inherit
     - :doc:`development/extending`
   * - Building the docs, and where each page lives
     - :doc:`development/documentation_build`
   * - What CI runs, and how to reproduce a failure locally
     - :doc:`development/ci`

Development Setup
-----------------

Clone the repository and install in editable mode with the development
extras:

.. code-block:: bash

   git clone https://github.com/earthai-tech/pycsamt.git
   cd pycsamt
   pip install -e ".[full,docs,dev]"

Verify the install with a fast subset of the suite:

.. code-block:: bash

   python -m pytest pycsamt/emtools/tests -q

The Contribution Workflow
-------------------------

#. **Branch** — create a feature branch from the active development branch
   rather than committing to it directly.

#. **Change with tests** — every subpackage keeps its tests beside the code
   (``pycsamt/<subpackage>/tests/``). A bug fix comes with a test that fails
   before the fix; a feature comes with tests that document its contract.
   Tests must not write inside the repository — use ``tmp_path`` or the system
   temp directory for generated files.

#. **Run the checks locally** — the same two commands CI runs first, so a
   failure costs you seconds instead of a round trip:

   .. code-block:: bash

      python -m ruff check pycsamt/
      python -m pytest pycsamt/ -v --tb=short

   While iterating, scope pytest to what you touched
   (``python -m pytest pycsamt/iot/tests -q``).
   :doc:`development/ci` covers the full matrix and how to reproduce any job.

#. **Build the docs if you changed them** — the layout, the build commands,
   and how to register a new gallery section are in
   :doc:`development/documentation_build`. Gallery examples execute at build
   time, so a new ``plot_*.py`` must be self-contained and run from bundled
   sample data or synthetics.

#. **Open a pull request** — small and focused beats large and mixed. Say what
   changed and why, name the affected subpackages, and include before/after
   screenshots for documentation or application UI changes.

Commit And Pull Request Conventions
-----------------------------------

* Imperative subject line (≤ 72 characters); use the body to explain *why*
  when it is not self-evident from the diff.
* Group related changes into one commit rather than many fragments.
* Co-author line:
  ``Co-Authored-By: earthai-tech <earthai-tech@users.noreply.github.com>``.

Reporting Issues
----------------

Use the `GitHub issue tracker
<https://github.com/earthai-tech/pycsamt/issues>`_. A report that can be acted
on quickly includes:

* the pyCSAMT version
  (``python -c "import pycsamt; print(pycsamt.__version__)"``), Python
  version, and operating system;
* the smallest command or snippet that reproduces the problem — ideally
  against the bundled sample data (``data/``);
* the full traceback or the wrong output, and what you expected instead;
* for data-format issues: a minimal anonymised sample file when sharing is
  possible.

Feature requests are welcome through the same tracker — describe the workflow
you are trying to achieve rather than only the API you expect.

License And Citation
--------------------

pyCSAMT is distributed under the LGPL-3.0 license; contributions are accepted
under the same terms. If pyCSAMT supports published work, please cite the
project — see :doc:`references` for the citable papers.
