.. _development-documentation-build:

Documentation Build
===================

This page explains how to build, inspect, and maintain the pyCSAMT v2
documentation.  It is written for contributors who need to edit the Sphinx
source files, generate the API reference from docstrings, diagnose warnings,
and prepare the documentation for continuous integration.

The documentation stack is intentionally close to the scientific Python
ecosystem:

* Sphinx for the documentation engine;
* pydata-sphinx-theme for the HTML theme;
* autodoc and autosummary for API reference pages;
* numpydoc and napoleon for NumPy-style docstrings;
* MyST for Markdown support;
* sphinx-copybutton and sphinx-design for reader-friendly pages.


Source layout
-------------

The documentation source lives under ``docs/source``.

.. code-block:: text
   :linenos:

   docs/
     Makefile
     requirements-docs.txt
     source/
       conf.py
       index.rst
       getting_started/
       tutorials/
       user_guide/
       agents/
       pipeline/
       cli/
       theory/
       api/
       development/
       release_notes/

Important files:

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - File
     - Purpose
   * - ``docs/source/conf.py``
     - Sphinx configuration, extensions, theme options, autodoc settings, and
       docs-build environment variables.
   * - ``docs/requirements-docs.txt``
     - Minimal documentation dependency list.
   * - ``pyproject.toml``
     - Project extras, including the ``docs`` optional dependency group.
   * - ``docs/Makefile``
     - Convenience commands for HTML, PDF, clean builds, and live rebuilds.
   * - ``docs/source/api/``
     - API reference pages and autosummary entry points.
   * - ``docs/source/development/``
     - Contributor rules for API policy, docstrings, and docs builds.


Install documentation dependencies
----------------------------------

From a fresh checkout, install the package and documentation extras.

.. code-block:: bash
   :linenos:

   cd /home/daniel/projects/pycsamt
   python -m pip install -e ".[docs]"

Alternative, using the docs requirements file:

.. code-block:: bash
   :linenos:

   cd /home/daniel/projects/pycsamt
   python -m pip install -e .
   python -m pip install -r docs/requirements-docs.txt

For broad local development, the full contributor setup may be useful:

.. code-block:: bash
   :linenos:

   python -m pip install -e ".[dev,docs]"

Optional features such as AI backends, GIS export, and app support may require
additional extras.  The docs must still import cleanly when optional runtime
packages are absent.


Build the HTML documentation
----------------------------

The normal local build is:

.. code-block:: bash
   :linenos:

   cd /home/daniel/projects/pycsamt/docs
   make html

The generated site is written to:

.. code-block:: text
   :linenos:

   docs/build/html/index.html

The equivalent direct Sphinx command is:

.. code-block:: bash
   :linenos:

   cd /home/daniel/projects/pycsamt/docs
   sphinx-build -b html source build/html


Clean builds
------------

Use a clean build when changing ``conf.py``, autosummary settings, toctrees,
theme assets, or generated API pages.

.. code-block:: bash
   :linenos:

   cd /home/daniel/projects/pycsamt/docs
   make clean
   make html

Direct command:

.. code-block:: bash
   :linenos:

   cd /home/daniel/projects/pycsamt/docs
   sphinx-build -E -a -b html source build/html

Flags:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Flag
     - Meaning
   * - ``-E``
     - Ignore the saved Sphinx environment and reread all source files.
   * - ``-a``
     - Write all output files, even if Sphinx thinks they are unchanged.
   * - ``-b html``
     - Build the HTML target.


Selected-page builds
--------------------

When editing one page, use a selected-page build for faster feedback.

.. code-block:: bash
   :linenos:

   cd /home/daniel/projects/pycsamt/docs
   sphinx-build -q -b html source /tmp/pycsamt-docs-html \
       source/development/documentation_build.rst

This still loads project configuration and may still trigger some autosummary
work, but it is faster than a full rebuild and good for checking reST syntax.


Long source-code examples
-------------------------

Keep short interactive examples visible with ``code-block:: pycon`` so readers
can see the input and captured output together. For a long, reusable function
stored in a project script, use the local ``code-dropdown`` directive instead
of exposing a full ``literalinclude`` unconditionally:

.. code-block:: rst

   .. code-dropdown:: ../../../scripts/generate_ai_inversion_figures.py
      :language: python
      :pyobject: make_losses_penalty_anatomy
      :linenos:
      :title: View penalty-anatomy source code

The HTML presentation is collapsed initially and uses native
``details``/``summary`` elements. It therefore remains keyboard accessible and
operates without JavaScript. The complete highlighted source stays in the page
for copying, searching, and printing; non-HTML builders display the source
normally rather than hiding it.

``code-dropdown`` delegates source handling to Sphinx's ``literalinclude``
implementation. It accepts the normal extraction options, including
``language``, ``pyobject``, ``lines``, ``start-after``, ``end-before``,
``linenos``, ``lineno-start``, ``emphasize-lines``, ``caption``, and ``name``.
It adds two presentation options:

``title``
   Replaces the default ``View source code`` label with a concise description.

``open``
   Expands the panel initially. Reserve this for code that must be visible on
   first reading; long supporting scripts should normally remain collapsed.

Do not use the dropdown to conceal required outputs, warnings, equations, or
the few lines a reader must understand before continuing. It is intended for
complete plotting functions and extended workflows that support the nearby
explanation but would otherwise interrupt its flow.

After a selected-page build, inspect the rendered file:

.. code-block:: bash
   :linenos:

   test -f /tmp/pycsamt-docs-html/development/documentation_build.html
   grep -n '</html>' /tmp/pycsamt-docs-html/development/documentation_build.html


Live rebuilds
-------------

For writing prose-heavy pages, use live rebuilds.

.. code-block:: bash
   :linenos:

   cd /home/daniel/projects/pycsamt/docs
   make livehtml

This uses ``sphinx-autobuild`` from the docs dependency set.  It watches source
files and serves the site locally.  Use it for authoring, not for CI, because
it is interactive and long-running.


Strict builds
-------------

A strict build treats warnings as failures.

.. code-block:: bash
   :linenos:

   cd /home/daniel/projects/pycsamt/docs
   sphinx-build -W --keep-going -b html source build/html

Use strict builds before releases and after large API documentation changes.
During the active v2 documentation migration, the project may still have known
global warnings.  In that phase, contributors should still run local builds and
avoid introducing new warnings in the files they touch.

Recommended release target:

.. code-block:: bash
   :linenos:

   cd /home/daniel/projects/pycsamt/docs
   sphinx-build -W --keep-going -E -a -b html source build/html


Docs build environment
----------------------

``docs/source/conf.py`` sets a small number of environment defaults:

.. code-block:: python
   :linenos:

   os.environ.setdefault("PYCSAMT_DOCS_BUILD", "1")
   os.environ.setdefault("MPLCONFIGDIR", "/tmp/pycsamt-matplotlib")

``PYCSAMT_DOCS_BUILD`` lets package code detect documentation builds and avoid
side effects that do not belong in API generation.  ``MPLCONFIGDIR`` keeps
matplotlib cache files out of user directories and project source folders.

Documentation imports should be side-effect-light:

* do not start web servers at import time;
* do not call LLM providers at import time;
* do not download model checkpoints at import time;
* do not require PyTorch, TensorFlow, GDAL, rasterio, or geopandas unless the
  documented object actually executes that feature;
* do not write logs or generated data outside the build tree or temporary
  directories.


Autosummary and generated API pages
-----------------------------------

pyCSAMT enables:

.. code-block:: python
   :linenos:

   autosummary_generate = True

This means Sphinx can generate stub pages for modules listed in autosummary
directives.  The API source pages live in ``docs/source/api``.  Generated pages
are commonly written under ``docs/source/api/generated``.

When adding a new public module:

1. Add the public object to the correct ``__all__``.
2. Add or improve the object docstring using :ref:`development-docstring-style`.
3. Add the module to the relevant API page in ``docs/source/api``.
4. Build the docs and inspect the generated page.
5. Fix warnings caused by the new module before merging.

Do not make private underscored modules the main documented entry point unless
the page is intentionally describing internal architecture for contributors.


Docstring warnings
------------------

Most API warnings come from docstrings parsed by numpydoc.  The most common
patterns are custom section headings such as ``Input Keys``, ``Output Data
Keys``, ``Resolution Rules``, or ``Recognised Variables``.

Use the guidance in :ref:`development-docstring-style`:

* keep recognized NumPy sections at the top level;
* move custom schema details into ``Notes``;
* document agent input and output mappings under ``Notes``;
* avoid fragile reST tables inside docstrings;
* keep section underline lengths correct.

Example fix:

.. code-block:: text
   :linenos:

   Notes
   -----
   Execution input mapping:

   ``sites`` : Sites
       Survey object to process.
   ``path`` : path-like
       Directory or file used when ``sites`` is not supplied.

This avoids a numpydoc ``Unknown section Input Keys`` warning while preserving
the useful information.


Intersphinx warnings
--------------------

When the machine has no network access, Sphinx may warn that external
inventories cannot be fetched:

.. code-block:: text
   :linenos:

   intersphinx inventory 'https://docs.python.org/3/objects.inv'
   not fetchable

This is expected in offline or restricted environments.  It does not usually
mean the local documentation page is broken.  Before release or CI, run the
build in an environment that can fetch inventories or configure cached
inventories.


GDAL and optional dependency warnings
-------------------------------------

During imports, optional geospatial support may emit a warning such as:

.. code-block:: text
   :linenos:

   GDAL_DATA is unavailable. Geospatial features will degrade.

This is acceptable for a base docs build if the affected geospatial features
are not executed.  API pages should still import and render.  Pages that
document GIS execution should mention the required GIS extras.

Optional dependency failures are not acceptable when simply importing a public
module for documentation.  Move heavy imports inside functions or methods when
needed.


Citation warnings
-----------------

Sphinx may warn about unreferenced citations if docstrings or generated API
pages contain bibliography entries that are not cited in prose.

Preferred policy:

* centralize project references in ``docs/source/references.rst``;
* cite references from theory, tutorials, or user guide pages;
* avoid long bibliography blocks inside API docstrings;
* keep API docstrings focused on behavior, inputs, outputs, and assumptions.


Build outputs and cleanup
-------------------------

Local build outputs:

.. code-block:: text
   :linenos:

   docs/build/html/
   docs/build/doctrees/

Temporary smoke-build outputs often use ``/tmp``:

.. code-block:: text
   :linenos:

   /tmp/pycsamt-docs-html
   /tmp/pycsamt-docs-html-api-policy
   /tmp/pycsamt-docs-html-docstring-style

Clean project build outputs with:

.. code-block:: bash
   :linenos:

   cd /home/daniel/projects/pycsamt/docs
   make clean


Recommended contributor workflow
--------------------------------

For most documentation edits:

.. code-block:: bash
   :linenos:

   cd /home/daniel/projects/pycsamt/docs
   sphinx-build -q -b html source /tmp/pycsamt-docs-html \
       source/path/to/page.rst

Then inspect the generated page and run a broader build before finishing a
large section:

.. code-block:: bash
   :linenos:

   make html

For API/docstring edits:

.. code-block:: bash
   :linenos:

   cd /home/daniel/projects/pycsamt/docs
   sphinx-build -q -b html source /tmp/pycsamt-docs-html-api \
       source/api/index.rst

For release preparation:

.. code-block:: bash
   :linenos:

   cd /home/daniel/projects/pycsamt/docs
   sphinx-build -W --keep-going -E -a -b html source build/html


Continuous integration
----------------------

The CI documentation job should eventually perform a strict clean build.

Recommended CI steps:

.. code-block:: yaml
   :linenos:

   - name: Install package with docs dependencies
     run: python -m pip install -e ".[docs]"

   - name: Build documentation
     working-directory: docs
     run: sphinx-build -W --keep-going -E -a -b html source build/html

During the v2 migration, CI may temporarily use a non-strict build while known
legacy warnings are cleaned:

.. code-block:: yaml
   :linenos:

   - name: Build documentation
     working-directory: docs
     run: sphinx-build --keep-going -E -a -b html source build/html

The long-term target should be warning-free strict builds.


Changelog workflow
-------------------

The changelog is split by major-version series under
``docs/source/changelog/`` (``v2.rst``, ``v1.rst``, …); ``changelog.rst``
is just an index page with a toctree over those files. During development,
contributors don't edit a series file directly -- each user-facing PR adds
a small fragment under ``docs/changelog.d/`` instead (see
``docs/changelog.d/README.rst`` for the exact filename convention and
type/badge/section table). This avoids merge conflicts on a shared file and
enforces the badge taxonomy by filename rather than editor discipline.

At release time:

.. code-block:: bash
   :linenos:

   python docs/scripts/changelog_release.py --version 2.2.0 --date 2026-08-15

This assembles every pending fragment into a new version block at the top
of the matching ``docs/source/changelog/vN.rst`` and deletes the consumed
fragments. Then, by hand:

1. Edit the ``TODO`` release-summary placeholder the script inserted.
2. Write the narrative ``docs/source/release_notes/vX.Y.Z.rst`` page.
3. If the script created a brand-new series file, add it to the toctree in
   ``docs/source/changelog.rst``.


Pre-release checklist
---------------------

Before publishing documentation:

.. code-block:: text
   :linenos:

   [ ] docs/source/index.rst links to every major section.
   [ ] All toctrees are reachable.
   [ ] API pages include all public packages intended for the release.
   [ ] New public APIs have NumPy-style docstrings.
   [ ] Agent pages link to relevant API pages.
   [ ] Pipeline pages link to step and preset APIs.
   [ ] Tutorials use public imports only.
   [ ] Optional dependency requirements are documented.
   [ ] No new numpydoc warnings are introduced.
   [ ] No unresolved internal references remain.
   [ ] Strict Sphinx build passes in CI.
   [ ] The generated homepage opens at docs/build/html/index.html.
   [ ] The newest docs/source/changelog/vN.rst version section has a
       "*Released YYYY-MM-DD.*" line (drives the site's "New" badge --
       an entry without one is treated as unreleased).
   [ ] docs/source/release_notes/vX.Y.Z.rst exists for that same version
       and is linked from release_notes/index.rst as "(latest)". The
       "New" badge's banner links directly to release_notes/vX.Y.Z.html
       regardless of release size -- a changelog entry without a matching
       page is a live 404, not just a missing narrative page.


Troubleshooting
---------------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Symptom
     - What to check
   * - ``Unknown section`` warnings
     - Move custom docstring headings under ``Notes``.  See
       :ref:`development-docstring-style`.
   * - Import fails during autodoc
     - Check optional dependencies and move heavy imports inside runtime
       functions.
   * - A page is missing from navigation
     - Add it to the nearest ``index.rst`` toctree.
   * - A generated API page is stale
     - Run a clean build with ``sphinx-build -E -a``.
   * - Cross-reference is unresolved
     - Add a label to the target page or use the fully qualified Python object
       path.
   * - Intersphinx fetch fails
     - Check network access or ignore during offline local builds.
   * - Matplotlib writes to user cache directories
     - Confirm ``MPLCONFIGDIR`` is set by ``docs/source/conf.py``.
   * - Theme assets do not update
     - Clean the build directory and rebuild.


In short
--------

Use selected-page builds while writing, full HTML builds before finishing a
section, and strict clean builds before release.  Keep imports lightweight,
move docstring schemas into recognized NumPy sections, and treat every new
warning as a sign that the documentation contract needs attention.
