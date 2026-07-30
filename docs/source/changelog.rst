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

The full log is split by major-version series so each page stays a
reasonable length. Entries accumulate during development as small files
under ``docs/changelog.d/`` and are assembled into the series page for each
release by ``docs/scripts/changelog_release.py`` — see
``docs/changelog.d/README.rst`` for the contributor workflow and the
"Changelog workflow" section of
:doc:`/development/documentation_build` for the release-time steps.

Latest release
--------------

.. _changelog-latest:

**2.2.0** — *2026-07-30* — a genuine 3-D Maxwell training-data
pipeline, an ``Inv3DAgent`` ``physics="mt3d"`` mode, a fixed
layered-earth benchmark for the research-only 3-D solver, and
``ModEm3DAdapter`` physics-validated against a real compiled ModEM
binary for the first time, plus cross-platform build tooling for the
external solvers.
:ref:`Full 2.2.0 entry <changelog-2-2-0>` · :ref:`Upgrade guidance
<release_v2_2_0>`.

.. toctree::
   :maxdepth: 1
   :caption: By series

   changelog/v2
   changelog/v1
