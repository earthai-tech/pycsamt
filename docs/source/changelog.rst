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

**2.3.0** — *2026-08-14* — adds a new :mod:`pycsamt.geology` package
(lithology classification, pluggable rock-property providers, and
structural-geology primitives), split out of :mod:`pycsamt.interp` with
backward-compatible re-exports. Deep-fills the lithology, petrophysics, and
monitoring pages of the interpretation user guide, fixing several stale
docstring examples and a real ``RockDatabase`` iteration bug along the way.
Running a real legacy CSAMT line end to end through AVG-to-EDI conversion
and the processing pipeline surfaced and fixed six further bugs, plus a
packaging gap that dropped the fallback EPSG table from installed wheels.
Also adds four method-aware pipeline presets with real near-field
correction and data-driven QC, a branded dashboard report, and data-driven
power-line-harmonic detection, alongside a full rewrite of the
config-driven-pipeline tutorial and CLI reference that surfaced three
further bugs. :ref:`Full 2.3.0 entry <changelog-2-3-0>` · :ref:`Release
notes <release_v2_3_0>`.

.. toctree::
   :maxdepth: 1
   :caption: By series

   changelog/v2
   changelog/v1
