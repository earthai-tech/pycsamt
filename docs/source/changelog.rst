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

**2.4.0** — *2026-08-18* — the station confidence profile
(:func:`~pycsamt.emtools.qc.plot_confidence_profile`,
:func:`~pycsamt.emtools.qc.station_confidence_table`,
:func:`~pycsamt.emtools.qc.frequency_confidence_table`) now computes real
inter-station distance from EDI coordinates via a dependency-free UTM
projection, instead of silently defaulting every station to a hardcoded
200 m spacing whenever east/north attributes were absent — which they
always are for ordinary EDI-backed ``Site`` objects. Adds an opt-in
``force_spacing`` parameter to bypass coordinates entirely for surveys
with unreliable positioning, and ``annotate_low_step`` to auto-declutter
low-confidence point labels on surveys where most stations are flagged.
Thanks to `@shahidalishah130-hub
<https://github.com/shahidalishah130-hub>`__ for the report
(`#76 <https://github.com/earthai-tech/pycsamt/issues/76>`__).
:ref:`Full 2.4.0 entry <changelog-2-4-0>` · :ref:`Release notes
<release_v2_4_0>`.

.. toctree::
   :maxdepth: 1
   :caption: By series

   changelog/v2
   changelog/v1
