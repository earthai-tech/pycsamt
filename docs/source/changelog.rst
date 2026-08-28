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

**2.5.0** — *2026-08-28* — a major interoperability release centered on two
complementary common-format families: PCSF/PCSM for subsurface models and
inversion results, and PCBH for spatial multi-borehole observations and
geology. Their explicit association supports reproducible model-plus-evidence
exchange and shared 3-D interpretation. The release also adds a survey-scale
confidence-ratio QC suite
with ten complementary plots, smart geographic contours, and CSV/Surfer
grid export. The release also fixes ``ModEmModel3D`` centre/rotation
retention, web 3-D topography uploads, and confidence-map/profile geometry,
with expanded PCSF and confidence-evaluation guides and reproducible
examples. Running real, full-scale Occam2D and ModEM 3-D inversions
against a real 128-station survey found and fixed four more real
solver-launch bugs, including a new ``ModEmForwardControl`` that finally
lets a 3-D covariance file reach Mod3DMT at all, and a Fortran
formatted-input decimal-point bug that had silently corrupted control
values (``target_rms`` included) in every prior 3-D run. The Map View's
resistivity modes were reworked for real inversion volumes: 3-D
resistivity/depth filters that mask instead of recolour, block and
iso-surface modes that render under topography at all, depth slices kept
inside the real model, a new inversion "Resistivity @ depth" plan-slice
over the basemap, opt-in volume smoothing, rotatable/thinnable station
labels, a Borehole Builder with shared web/MapView 3-D rendering, and a
Windows ``pycsamt-mapview`` logging-flood fix.
:ref:`Full 2.5.0 entry <changelog-2-5-0>` · :ref:`Release notes
<release_v2_5_0>`.

.. toctree::
   :maxdepth: 1
   :caption: By series

   changelog/v2
   changelog/v1
