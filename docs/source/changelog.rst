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

**2.5.2** — *2026-09-01* — a patch release. Converting a real,
independent third-party ModEM 3-D inversion to PCSF and exercising every
Map View mode against it exposed five pre-existing defects in the native
ModEM 3-D import path: above-topography air fill drove the automatic
colour scale; the positive-down ``.dat`` ``Z`` column was read as an
elevation, draping topography upside-down; depth slices were referenced
to the model top, not the ground surface; the block / iso-surface volume
came back empty when real multi-line station geometry produced duplicate
along-strike coordinates; and profile lines zig-zagged when
``InversionResult`` picked a 3-decimal ModEM ``-R`` coordinate echo over
the real input. ``modem3d_to_pcsf`` gains ``air_threshold_ohm_m`` and
``station_z_convention`` (with ``pycsamt format convert --air-threshold``
/ ``--station-z``) plus a quantized-coordinate warning; Map View gains
``crange_percentile`` / ``rho_display_max`` and their inspector controls,
re-references ``grid3d`` curtains to each station's own surface, trims the
deep boundary-condition padding, and de-duplicates the block/iso
interpolation axes; ``InversionResult`` now prefers a full-precision
input ``.dat`` over a rewrite echo.
:ref:`Full 2.5.2 entry <changelog-2-5-2>` · :ref:`Release notes
<release_v2_5_2>`.

**2.5.0** — *2026-08-28* — a major interoperability release centered on two
complementary common-format families: PCSF/PCSM for subsurface models and
inversion results, and PCBH for spatial multi-borehole observations and
geology. See the :ref:`full 2.5.0 entry <changelog-2-5-0>` and
:ref:`release notes <release_v2_5_0>`.

.. toctree::
   :maxdepth: 1
   :caption: By series

   changelog/v2
   changelog/v1
