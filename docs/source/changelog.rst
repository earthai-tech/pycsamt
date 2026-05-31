.. _changelog:

Changelog
=========

v2.0.0 (in development)
------------------------

This is a major rewrite.  The public API is **not** backwards-compatible
with v1.  See :ref:`migration` below for a mapping from v1 to v2 names.

New features
~~~~~~~~~~~~

* **Multi-backend AI engine** — all neural-network models (FCN, ResNet,
  CNN1D, UNet2D, DRCNN) dispatch transparently to PyTorch **or**
  TensorFlow via a unified backend abstraction layer
  (:mod:`pycsamt.backends`).
* **Occam2D wrapper** — full Python interface to the OCCAM2D Fortran
  solver: data builder, mesh generator, runner, result loader, and
  four plot classes (:mod:`pycsamt.models.occam2d`).
* **ModEM wrapper** — 2-D and 3-D ModEM Python interface
  (:mod:`pycsamt.models.modem`).
* **pycsamt.interp** — new geological interpretation package that
  supersedes ``geodrill``: method-agnostic model calibration, 25-rock
  EM database, Oasis Montaj / LAS 2.0 / VTK export
  (:mod:`pycsamt.interp`).
* **AI inversion** — 1-D (EMInverter1D), 2-D (EMInverter2D), and
  multi-modal joint inversion (JointInverter)
  (:mod:`pycsamt.ai.inversion`).
* **EDI as source of truth** — the :class:`~pycsamt.z.z.Z` object and
  :class:`~pycsamt.io.EDICollection` now drive all processing and
  inversion pipelines.
* **PyPI-installable** — pure ``pyproject.toml``; ``setup.py`` removed.

Breaking changes
~~~~~~~~~~~~~~~~

* ``pycsamt.geodrill`` — removed; use :mod:`pycsamt.interp`.
* ``pycsamt.ff`` — removed; use :mod:`pycsamt.emtools`.
* ``pycsamt.viewer`` — removed; use module-level ``.plot()`` methods.
* ``pycsamt.modeling`` — removed; use :mod:`pycsamt.models`.
* All v1 CLI entry points — removed.

.. _migration:

Migration from v1
-----------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - v1
     - v2 equivalent
   * - ``pycsamt.geodrill.geocore.GeoStratigraphy``
     - :class:`pycsamt.interp.ModelCalibrator`
   * - ``pycsamt.geodrill.geocore.Geodrill``
     - :class:`pycsamt.interp.ResistivityModel` + :class:`~pycsamt.interp.ModelCalibrator`
   * - ``pycsamt.geodrill.geodatabase.GeoDataBase``
     - :class:`pycsamt.interp.RockDatabase`
   * - ``pycsamt.modeling.occam2d``
     - :mod:`pycsamt.models.occam2d`
   * - ``pycsamt.viewer.plot``
     - Module-level ``plot_*()`` methods on result objects
   * - ``pycsamt.ff.processing``
     - :mod:`pycsamt.emtools`

v1.2.1 (2023)
-------------

Last stable v1 release.  See ``README_v1.md`` for the v1 changelog.
