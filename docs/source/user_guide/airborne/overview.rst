.. _user_guide_airborne_overview:

Airborne Data Model Overview
=============================

The core containers are technology-neutral:
:class:`~pycsamt.airborne.NavigationTrack` for one flight line's
sample-aligned position, attitude, and elevation;
:class:`~pycsamt.airborne.AirborneEMRecord` for one sample's
:class:`~pycsamt.emtf.EMTF` transfer-function payload plus any
auxiliary decoded fields; and
:class:`~pycsamt.airborne.AirborneEMLine`/
:class:`~pycsamt.airborne.AirborneEMDataset` for a line's, or a
survey's, records keyed by sample ID. See :doc:`data_model` for how
those four pieces fit together.
:class:`~pycsamt.airborne.site.AirborneSite`/
:class:`~pycsamt.airborne.site.AirborneSites` give a flat, ``Sites``-
shaped view over the same records, read straight from a directory of
EMTF-XML files via :func:`~pycsamt.airborne.site.ensure_asites` --
the airborne counterpart of :class:`~pycsamt.site.base.Site`/
:class:`~pycsamt.site.base.Sites`; see :doc:`site`. Around that data
model, :mod:`~pycsamt.airborne.registry` and :mod:`~pycsamt.airborne.io`
govern how a technology or a native file format is recognized and
dispatched (:doc:`registry_and_io`), and :mod:`~pycsamt.airborne.qc`
assesses structural completeness and metadata consistency without
inventing technology-specific signal thresholds (:doc:`quality_control`).

Three technology subpackages map decoded scientific arrays onto this
model, each with its own ``build_*_line``/``build_*_dataset``
constructors and a ``*SystemSpec`` describing the real instrument's
published characteristics:

* :mod:`pycsamt.airborne.afmag` -- AFMAG (comparator + AirMt), a
  scalar tilt angle or a :math:`(n_f, 3, 2)` interstation tensor.
* :mod:`pycsamt.airborne.ztem` -- ZTEM, a tipper,
  :math:`H_z = T_{zx}H_x + T_{zy}H_y`.
* :mod:`pycsamt.airborne.mobilemt` -- MobileMT, a :math:`(n_f, 3, 2)`
  admittance tensor, :math:`H = Y E`.

No proprietary vendor archive format is parsed anywhere in pyCSAMT --
every subpackage above only maps already-decoded arrays onto the
common model. Synthetic sample surveys for all three are committed
under ``data/ZTEM/``, ``data/AFMAG/``, and ``data/mobileMT/``, and are
used throughout this section's pages.

This section only covers the data model itself. For the science --
reading the sample surveys, computing tilt/divergence/admittance
diagnostics, and interpreting the resulting figures -- continue with
:ref:`emtools_afmag`, :ref:`emtools_ztem`, and :ref:`emtools_mobilemt`.
For the shared physics linking all three technologies, including why
MobileMT's admittance tensor is a genuinely different object from the
other two's tipper/interstation tensor, see :ref:`airborne_theory`.
For the complete callable reference, see :doc:`../../api/airborne`.
