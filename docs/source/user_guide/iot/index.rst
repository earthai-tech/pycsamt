.. _user_guide_iot:

IoT-Enabled Field Acquisition
=============================

The :mod:`pycsamt.iot` subpackage implements the :term:`IoT` layer that
documents and audits :term:`IoT`-enabled AMT, MT, CSAMT, and CSEM
acquisition. It does not replace the normal impedance, quality-control,
dimensionality, or inversion workflow. Instead, it records the
operational layer around those data: devices, stations, packets, edge
QC decisions, power, clock synchronisation, security configuration, and
provenance. It also bridges back into the science API, turning edge
impedance estimates into a :class:`pycsamt.z.z.Z` and a preliminary EDI,
and seeding a re-occupation session from an existing EDI survey.

Use this guide to configure and audit the operational layer of an
:term:`IoT`-enabled field survey: telemetry, edge QC, monitoring, power,
synchronisation, security, and provenance.

.. toctree::
   :maxdepth: 3
   :class: pycsamt-guide-toc

   basic_session
   edge_qc
   amt_diagnostics
   controlled_source
   method_profiles
   data_bridge
   telemetry
   monitoring
   visualization
   power_management
   clock_sync
   provenance
   security
   simulation
