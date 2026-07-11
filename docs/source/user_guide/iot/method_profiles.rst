.. _user_guide_iot_method_profiles:

Method-Aware QC
===============

An :class:`~pycsamt.iot.monitoring.EMMethod` names the survey method, but on
its own it changes no threshold. The :mod:`pycsamt.iot.methods` module
carries the canonical, per-method knowledge that turns a method name into
concrete QC defaults, so the same monitor behaves correctly for AMT, MT,
CSAMT, CSEM, and TDEM/TEM without hand-tuning.

Method profiles
---------------

A :class:`~pycsamt.iot.MethodProfile` records the typical acquisition band,
the channels a valid record must carry, a nominal sample rate, and whether
the method uses a controlled source or is sensitive to powerline noise:

.. code-block:: python

   from pycsamt.iot import method_profile

   p = method_profile("csamt")
   p.frequency_band_hz      # (0.125, 8192.0)
   p.required_channels      # ('ex', 'ey', 'hx', 'hy')
   p.controlled_source      # True
   p.powerline_sensitive    # True

.. list-table::
   :header-rows: 1
   :widths: 14 22 26 18 20

   * - Method
     - Band (Hz)
     - Required channels
     - Source
     - Powerline
   * - AMT
     - 1 – 10 000
     - ex, ey, hx, hy
     - natural
     - sensitive
   * - MT
     - 1e-4 – 1000
     - ex, ey, hx, hy, hz
     - natural
     - sensitive
   * - CSAMT
     - 0.125 – 8192
     - ex, ey, hx, hy
     - controlled
     - sensitive
   * - CSEM
     - 0.01 – 100
     - ex, ey, hx, hy
     - controlled
     - sensitive
   * - TDEM / TEM
     - (time-domain)
     - hz
     - controlled
     - not sensitive

The band edges are representative defaults, not hard standards; override
them per survey when acquisition parameters differ.

Seeding a monitor
-----------------

:meth:`~pycsamt.iot.MonitoringConfig.for_method` seeds the expected band
and required channels from a profile, so the method-mismatch,
missing-channel, and out-of-band checks in
:class:`~pycsamt.iot.TelemetryMonitor` become method-aware with no manual
tuning. Keyword overrides win over the profile defaults:

.. code-block:: python

   from pycsamt.iot import MonitoringConfig, assess_telemetry

   cfg = MonitoringConfig.for_method("csamt", min_battery_v=11.5)
   status = assess_telemetry(packets, config=cfg)
   # an MT packet, or one missing hx/hy, now raises method_mismatch /
   # required_channels_missing against the CSAMT expectations.

Method-aware edge QC
--------------------

Passing ``method`` to :func:`~pycsamt.iot.amt_edge_report` gates the
diagnostics by method. Powerline-harmonic detection runs only for
powerline-sensitive methods -- it is skipped for TDEM/TEM, where the report
carries ``powerline_applicable=False`` instead of a misleading "clean"
verdict -- and frequency coverage is scored against the method's target
bands:

.. code-block:: python

   from pycsamt.iot import amt_edge_report, target_bands_for_method

   report = amt_edge_report(window, sample_rate, method="amt")
   report["powerline_applicable"]                     # True for AMT
   report["frequency_coverage"]["coverage_fraction"]  # scored vs AMT bands

   target_bands_for_method("csamt")
   # [(0.125, 1.0), (1.0, 10.0), (10.0, 100.0), (100.0, 1000.0), (1000.0, 8192.0)]

Passing no method preserves the original, method-agnostic behaviour.
