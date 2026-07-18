.. _user_guide_iot_controlled_source:

Controlled-Source Edge QC (CSAMT / CSEM)
========================================

Natural-source :term:`AMT`/:term:`MT` diagnostics (:doc:`amt_diagnostics`)
assume a :term:`plane-wave field` with no operator-controlled
transmitter. :term:`Controlled-source` methods, especially
:term:`CSAMT` and :term:`CSEM`, add a :term:`grounded dipole transmitter`.
That changes the field questions at the edge: is the receiver far enough
from the source for a plane-wave interpretation, is there resolvable
energy at every transmitted line, and is the :term:`source current`
steady enough to trust the sounding? The
:mod:`pycsamt.iot.edge_csamt` and :mod:`pycsamt.iot.edge_csem` modules
answer these questions with numpy-only checks on short edge windows.

A ``transmitter`` device role and a ``source`` :term:`telemetry packet`
let a transmitter node report its state alongside the receivers. The
examples below use synthetic data so the same outputs can be reproduced
without an EDI file or a field logger export.

Build A Controlled-Source Window
--------------------------------

This first block creates one CSAMT-like receiver window and matching
transmitter telemetry. The expected :term:`transmitter frequency comb` is
8, 32, 128, and 512 Hz, but the synthetic signal only contains the first
three lines. That deliberate missing line makes the frequency-comb QC
visible. The CSEM amplitude curve also contains a small bump at 6000 m so
the monotonic-decay check has something real to flag.

.. code-block:: python
   :linenos:

   import numpy as np

   rng = np.random.default_rng(51)
   sample_rate = 2048.0
   n_samples = 8192
   t = np.arange(n_samples) / sample_rate

   tx_frequencies = np.array([8.0, 32.0, 128.0, 512.0])
   window = (
       0.18 * np.sin(2 * np.pi * 8.0 * t)
       + 0.32 * np.sin(2 * np.pi * 32.0 * t)
       + 0.24 * np.sin(2 * np.pi * 128.0 * t)
       + 0.035 * rng.standard_normal(n_samples)
   )

   tx_current = (
       9.8
       + 0.08 * np.sin(2 * np.pi * 0.6 * t)
       + 0.03 * rng.standard_normal(n_samples)
   )
   tx_voltage = 248.0 + 1.2 * rng.standard_normal(n_samples)

   offsets_m = np.array([1000, 2000, 4000, 6000, 8000, 10000],
                        dtype=float)
   amplitudes = np.array([8.0e-12, 4.1e-12, 2.1e-12,
                          2.7e-12, 7.0e-13, 2.0e-13])
   phases = np.array([-8.0, -16.0, -29.0, -44.0, -61.0, -79.0])

Field Zones From Skin Depth
---------------------------

The controlling CSAMT quantity is not offset alone, but the
:term:`transmitter-receiver offset` expressed in :term:`skin depth`
units. For resistivity :math:`\rho` in :math:`\Omega\,m` and frequency
:math:`f` in Hz, pyCSAMT uses

.. math::

   \delta(f) = 503.29\sqrt{\frac{\rho}{f}},
   \qquad q(f)=\frac{r}{\delta(f)}.

Here :math:`r` is the transmitter-receiver separation and :math:`q` is
the offset ratio. A frequency is labelled ``near`` when
:math:`q \le 1`, ``far`` when :math:`q \ge 3`, and ``transition`` between
those limits. High frequencies have smaller skin depth, so they usually
reach the far field first.

.. code-block:: python
   :linenos:

   from pycsamt.iot import classify_field_zones

   freqs = np.array([4096.0, 1024.0, 256.0, 64.0, 16.0, 4.0, 1.0])
   cov = classify_field_zones(
       freqs,
       resistivity=100.0,
       offset_m=5000.0,
   )

   for f, d, ratio, zone in zip(
       cov.freq_hz, cov.skin_depth_m, cov.offset_ratio, cov.zones
   ):
       print(
           f"{f:7.1f} Hz  skin={d:8.1f} m  "
           f"r/delta={ratio:5.2f}  zone={zone}"
       )
   print(f"all_far_field: {cov.all_far_field}")
   print(f"correction_recommended: {cov.correction_recommended}")
   print(f"first_far_field_hz: {cov.first_far_field_hz():.1f}")

Output:

.. code-block:: text

    4096.0 Hz  skin=    78.6 m  r/delta=63.58  zone=far
    1024.0 Hz  skin=   157.3 m  r/delta=31.79  zone=far
     256.0 Hz  skin=   314.6 m  r/delta=15.90  zone=far
      64.0 Hz  skin=   629.1 m  r/delta= 7.95  zone=far
      16.0 Hz  skin=  1258.2 m  r/delta= 3.97  zone=far
       4.0 Hz  skin=  2516.5 m  r/delta= 1.99  zone=transition
       1.0 Hz  skin=  5032.9 m  r/delta= 0.99  zone=near
   all_far_field: False
   correction_recommended: True
   first_far_field_hz: 16.0

The low-frequency rows are the ones to treat carefully. A
:term:`near-field correction` is recommended because the 4 Hz row is in
transition and the 1 Hz row is near field.

Transmitter Frequency Comb
--------------------------

CSAMT transmits a discrete frequency comb rather than the continuous
natural spectrum used by passive AMT/MT. For each expected line
:math:`f_k`, pyCSAMT compares the peak :term:`power spectral density`
inside a narrow band around :math:`f_k` with the median positive PSD
floor :math:`P_{50}`:

.. math::

   \operatorname{SNR}_k =
   10\log_{10}\left(\frac{\max P(f \approx f_k)}{P_{50}}\right).

A line is detected when that SNR reaches the configured threshold and the
line is below Nyquist.

.. code-block:: python
   :linenos:

   from pycsamt.iot import detect_transmitter_frequencies

   comb = detect_transmitter_frequencies(
       window,
       sample_rate=sample_rate,
       tx_frequencies=tx_frequencies,
       snr_threshold_db=8.0,
   )

   print(f"detected lines: {comb.n_detected} of {comb.n_expected}")
   print(f"missing lines: {comb.missing()}")
   for line in comb.lines:
       print(
           f"{line.frequency_hz:6.1f} Hz  "
           f"snr={line.snr_db:6.2f} dB  "
           f"detected={line.detected}"
       )

Output:

.. code-block:: text

   detected lines: 3 of 4
   missing lines: [512.0]
      8.0 Hz  snr= 25.89 dB  detected=True
     32.0 Hz  snr= 35.50 dB  detected=True
    128.0 Hz  snr= 32.98 dB  detected=True
    512.0 Hz  snr=  0.87 dB  detected=False

Source-Signal Stability
-----------------------

The transmitter current sets the signal level of every sounding, so its
steadiness bounds the receiver-side data quality. The source-stability
check first identifies the on-state samples as those above a fraction of
the peak absolute current. On those samples it reports the mean current
:math:`\bar I`, the :term:`coefficient of variation`

.. math::

   \mathrm{CV}_I = \frac{\operatorname{std}(I_\mathrm{on})}
                        {|\operatorname{mean}(I_\mathrm{on})|},

and a simple drift amplitude, estimated as the absolute slope of
on-state current versus sample index multiplied by the number of on-state
samples. A source is unstable when the current CV exceeds ``max_cv`` or
when no finite current is present.

.. code-block:: python
   :linenos:

   from pycsamt.iot import assess_source_stability

   status = assess_source_stability(tx_current, tx_voltage=tx_voltage)
   print(f"stable: {status.stable}")
   print(f"current_mean_a: {status.current_mean_a:.3f}")
   print(f"current_cv: {status.current_cv:.4f}")
   print(f"current_drift_a: {status.current_drift_a:.3f}")
   print(f"on_fraction: {status.on_fraction:.2f}")
   print(f"voltage_mean_v: {status.voltage_mean_v:.2f}")

Output:

.. code-block:: text

   stable: True
   current_mean_a: 9.810
   current_cv: 0.0066
   current_drift_a: 0.005
   on_fraction: 1.00
   voltage_mean_v: 247.98

CSEM: Magnitude/Phase Versus Offset
-----------------------------------

:term:`CSEM` records a dipole source with a :term:`receiver array`, and
its signature edge product is the response as a function of offset at
each frequency. :func:`~pycsamt.iot.field_vs_offset` builds the
:term:`magnitude-versus-offset` and :term:`phase-versus-offset` curve,
finds the :term:`detectability limit`, and checks that detectable
amplitude decays monotonically with offset.

For detectable amplitudes :math:`A_i`, the dynamic range is

.. math::

   D_\mathrm{dB} = 20\log_{10}
   \left(\frac{\max(A_i)}{\min(A_i)}\right).

The monotonic check allows a small tolerance :math:`\tau` and flags a
reading when :math:`A_{i+1} > A_i(1+\tau)`. A bump can mean a bad
receiver, a geometry error, or real 3-D structure worth revisiting.

.. code-block:: python
   :linenos:

   from pycsamt.iot import field_vs_offset

   resp = field_vs_offset(
       offsets_m=offsets_m,
       amplitudes=amplitudes,
       phases_deg=phases,
       noise_floor=1e-13,
       frequency_hz=1.0,
   )
   print(f"n_detectable: {resp.n_detectable} of {resp.n_offsets}")
   print(f"max_detectable_offset_m: {resp.max_detectable_offset_m:.0f}")
   print(f"monotonic_decay: {resp.monotonic_decay}")
   print(f"dynamic_range_db: {resp.dynamic_range_db:.2f}")
   print(f"above_noise: {resp.above_noise}")

Output:

.. code-block:: text

   n_detectable: 6 of 6
   max_detectable_offset_m: 10000
   monotonic_decay: False
   dynamic_range_db: 32.04
   above_noise: [True, True, True, True, True, True]

Transmitter Telemetry
---------------------

A transmitter node reports its state as a ``source`` packet, parsed by the
``SourcePayload`` schema with the same tolerant alias folding and range
validation as the other payloads. In the example below, ``current``,
``tx_voltage``, ``frequency``, ``ab_m``, and ``tx_rx_offset`` are all
normalised to canonical payload fields.

.. code-block:: python
   :linenos:

   from pycsamt.iot import DeviceConfig, FieldSession, parse_payload

   tx = DeviceConfig("tx-1", role="transmitter")
   payload = parse_payload("source", {
       "tx_id": "TX1",
       "current": 9.8,
       "tx_voltage": 250.0,
       "frequency": 32.0,
       "ab_m": 100.0,
       "tx_rx_offset": 5000.0,
   })
   session = FieldSession("CS1", method="csamt", devices=[tx])
   session.add_packet({
       "device_id": "tx-1",
       "timestamp": 10.0,
       "topic": tx.topic("source"),
       "kind": "source",
       "payload": payload.as_dict(),
   })

   print(
       "payload current/frequency/offset: "
       f"{payload.tx_current_a:.1f} A, "
       f"{payload.tx_frequency_hz:.1f} Hz, "
       f"{payload.offset_m:.0f} m"
   )
   print(f"session packets: {len(session.packets)}")

Output:

.. code-block:: text

   payload current/frequency/offset: 9.8 A, 32.0 Hz, 5000 m
   session packets: 1

Static Shift And Transmitter Timing Lock
----------------------------------------

Two further controlled-source concerns are worth keeping next to the
edge report. A galvanic :term:`static shift` multiplies
:term:`apparent resistivity` by a frequency-independent factor while
leaving :term:`phase` nearly unchanged. pyCSAMT estimates the per-period
log split between the ``xy`` and ``yx`` modes,

.. math::

   s_i = \log_{10}\rho_{xy,i} - \log_{10}\rho_{yx,i},
   \qquad G = 10^{\operatorname{median}(s_i)}.

The result is a static-shift candidate when the median split is large
enough, the standard deviation of :math:`s_i` is small enough, and the
mean phase difference stays below the phase threshold.

.. code-block:: python
   :linenos:

   from pycsamt.iot import estimate_static_shift

   periods = np.logspace(-3, 1, 12)
   res_yx = 80.0 * (1 + 0.08 * np.sin(np.linspace(0, np.pi,
                                                 periods.size)))
   res_xy = 3.0 * res_yx
   phi_yx = 42.0 + 1.5 * np.sin(np.linspace(0, 2 * np.pi,
                                           periods.size))
   phi_xy = phi_yx + 1.2

   ss = estimate_static_shift(
       res_xy,
       res_yx,
       phase_xy=phi_xy,
       phase_yx=phi_yx,
   )
   print(f"static_shift: {ss.static_shift}")
   print(f"shift_factor: {ss.shift_factor:.2f}")
   print(f"split_decades: {ss.split_decades:.3f}")
   print(f"consistency_std: {ss.consistency_std:.3f}")
   print(f"phase_diff_deg: {ss.phase_diff_deg:.2f}")

Output:

.. code-block:: text

   static_shift: True
   shift_factor: 3.00
   split_decades: 0.477
   consistency_std: 0.000
   phase_diff_deg: 1.20

A CSAMT/CSEM receiver can also report its :term:`transmitter timing lock`
alongside normal clock synchronisation, using the ``tx_locked``,
``tx_sync_offset_ms``, and ``tx_id`` fields of the ``sync`` payload.

.. code-block:: python
   :linenos:

   sync = parse_payload("sync", {
       "offset_ms": 0.4,
       "transmitter_locked": True,
       "tx_offset_ms": 0.2,
       "tx_id": "TX1",
   })
   print(f"tx_locked: {sync.tx_locked}")
   print(f"tx_sync_offset_ms: {sync.tx_sync_offset_ms:.1f}")

Output:

.. code-block:: text

   tx_locked: True
   tx_sync_offset_ms: 0.2

Plot The Controlled-Source Checks
---------------------------------

The same synthetic objects can be rendered as a compact four-panel field
summary: field-zone ratio, transmitter-comb SNR, source-current
stability, and the CSEM offset response. The figure is generated outside
the documentation build and stored in the normal IoT image folder.

.. code-block:: python
   :linenos:

   from pathlib import Path

   import matplotlib.pyplot as plt

   out_dir = Path("docs/source/images/user_guide/iot")
   out_dir.mkdir(parents=True, exist_ok=True)

   fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.2),
                            constrained_layout=True)

   zone_colors = {
       "near": "tab:red",
       "transition": "tab:orange",
       "far": "tab:green",
   }
   axes[0, 0].scatter(
       cov.freq_hz,
       cov.offset_ratio,
       c=[zone_colors[z] for z in cov.zones],
       s=56,
   )
   axes[0, 0].axhline(1.0, color="0.35", ls=":", lw=1.0)
   axes[0, 0].axhline(3.0, color="0.35", ls="--", lw=1.0)
   axes[0, 0].set_xscale("log")
   axes[0, 0].invert_xaxis()
   axes[0, 0].set_xlabel("Frequency (Hz)")
   axes[0, 0].set_ylabel("r / skin depth")
   axes[0, 0].set_title("CSAMT field zones")
   axes[0, 0].grid(alpha=0.25, which="both")

   axes[0, 1].bar(
       [line.frequency_hz for line in comb.lines],
       [line.snr_db for line in comb.lines],
       width=[line.frequency_hz * 0.18 for line in comb.lines],
       color=["tab:blue" if line.detected else "tab:red"
              for line in comb.lines],
   )
   axes[0, 1].axhline(8.0, color="0.25", ls="--", lw=1.0)
   axes[0, 1].set_xscale("log")
   axes[0, 1].set_xlabel("Expected line (Hz)")
   axes[0, 1].set_ylabel("SNR (dB)")
   axes[0, 1].set_title("Transmitter comb")
   axes[0, 1].grid(axis="y", alpha=0.25)

   axes[1, 0].plot(t[:1200], tx_current[:1200], color="tab:purple",
                   lw=1.0)
   axes[1, 0].axhline(status.current_mean_a, color="0.25", ls="--",
                      lw=1.0)
   axes[1, 0].set_xlabel("Time (s)")
   axes[1, 0].set_ylabel("Current (A)")
   axes[1, 0].set_title("Source-current stability")
   axes[1, 0].grid(alpha=0.25)

   axes[1, 1].semilogy(offsets_m, amplitudes, marker="o",
                       color="tab:brown")
   axes[1, 1].axhline(1e-13, color="0.25", ls="--", lw=1.0)
   axes[1, 1].set_xlabel("Offset (m)")
   axes[1, 1].set_ylabel("Amplitude")
   axes[1, 1].set_title("CSEM MVO with suspect bump")
   axes[1, 1].grid(alpha=0.25, which="both")

   fig.savefig(
       out_dir / "user-guide-iot-controlled-source-01.png",
       dpi=180,
   )

.. image:: ../../images/user_guide/iot/user-guide-iot-controlled-source-01.png
   :width: 100%

Aggregation
-----------

:func:`~pycsamt.iot.csamt_edge_report` and
:func:`~pycsamt.iot.csem_edge_report` collate the per-channel diagnostics,
and :func:`~pycsamt.iot.csamt_edge_table` /
:func:`~pycsamt.iot.csem_edge_table` flatten them into pyCSAMT tables for
reporting. In practice, keep the source packet, comb result, field-zone
coverage, and CSEM offset response together in the field manifest so the
later interpretation can see both the measured response and the source
state that produced it.
