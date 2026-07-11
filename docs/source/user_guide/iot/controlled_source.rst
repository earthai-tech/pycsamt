.. _user_guide_iot_controlled_source:

Controlled-Source Edge QC (CSAMT / CSEM)
========================================

Natural-source AMT/MT diagnostics (:doc:`amt_diagnostics`) assume a
plane-wave field with no operator-controlled transmitter. Controlled-source
methods -- CSAMT and CSEM -- add a grounded dipole transmitter, and with it
a distinct set of field questions: is the receiver in the far field where
the plane-wave assumption holds, is there energy at every transmitted
frequency, and is the source current steady? The
:mod:`pycsamt.iot.edge_csamt` and :mod:`pycsamt.iot.edge_csem` modules
answer these, numpy-only, on short edge windows.

A ``transmitter`` device role and a ``source`` telemetry packet let a
transmitter node report its state alongside the receivers.

Field zones from skin depth
---------------------------

The controlling quantity is the transmitter-receiver offset expressed in
skin depths, :math:`r/\delta`. A larger ratio means the receiver is
electrically farther from the source and the plane-wave assumption is
safer. :func:`~pycsamt.iot.classify_field_zones` labels each frequency
``near``, ``transition``, or ``far``:

.. code-block:: python

   import numpy as np
   from pycsamt.iot import classify_field_zones, skin_depth_m

   freqs = np.array([4096.0, 1024.0, 256.0, 64.0, 16.0, 4.0, 1.0])
   cov = classify_field_zones(freqs, resistivity=100.0, offset_m=5000.0)

   cov.all_far_field            # False -- some low frequencies roll off
   cov.correction_recommended   # True  -- a near-field correction is due
   cov.first_far_field_hz()     # lowest frequency safe for plane-wave MT

High frequencies (shallow skin depth) reach the far field first; the lowest
frequencies fall into the transition or near field, where CSAMT apparent
resistivities roll off and require a near-field correction.

Transmitter frequency comb
--------------------------

CSAMT transmits a *discrete* set of frequencies. QC therefore checks that
resolvable energy is actually present at each expected line, using a robust
median noise floor so closely spaced lines never mask one another:

.. code-block:: python

   from pycsamt.iot import detect_transmitter_frequencies

   comb = detect_transmitter_frequencies(
       window, sample_rate=2048.0,
       tx_frequencies=[8.0, 32.0, 128.0, 512.0],
   )
   comb.n_detected, comb.n_expected   # e.g. (3, 4)
   comb.missing()                     # frequencies with no resolvable energy

Source-signal stability
------------------------

The transmitter current sets the signal level of every sounding, so its
steadiness bounds data quality. :func:`~pycsamt.iot.assess_source_stability`
reports the on-state current coefficient of variation, drift, and the
on/off keying fraction:

.. code-block:: python

   from pycsamt.iot import assess_source_stability

   status = assess_source_stability(tx_current, tx_voltage=tx_voltage)
   status.stable, status.current_cv, status.on_fraction

CSEM: magnitude/phase versus offset
-----------------------------------

CSEM records a dipole source with a *receiver array*, and its signature
data product is the response as a function of source-receiver **offset** at
each frequency. :func:`~pycsamt.iot.field_vs_offset` builds the
magnitude-versus-offset (MVO) and phase-versus-offset (PVO) curve, finds
where the signal crosses the noise floor, and checks that amplitude decays
monotonically -- a bump usually means a bad receiver, a geometry error, or
genuine 3-D structure worth a second look:

.. code-block:: python

   from pycsamt.iot import field_vs_offset

   resp = field_vs_offset(
       offsets_m=[1000, 2000, 4000, 6000, 8000, 10000],
       amplitudes=amplitudes, phases_deg=phases,
       noise_floor=1e-13, frequency_hz=1.0,
   )
   resp.max_detectable_offset_m   # detectability limit
   resp.monotonic_decay           # False flags a suspect reading
   resp.dynamic_range_db

Transmitter telemetry
---------------------

A transmitter node reports its state as a ``source`` packet, parsed by the
``SourcePayload`` schema with the same tolerant alias folding and range
validation as the other payloads:

.. code-block:: python

   from pycsamt.iot import DeviceConfig, FieldSession, parse_payload

   tx = DeviceConfig("tx-1", role="transmitter")
   payload = parse_payload("source", {
       "tx_id": "TX1", "current": 9.8, "tx_voltage": 250.0,
       "frequency": 32.0, "ab_m": 100.0, "tx_rx_offset": 5000.0,
   })
   payload.tx_current_a, payload.tx_frequency_hz, payload.offset_m

   session = FieldSession("CS1", method="csamt", devices=[tx])
   session.add_packet({"device_id": "tx-1", "timestamp": 10.0,
                       "topic": tx.topic("source"), "kind": "source",
                       "payload": payload.as_dict()})

Static shift and transmitter timing lock
----------------------------------------

Two further controlled-source concerns. A galvanic **static shift**
multiplies apparent resistivity by a frequency-independent factor while
leaving phase unchanged; :func:`~pycsamt.iot.estimate_static_shift`
detects it as a persistent, phase-neutral split between the ``xy`` and
``yx`` resistivity modes and separates it from anisotropy:

.. code-block:: python

   from pycsamt.iot import estimate_static_shift

   ss = estimate_static_shift(res_xy, res_yx, phase_xy=phi_xy, phase_yx=phi_yx)
   ss.static_shift, ss.shift_factor   # e.g. (True, 3.0)

And a CSAMT/CSEM receiver can report its **transmitter timing lock**
alongside the clock sync, using the ``tx_locked``, ``tx_sync_offset_ms``,
and ``tx_id`` fields of the ``sync`` payload:

.. code-block:: python

   from pycsamt.iot import parse_payload

   sync = parse_payload("sync", {"offset_ms": 0.4, "transmitter_locked": True,
                                 "tx_offset_ms": 0.2, "tx_id": "TX1"})
   sync.tx_locked, sync.tx_sync_offset_ms

Aggregation
-----------

:func:`~pycsamt.iot.csamt_edge_report` and
:func:`~pycsamt.iot.csem_edge_report` collate the per-channel diagnostics,
and :func:`~pycsamt.iot.csamt_edge_table` /
:func:`~pycsamt.iot.csem_edge_table` flatten them into pyCSAMT tables for
reporting.
