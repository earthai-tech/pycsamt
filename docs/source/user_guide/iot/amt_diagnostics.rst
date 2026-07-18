.. _user_guide_iot_amt_diagnostics:

AMT/CSAMT Edge Diagnostics
==========================

The :mod:`pycsamt.iot.edge_amt` module adds :term:`edge diagnostics` that
are specific to :term:`AMT`, :term:`MT`, and :term:`CSAMT` field
telemetry. These routines are intended for the edge side of a survey:
they run on short time windows, produce compact metrics, and help decide
whether a packet should be accepted, warned, or rejected before
downstream :term:`impedance tensor` processing.

The examples below use synthetic data. No EDI or field logger file is
required. The synthetic window is deliberately built with useful AMT-like
energy, 50 Hz and 100 Hz powerline contamination, a clipped channel, a
flat dropout interval, NaNs, and two sets of synthetic impedance windows.
That makes the failure modes explicit and keeps the example reproducible.

What The Diagnostics Check
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 28 44 28

   * - Function
     - Field question
     - Typical action
   * - ``detect_powerline_harmonics``
     - Are 50/60 Hz mains harmonics dominating this window?
     - Warn or reject the packet, then inspect grounding and nearby power
       infrastructure.
   * - ``estimate_channel_snr``
     - Is useful AMT-band energy above the noise floor?
     - Prefer higher-SNR windows for transfer-function estimates.
   * - ``check_channel_saturation``
     - Is the ADC or logger rail clipping?
     - Lower gain, check coil/electrode wiring, or discard saturated
       intervals.
   * - ``check_contact_resistance``
     - Do electric channels show high drift or high noise proxies?
     - Re-seat electrodes, wet contacts, or review dipole wiring.
   * - ``estimate_frequency_coverage``
     - Which target frequency bands are actually resolved?
     - Mark missing bands before inversion or reporting.
   * - ``assess_impedance_stability``
     - Are per-window impedance estimates repeatable?
     - Keep stable windows and down-weight unstable windows.
   * - ``detect_sensor_dropout``
     - Did a channel go flat or contain gaps?
     - Reject the affected interval and inspect sensor/logger health.

Build A Synthetic Edge Window
-----------------------------

This first block creates a deterministic synthetic window. ``Ex`` contains
useful low-frequency energy plus 50 Hz and 100 Hz
:term:`powerline harmonics`. ``Ey`` is also clipped to mimic an overloaded
:term:`ADC` channel. ``Ex`` then receives a short flatline and two NaNs to
mimic a :term:`sensor dropout`.

.. code-block:: python
   :linenos:

   import numpy as np

   sample_rate = 512.0
   n_samples = 8192
   t = np.arange(n_samples) / sample_rate
   rng = np.random.default_rng(42)

   ex = (
       0.80 * np.sin(2 * np.pi * 7.5 * t)
       + 0.35 * np.sin(2 * np.pi * 23.0 * t)
       + 0.45 * np.sin(2 * np.pi * 50.0 * t)
       + 0.18 * np.sin(2 * np.pi * 100.0 * t)
       + 0.08 * rng.standard_normal(n_samples)
   )
   ey = (
       0.65 * np.sin(2 * np.pi * 9.0 * t + 0.35)
       + 0.22 * np.sin(2 * np.pi * 31.0 * t)
       + 0.06 * rng.standard_normal(n_samples)
   )
   hx = (
       0.32 * np.sin(2 * np.pi * 7.5 * t + 0.9)
       + 0.05 * rng.standard_normal(n_samples)
   )

   ex_bad = ex.copy()
   ex_bad[2500:2560] = ex_bad[2499]
   ex_bad[6000] = np.nan
   ex_bad[6001] = np.nan

   ey_clipped = np.clip(2.6 * ey, -1.0, 1.0)

Before sending that window into the checks, it helps to know what each
number means. The spectral routines work on a finite channel
:math:`x_n`, sampled at rate :math:`f_s`; non-finite values are replaced
by the finite channel mean before the spectrum is estimated. The
one-sided :term:`power spectral density` is computed with Welch averaging
and is written below as :math:`P(f)`.

The mains check then integrates the spectrum around each harmonic
:math:`f_k = k f_m`, where :math:`f_m` is 50 or 60 Hz:

.. math::

   r_k =
   \frac{\int_{f_k-\Delta f}^{f_k+\Delta f} P(f)\,df}
        {\int_0^{f_s/2} P(f)\,df}.

The total contamination ratio is :math:`R = \sum_k r_k`. A harmonic is
flagged when :math:`r_k` is greater than or equal to the configured
``threshold_ratio``. In this example, :math:`f_m=50` Hz and the 50 Hz
line is expected to be the dominant flagged harmonic. The same spectrum
also gives the channel :term:`SNR`: useful in-band power is compared with
the remaining spectral power,

.. math::

   \operatorname{SNR}_{dB} =
   10 \log_{10}
   \left(
   \frac{\int_{f_1}^{f_2} P(f)\,df}
        {\int_0^{f_s/2} P(f)\,df - \int_{f_1}^{f_2} P(f)\,df}
   \right).

where :math:`f_1=4` Hz and :math:`f_2=40` Hz for the AMT band used below.
The saturation check is more direct. With full-scale limit :math:`L`, a
finite sample is clipped when :math:`|x_n| \ge L`, and the clip fraction
is

.. math::

   c = \frac{N(|x_n| \ge L)}{N_\mathrm{finite}}.

The channel is marked saturated when :math:`c` is greater than
``max_clip_fraction``. The contact check should be read in the same
field-practical spirit: ``check_contact_resistance`` is a
:term:`contact resistance proxy`, not a direct injected-current
measurement, and it uses robust electric-channel noise/drift symptoms to
flag channels that deserve field inspection.

For frequency coverage, pyCSAMT takes the median positive PSD value as a
noise floor :math:`P_{50}`. A frequency is resolved when

.. math::

   10 \log_{10}\left(\frac{P(f)}{P_{50}}\right) \ge S_\mathrm{floor}.

The lowest and highest resolved frequencies define
:math:`[f_\mathrm{low}, f_\mathrm{high}]`; the number of covered decades
is :math:`\log_{10}(f_\mathrm{high}/f_\mathrm{low})`. A target band is
counted as covered only when its lower and upper bounds both lie inside
that resolved interval. Finally, impedance stability is checked on a
stack of complex :term:`impedance tensor` estimates :math:`Z_{ij}` from
repeated windows. For each frequency or component, pyCSAMT computes the
magnitude :term:`coefficient of variation`

.. math::

   \mathrm{CV}(|Z|) =
   \frac{\operatorname{std}(|Z|)}{\operatorname{mean}(|Z|)}

and the standard deviation of the impedance phase
:math:`\arg(Z)` in degrees. The returned summary is the average CV and
average phase scatter over the supplied impedance columns. A window set
is stable only when both metrics stay below ``max_cv`` and
``max_phase_std_deg``.

Finally, ``detect_sensor_dropout`` scans the original sample sequence.
It reports the number of NaNs and the longest run of consecutive finite
samples satisfying :math:`|x_n - x_{n-1}| \le \epsilon`. A channel is a
dropout candidate if any NaN is present or if a flat run reaches
``min_flat_run`` samples.

Run The Edge Diagnostics
------------------------

The same window can be passed through the AMT-specific checks. The
``target_bands`` argument says which bands the deployment cares about,
while ``signal_band_hz`` controls the :term:`SNR` estimate.

.. code-block:: python
   :linenos:

   from pycsamt.iot import (
       assess_impedance_stability,
       check_channel_saturation,
       check_contact_resistance,
       detect_powerline_harmonics,
       detect_sensor_dropout,
       estimate_channel_snr,
       estimate_frequency_coverage,
   )

   harmonics = detect_powerline_harmonics(
       ex_bad,
       sample_rate,
       mains_hz=50.0,
       n_harmonics=5,
       threshold_ratio=0.02,
   )
   snr_ex = estimate_channel_snr(
       ex_bad,
       sample_rate,
       signal_band_hz=(4.0, 40.0),
   )
   snr_hx = estimate_channel_snr(
       hx,
       sample_rate,
       signal_band_hz=(4.0, 40.0),
   )
   saturation = check_channel_saturation(
       ey_clipped,
       limit=1.0,
       max_clip_fraction=0.01,
   )
   contact = check_contact_resistance(
       ex_bad,
       ey,
       sample_rate=sample_rate,
       noise_rms_threshold=0.45,
   )
   coverage = estimate_frequency_coverage(
       ex_bad,
       sample_rate,
       target_bands=[(4.0, 40.0), (40.0, 120.0), (120.0, 240.0)],
       snr_floor_db=6.0,
   )
   dropout = detect_sensor_dropout(ex_bad, min_flat_run=16)

   z_stable = (1.0 + 0.45j) * (
       1.0 + 0.03 * rng.standard_normal((16, 4))
   ) * np.exp(1j * np.deg2rad(1.5 * rng.standard_normal((16, 4))))

   z_unstable = (1.0 + 0.45j) * (
       1.0 + 0.35 * rng.standard_normal((16, 4))
   ) * np.exp(1j * np.deg2rad(18.0 * rng.standard_normal((16, 4))))

   stable = assess_impedance_stability(z_stable)
   unstable = assess_impedance_stability(z_unstable)

   dominant = harmonics.dominant
   print("Synthetic AMT edge-diagnostic output")
   print(f"Powerline contaminated: {harmonics.contaminated}")
   print(f"Dominant harmonic: {dominant.frequency_hz:.1f} Hz")
   print(f"Total harmonic power ratio: {harmonics.total_ratio:.3f}")
   print(f"Ex SNR in 4-40 Hz band: {snr_ex:.2f} dB")
   print(f"Hx SNR in 4-40 Hz band: {snr_hx:.2f} dB")
   print(
       "Ey clipped samples: "
       f"{saturation['n_clipped']} of {saturation['n_samples']} "
       f"({100.0 * saturation['clip_fraction']:.1f}%)"
   )
   print(f"Ey saturated: {saturation['saturated']}")
   print(f"Contact proxy ok: {contact['ok']}")
   print(f"Contact flags: {', '.join(contact['flags'])}")
   print(
       "Frequency coverage: "
       f"{coverage.f_low_hz:.2f}-{coverage.f_high_hz:.2f} Hz "
       f"({coverage.n_decades:.2f} decades)"
   )
   print(f"Target-band coverage fraction: {coverage.coverage_fraction:.2f}")
   print(f"Missing target bands: {coverage.missing_bands}")
   print(
       "Dropout: "
       f"{dropout['dropout']}, "
       f"NaNs={dropout['n_nan']}, "
       f"longest flat run={dropout['longest_flat_run']} samples"
   )
   print(
       "Stable impedance windows: "
       f"{stable.stable} "
       f"(CV={stable.cv_magnitude:.3f}, "
       f"phase std={stable.phase_std_deg:.2f} deg)"
   )
   print(
       "Unstable impedance windows: "
       f"{unstable.stable} "
       f"(CV={unstable.cv_magnitude:.3f}, "
       f"phase std={unstable.phase_std_deg:.2f} deg)"
   )

Output:

.. code-block:: text

   Synthetic AMT edge-diagnostic output
   Powerline contaminated: True
   Dominant harmonic: 50.0 Hz
   Total harmonic power ratio: 0.232
   Ex SNR in 4-40 Hz band: 4.85 dB
   Hx SNR in 4-40 Hz band: 13.30 dB
   Ey clipped samples: 4601 of 8192 (56.2%)
   Ey saturated: True
   Contact proxy ok: False
   Contact flags: ex_noise_above_threshold, ey_noise_above_threshold
   Frequency coverage: 2.00-102.00 Hz (1.71 decades)
   Target-band coverage fraction: 0.33
   Missing target bands: [(40.0, 120.0), (120.0, 240.0)]
   Dropout: True, NaNs=2, longest flat run=61 samples
   Stable impedance windows: True (CV=0.026, phase std=1.60 deg)
   Unstable impedance windows: False (CV=0.308, phase std=18.63 deg)

Plot The Diagnostics
--------------------

The same variables can be plotted for review or embedded into an edge
report. The first figure shows the synthetic window and the live
:term:`power spectral density`. The second figure condenses the key
:term:`quality control` fractions and compares stable versus unstable
impedance windows.

.. code-block:: python
   :linenos:

   from pathlib import Path

   import matplotlib.pyplot as plt

   from pycsamt.iot import compute_live_spectra

   out_dir = Path("docs/source/images/user_guide/iot")
   out_dir.mkdir(parents=True, exist_ok=True)

   spectra = compute_live_spectra(ex_bad, sample_rate)
   freq = spectra["frequency_hz"]
   psd = spectra["psd"]

   fig, axes = plt.subplots(2, 1, figsize=(8.2, 5.8),
                            constrained_layout=True)
   axes[0].plot(t[:1600], ex_bad[:1600], lw=1.0, label="Ex synthetic")
   axes[0].plot(t[:1600], ey[:1600], lw=0.9, alpha=0.75,
                label="Ey synthetic")
   axes[0].set_xlabel("Time (s)")
   axes[0].set_ylabel("Amplitude")
   axes[0].set_title("Synthetic AMT edge window")
   axes[0].legend(loc="upper right")
   axes[0].grid(alpha=0.25)

   axes[1].semilogy(freq[1:], psd[1:], color="tab:blue")
   for peak in harmonics.peaks:
       if peak.flagged:
           axes[1].axvline(peak.frequency_hz, color="tab:red",
                           ls="--", lw=1.0, alpha=0.75)
   axes[1].set_xlim(1, 180)
   axes[1].set_xlabel("Frequency (Hz)")
   axes[1].set_ylabel("PSD")
   axes[1].set_title("Live spectrum with flagged mains harmonics")
   axes[1].grid(alpha=0.25, which="both")
   fig.savefig(
       out_dir / "user-guide-iot-amt-diagnostics-01.png",
       dpi=180,
   )

   fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.6),
                            constrained_layout=True)
   axes[0].bar(
       ["Harmonic\nratio", "Ey clip\nfraction", "Missing\nband frac"],
       [
           harmonics.total_ratio,
           saturation["clip_fraction"],
           1.0 - coverage.coverage_fraction,
       ],
       color=["tab:red", "tab:orange", "tab:purple"],
   )
   axes[0].set_ylim(0, 1)
   axes[0].set_ylabel("Fraction")
   axes[0].set_title("Edge QC fractions")
   axes[0].grid(axis="y", alpha=0.25)

   axes[1].bar(
       ["Stable Z", "Unstable Z"],
       [stable.cv_magnitude, unstable.cv_magnitude],
       color=["tab:green", "tab:red"],
   )
   axes[1].axhline(0.15, color="0.25", ls="--", lw=1.0,
                   label="max_cv")
   axes[1].set_ylabel("Coefficient of variation")
   axes[1].set_title("Impedance-window stability")
   axes[1].legend(loc="upper left")
   axes[1].grid(axis="y", alpha=0.25)
   fig.savefig(
       out_dir / "user-guide-iot-amt-diagnostics-02.png",
       dpi=180,
   )

.. grid:: 1 1 2 2

   .. grid-item::

      .. image:: ../../images/user_guide/iot/user-guide-iot-amt-diagnostics-01.png
         :width: 100%

   .. grid-item::

      .. image:: ../../images/user_guide/iot/user-guide-iot-amt-diagnostics-02.png
         :width: 100%

Interpret The Result
--------------------

The synthetic packet should not be treated as a clean acquisition window.
The mains detector flags contamination because the 50 Hz and 100 Hz bands
carry a large fraction of the spectral power. The clipped ``Ey`` channel
is also unusable for impedance estimation because more than half of its
samples are at the :term:`ADC` limit. The :term:`contact resistance proxy`
is a warning rather than a true resistance measurement; passive
:term:`AMT` data cannot measure contact resistance directly, but high
drift and high channel noise are useful field-side symptoms.

The :term:`frequency coverage` result says that only one of the requested
target bands is fully represented by this short window. That does not
mean the survey failed; it means this packet alone should not be used to
support the missing bands. The :term:`impedance stability` check
demonstrates the same principle at :term:`transfer function` level: low
:term:`coefficient of variation` and low phase scatter mark stable
windows, while high magnitude and phase scatter mark windows that should
be down-weighted or rejected.

Using These Metrics In Telemetry
--------------------------------

In an IoT acquisition workflow, these diagnostics are usually attached to
``qc`` packets or folded into the edge-processing result. A practical
edge policy might accept packets with finite data and stable impedance,
warn on mild mains contamination, and reject packets with saturation,
dropout, or severe missing-band coverage. The exact thresholds should be
survey-specific and should be recorded in the deployment manifest so that
field decisions remain reproducible.
