.. _tdem_basics:

TDEM Basics
===========

Time-domain electromagnetic methods, often called TEM or TDEM, measure the
earth's transient response after a controlled transmitter current changes
with time. In contrast to MT, AMT, and CSAMT, where the primary observations
are frequency-domain impedances, TDEM records decay curves at discrete time
gates. Early gates are mainly sensitive to shallow structure; later gates are
influenced by progressively deeper currents.

pyCSAMT uses TDEM in three related ways:

* reading and plotting field TDEM products such as TEMAVG and ZPLOT files;
* converting TDEM soundings to MT-like impedance/EDI products when a
  workflow needs frequency-domain compatibility;
* inverting TDEM decay data directly through the inversion API when the data
  vector contains ``times`` and ``values``.

This page explains the physical concepts behind those workflows and the
practical choices users must document.

Transient Diffusion
-------------------

In a common ground TDEM survey, a loop transmitter carries a steady current.
At switch-off, the primary magnetic field collapses. Faraday induction
creates eddy currents in the ground, and those currents diffuse outward and
downward while decaying with time. The receiver measures a secondary-field
quantity, commonly voltage, :math:`\partial B_z/\partial t`, or
:math:`\partial H_z/\partial t`.

The diffusion character of TDEM is central. A useful dimensional estimate is

.. math::

   z(t) \propto \sqrt{\frac{\rho t}{\mu_0}},

where :math:`z` is a characteristic investigation depth, :math:`\rho` is
resistivity, :math:`t` is time after switch-off, and :math:`\mu_0` is the
magnetic permeability of free space. This is only a scale estimate, not a
resolution guarantee. It says that later times and more resistive earth tend
to sample larger depths.

The TDEM response is therefore not a direct depth log. It is a time series
controlled by transmitter moment, receiver geometry, waveform, earth
conductivity, and data processing.

What Is Measured
----------------

TDEM instruments usually record receiver voltage over a sequence of off-time
windows. The data are then processed into gate-centred values. In pyCSAMT,
one sounding is represented by :class:`pycsamt.tdem.TEMSounding`, whose main
fields are:

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Field
     - Meaning
   * - ``time_gates``
     - Centre times of the off-time measurement windows in seconds.
   * - ``data``
     - Decay value at each gate.
   * - ``data_type``
     - Units or normalization of the decay values: ``"dBdt"``,
       ``"dHdt"``, ``"voltage"``, or ``"normalized_voltage"``.
   * - ``current``
     - Transmitter current in amperes.
   * - ``tx_area`` and ``tx_turns``
     - Transmitter loop area and turns, used to compute magnetic moment.
   * - ``rx_area`` and ``rx_turns``
     - Receiver coil geometry, needed for raw-voltage normalization.
   * - ``offset`` or ``rx_position``
     - Receiver position relative to the transmitter loop.
   * - ``waveform``
     - Transmitter-current waveform used for transform corrections.
   * - ``error``
     - Optional absolute uncertainty at each time gate.

The transmitter magnetic moment is

.. math::

   M = I\,n_{\mathrm{tx}}\,A_{\mathrm{tx}},

where :math:`I` is current, :math:`n_{\mathrm{tx}}` is transmitter turns, and
:math:`A_{\mathrm{tx}}` is transmitter area. A larger moment increases signal
amplitude and may improve late-time signal-to-noise, but it does not remove
the need for correct geometry and noise handling.

Decay Curves
------------

A TDEM sounding is normally inspected as a decay curve: response amplitude
versus time. Depending on the instrument and processing, the plotted value
may be signed or absolute, linear or logarithmic, voltage or normalized
field derivative.

For a simple central-loop late-time response over a uniform halfspace,
:math:`|\partial B_z/\partial t|` approximately follows a power-law decay:

.. math::

   \left|\frac{\partial B_z}{\partial t}\right|
   \propto
   M\,\rho^{-3/2}\,t^{-5/2}.

This motivates the late-time apparent-resistivity estimate used by many TEM
workflows:

.. math::

   \rho_a(t)
   =
   \frac{\mu_0}{\pi}
   \left(
      \frac{\mu_0^2 M}
           {20\,|\partial B_z/\partial t|\,t^{5/2}}
   \right)^{2/3}.

This formula is built for an idealized central-loop configuration and late
times. It is useful for quality control and approximate transformation, but
it should not be mistaken for a full inversion.

Time Gates
----------

TDEM data are not sampled at every instant. Measurements are integrated over
time windows, then represented by gate centres. Gate design controls what
the survey can see.

Early gates:

* are most sensitive to shallow structure;
* can be affected by transmitter ramp, receiver saturation, and system
  response;
* may contain strong cultural or capacitive noise;
* require careful waveform correction when the ramp is not negligible.

Late gates:

* carry deeper information;
* often have lower signal-to-noise ratio;
* are sensitive to stacking quality and background noise;
* may become negative, unstable, or clipped after processing.

Do not invert gates blindly. Inspect the decay curve and mask gates that are
outside the reliable signal band.

Waveforms
---------

The transmitter current does not always switch off instantaneously. pyCSAMT
models several waveform types in :mod:`pycsamt.tdem.waveform`:

.. list-table::
   :header-rows: 1
   :widths: 24 36 40

   * - Waveform
     - Meaning
     - When to use
   * - ``SquareWaveform``
     - Ideal square-wave current with zero ramp time.
     - Late-time processing when ramp effects are negligible.
   * - ``RampWaveform``
     - Linear current turn-off over ``ramp_off`` seconds.
     - Early-time processing where finite ramp matters.
   * - ``HalfSineWaveform``
     - Half-sine current envelope.
     - Systems whose transmitted pulse follows a sinusoidal envelope.
   * - ``CustomWaveform``
     - User-supplied time/current samples.
     - Instrument-specific waveform corrections.

The switch-off waveform matters most when the first gates are close to the
ramp duration. If a gate centre :math:`t` is much larger than the ramp time,
the ideal step approximation may be adequate. If :math:`t` is comparable to
the ramp, ignoring waveform shape can bias shallow apparent resistivity and
transformed frequency-domain products.

Central-Loop and Offset Geometry
--------------------------------

Many ground TDEM surveys use a central-loop or coincident-loop geometry,
where the receiver is near the centre of the transmitter loop. Other surveys
use offset-loop, in-loop, slingram, or moving-loop layouts. Geometry changes
the amplitude and spatial sensitivity of the response.

pyCSAMT stores:

* transmitter area and turns;
* current;
* receiver coil area and turns;
* loop shape and dimensions;
* receiver offset or 2-D receiver position.

For non-central-loop configurations, a geometry correction may be required.
The CLI conversion command exposes this as ``--no-geometry-correction`` to
disable the default in-loop correction when the user intentionally wants raw
or externally corrected values.

Apparent Resistivity
--------------------

TDEM apparent resistivity is a derived diagnostic, not a direct measurement.
It asks: "What uniform halfspace would produce this gate value under the
assumed geometry?" The answer varies with time because the real earth is not
normally a uniform halfspace.

When plotted against time or pseudo-depth, apparent resistivity helps reveal:

* conductive overburden;
* resistive basement;
* late-time conductors;
* noisy gates;
* inconsistent station responses;
* bad normalization or geometry metadata.

But apparent resistivity does not uniquely identify geology. A conductor can
come from clay, saline water, sulfides, graphite, cultural infrastructure, or
some combination. Interpretation must use geology and independent constraints.

Pseudo-Frequency
----------------

Some pyCSAMT workflows convert TDEM data to MT-like EDI files. This requires
assigning an equivalent frequency to each time gate. A common convention is

.. math::

   f_{\mathrm{eq}}(t) =
   \frac{1}{2\pi t}.

This relation is a convention, not a physical statement that the time-domain
measurement is identical to a sinusoidal MT measurement at one frequency.
It is useful for approximate comparison, plotting, and frequency-domain
pipeline compatibility.

pyCSAMT exposes pseudo-frequency choices through TDEM transforms and the CLI
option ``--freq-convention``. Users should report the convention used when
TDEM data are converted to EDI or merged with frequency-domain products.

Late-Time Transform
-------------------

The late-time transform is fast and practical. It estimates apparent
resistivity from the decay curve and then synthesizes an impedance magnitude
at pseudo-frequency:

.. math::

   |Z(f_{\mathrm{eq}})|
   =
   \sqrt{\mu_0 \omega \rho_a}.

A phase must also be assigned to build a complex impedance. pyCSAMT exposes
phase modes such as ``"homogeneous"`` and ``"weidelt"`` in the TDEM
conversion workflow.

Use the late-time transform when:

* data quality is sufficient and the target use is approximate
  frequency-domain compatibility;
* the time gates are late enough for the approximation;
* a fast survey-level conversion is more important than waveform-rigorous
  modelling.

Do not use it as proof that TDEM and MT data are interchangeable. The
transformed product is an interpretation bridge.

Fourier Transform
-----------------

The more rigorous transform uses numerical Fourier integration to connect
the time-domain response to frequency-domain quantities. It can include
waveform effects and is better suited when early gates or non-ideal
transmitter currents matter.

Use the Fourier transform when:

* waveform metadata are reliable;
* early-time information is important;
* the transmitter waveform departs from an ideal step;
* the transformed impedance will be used in a sensitive inversion or joint
  interpretation.

The cost is that it requires more complete acquisition metadata and is more
computationally involved.

pyCSAMT Data Flow
-----------------

The TDEM subpackage follows a practical data flow:

#. read TEMAVG, ZPLOT, coordinate, or log files;
#. build :class:`pycsamt.tdem.TEMSounding` objects;
#. attach geometry, waveform, station names, and coordinates;
#. inspect decay curves, sections, maps, and gate profiles;
#. optionally transform soundings to EDI collections;
#. pass transformed or direct TDEM data into inversion and interpretation
   workflows.

Minimal Python example:

.. code-block:: python
   :linenos:

   import numpy as np
   from pycsamt.tdem import TEMSounding, TEMtoEDI

   times = np.logspace(-5, -2, 30)
   dbdt = 5e-5 * times ** (-2.5)

   sounding = TEMSounding(
       times,
       dbdt,
       current=8.0,
       tx_area=100.0 ** 2,
       data_type="dBdt",
       station_name="S01",
       x=1000.0,
       y=500.0,
   )

   converter = TEMtoEDI(method="late_time", phase_mode="weidelt")
   collection = converter.transform(sounding)

Survey-folder conversion:

.. code-block:: python
   :linenos:

   from pycsamt.tdem import read_temavg_soundings, transform_temavg_survey

   soundings = read_temavg_soundings(
       "data/TEMAVG/JIANGSU",
       component="Hz",
       pattern="*.AVG",
   )

   result = transform_temavg_survey(
       "data/TEMAVG/JIANGSU",
       component="Hz",
       method="late_time",
       freq_convention="skin_depth",
       phase_mode="homogeneous",
       savepath="outputs/tdem_edi",
   )

Command-line conversion:

.. code-block:: console
   :linenos:

   pycsamt tdem info data/TEMAVG/JIANGSU
   pycsamt tdem convert data/TEMAVG/JIANGSU --dry-run
   pycsamt tdem convert data/TEMAVG/JIANGSU \
       --output-dir outputs/tdem_edi \
       --method late_time \
       --component Hz

Plotting and QC
---------------

Before conversion or inversion, plot the TDEM data. pyCSAMT provides plot
helpers for:

* decay curves;
* transformed apparent resistivity;
* TEMAVG sections;
* ZPLOT sections;
* station maps;
* elevation profiles;
* gate profiles;
* compact dashboards.

Example:

.. code-block:: python
   :linenos:

   from pycsamt.tdem import read_temavg_soundings, plot_decay

   soundings = read_temavg_soundings("data/TEMAVG/JIANGSU")
   ax = plot_decay(soundings)

When reviewing plots, look for:

* sign reversals not explained by expected physics;
* late-time noise floors;
* early-time saturation or ramp artifacts;
* stations with inconsistent amplitude after geometry normalization;
* gates missing across many stations;
* profile-aligned cultural features;
* elevation or topography effects near station clusters.

TDEM Inversion
--------------

TDEM can be inverted directly as a decay curve. In pyCSAMT's inversion
configuration, this means ``method="tdem"`` and data containing time gates
and response values. The shared error model treats TDEM values with relative
and absolute floors, using configuration options such as ``error_floor`` and
backend-specific ``tdem_relative`` or ``tdem_absolute``.

Minimal inversion sketch:

.. code-block:: python
   :linenos:

   from pycsamt.inversion import InversionConfig, InversionWorkflow

   cfg = InversionConfig(
       method="tdem",
       dimension="1d",
       backend="builtin",
       data={
           "times": [1e-5, 3e-5, 1e-4, 3e-4, 1e-3],
           "values": [2.0e-7, 5.0e-8, 1.2e-8, 2.5e-9, 6.0e-10],
       },
       n_layers=5,
       error_floor=0.05,
       regularization="smooth",
   )

   result = InversionWorkflow(cfg).run()
   print(result.rms)

Direct TDEM inversion avoids the interpretive step of synthesizing EDI files.
Use it when the inversion backend supports the required TDEM physics and the
survey geometry is represented adequately.

Combining TDEM With AMT, MT, or CSAMT
-------------------------------------

TDEM is often combined with frequency-domain methods because they complement
each other:

* TDEM is strong for shallow conductive structure and static-shift control;
* AMT/CSAMT may provide better frequency-domain continuity along profiles;
* MT can extend sensitivity to greater depths;
* transformed TDEM can help initialize or constrain frequency-domain
  interpretation.

However, the data are not identical measurements. When combining them:

* keep original TDEM decay data;
* record transform method and pseudo-frequency convention;
* avoid over-weighting transformed TDEM points as if they were measured MT
  impedance;
* check overlap consistency rather than forcing agreement;
* document whether TDEM is used as direct inversion data, transformed EDI,
  or external constraint.

Noise and Error Floors
----------------------

TDEM noise is time dependent. Early gates can be affected by transmitter and
receiver system response; late gates can approach the noise floor. A single
percentage error is rarely perfect, but it is better than pretending all gates
are exact.

The normalized residual for inversion is

.. math::

   r_i =
   \frac{d_{\mathrm{obs},i} - d_{\mathrm{pred},i}}{\sigma_i}.

For TDEM data, a practical error model combines relative and absolute floors:

.. math::

   \sigma_i =
   \max(
      f_{\mathrm{tdem}} |d_i|,
      \sigma_{\mathrm{abs}}
   ).

The absolute floor prevents tiny late-time values from receiving excessive
weight. The relative floor prevents high-amplitude gates from dominating only
because their values are large.

Depth of Investigation
----------------------

TDEM depth of investigation is not a single number. It depends on:

* transmitter moment;
* receiver noise floor;
* loop size and geometry;
* earth resistivity;
* gate range;
* stacking and processing quality;
* target geometry;
* inversion regularization.

A late-time gate may not contain useful depth information if it is dominated
by noise. Conversely, a shallow conductive layer can screen deeper structure,
reducing the ability of later gates to resolve resistive basement.

A good report should state the reliable time range and avoid interpreting
features below the supported depth range.

Common Mistakes
---------------

Avoid these mistakes:

* mixing milliseconds and seconds in ``time_gates``;
* using raw voltage without receiver area and turns;
* forgetting transmitter current, area, or turns;
* treating pseudo-frequency as measured MT frequency;
* applying a late-time transform to early-time ramp-contaminated data;
* interpreting noisy late gates as deep conductors;
* dropping sign information without documenting why;
* converting to EDI without recording transform method and phase mode;
* combining TDEM and AMT/MT data with incompatible weights;
* reporting a smooth apparent-resistivity curve as a unique geological model.

Recommended Workflow
--------------------

A conservative TDEM workflow is:

#. Load survey files and coordinates.
#. Confirm units, current, loop area, turns, receiver geometry, and time
   units.
#. Plot decay curves and gate profiles.
#. Mask saturated, noisy, missing, or sign-unstable gates.
#. Choose whether the task needs direct TDEM inversion or EDI conversion.
#. If converting, select late-time or Fourier transform and record waveform
   assumptions.
#. Compare transformed products with AMT/MT/CSAMT only in compatible bands.
#. Invert with realistic relative and absolute error floors.
#. Interpret resistivity with geological constraints.
#. Archive original data, processed data, transform settings, and plots.

Reporting Checklist
-------------------

Include the following in reports or project metadata:

* instrument and file type, such as TEMAVG or ZPLOT;
* transmitter current, loop size, loop turns, and receiver geometry;
* waveform type, base frequency, ramp time, or custom waveform file;
* time-gate range and gates removed;
* data type and units;
* coordinate file and station naming;
* transform method, pseudo-frequency convention, and phase mode;
* error floors used for inversion;
* plots of decay curves and gate/profile sections;
* limitations on depth of investigation and geological interpretation.

Next Steps
----------

For implementation details, see:

* :doc:`../api/tdem` for the generated TDEM API;
* :doc:`../cli/tdem` for command-line conversion and plotting;
* :ref:`inversion_concepts` for inversion objective functions and errors;
* :doc:`csamt_amt_mt_overview` for method comparison with MT, AMT, and CSAMT.

References
----------

The transient diffusion and apparent-resistivity concepts follow standard
TDEM references [NabighianMacnae1991]_ and [Christiansen2009]_. General EM
field theory follows [WardHohmann1988]_.
