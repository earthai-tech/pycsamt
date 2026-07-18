.. _glossary:

Glossary
========

Core terms used throughout the pyCSAMT documentation and codebase. Most
entries are cross-referenced elsewhere with ``:term:`` roles, so the
definitions here are the single source of truth.

.. glossary::
   :sorted:

   MT
   Magnetotellurics
      A passive electromagnetic method that images subsurface electrical
      resistivity from natural time variations of the Earth's electric and
      magnetic fields. The ratio of horizontal electric to magnetic field
      components yields the :term:`impedance tensor`.

   AMT
   Audio-frequency magnetotellurics
      Magnetotellurics in the audio-frequency band (roughly 1 Hz to 10 kHz),
      giving shallow-to-intermediate depth resolution. The natural source is
      weak in the AMT "dead band" near 1–5 kHz.

   CSAMT
   Controlled-source audio-frequency magnetotellurics
      An active-source variant of :term:`AMT` that uses a grounded electric
      dipole or magnetic loop transmitter to overcome the weak natural signal,
      at the cost of near-field and source-overprint corrections.

   CSUMT
   Controlled-source ultra-audio magnetotellurics
      A controlled-source magnetotelluric method operating above the usual
      audio-frequency AMT band, commonly in the kilohertz to hundreds of
      kilohertz range. In pyCSAMT, CSUMT tools are used for shallow target-depth
      planning, Bostick-depth estimates, and transmitter frequency scheduling.

   CSEM
      Controlled-source electromagnetics — the broader family of active-source
      EM methods to which :term:`CSAMT` belongs.

   TDEM
   TEM
   Time-domain electromagnetics
      A transient EM method: a controlled-source current is switched off
      abruptly and the induced secondary-field decay is recorded over time,
      rather than the continuous-wave :term:`impedance tensor` that
      :term:`AMT`/:term:`MT`/:term:`CSAMT` estimate. There is no steady
      spectrum to test for mains contamination, so pyCSAMT's method-aware
      :term:`edge diagnostics` skip powerline-harmonic detection for
      TDEM/TEM streams.

   Controlled-source
      An electromagnetic acquisition mode in which the source field is generated
      by field equipment rather than by natural variations. CSAMT and CSEM are
      controlled-source methods.

   Plane-wave field
      An electromagnetic field approximation in which the wavefront is treated
      as locally planar at the receiver. Standard MT-style interpretation is
      most defensible when the transmitter is far enough away for this
      approximation to hold.

   Grounded dipole transmitter
      A controlled-source transmitter that injects current between two grounded
      electrodes. Its length, current, frequency, and receiver offset are part
      of CSAMT/CSEM acquisition metadata.

   Transmitter-receiver offset
      The separation between the controlled-source transmitter and a receiver.
      In CSAMT field-zone checks it is compared with skin depth to classify
      near, transition, and far field behaviour.

   Tensor
      In magnetotellurics, a frequency-dependent 2×2 matrix relating two
      vector field quantities — most importantly the :term:`impedance tensor`
      and the :term:`phase tensor`.

   Impedance tensor
   Z
      The frequency-dependent 2×2 complex tensor :math:`\mathbf{Z}` relating the
      horizontal electric and magnetic fields, :math:`\mathbf{E} =
      \mathbf{Z}\,\mathbf{H}`. Its elements (:math:`Z_{xx}, Z_{xy}, Z_{yx},
      Z_{yy}`) are the primary MT observable, from which
      :term:`apparent resistivity` and :term:`phase` are derived.

   Phasor
      A complex value represented by magnitude and phase, often drawn as a point
      or vector in the complex plane. Impedance phasor plots show how
      :term:`impedance tensor` components move with period before they are
      reduced to apparent resistivity and phase.

   Determinant response
      A compact impedance summary based on :math:`\det(\mathbf{Z})`. Because the
      determinant is unchanged by horizontal coordinate rotation, it is useful
      for station-level checks that should not depend on a chosen strike angle.

   Off-diagonal antisymmetry
      The ideal 1-D/2-D impedance relation :math:`Z_{xy}\approx -Z_{yx}`. A large
      antisymmetry residual means the off-diagonal modes no longer cancel as a
      simple 1-D/2-D response would suggest.

   Apparent resistivity
      The resistivity of an equivalent uniform half-space that would produce the
      observed impedance at a given frequency, :math:`\rho_a = 0.2\,|Z|^2 / f`
      (with :math:`Z` in field units and :math:`f` in Hz). It is *apparent*
      because the real earth is layered/heterogeneous.

   Anisotropy
      Direction-dependent electrical resistivity. In CSAMT/AMT diagnostics it
      usually means that responses measured in two horizontal directions are not
      equivalent, so the :math:`Z_{xy}` and :math:`Z_{yx}` modes can imply
      different :term:`apparent resistivity` or phase behavior.

   Axial anisotropy
      A simplified anisotropic-earth case where electrical properties have
      preferred horizontal axes. It can produce systematic differences between
      off-diagonal impedance modes without requiring every response to be fully
      3-D.

   Phase
      The phase angle of an :term:`impedance tensor` element versus period; for a
      layered earth it tracks whether resistivity increases or decreases with
      depth. Reported in degrees.

   Phase tensor
      A real 2×2 tensor derived from the real and imaginary parts of
      :math:`\mathbf{Z}` that is provably immune to galvanic
      :term:`static shift`. Its ellipse and skew angle summarise dimensionality
      and strike without distortion.

   Static shift
      A frequency-independent vertical shift of the :term:`apparent resistivity`
      curve caused by galvanic charges on small near-surface heterogeneities.
      Left uncorrected it biases inverted depths and resistivities.

   Galvanic distortion
      Frequency-independent distortion of the measured electric field caused by
      near-surface conductivity heterogeneities. It mixes and scales impedance
      tensor components without carrying the same inductive depth information as
      the regional earth response.

   Groom-Bailey decomposition
   Groom-Bailey
      A galvanic-distortion model that represents the observed impedance as a
      real distortion matrix multiplying an underlying regional tensor. It is
      commonly summarized by gain, twist, shear, and anisotropy-style
      parameters.

   Distortion matrix
      The real 2 x 2 matrix used in a :term:`Groom-Bailey decomposition` to
      describe local, frequency-independent mixing and scaling of the electric
      field before the regional impedance is observed.

   Regional tensor
      The impedance tensor that would be measured without local galvanic
      distortion. In 2-D Groom-Bailey workflows it is usually approximated as an
      anti-diagonal tensor after rotation to geoelectric strike.

   Twist
      A Groom-Bailey parameter describing rotational mixing of the two horizontal
      electric-field components by local galvanic distortion. It is reported as
      an angle, but should not be interpreted as geoelectric strike.

   Shear
      A Groom-Bailey parameter describing non-orthogonal mixing of horizontal
      electric-field components. Large absolute shear can indicate strong local
      distortion or a poor 2-D assumption.

   Source overprint
      A controlled-source effect where the finite transmitter field influences
      the measured response strongly enough that a plane-wave MT interpretation
      becomes unreliable.

   Near field
      The source-proximal regime where transmitter geometry strongly affects the
      measured electromagnetic field. CSAMT near-field rows usually require
      explicit review before plane-wave interpretation.

   Near-field correction
      A correction or review step applied when controlled-source measurements
      are not far enough from the transmitter for a plane-wave approximation.
      It is commonly triggered by near-field or transition-field CSAMT rows.

   Far field
      The source-distant regime where the transmitter field is close enough to a
      plane wave for standard MT-style interpretation to be more defensible.

   Transition field
      The intermediate regime between near field and far field, where source
      effects may be present but are not as dominant as in the near field.

   Phased-array source
   PAS
      A controlled-source transmitter layout made from multiple source dipoles
      with controlled spacing and phase shifts, used to steer or narrow the
      transmitted field pattern.

   Single-dipole antenna source
   SDAS
      A single controlled-source dipole transmitter element. A phased-array
      source combines several SDAS elements.

   Array factor
      The interference pattern produced by the spacing and phase shifts of
      multiple source elements, independent of the radiation pattern of one
      element.

   Grating lobe
      An unintended strong beam direction produced when array spacing is too
      large relative to wavelength or when steering creates additional valid
      main-lobe solutions.

   Directivity
      The ratio between peak radiation intensity and average radiation intensity
      over all directions. Higher directivity means energy is concentrated into
      a narrower angular region.

   Cross-power spectrum
   Cross-spectra
      A complex frequency-domain estimate of how two channels vary together. The
      diagonal entries are auto-power spectra; off-diagonal entries carry
      amplitude and phase relationships between channels.

   Power spectral density
   PSD
      The auto-power spectrum of one channel as a function of frequency. In a
      cross-power matrix it is stored on the diagonal.

   Powerline harmonics
      Spectral peaks at integer multiples of the local mains frequency, usually
      50 or 60 Hz. In AMT/CSAMT edge diagnostics they are treated as cultural
      noise because they can dominate natural-field energy in short packets.

   Noise floor
      A background signal level used as the detection threshold for spectra or
      amplitudes. pyCSAMT often uses robust median estimates so isolated peaks
      do not define the floor.

   Transmitter frequency comb
      The discrete set of frequencies emitted by a controlled-source
      transmitter. Edge QC checks whether each expected line has resolvable
      energy in the receiver window.

   Source current
      The current injected by a controlled-source transmitter. Its mean,
      variation, and drift constrain the reliability of receiver-side
      amplitudes.

   Source stability
      A transmitter QC summary based on the on-state source current and optional
      voltage. Stable source records have low current variation and no missing
      current.

   SNR
   Signal-to-noise ratio
      The ratio between useful signal power and estimated noise power, commonly
      reported in decibels as :math:`10\log_{10}(P_\mathrm{signal}/P_\mathrm{noise})`.
      pyCSAMT edge diagnostics can estimate it from in-band versus out-of-band
      spectral power.

   Coherence
      A normalized measure of the linear relationship between two channels at a
      frequency. Squared coherence ranges from 0 to 1 and is often used as a
      frequency-band quality-control metric.

   Frequency coverage
      The frequency interval in which a packet has spectral power above the
      configured noise floor. It is often reported as a low/high frequency,
      covered decades, and the fraction of survey target bands represented.

   Method profile
      The canonical per-method acquisition characteristics used by
      pyCSAMT's IoT layer: typical frequency band, required channels, a
      nominal sample rate, and whether the method is controlled-source and
      powerline-sensitive. It turns an :term:`AMT`/:term:`MT`/:term:`CSAMT`/
      :term:`CSEM`/:term:`TDEM` label into concrete :term:`quality control`
      defaults instead of leaving every threshold to be set by hand.

   Quality control
   QC
      The set of checks used to decide whether data are acceptable, need review,
      or should be rejected before interpretation or inversion.

   Timestamp
      A numeric time label attached to a packet or sample. In IoT acquisition it
      may be an epoch time or a relative survey time, but it must be finite and
      non-negative so clocks can be compared.

   Reference clock
      The timing source treated as correct when auditing field-node clocks, such
      as GPS time, a disciplined base-station clock, or a laboratory timing
      standard.

   GPS
      Global Positioning System. In field acquisition it is commonly used as a
      timing reference as well as a positioning system.

   GPS lock
      The state in which a receiver reports that it is actively disciplined by
      GPS or an equivalent reference. Loss of lock means the node may be
      free-running even if its current offset is still small.

   Clock offset
      The local-minus-reference timestamp difference at a sample time. pyCSAMT
      reports the median offset in milliseconds for clock-sync status.

   Clock drift
      The rate at which clock offset changes with time, estimated as the slope
      of local-minus-reference error versus reference time.

   Parts per million
   PPM
      A relative rate unit equal to one part in :math:`10^6`. For clock drift,
      1 ppm means a clock gains or loses about one microsecond per second.

   Timing jitter
      Short-timescale timing scatter after the best linear clock drift has been
      removed. pyCSAMT reports it as the standard deviation of residual offset in
      milliseconds.

   Synchronisation quality
      A compact grade summarising clock offset, drift, jitter, and reference-lock
      state for one field node.

   Remote reference
      A processing workflow that uses a separate, synchronised station as a noise
      reference for transfer-function estimation. It is sensitive to clock errors
      between stations.

   Ridge regularization
      A small positive value added to a matrix diagonal before inversion to
      stabilize poorly conditioned least-squares estimates.

   ADC
   Analog-to-digital converter
      The logger component that converts a continuous voltage into discrete
      digital samples. If the input exceeds the configured full-scale range, the
      samples clip at the rail and the channel is saturated.

   Edge diagnostics
      Lightweight checks run close to acquisition, often on the logger or an
      edge gateway, before full transfer-function processing. They summarize
      packet health, spectral content, channel faults, and field conditions.

   Decimation
      Reduction of a time series by keeping every nth sample or otherwise
      lowering the sample rate. In generic edge QC, decimation controls how many
      samples are emitted in the compact payload.

   Finite coverage
      The fraction of samples in a window or channel that are finite numbers
      rather than NaN or infinity. Low finite coverage usually indicates gaps,
      logger faults, or corrupted packets.

   Edge decision
      The compact accept/warning/reject state assigned by edge-side quality
      control. It travels with QC telemetry so downstream monitoring can audit
      what the field node decided.

   Channel summary
      A per-channel edge-QC record containing finite coverage, RMS, basic
      statistics, spike fraction, acceptance state, and rejection reasons.

   Robust spike fraction
      The fraction of finite samples that exceed a robust median/MAD threshold.
      It is less sensitive to a few extreme samples than a mean/std-only spike
      detector.

   Median absolute deviation
   MAD
      The median of absolute deviations from the median. Multiplying MAD by
      1.4826 gives a robust estimate of standard deviation for normally
      distributed noise.

   Time series
      A sequence of measurements ordered by time, such as live Ex, Ey, Hx, and
      Hy samples recorded by an edge device.

   Receiver array
      A set of receivers deployed at multiple offsets or stations to sample the
      controlled-source field. In CSEM it is used to build magnitude- and
      phase-versus-offset curves.

   Magnitude-versus-offset
   MVO
      A CSEM curve showing response amplitude as a function of
      transmitter-receiver offset at one frequency.

   Phase-versus-offset
   PVO
      A CSEM curve showing response phase as a function of
      transmitter-receiver offset at one frequency.

   Detectability limit
      The farthest offset or highest/lowest frequency where a signal remains
      above the configured noise floor.

   Dynamic range
      The ratio between the largest and smallest usable signal amplitudes,
      commonly reported in decibels.

   Transmitter timing lock
      A receiver-side sync status indicating that the receiver timing is locked
      to, or explicitly checked against, the controlled-source transmitter
      timing.

   Telemetry packet
      A single timestamped message reported by an IoT field device — data,
      :term:`QC`, heartbeat, or event — carrying a JSON-like payload.
      pyCSAMT represents it as a ``TelemetryPacket`` and aggregates the
      stream into station tables and a :term:`monitoring status`.

   Store-and-forward
      A telemetry delivery pattern in which a client queues a
      :term:`telemetry packet` instead of dropping it when the transport is
      unavailable, then drains the queue in order once connectivity
      returns. pyCSAMT's ``StoreAndForwardClient`` wraps any transport this
      way, with an optional spool file so the backlog survives a restart.

   Field session
      The operational record of one IoT-enabled survey, implemented as
      ``pycsamt.iot.FieldSession``. It groups device and station inventory
      with the accumulated :term:`telemetry packet` stream and can produce a
      :term:`monitoring status` and a :term:`pipeline hand-off` for
      downstream processing.

   Packet success rate
      The fraction of collected :term:`telemetry packet`\ s whose transport
      acknowledgement (``ack_ok``) is true. It reflects link/transport
      reliability, independent of whether the payload itself was judged
      acceptable by :term:`edge diagnostics`.

   Packet acknowledgement
      A transport or storage confirmation showing whether a packet was
      successfully received, written, or otherwise accepted by the next system.
      In monitoring payloads this is usually the ``ack_ok`` field.

   Edge acceptance rate
      The fraction of :term:`telemetry packet`\ s whose edge :term:`quality
      control` decision is "accept", counted only among packets that carry
      an edge decision. Packets without one are excluded rather than
      treated as rejected.

   Monitoring status
      The per-stream health summary produced by assessing a
      :term:`telemetry packet` stream against a monitoring configuration:
      an overall level (ok/warning/critical/no_data), :term:`packet success
      rate`, :term:`edge acceptance rate`, minimum battery voltage, maximum
      clock offset, maximum packet gap, and the list of threshold issues
      that were triggered.

   Latency
      The delay between packet acquisition and packet arrival or assessment.
      pyCSAMT uses a payload ``latency_s`` value when present, otherwise it can
      estimate latency from ``now - timestamp``.

   Packet gap
      The elapsed time between consecutive packet timestamps after sorting the
      telemetry stream. Large gaps usually indicate dropouts, buffering, or
      communication loss.

   Battery voltage
      The device supply voltage reported in telemetry. Monitoring compares its
      minimum observed value with the configured low-voltage threshold.

   Battery capacity
      The stored energy available from a battery, commonly expressed in
      watt-hours. pyCSAMT power budgets treat this as the starting energy before
      reserve is held back.

   Energy reserve
      The fraction or amount of battery energy intentionally kept unused so a
      node is not planned down to complete depletion.

   Active/sleep duty cycle
      The fraction of time a node spends in its active acquisition state versus
      its lower-power sleep state.

   Regulator efficiency
      The fraction of battery energy delivered to useful electronics after DC/DC
      conversion losses. Lower efficiency increases daily load.

   Telemetry window
      The daily time spent transmitting or receiving telemetry. Its energy is
      transmitter power multiplied by telemetry seconds per day.

   Edge-processing overhead
      Extra daily energy consumed by local processing on the node, separate from
      baseline acquisition and radio telemetry.

   Auxiliary load
      Any additional daily energy draw not captured by recorder, telemetry, or
      edge-processing terms, such as heaters, relays, or external sensors.

   Solar harvesting
      Energy collected from a solar panel or similar harvester. pyCSAMT applies
      charge efficiency before comparing harvested energy with daily load.

   Power budget
      A deployment estimate balancing battery capacity, daily load, harvest,
      reserve, and minimum runtime requirements.

   No-harvest autonomy
      The runtime available from usable battery energy if no daily harvest is
      available.

   Power state
      The compact power-budget classification: sustaining, ok, warning, or
      critical.

   Daily energy deficit
      A positive net daily draw after harvest is subtracted from daily load. A
      sustained deficit means the battery will eventually be depleted.

   Synthetic data
      Data generated by a model or simulator rather than recorded by field
      hardware. In the IoT guide it is used for reproducible examples, tests,
      and demonstrations, and should be labelled clearly when mixed with real
      survey workflows.

   Random seed
      The initial value passed to a pseudorandom number generator so a sequence
      of random-looking draws can be reproduced exactly.

   Gaussian noise
      Random noise whose samples follow a normal distribution. pyCSAMT
      simulation examples use it for background channel noise, clock jitter, and
      small battery-voltage perturbations.

   Dropout gap
      A contiguous interval of missing samples inserted into a simulated or
      observed time series. In arrays it is commonly represented by ``NaN``
      values and lowers finite coverage.

   Packet loss
      The removal or non-arrival of telemetry packets during transport. The
      simulator reproduces it by drawing a seeded random keep/drop mask for the
      packet stream.

   Battery decay
      A decreasing battery-voltage curve used to mimic discharge through time.
      pyCSAMT simulation uses an exponential sag with small seeded voltage
      noise.

   Configured frequency band
      The frequency interval expected for a method or deployment. Monitoring
      flags packets whose reported frequency band extends outside this interval.

   Pipeline hand-off
      The compact per-station summary produced by
      ``FieldSession.to_pipeline_input``, carrying packet counts, edge
      acceptance rate, and the accepted frequency band forward into
      processing code that does not have direct access to the raw
      :term:`telemetry packet` stream.

   Provenance manifest
      A reproducibility record for a field session or processing run. It keeps
      acquisition metadata, thresholds, QC decisions, accepted bands, rejected
      windows, and other audit information together with the data products.

   Transport security
      The protection applied while telemetry moves between a field node and a
      receiver, gateway, broker, or server. In pyCSAMT IoT examples this is
      configured through TLS options and then enforced by the concrete
      transport implementation.

   TLS
   Transport Layer Security
      The standard encrypted transport protocol used by HTTPS, secure MQTT, and
      similar clients. pyCSAMT stores TLS settings such as certificate paths,
      verification mode, and minimum version, but the cryptographic handshake is
      performed by the underlying transport library.

   TLS material
      The certificate and key configuration needed to initialise a TLS-enabled
      client, including CA certificates, optional client certificate files, and
      optional private-key files.

   Authentication
      The process of proving the identity or authorisation of a telemetry
      client before a server accepts its packets. pyCSAMT represents this with a
      ``Credential`` using bearer, API-key, basic, or no-auth schemes.

   Authentication header
      A protocol header carrying credential material, such as an HTTP
      ``Authorization`` header or a deployment-specific API-key header.

   Secret
      A value that should not appear in logs, notebooks, manifests, or version
      control, such as a bearer token, password, API key, or private key.

   Redaction
      Replacement of a non-empty secret value with a fixed placeholder before a
      configuration is printed or serialised. pyCSAMT uses redaction for
      credential summaries while preserving non-secret fields such as
      certificate paths.

   Telemetry protocol
      A transport family used to move field telemetry, such as HTTPS, MQTT,
      WebSocket, serial, or file-backed replay.

   Bearer token
      A credential string sent with a ``Bearer`` authentication scheme. Whoever
      presents the token is treated as authorised until the server rejects or
      expires it.

   API key
      A credential string sent in a configured header field to identify or
      authorise a client. API keys are secrets and should be redacted outside
      the live transport boundary.

   Basic authentication
      A username/password authentication scheme whose header contains the
      Base64 encoding of ``username:password``. The encoding is not encryption,
      so basic authentication should be used only over TLS.

   Environment variable
      A process-level key/value setting used to inject deployment-specific
      configuration at runtime. pyCSAMT security helpers can read
      ``PYCSAMT_IOT_*`` variables for credentials and TLS paths.

   Protocol policy
      A local allow-list that decides whether a requested telemetry protocol is
      permitted before a client is built.

   Short-lived credential
      A token, API key, or password issued with a limited validity window. Short
      lifetimes reduce the impact of accidental exposure because the credential
      expires without needing to reproduce the full acquisition.

   Operational acquisition plot
      A figure that explains field-system status rather than subsurface
      response. IoT operational plots show telemetry, QC, power, timing, and
      station health before geophysical transfer functions or inversions are
      interpreted.

   Matplotlib figure
      The Python object returned by Matplotlib plotting calls. pyCSAMT IoT
      plotting helpers return this object and attach the rows used to draw it so
      report figures remain auditable.

   Normalised plot data
      The table-like dictionaries derived from sessions, packets, or result
      objects before drawing a figure. They remove input-format differences so
      the same plotting code can handle live objects and serialised mappings.

   Field dashboard
      A compact four-panel overview of an IoT field session, combining station
      health, edge-QC acceptance, power or synchronisation state, and packet
      timing.

   Energy estimate
      The computed result of a power-budget calculation, including daily load,
      harvest, runtime, autonomy, power state, and triggered issues.

   Synchronisation status
      The per-device result of a clock-sync assessment, including offset, drift,
      jitter, reference support, GPS lock, and an overall quality grade.

   Content hash
      A SHA-256 digest of a :term:`provenance manifest` or audit record's
      canonical JSON encoding -- keys sorted, no whitespace -- computed over
      every field except itself, so the same content always hashes the same
      way regardless of field-insertion order. Recomputing it independently
      is how a reviewer checks that an exported manifest was not altered.

   Hash chain
      A tamper-evident sequence of records where each entry's hash folds in
      the previous entry's hash, so altering, inserting, or reordering any
      entry breaks every hash from that point forward. pyCSAMT chains QC
      decision logs this way so silent edits to the audit trail become
      detectable.

   Manifest signature
      An HMAC-SHA256 signature over a :term:`provenance manifest`'s
      canonical JSON, computed with a shared key. Unlike a :term:`content
      hash`, which anyone can recompute, a valid signature also proves the
      manifest was produced (or re-signed) by a party holding that key.

   Sensor dropout
      A missing or stuck sensor interval, visible as NaN gaps or unusually long
      flat runs in a time series.

   Contact resistance proxy
      An indirect field-side indicator for electric electrode contact quality.
      Passive AMT cannot measure true contact resistance without an injected
      test current, so pyCSAMT uses drift and noise symptoms as a proxy warning.

   Impedance stability
      A repeatability measure for windowed complex impedance estimates. Stable
      windows have low variation in impedance magnitude and low phase scatter
      across repeated estimates.

   Coefficient of variation
   CV
      A unitless relative-spread statistic defined as standard deviation divided
      by the mean. It is useful for comparing impedance-magnitude variability
      across frequencies with different absolute amplitudes.

   Transfer function
      A frequency-domain relation mapping input field components to output field
      components. The MT/AMT impedance tensor and tipper are transfer functions.

   Tipper
      The complex vertical magnetic transfer function relating the vertical to
      the horizontal magnetic field. Its induction arrows point toward (or away
      from) lateral conductivity contrasts.

   Strike
   Geoelectric strike
      The azimuth of the principal geoelectric direction, estimated from the
      impedance or :term:`phase tensor`. Rotating data to strike separates the
      TE and TM modes for 2-D interpretation.

   Dimensionality
      A classification of the subsurface electrical structure sensed by a site as
      1-D (layered), 2-D (strike-invariant), or 3-D, typically assessed from
      :term:`phase tensor` skew and impedance invariants.

   Skin depth
      The depth at which an EM field attenuates to :math:`1/e` of its surface
      amplitude, :math:`\delta \approx 503\,\sqrt{\rho / f}` metres. It sets the
      depth of investigation for a given period and resistivity.

   TE mode
   TM mode
      The transverse-electric (electric field along strike) and transverse-magnetic
      (magnetic field along strike) polarisations into which 2-D MT data separate
      after rotation to :term:`strike`.

   EDI
      The Electrical Data Interchange file format — the SEG standard text format
      for storing MT/AMT impedances, :term:`tipper`, and metadata per station.

   AVG file
      Zonge instrument-averaged CSAMT/AMT export format; pyCSAMT reads it and can
      transform it to :term:`EDI`.

   J-file
      The A.G. Jones ("BIRRP") ASCII format for MT transfer functions, readable by
      pyCSAMT and convertible to :term:`EDI`.

   Occam2D
      A 2-D regularised (smooth) inversion scheme and file format for MT data;
      pyCSAMT can write its data, mesh, and startup inputs.

   ModEM
      A widely used 3-D MT inversion code; pyCSAMT can prepare its data files.

   MARE2DEM
      A 2-D/2.5-D goal-oriented adaptive finite-element inversion code for MT and
      CSEM; pyCSAMT can build its input files.

   EMAP
      Electromagnetic array profiling — a spatially continuous acquisition and
      filtering approach used to suppress :term:`static shift`.

   RMS
   RMS misfit
      Root-mean-square misfit between observed and predicted responses, normalised
      by data error; the primary goodness-of-fit measure for an inversion.

   MapData
      The normalized survey container used by ``pycsamt.map``. It stores the
      loaded site objects, one :term:`station record` per station, one
      :term:`profile line` per survey line, and loader metadata so different
      map renderers use the same station order and grouping.

   Station record
      A single normalized station row in :term:`MapData`, containing the station
      identifier, latitude, longitude, elevation, profile-line name, zero-based
      station index, and the original EDI-like source object.

   Profile line
      A named ordered group of station records. Map loaders use profile lines to
      keep stations from different survey traverses separate while still
      allowing combined 2-D and 3-D views.

   Pseudosection
      A profile plot where stations or along-line distance form the horizontal
      axis and period or frequency forms the vertical axis. Values such as
      apparent resistivity or phase are sampled from each station response and
      gridded for visual continuity; the result is a display aid, not a true
      depth section.

   Frequency grid
      The ordered set of frequencies available in one station response. Nearby
      stations may not share exactly the same grid, so pyCSAMT records the
      requested frequency, selected frequency, absolute difference, and relative
      difference when extracting map values.

   Coordinate reference system
   CRS
      The coordinate definition used to interpret map coordinates, including
      datum, projection, units, and axis order. A CRS transform makes station
      coordinates comparable when one source is projected in metres and another
      expects WGS84 longitude and latitude.

   Basemap
      The geographic tile layer and map-camera settings used behind station
      traces, labels, contours, and profile lines. In pyCSAMT map helpers the
      basemap configuration stores style, center, zoom, and bearing separately
      from the data traces.

   Contour overlay
      A map layer made by interpolating scattered station values onto a regular
      grid and drawing filled bands, contour lines, or both. It is useful for
      visual continuity between stations, but it should be interpreted as an
      interpolation of sampled values rather than a measured continuous field.

   Topography overlay
      A 3-D terrain layer built from scattered station elevations or a regular
      elevation grid. It provides surface context for map and volume views
      without changing the underlying electromagnetic response values.

   Station map
      A 2-D map view that shows station positions or station order, optional
      profile-line traces and station labels, and one scalar overlay value per
      station. It is usually the first spatial quality-control view for a
      loaded survey.

   Scalar overlay
      A single numeric value assigned to each station for coloring a map marker
      or building an interpolated layer. Station maps commonly use station
      index, elevation, apparent resistivity, phase, or skin depth as scalar
      overlays.

   Density layer
      A Plotly map layer that displays the spatial concentration or intensity of
      finite station values beneath the station markers. In station maps it is
      used as a quick visual trend layer, not as a replacement for measured
      station values.

   Figure export
      The process of writing a map figure to a persistent artifact such as HTML,
      PNG, SVG, PDF, JSON, or a dictionary-style figure specification. pyCSAMT
      map exports return the final path so workflows can record exactly which
      artifact was produced.

   Figure specification
      The serialized structure of a Plotly figure, including data traces,
      layout, color scales, and map settings. It is useful for testing and
      audit workflows because it can be compared without rendering pixels.

   3-D quick-look map
      A non-inversion 3-D visualization built from station or pseudosection
      values. In pyCSAMT volume views, apparent resistivity and period are used
      to estimate pseudo-depth so survey trends can be inspected before a
      constrained inversion model is available.

   Pseudo-depth
      A depth-like plotting coordinate estimated from electromagnetic sampling
      scale rather than recovered by inversion. For pyCSAMT EDI volume maps it
      is computed from apparent resistivity and period using the skin-depth
      relation :math:`z \approx 503\sqrt{\rho_a T}`.

   Depth slice
      A horizontal 3-D quick-look surface sampled at one pseudo-depth. Values
      are interpolated from the station/profile grids and should be interpreted
      as a visualization of the pseudosection-derived point cloud, not a
      geological layer boundary.

   Block volume
      A sparse 3-D volume rendering built from the finite pseudo-depth samples
      across all survey lines. It gives a compact impression of the full survey
      volume, but it can hide individual line structure when station spacing is
      sparse.

   Isosurface
      A 3-D surface connecting points with the same plotted value. In pyCSAMT
      volume maps it is built from a dense interpolation of the finite
      pseudo-depth point cloud and is controlled by ``iso_range`` and
      ``surface_count``.

   Static image export
      Writing a non-interactive image such as PNG, SVG, PDF, or WebP from a
      figure. Matplotlib static image export uses ``savefig`` directly, while
      Plotly static image export requires an image engine such as Kaleido.

   Kaleido
      Plotly's static-image export engine. It converts Plotly figures to image
      files such as PNG, SVG, PDF, and WebP when interactive HTML is not the
      desired artifact.

   Agent
      In :mod:`pycsamt.agents`, a small composable unit that performs one step of
      a workflow (routing, parsing, QC, forward modelling, inversion, reporting)
      and returns a standardised ``AgentResult``.

   MapView session
      The in-memory, code-first survey handle created by
      ``pycsamt.map.MapView``. It wraps one normalized
      :term:`MapData` object, possibly spanning several survey lines,
      so that station maps, :term:`pseudosection`\ s, 3-D
      :term:`fence view`\ s, and exports all read from the same loaded
      data instead of re-parsing EDI files for every figure.

   Fence view
      A 3-D "fence diagram" that renders each survey line as its own
      vertical resistivity curtain, positioned in 3-D by station offset
      and line spacing, with the vertical axis converted from period to
      a pseudo-depth via :math:`\delta \approx 503\sqrt{\bar\rho\,T}` --
      the same :term:`skin depth` relation used elsewhere, evaluated
      with the per-period median :term:`apparent resistivity`
      :math:`\bar\rho` across the line and period :math:`T`. It is one
      of the modes built by :mod:`pycsamt.map.volume`, alongside block,
      depth-slice, and surface modes, and remains a pseudo-depth
      visualization rather than a constrained inversion model.

   Station distance
      The cumulative along-line separation between stations, in
      kilometres, used as the ``x_axis="distance"`` option on
      :term:`profile line` and :term:`pseudosection` views. pyCSAMT
      projects latitude/longitude to a local equirectangular frame,
      :math:`x=\lambda\cdot111.320\cos\bar\phi` and
      :math:`y=\phi\cdot110.574`, with :math:`\bar\phi` the mean
      station latitude, then sums consecutive point separations
      :math:`\sqrt{\Delta x^2+\Delta y^2}`. It falls back to plain
      station order when any station is missing a finite coordinate,
      so a distance axis is never silently wrong -- only degraded to
      an index.
