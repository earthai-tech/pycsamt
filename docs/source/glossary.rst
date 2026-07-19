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

   1D
      A one-dimensional earth representation in which resistivity varies only
      with depth. In magnetotelluric and audio-magnetotelluric inversion, a 1D
      model is commonly represented as horizontal layers with resistivity and
      thickness parameters.

   2-D
   2D
      A two-dimensional earth representation in which resistivity varies along a
      profile and with depth, :math:`\rho=\rho(x,z)`, while being assumed
      approximately constant along geoelectric strike. It is the usual
      assumption behind profile inversions such as :term:`Occam2D`.

   3D
      A three-dimensional earth representation in which resistivity can vary
      with horizontal position and depth. A 3D interpretation is needed when
      strike-invariant or layered assumptions cannot explain the measured
      electromagnetic response.

   Backend
      The computational route used to prepare, solve, or load an inversion
      problem. In pyCSAMT a backend may be pure Python, an optional scientific
      library, or an adapter around an external native engine.

   Inversion backend
      The backend selected by ``InversionConfig.backend`` in
      :mod:`pycsamt.inversion`. It gives the backend-neutral workflow a named
      solver route while preserving a common configuration and result interface.

   Inversion model
      The resistivity or conductivity distribution recovered by an inversion
      workflow from observed electromagnetic data and modelling assumptions.
      In reporting it is the reviewed source object from which interpreted,
      calibrated, gridded, and plotted products must be traced.

   Calculated resistivity model
   CRM
      The source resistivity model calculated by inversion before
      interpretation calibration is applied. In pyCSAMT interpretation reports
      it is usually stored as :math:`m=\log_{10}\rho`, with :math:`\rho` in
      ohm metres, and is kept distinct from field observations and from the
      calibrated model used for delivery products.

   Calibrated new model
   NM
      The resistivity model after applying documented calibration constraints,
      borehole evidence, or interpretation parameters to the calculated
      resistivity model. It remains model-derived evidence, not a direct
      geological observation.

   Archie's law
      The petrophysical relation :math:`\rho = a\,\rho_w\,\phi^{-m}\,S_w^{-n}`
      linking formation resistivity to porosity :math:`\phi`, water saturation
      :math:`S_w`, and pore-water resistivity :math:`\rho_w`, implemented as
      :class:`pycsamt.interp.petrophysics.ArchieModel`. Most defensible in
      relatively clean, clay-poor material; see :term:`Waxman-Smits model` for
      clay-bearing formations.

   Waxman-Smits model
      A petrophysical model extending :term:`Archie's law` with a surface
      (clay) conductivity term, implemented as
      :class:`pycsamt.interp.petrophysics.WaxmanSmitsModel`. In the current
      :class:`~pycsamt.interp.hydromodel.EMHydroModel`, its parameters are
      converted to an Archie-form approximation for water-table detection and
      cell-wise inversion; ``sigma_s`` is not yet propagated through those
      inverse steps.

   Water table (hydrogeophysical)
      In :mod:`pycsamt.interp`, the shallowest depth per resistivity-model
      column where Archie-inverse water saturation first reaches
      ``Sw_water_table_threshold``. It is an operational detection threshold
      derived from resistivity and configured parameters, not a directly
      measured phreatic surface; ``nan`` marks a column where no qualifying
      transition was found.

   Transmissivity
      Aquifer transmissivity :math:`T=\int K\,dz` integrated over the
      saturated interval represented by a resistivity model column, returned
      by :class:`pycsamt.interp.hydromodel.EMHydroModel` and propagated by
      :class:`~pycsamt.interp.uncertainty.MonteCarloHydro`. It inherits all
      uncertainty in water-table detection, hydraulic conductivity, and the
      model's represented depth range.

   Dar-Zarrouk parameters
      Transverse resistance :math:`TR=\sum \rho_i h_i` and longitudinal
      conductance :math:`S=\sum h_i/\rho_i`, integrated over a resistivity
      model column's complete depth range. They summarize the whole column
      for comparative screening but do not by themselves determine aquifer
      productivity.

   Model integration
      A direct interface under :mod:`pycsamt.models` for working with a native
      inversion engine's files, runner, and result objects. It is used when the
      native project itself is part of the reproducible deliverable.

   Mesh
      The discretised numerical domain on which a forward or inverse problem is
      solved. A mesh may consist of 1-D layers, a 2-D profile grid, triangular
      finite elements, or a 3-D volume of cells, depending on the backend.

   Objective function
      The scalar quantity minimized during inversion, usually combining a data
      misfit term with a regularization term so the recovered model fits the
      observations without becoming unnecessarily rough or unstable.

   Regularization
      Constraints or penalties added to an inverse problem to stabilize
      non-unique solutions. Common examples penalize roughness, departure from a
      reference model, or implausible parameter values.

   Error floor
      A minimum uncertainty assigned to a datum or data component before
      inversion. It prevents unrealistically small formal errors from forcing
      the inversion to fit noise, processing artefacts, or modelling
      assumptions too strongly.

   Starting model
      The initial resistivity or conductivity model supplied to an inversion
      before iterations begin. It can influence convergence and should be
      archived with the run because a different starting model may reach a
      different acceptable solution.

   Residual
      The difference between observed and predicted data,
      :math:`r_i=d_{obs,i}-d_{pred,i}`, usually inspected after weighting by
      the assigned data uncertainty. Residual patterns by station, component,
      or frequency are often more informative than a single global RMS value.

   Covariance
      In inversion file formats such as ModEM, a description of how model
      cells are smoothed, linked, masked, or otherwise regularized. It is part
      of the model prior, not merely an output uncertainty table.

   Topography
      The elevation of the ground surface represented in a model or plotting
      workflow. In inversion, topography can affect mesh geometry, air cells,
      receiver elevations, and the interpretation of near-surface structure.

   Sensitivity
      The degree to which a change in a model parameter affects the predicted
      data. Low-sensitivity regions may be displayed in a model mesh but should
      not be interpreted with the same confidence as well-sampled regions.

   Resolution
      The ability of an acquisition and inversion setup to distinguish one
      subsurface feature from another. Resolution depends on frequency or time
      coverage, survey geometry, errors, regularization, and the physics of the
      selected backend.

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

   Forward operator
      The function :math:`F` mapping a resistivity model :math:`m` to
      predicted data, :math:`d_{\mathrm{pred}} = F(m)`. In :mod:`pycsamt.forward`
      it is implemented by the 1-D, 2-D, and :term:`quasi-3-D` solvers, and in
      :mod:`pycsamt.models` by external engines such as :term:`Occam2D`,
      :term:`ModEM`, and :term:`MARE2DEM`. The same operator sits inside the
      inversion objective function, so a forward assumption that does not
      match the true physics can bias the recovered model even while the data
      misfit looks small.

   Forward model
      A concrete earth model and solver setup used to compute synthetic
      electromagnetic data from assumed subsurface parameters. In validation,
      related noise realizations generated from the same forward model should
      stay in the same data partition.

   Half-space
      The lowermost, infinitely thick layer of a :term:`layered model`,
      assigned a resistivity but no thickness. It represents the electrical
      properties below the deepest resolved interface and anchors the
      long-period asymptote of the :term:`apparent resistivity` curve.

   Layered model
      A 1-D resistivity model built from horizontal layers,
      :math:`\rho(\mathbf{x}) = \rho(z)`, each carrying a resistivity and a
      thickness except the terminal :term:`half-space`. pyCSAMT's
      ``LayeredModel`` builds one from explicit values or from ``random``,
      ``blocky``, ``smooth``, and ``from_geology`` priors.

   Quasi-3-D
      A forward-modelling approximation that assembles an approximate 3-D
      tensor response from orthogonal 2-D slices through a 3-D grid, instead
      of solving the full 3-D Maxwell equations. pyCSAMT's ``MT3DForward``
      uses it for survey-scale synthetic experiments where a full production
      3-D solver would be too costly. It is a distinct idea from the
      site-level :term:`dimensionality` classification used when interpreting
      recorded MT data.

   Feature array
      A flattened numeric array built from a forward response's apparent
      resistivity and phase, produced by a response object's ``to_array`` or
      ``to_feature_array`` method. It is the data-vector shape expected by AI
      training and inversion code, as distinct from the physical per-frequency
      or per-station arrays a solver returns directly.

   Synthetic recovery
      A validation workflow that forward-models a known model, adds
      controlled noise, inverts the resulting response, and compares the
      recovered model with the known one. A successful synthetic recovery
      test shows that an inversion workflow can recover a known model under
      controlled assumptions; it does not, by itself, prove that a field
      inversion of unknown structure is correct.

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

   Forward modelling
      The calculation of a synthetic electromagnetic response from a prescribed
      earth model, survey geometry, source description, and solver setup. It is
      the reproducible "given the model, predict the data" counterpart to
      inversion, which tries to recover a model from observed data.

   Forward response
      The predicted data produced by a forward solver for a specified model and
      survey setup. Depending on method and dimensionality it may contain
      impedance, apparent resistivity, phase, transient decay values, station
      positions, and tensor components.

   Synthetic dataset
      A collection of model parameters, computed forward responses, metadata,
      and optional train/validation/test splits generated from known inputs
      rather than acquired in the field. It is useful for algorithm development
      because the target model is known.

   Feature vector
      The numeric input row passed to a learning algorithm or diagnostic plot.
      In pyCSAMT forward datasets it is usually built from transformed response
      quantities, such as log apparent resistivity followed by phase.

   Target vector
      The numeric output row that a supervised learning model is expected to
      predict. For layered-earth forward datasets it contains log-resistivities
      and layer thicknesses, with NaN padding when different samples have
      different layer counts.

   Dataset split
      A deterministic partition of a dataset into training, validation, and
      test subsets. Keeping the split seed fixed makes model-performance
      comparisons reproducible.

   Layered earth
      A one-dimensional earth model in which electrical resistivity changes only
      with depth. Each layer has a resistivity and, except for the bottom
      half-space, a thickness.

   Survey geometry
      The spatial and source-receiver arrangement used by a simulation or field
      acquisition, including station positions, profile layout, transmitter
      geometry, offsets, and dimensionality.

   Configuration file
      A persistent text file that records the parameters used by a run, such as
      solver type, frequency or time sampling, model bounds, station layout,
      noise settings, random seed, and output paths. In pyCSAMT forward
      workflows it is treated as the source of truth for rebuilding a synthetic
      dataset or response.

   Model prior
      The assumptions used before simulation or inversion to restrict plausible
      earth models. A forward-model prior may define layer-count limits,
      resistivity bounds, anomaly geometry, geological class, or spatial
      correlation length.

   AI inversion
      An inversion workflow that uses a learned or differentiable model to map
      electromagnetic observations to candidate earth-model parameters. Its
      result is conditional on the training distribution, feature contract,
      forward operator, regularization, and validation evidence.

   Supervised AI inversion
      A learned inverse mapping trained from response-model pairs, usually
      synthetic examples where both the forward response and target vector are
      known. It minimizes target error on the sampled examples and must still
      be checked in response space.

   Physics-informed inversion
      An inversion workflow that includes a differentiable physics residual in
      the optimization objective, commonly combining data fit with model
      regularization. The word physics-informed does not imply exact physics,
      uniqueness, or automatic field validity.

   Hybrid inversion
      A workflow that combines an AI estimate with physics-based refinement,
      for example by using the learned prediction as a starting model, prior,
      or proposal before iterative inversion.

   Feature contract
      The complete agreement between training and inference arrays: feature
      names, order, units, transformations, frequency or period grid, masking,
      interpolation, normalization statistics, component convention, station
      order, and padding.

   Training distribution
      The probability distribution that generated the training examples. In AI
      inversion it acts as a model prior because the network is mainly tested
      on structures, responses, noise, and nuisance effects sampled from it.

   Domain gap
      The mismatch between examples used to train or validate a model and the
      field observations where the model is applied. It can arise from physics,
      noise, survey geometry, processing, dimensionality, or geology outside
      the synthetic prior.

   Non-uniqueness
      The property that more than one earth model can fit the same
      electromagnetic observations within uncertainty. It is a physical
      limitation of the inverse problem, not a limitation that disappears
      because a neural network predicts one model quickly.

   Calibration set
      Held-out data used to calibrate uncertainty or interval coverage after a
      base model has been trained. It should not be reused for fitting network
      weights or selecting the final test result.

   Exchangeability
      The assumption that calibration examples and future examples are drawn in
      a comparable way, so their residuals can be treated as interchangeable
      for coverage calculations such as conformal prediction.

   Aleatoric uncertainty
      Uncertainty associated with observation noise, incomplete measurements,
      or irreducible ambiguity in the data.

   Epistemic uncertainty
      Uncertainty in learned model parameters or predictions caused by limited
      or incomplete training evidence.

   Distributional uncertainty
      Uncertainty caused by applying a model to inputs that differ from the
      training or calibration distribution.

   Structural uncertainty
      Uncertainty caused by an incorrect dimensionality, parameterization,
      forward physics, architecture, or geological assumption.

   Model-space metric
      A validation metric computed directly on earth-model parameters, such as
      log-resistivity error, thickness error, interface depth error, or section
      similarity when synthetic truth is known.

   Response-space metric
      A validation metric computed after forwarding the predicted model and
      comparing the reconstructed response with observed or synthetic data.

   Out-of-distribution diagnostic
      A check that estimates whether an input lies outside the distribution
      represented by training, validation, or calibration examples.

   Dataset card
      A structured documentation record describing a dataset's purpose,
      generation process, field sources, feature and target contracts, splits,
      known gaps, limitations, and intended use.

   Model card
      A structured documentation record describing a model's identity,
      intended use, training data, architecture, evaluation, uncertainty,
      limitations, and operational constraints.

   Checkpoint
      A saved model artifact containing learned parameters and enough
      architecture/backend information to restore inference. It is not
      scientifically interpretable without its feature contract, dataset
      record, checksum, and validation evidence.

   AgentResult
      The standardised return value of every :term:`agent` in
      :mod:`pycsamt.agents`: an execution ``status`` (success, failed, or
      needs_review), a human-readable summary, agent-specific arrays and
      figures under ``data``, non-fatal ``warnings``, and elapsed time and
      cost fields. Its uniform shape lets orchestration code branch on
      outcome without knowing which agent produced it, but a success status
      only means the programmed workflow returned, not that the result
      meets scientific acceptance criteria.

   Agent coordinator
      The explicit workflow runner in :mod:`pycsamt.agents` that executes a
      named sequence of registered agent steps, passes selected outputs from
      previous :term:`AgentResult`\ s into later inputs, writes resumable
      workflow checkpoints, and returns one workflow-level
      :term:`AgentResult`. It is used when the step order is already known and
      should be reproducible.

   Workflow step
      One named unit inside an :term:`agent coordinator` workflow. A step
      stores an agent instance, a stable step name, an optional input-mapping
      callback, a human-readable description, and a required/optional flag
      that controls whether failure aborts the workflow.

   Input mapping
      The explicit conversion from accumulated previous step results to the
      input dictionary expected by the next agent. In
      ``AgentCoordinator.add_step`` this is the ``input_fn`` callback. It is a
      reproducibility boundary because it documents exactly which upstream
      output keys are consumed by each downstream step.

   Workflow checkpoint
      A serialized record written after an :term:`agent coordinator` step so a
      later run can resume without recomputing completed work. pyCSAMT writes a
      pickled, figure-stripped ``AgentResult`` for execution and a JSON sidecar
      for human inspection; checkpoints accelerate resume but do not replace
      archived scientific deliverables.

   Dry run
      A non-executing preview of a command or workflow. In the pyCSAMT agent
      layer, ``dry_run=True`` returns an :term:`AgentResult` containing the
      planned step order, agent classes, required/optional flags, and LLM
      configuration without calling the registered agents' ``execute`` methods
      or writing workflow outputs. An :term:`agent coordinator` returns the
      formatted preview under ``data["plan"]``, while
      ``WorkflowOrchestratorAgent`` -- which builds a coordinator internally
      after classifying a natural-language request -- returns it under
      ``data["workflow_plan"]``; both expose the same structured
      ``data["steps"]`` list.

   Ensemble inversion
   Deep ensemble
      An uncertainty-aware :term:`supervised AI inversion` that trains
      several independent estimators from the same architecture and
      synthetic dataset with different random seeds, then reports the
      spread across members as one uncertainty source. It captures
      training-driven variability but not the :term:`calibration set`'s
      domain limits, so its empirical coverage still needs to be checked
      before it is read as a field-valid confidence interval.

   Conformal prediction
      A distribution-free calibration method that turns point predictions
      into prediction intervals with a guaranteed marginal coverage on
      held-out data, provided the :term:`calibration set` and future
      inputs are :term:`exchangeability`-compatible. It does not certify
      per-station coverage, and its guarantee degrades under
      :term:`domain gap`.

   Monte Carlo dropout
      A stochastic-inference technique that keeps a network's dropout
      layers active at prediction time and repeats the forward pass to
      obtain a spread of outputs, treated as one epistemic-uncertainty
      estimate. It is one contributor to :term:`epistemic uncertainty`
      alongside :term:`ensemble inversion`, not a substitute for
      :term:`aleatoric uncertainty` or :term:`distributional uncertainty`
      sources it does not model.

   Model zoo
      A registry of named, versioned pre-trained :term:`checkpoint`\ s
      with recorded architecture, layer count, and solver metadata, so a
      released model can be listed, downloaded, and applied without
      repeating training. Using a zoo entry still requires checking that
      its :term:`feature contract` and :term:`training distribution`
      match the field survey.

   Report package
      The versioned collection of structured outputs, cards, manifests,
      figures, metrics, predictions, review records, and rendered narrative
      used to release an AI inversion result for a declared purpose.

   Response reconstruction
      The diagnostic step that forwards a predicted earth model and compares
      the synthetic response with observed data in a declared residual space.

   Station status
      The per-station or per-profile decision attached to an inference result,
      commonly accepted, needs review, or rejected, with a reason and domain
      evidence.

   Validation leakage
      Any path by which information from validation, calibration, test, field,
      or challenge data influences model fitting, preprocessing, model
      selection, thresholds, or interpretation before those data are formally
      evaluated.

   Challenge set
      A held-out validation set deliberately shifted toward difficult or
      boundary cases, used to map failure modes and operating limits rather
      than to tune the model.

   Operating envelope
      The documented range of methods, components, frequencies, geometry,
      geology, noise, missingness, and quality conditions under which a model
      may be used with its stated acceptance evidence.

   Baseline model
      A simpler or established method used as a comparison point, such as a
      median target predictor, nearest-neighbour regression, or classical
      inversion workflow evaluated on the same held-out cases.

   Ablation study
      A controlled comparison in which one input, component, architecture
      feature, loss term, or augmentation is removed or changed to test whether
      it materially contributes to validation performance.

   Geological prior
      A named model prior tied to an expected geological setting, such as a
      geothermal, marine, sedimentary, crystalline, or permafrost target. It
      narrows synthetic-model sampling so generated examples resemble the
      intended application instead of arbitrary resistivity variation.

   Noise model
      The rule used to perturb a synthetic response so it resembles measured
      data. It defines the error distribution, scale, and sometimes
      field-style behaviour applied after the noise-free forward response is
      computed.

   Station layout
      The receiver positions used to sample a modelled response. For profile
      simulations this is usually an along-line station count and spacing; for
      map or quasi-3-D simulations it may be a two-dimensional receiver grid.

   Finite-difference grid
      A discretised numerical mesh on which derivatives in the governing
      electromagnetic equations are approximated by differences between
      neighbouring cells. Cell size, padding, and model extent control both
      numerical accuracy and boundary effects.

   Model container
      A Python object that stores the earth model and survey geometry needed by
      a forward solver, without itself solving the electromagnetic equations.
      Examples include ``LayeredModel``, ``Grid2D``, and ``Grid3D``.

   Response container
      A Python object that stores predicted fields and derived quantities from a
      forward run. Response containers keep physical arrays, coordinates, and
      feature-array helpers together so plotting, inversion handoff, and machine
      learning use the same computed result.

   Halfspace
      A uniform earth model that extends infinitely downward. In layered-earth
      notation it is the final layer, which has resistivity but no finite
      thickness.

   Padding cells
      Numerical buffer cells added outside the scientific core of a finite
      difference grid. They reduce boundary effects but should not be interpreted
      as part of the target model.

   Time gate
      One sample time in a transient electromagnetic decay curve. A TEM
      configuration uses a sequence of time gates rather than a frequency grid.

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

   Multiplicative noise
      A noise model that perturbs a response in log-space, so the added
      scatter scales with signal magnitude instead of being a fixed absolute
      value. It suits responses such as :term:`apparent resistivity` that
      span several orders of magnitude better than :term:`Gaussian noise`
      alone.

   Field-realistic noise
      A noise model that layers frequency-dependent uncertainty, AMT
      dead-band-style degradation, and :term:`powerline harmonics`-like
      contamination onto a clean synthetic response, approximating field data
      quality more closely than :term:`Gaussian noise` or
      :term:`Multiplicative noise` alone.

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

   Processing pipeline
      The ordered pyCSAMT workflow that loads a site collection, applies one or
      more configured processing steps, records per-step results, and writes
      optional processed EDI files, plots, reports, and a reproduced
      ``pipeline.yaml``.

   Pipeline configuration file
      A YAML, JSON, or trusted Python file that serialises the name, output
      directory, optional preset, ordered step list, and parameter overrides
      for a pyCSAMT processing pipeline.

   Pipeline output directory
      The root folder created for one output-enabled pipeline run. It contains
      the reproduced ``pipeline.yaml``, processed EDI files, QC figures,
      reports, and any optional intermediate snapshots for that run.

   Output artifact
      A file or in-memory record produced by a pipeline run, such as a
      processed EDI file, QC figure, report, step result, or pipeline snapshot.

   Site collection
      The ordered in-memory group of EDI-like station objects passed between
      pyCSAMT site, emtools, and pipeline functions. Pipeline steps treat this
      collection as the current survey state.

   Pipeline preset
      A named built-in processing workflow, such as ``basic_qc`` or
      ``publication_ready``, that expands to an ordered list of registered
      pipeline steps with default parameters.

   Preset expansion
      The act of turning a named pipeline preset into its explicit ordered
      list of step labels, registry codes, and parameters. Expansion is useful
      before review because it exposes the exact recipe that will run.

   Preset comparison run
      A controlled set of pipeline runs where the same input site collection is
      processed with different presets and separate output roots so reports,
      plots, processed EDIs, and pipeline snapshots can be compared.

   Pipeline step code
      A stable short identifier for a registered processing operation, such as
      ``NR001`` or ``QC001``. Step codes are accepted by the CLI and config
      files so workflows can refer to operations without depending on display
      labels.

   Step label
      The user-facing name assigned to one occurrence of a pipeline step inside
      a configuration file or constructed pipeline. It can describe the role of
      the occurrence, for example ``trim_to_amt_band``, even when the registry
      operation is still ``FREQ001``.

   Parameter override
      A value supplied in a pipeline configuration that replaces the registered
      default for the same step parameter while leaving unspecified defaults
      unchanged.

   Step registry
      The catalogue of processing operations known to the pipeline engine. Each
      registry entry defines a code, name, category, default parameters,
      callable transform, optional plot functions, and whether the operation
      returns a modified site collection.

   StepSpec
      The registry metadata object for one pipeline operation. It binds the
      stable step code, snake-case name, category, transform function, default
      parameters, QC plot functions, and diagnostic-vs-transform behavior.

   Transform step
      A pipeline step whose function returns a modified site collection that
      becomes the input to the next step.

   Diagnostic step
      A pipeline step used for QC, figures, or summaries that intentionally
      passes the site collection through unchanged.

   QC plot function
      A plotting callable attached to a step registry entry. After a successful
      step transform, the pipeline calls each registered QC plot function and
      saves any returned Matplotlib figures when plot output is enabled.

   Step ordering
      The scientific and computational sequence in which pipeline steps are
      applied. Ordering matters because every transform receives the survey
      state produced by the previous transform.

   Canonical pipeline snapshot
      The normalised YAML, JSON, or Python representation exported from a
      constructed pipeline. It records the resolved step order and parameters
      so a reviewer can compare what was intended with what was run.

   Processed EDI
      An EDI file written after a processing pipeline has transformed a site
      collection. It should be stored separately from raw EDI files because it
      reflects filtering, editing, correction, or QC decisions.

   QC figure
      A diagnostic Matplotlib figure produced by a registered pipeline step to
      show whether the transformed data remain physically plausible and useful
      for the next processing or inversion stage.

   PipelineResult
      The Python object returned by ``Pipeline.run``. It carries the input and
      output site collections, per-step results, saved paths, output root,
      runtime, and overall success state.

   StepResult
      The per-step record stored inside a ``PipelineResult``. It records the
      step index, user label, registry code, parameters, elapsed time, site
      counts, saved plot paths, and any captured error.

   Active survey context
      The survey path remembered by the survey CLI for the current project or
      session. It is convenient for interactive command sequences, but an
      explicit positional path or ``--survey`` option is clearer for
      reproducible scripts.

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

   EDI-like object
      Any Python object that behaves like an EDI station record for pyCSAMT site
      tools. At minimum, computed diagnostics expect a ``get_section`` method
      and a ``Z`` section exposing frequency and impedance arrays; tipper
      diagnostics also look for ``Tip``, ``TIP``, ``T``, or ``Z``-attached
      tipper arrays.

   Off-diagonal component
      One of the cross-coupled impedance tensor elements :math:`Z_{xy}` or
      :math:`Z_{yx}`. These components usually carry the primary TE/TM
      information for 1-D and 2-D MT-style interpretation, while diagonal
      components are expected to be small after rotation to strike.

   Native frequency
      A frequency value that is already present in a station's recorded
      frequency grid, before interpolation or resampling to a requested common
      comparison frequency.

   Frequency decade
      A factor-of-ten interval in frequency. A slope reported in degrees per
      decade means the fitted phase change for each unit increase in
      :math:`\log_{10}(f)`.

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
      Root-mean-square misfit between observed and predicted responses,
      normalised by data error. Writing the :term:`residual` at datum
      :math:`i` as :math:`r_i` and its assigned uncertainty as
      :math:`\sigma_i`,

      .. math::

         \mathrm{RMS} = \sqrt{\frac{1}{N}\sum_{i=1}^{N}
             \left(\frac{r_i}{\sigma_i}\right)^2},

      over the :math:`N` fitted data. It is the primary goodness-of-fit
      measure for an inversion; a value near :math:`1` means the fit is
      consistent with the assigned :term:`error floor` and uncertainties,
      while a much larger value points to underestimated errors, a poor
      starting model, or physics the :term:`forward operator` cannot
      represent.

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

   RAG
   Retrieval-augmented generation
      Grounding a model's answer in evidence retrieved from a corpus at
      query time, rather than relying only on parameters learned during
      training. :mod:`pycsamt.assistant.rag` implements the pyCSAMT
      instance: an offline, deterministic :term:`BM25` retriever over an
      indexed :term:`chunk` corpus, with dense embeddings as an optional
      addition rather than a requirement.

   Chunk
      One retrievable unit of the RAG corpus, represented by
      :class:`pycsamt.assistant.rag.schemas.RAGChunk`. A chunk carries
      searchable text plus metadata -- source path, ``kind`` (for example
      ``python_symbol``, ``doc_section``, or ``recipe``), an optional
      workflow tag, a static :term:`priority`, and a stable id derived from
      its source location -- so a retrieval result can be ranked, filtered,
      and cited back to a real file and line range.

   Priority
      A static importance score attached to a :term:`chunk` from its
      source path alone (``pycsamt/agents``, ``pycsamt/emtools``, and
      similar high-value implementation paths score highest, general
      documentation next, everything else the baseline). It is a ranking
      boost applied at retrieval time, not a measure of a chunk's textual
      relevance to any specific query.

   BM25
   Okapi BM25
      A term-frequency ranking function that scores how well a document
      matches a query without semantic embeddings, used as the default
      lexical ranker in :mod:`pycsamt.assistant.rag`. For a document
      (:term:`chunk`) :math:`c` and query terms :math:`t`, pyCSAMT's
      implementation scores

      .. math::

         \mathrm{BM25}(c, q) = \sum_{t \in q}
         \mathrm{idf}(t)\,
         \frac{f_{t,c}\,(k_1 + 1)}
              {f_{t,c} + k_1\left(1 - b + b\,\dfrac{|c|}{\mathrm{avgdl}}\right)},
         \qquad
         \mathrm{idf}(t) = \ln\!\left(1 + \frac{N - n_t + 0.5}{n_t + 0.5}\right),

      where :math:`f_{t,c}` is how often :math:`t` occurs in :math:`c`,
      :math:`|c|` is the chunk's token count, :math:`\mathrm{avgdl}` is the
      corpus's average chunk length, :math:`N` is the corpus size, and
      :math:`n_t` is the number of chunks containing :math:`t`. pyCSAMT
      uses the conventional :math:`k_1=1.5`, :math:`b=0.75`. A chunk that
      shares no term with the query scores exactly ``0`` and is dropped
      before any boost is applied.

   Query expansion
      A deterministic, zero-dependency substitute for semantic matching in
      :mod:`pycsamt.assistant.rag`: a fixed table maps specific
      natural-language trigger phrases (for example "vertical offset") to
      extra domain terms that are known to exist in the corpus vocabulary
      (``static``, ``shift``, ``galvanic``). Expansion terms are scored by
      the same :term:`BM25` function but added at a fixed ``0.35`` weight,
      so they can surface the right :term:`chunk` when a user's wording and
      the code's vocabulary do not overlap, without ever outweighing the
      user's own words.

   Reciprocal Rank Fusion
   RRF
      A rank-based (not score-based) method for combining two or more
      ranked lists into one, used by :mod:`pycsamt.assistant.rag` to blend
      lexical (:term:`BM25`) and optional dense retrieval. Each ranking
      :math:`r` contributes

      .. math::

         \mathrm{RRF}(i) = \sum_r \frac{w_r}{k + \mathrm{rank}_r(i)},

      summed over the rankings :math:`i` appears in (0-based rank; absence
      from a ranking contributes nothing), with pyCSAMT's default
      :math:`k=60` and equal weights :math:`w_r=1`. Rank-based fusion sidesteps
      the problem of combining BM25 scores and cosine similarities directly:
      the two live on incomparable numeric scales, while ranks do not.

   Confidence proxy
      A single numeric stand-in for "does the assistant actually know the
      answer," computed as the top :term:`chunk`'s post-boost retrieval
      score (``RetrievedContext.top_score``) with no query-independent
      calibration. :class:`~pycsamt.assistant.rag.context_builder.AssembledContext`
      is considered confident when project context was resolved or this
      score reaches a fixed floor (``25.0``); below that, with no session
      context to fall back on, the assistant asks a clarifying question
      instead of guessing.

   Hallucination guard
      An evaluation check that flags a forbidden string -- typically an
      invented import or API name -- appearing anywhere in the text of the
      :term:`chunk`\ s retrieved for a query. It is declared per-record via
      a suite's ``must_not_contain`` list and reported as a hard violation
      count by :mod:`pycsamt.assistant.evals`, distinct from the
      generated-code validator, which checks a specific script's imports
      rather than the retrieved evidence itself.

   Corpus fingerprint
      A SHA-256 digest of the indexed source tree used to detect a stale
      persisted RAG index: every indexable file is hashed individually,
      and the fingerprint folds each file's repository-relative path and
      content hash into one running digest over paths sorted
      deterministically. Because it hashes content rather than
      modification time, a tool that rewrites a file byte-for-byte does
      not trigger a false "stale index" warning; any real change to
      indexed source, documentation, examples, or recipes does.

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

   Station identity
      The normalized name pyCSAMT assigns to one site container. It is
      resolved from EDI ``HEAD`` fields in a fixed order --
      ``dataid``, ``station``, ``sitename``, ``name``, ``STATION``,
      falling back to the file stem when none are present -- and, once
      resolved, is written back into ``dataid`` (and ``station`` when
      absent) so that later name-based lookups and joins see one
      consistent label per station.

   Geodetic distance
      The great-circle separation between two ``(lat, lon)`` points on
      a spherical Earth of radius :math:`R`, used by station
      nearest-neighbour search. pyCSAMT evaluates the haversine formula

      .. math::

         d = 2R \arcsin\left(\sqrt{\sin^2\tfrac{\Delta\phi}{2}
             + \cos\phi_1\cos\phi_2\sin^2\tfrac{\Delta\lambda}{2}}\right)

      with latitudes and longitudes in radians, returning a distance in
      metres.

   Chainage
      The signed distance of a station from a chosen origin, measured
      along a survey line's azimuth :math:`A` rather than along
      latitude/longitude directly. Writing local east/north offsets
      from the origin as :math:`dx` and :math:`dy`, chainage is

      .. math::

         s = dx\sin A + dy\cos A,

      so stations ahead of the origin along the line get positive
      values and stations behind it get negative values. It is the
      coordinate that :term:`profile line` construction sorts stations
      by, and it is distinct from :term:`station distance`, which
      accumulates separation without reference to a single azimuth.

   Native file
      A file format read or written directly by an external modelling or
      inversion engine, such as an :term:`Occam2D` mesh or a :term:`ModEM`
      model file. It differs from a pyCSAMT :term:`configuration file` in
      that it is the exact record the engine's binary consumed or produced,
      not the editable parameters used to build it.

   Interpretation package
      A controlled set of interpretation deliverables, usually including source
      run identifiers, configuration, evidence tables, exported grids or logs,
      figures, narrative report text, review status, checksums, and a
      :term:`provenance manifest`. It is the unit a reviewer or client can
      audit and, when approved, archive.

   LAS
      Log ASCII Standard, a text format widely used for well-log curves. In
      pyCSAMT interpretation exports, LAS files can carry EM-derived station
      depth curves, but they should not be described as drilled well logs unless
      the station is actually tied to a borehole and validated as such.

   VTK
      Visualization Toolkit file format family. pyCSAMT interpretation exports
      use an ASCII rectilinear-grid VTK file for model resistivity values so the
      grid can be opened in tools such as ParaView or GIS viewers that support
      VTK.

   Startup file
      The :term:`Occam2D` control file that names the paired data, mesh, and
      model files and sets the roughness type, initial resistivity, and
      iteration limits for the first solver call. pyCSAMT's
      ``OccamStartup.from_model`` writes it as the last step of
      ``InputBuilder.build``, after the data, mesh, and model files.

   Control file
      The :term:`ModEM` settings file that schedules the non-linear
      conjugate-gradient search: starting and exit values for the trade-off
      parameter :math:`\lambda`, the step-size scale :math:`\alpha`, the
      smoothing-iteration count, and the maximum iterations and RMS
      tolerance that stop the run.

   Response file
      The engine-written file holding predicted data for the current model,
      such as an Occam2D ``.resp`` file or a MARE2DEM ``.EMResp`` file. It is
      compared against the observed data file to obtain the :term:`RMS
      misfit`; loading it without also inspecting the misfit and residual
      pattern can hide a model that fits poorly in specific components.

   Iteration file
      A per-iteration engine record, such as an Occam2D ``.iter`` file,
      pairing one candidate model with its :term:`RMS misfit` and
      roughness at that step. Keeping the full sequence rather than only the
      final iteration makes it possible to check convergence behaviour and
      to roll back to an earlier, better-regularized model.

   Roughness
      A scalar penalty on how much a model changes between adjacent cells,
      minimized alongside the :term:`RMS misfit` in an
      :term:`objective function`. In :term:`Occam2D`, roughness typically
      *rises* across iterations while RMS falls, because the smoothest
      model that still fits the target misfit is usually rougher than the
      near-uniform starting half-space; a run whose roughness keeps
      climbing after RMS has already reached target is a sign to stop
      rather than keep iterating.

   Lagrange multiplier
      The trade-off weight between data misfit and :term:`roughness` at one
      inversion step, analogous to the ModEM :term:`control file`'s
      :math:`\lambda`. Occam reduces it during a line search --
      ``stepsize_cut_count`` bounds how many reductions are tried -- until
      the resulting model improves the objective function.

   Type code
      The integer tag on an :term:`Occam2D` data row identifying which
      response component it carries: ``1``/``2`` for TE apparent
      resistivity/phase and ``5``/``6`` for TM apparent resistivity/phase.
      Plot and diagnostic helpers select rows by type code rather than by
      column position.

   NLCG
      Nonlinear conjugate-gradient, the iterative search algorithm
      :term:`ModEM` uses to minimize its :term:`objective function`. Each
      step moves along a search direction built from the current gradient
      and the previous step, scaled by the trade-off parameter tracked in
      the :term:`control file`; ``Modular_NLCG`` in ModEM's log and output
      filenames names this algorithm, not the survey or model.

   PSLG
      Planar straight-line graph: the boundary representation Triangle-based
      finite-element :term:`mesh` generators build from, given as nodes,
      connecting segments, and interior region seed points. :term:`MARE2DEM`
      reads and writes this geometry as a ``.poly`` file before it (or an
      external Triangle call) refines it into the ``.node``/``.ele`` files
      an inversion actually solves on.

   Free parameter
      A resistivity region in a :term:`MARE2DEM` model that the inversion is
      allowed to adjust, as opposed to a fixed region such as air, ocean, or
      a boundary padding cell held at a reference value. A region file lists
      every region regardless of type, so filtering to free parameters only
      is often necessary before a resistivity histogram or summary statistic
      is meaningful.
