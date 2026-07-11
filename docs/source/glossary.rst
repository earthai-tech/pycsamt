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

   CSEM
      Controlled-source electromagnetics — the broader family of active-source
      EM methods to which :term:`CSAMT` belongs.

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

   Apparent resistivity
      The resistivity of an equivalent uniform half-space that would produce the
      observed impedance at a given frequency, :math:`\rho_a = 0.2\,|Z|^2 / f`
      (with :math:`Z` in field units and :math:`f` in Hz). It is *apparent*
      because the real earth is layered/heterogeneous.

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

   Agent
      In :mod:`pycsamt.agents`, a small composable unit that performs one step of
      a workflow (routing, parsing, QC, forward modelling, inversion, reporting)
      and returns a standardised ``AgentResult``.
