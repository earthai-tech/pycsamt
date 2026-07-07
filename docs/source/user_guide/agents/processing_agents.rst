Processing And Diagnostics Agents
=================================

These agents operate after data loading and before modelling or reporting.
They clean, diagnose, transform, and prepare survey data for interpretation or
inversion.

.. _agent-data-qc:

DataQCAgent
-----------

``DataQCAgent`` evaluates survey quality after loading.  It is the normal
first processing step after ``MTLoaderAgent`` and before correction or
analysis.

Use it for station scoring, frequency coverage checks, SNR-like indicators,
dead-band detection, outlier flags, and QC summaries.

.. code-block:: python
   :linenos:

   from pycsamt.agents import MTLoaderAgent, DataQCAgent

   loaded = MTLoaderAgent().execute({"path": "/data/WILLY_EDIs"})
   qc = DataQCAgent().execute({
       "sites": loaded["sites"],
       "output_dir": "/out/willy/qc",
   })

   print(qc.summary)
   print(qc.get("qc_table"))

.. _agent-static-shift:

StaticShiftAgent
----------------

``StaticShiftAgent`` detects and corrects static-shift effects in apparent
resistivity responses.  It belongs after QC and before phase analysis or
inversion preparation.

Use it when station responses show amplitude offsets that may be caused by
near-surface conductivity heterogeneity.

.. code-block:: python
   :linenos:

   corrected = StaticShiftAgent().execute({
       "sites": qc["sites"],
       "output_dir": "/out/willy/static_shift",
   })

   print(corrected.summary)

.. _agent-phase-analysis:

PhaseAnalysisAgent
------------------

``PhaseAnalysisAgent`` runs phase-tensor and dimensionality diagnostics.  It
is useful for deciding whether a 1-D, 2-D, or 3-D interpretation is justified
and for preparing inversion orientation choices.

Typical products include phase-tensor summaries, strike estimates, skew or
dimensionality indicators, Mohr-style diagnostics, and Argand-style
representations.

.. code-block:: python
   :linenos:

   phase = PhaseAnalysisAgent().execute({
       "sites": corrected.get("corrected_sites", corrected.get("sites")),
       "output_dir": "/out/willy/phase",
   })

   print(phase.summary)
   print(phase.llm_interpretation)

.. _agent-tensor-rotation:

TensorRotationAgent
-------------------

``TensorRotationAgent`` rotates impedance tensors into a selected coordinate
frame.  It is useful after strike analysis, before export, or before inversion
preparation when a consistent orientation is required.

Use it when a workflow has an estimated strike angle or a target coordinate
frame that should be applied consistently across stations.

.. code-block:: python
   :linenos:

   rotated = TensorRotationAgent().execute({
       "sites": loaded["sites"],
       "angle_deg": 35.0,
       "output_dir": "/out/willy/rotated",
   })

.. _agent-tipper-analysis:

TipperAnalysisAgent
-------------------

``TipperAnalysisAgent`` analyzes vertical magnetic transfer functions where
tipper data are available.  Use it for induction-arrow products and
tipper-amplitude or phase diagnostics that complement impedance analysis.

.. code-block:: python
   :linenos:

   tipper = TipperAnalysisAgent().execute({
       "sites": loaded["sites"],
       "periods": [0.01, 0.1, 1.0, 10.0],
       "output_dir": "/out/willy/tipper",
   })

.. _agent-frequency-decimation:

FrequencyDecimationAgent
------------------------

``FrequencyDecimationAgent`` selects a stable subset of frequencies or
periods.  It is useful before inversion when data are oversampled, unevenly
sampled, or affected by dead bands.

Use it to produce a compact period set for Occam2D, ModEM, AI inversion, or
publication plots.

.. code-block:: python
   :linenos:

   decimated = FrequencyDecimationAgent().execute({
       "sites": loaded["sites"],
       "n_periods": 32,
       "output_dir": "/out/willy/decimated",
   })

.. _agent-denoising:

DenoisingAgent
--------------

``DenoisingAgent`` removes or suppresses noise before analysis or inversion.
It can be used in classical workflows before phase analysis and in AI
workflows before neural inversion.

Use it when QC indicates noisy bands, outlier stations, or responses that need
robust preprocessing.

.. code-block:: python
   :linenos:

   denoised = DenoisingAgent().execute({
       "sites": loaded["sites"],
       "method": "rpca",
       "output_dir": "/out/willy/denoised",
   })

Processing Order
----------------

A common processing chain is:

.. code-block:: text

   MTLoaderAgent -> DataQCAgent -> StaticShiftAgent
   -> PhaseAnalysisAgent -> FrequencyDecimationAgent

Add ``TensorRotationAgent`` after strike estimation when inversion or export
requires a consistent coordinate frame.  Add ``DenoisingAgent`` before
analysis or inversion when QC flags noisy data.
