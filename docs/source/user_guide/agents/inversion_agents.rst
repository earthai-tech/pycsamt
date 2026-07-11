Forward And Inversion Workflow Agents
=====================================

These agents prepare, run, evaluate, or compare physics-based modelling and
inversion workflows.

.. _agent-forward-model:

ForwardModelAgent
-----------------

``ForwardModelAgent`` runs MT forward modelling from resistivity models.  Use
it for synthetic checks, sensitivity experiments, or validating how a proposed
model should appear in impedance or apparent-resistivity space.

.. code-block:: python
   :linenos:

   from pycsamt.agents import ForwardModelAgent

   result = ForwardModelAgent().execute({
       "model": layered_model,
       "periods": [0.01, 0.1, 1.0, 10.0],
       "output_dir": "/out/forward",
   })

.. _agent-inversion-prep:

InversionPrepAgent
------------------

``InversionPrepAgent`` is the general pre-inversion interface.  It prepares
MT data products for inversion codes and is useful when a workflow should be
backend-agnostic or when the exact inversion target is decided later.

.. code-block:: python
   :linenos:

   result = InversionPrepAgent().execute({
       "sites": processed_sites,
       "output_dir": "/out/inversion_prep",
   })

.. _agent-occam2d:

Occam2DAgent
------------

``Occam2DAgent`` is the specialized Occam2D preparation agent.  It writes the
data, mesh, startup, and related files needed to start an Occam2D workflow.

Use it after QC, correction, optional rotation, and frequency selection.

.. code-block:: python
   :linenos:

   occam = Occam2DAgent().execute({
       "sites": processed_sites,
       "output_dir": "/out/willy_occam2d",
       "run_external": False,
   })

.. _agent-modem:

ModEmAgent
----------

``ModEmAgent`` prepares ModEM 3-D impedance data files.  Use it when the
survey and interpretation target require 3-D inversion products.

.. code-block:: python
   :linenos:

   modem = ModEmAgent().execute({
       "sites": processed_sites,
       "output_dir": "/out/willy_modem",
   })

.. _agent-mare2dem:

Mare2DEMAgent
-------------

``Mare2DEMAgent`` orchestrates MARE2DEM 2.5-D EM inversion preparation and,
when a compiled binary is available, execution.  It is the specialized agent
to use when the target workflow needs MARE2DEM ``.emdata``, ``.resistivity``,
``.settings``, and run-directory products rather than Occam2D or ModEM files.

The agent can work in three modes:

``"prepare"``
    Write input files only.  This is the safest default and is appropriate
    when the run will be launched manually on a workstation or cluster.

``"run"``
    Write inputs and launch the configured MARE2DEM binary.  Use this only
    after confirming the binary, MPI configuration, and source tree.

``"report"``
    Inspect an existing run directory and summarize result state without
    writing new inputs or launching the binary.

Typical input keys include ``emdata`` for an existing data file, ``mt`` or
``csem`` dictionaries for building survey data, optional ``topo``, optional
``resistivity``, ``output_dir``, ``source_dir``, ``n_procs``,
``max_iterations``, ``target_rms``, and ``initial_rho``.

.. code-block:: python
   :linenos:

   from pycsamt.agents import Mare2DEMAgent

   prepared = Mare2DEMAgent(n_procs=8).execute({
       "emdata": "/data/willy.emdata",
       "output_dir": "/out/willy_mare2dem",
       "mode": "prepare",
       "initial_rho": 100.0,
       "target_rms": 1.0,
   })

   print(prepared["data_path"])
   print(prepared["settings_path"])
   print(prepared["binary_found"])

When building from survey parameters, pass the MARE2DEM survey configuration
through ``mt`` or ``csem``:

.. code-block:: python
   :linenos:

   import numpy as np
   from pycsamt.agents import Mare2DEMAgent

   result = Mare2DEMAgent(n_procs=16).execute({
       "mt": {
           "frequencies": list(np.logspace(-3, 3, 20)),
           "rx_y": list(np.linspace(-5000, 5000, 20)),
           "rx_type": "land",
           "lTE": True,
           "lTM": True,
       },
       "topo": 0.0,
       "output_dir": "/out/synthetic_mare2dem",
       "mode": "prepare",
   })

Use ``Mare2DEMAgent`` after QC, static-shift correction, tensor rotation, and
frequency selection.  Use ``InversionEvaluationAgent`` afterwards when a run
directory contains results that need RMS, convergence, or residual summaries.

.. _agent-inversion-backend:

InversionBackendAgent
---------------------

``InversionBackendAgent`` connects agent workflows to pyCSAMT inversion
backends.  It is the right place for workflows that should select or execute a
configured backend rather than only write preparation files.

.. code-block:: python
   :linenos:

   backend = InversionBackendAgent().execute({
       "sites": processed_sites,
       "backend": "occam2d",
       "output_dir": "/out/backend_run",
   })

.. _agent-inversion-evaluation:

InversionEvaluationAgent
------------------------

``InversionEvaluationAgent`` evaluates inversion results.  Use it after an
Occam2D, ModEM, or backend run to summarize RMS, residuals, misfit sections,
and quality indicators.

.. code-block:: python
   :linenos:

   evaluation = InversionEvaluationAgent().execute({
       "result_path": "/out/willy_occam2d",
       "observed_sites": processed_sites,
       "output_dir": "/out/willy_eval",
   })

.. _agent-inversion-comparison:

InversionComparisonAgent
------------------------

``InversionComparisonAgent`` compares inversion sections or model outputs.
Use it for before/after studies, backend comparisons, parameter sweeps, or
checking whether correction changed the interpreted structure.

.. code-block:: python
   :linenos:

   comparison = InversionComparisonAgent().execute({
       "reference": "/out/run_a/model.dat",
       "candidate": "/out/run_b/model.dat",
       "output_dir": "/out/comparison",
   })

Typical Physics-Based Chain
---------------------------

.. code-block:: text

   MTLoaderAgent -> DataQCAgent -> StaticShiftAgent
   -> PhaseAnalysisAgent -> FrequencyDecimationAgent
   -> Occam2DAgent or ModEmAgent or Mare2DEMAgent
   -> InversionEvaluationAgent

Use ``InversionComparisonAgent`` when multiple model outputs need to be
compared and documented.
:html_theme.sidebar_secondary.remove:
