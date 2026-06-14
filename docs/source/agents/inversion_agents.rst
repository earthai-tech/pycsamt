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
   -> Occam2DAgent or ModEmAgent
   -> InversionEvaluationAgent

Use ``InversionComparisonAgent`` when multiple model outputs need to be
compared and documented.
