AI And Model-Zoo Agents
=======================

These agents wrap neural inversion, uncertainty, anomaly detection, and
pre-trained model access.  They may require optional AI dependencies and model
checkpoints.

.. _agent-ai-inversion:

AIInversionAgent
----------------

``AIInversionAgent`` runs end-to-end 1-D AI inversion.  It can be
instantiated directly or created from a model-zoo checkpoint when available.

.. code-block:: python
   :linenos:

   from pycsamt.agents import AIInversionAgent

   agent = AIInversionAgent.from_pretrained("mt1d-resnet-5layer-v1")
   result = agent.execute({
       "path": "/data/WILLY_EDIs",
       "output_dir": "/out/willy_ai1d",
   })

   print(result.get("rms_global"))

.. _agent-inv2d:

Inv2DAgent
----------

``Inv2DAgent`` performs 2-D profile inversion using U-Net style models.  Use
it when lateral continuity along a profile is important and a compatible model
or training setup is available.

.. code-block:: python
   :linenos:

   result = Inv2DAgent().execute({
       "path": "/data/WILLY_profile",
       "output_dir": "/out/willy_inv2d",
   })

.. _agent-inv3d:

Inv3DAgent
----------

``Inv3DAgent`` performs 3-D spatial AI inversion with graph-based models.
Use it when inter-station geometry and spatial relationships are part of the
inversion target.

.. code-block:: python
   :linenos:

   result = Inv3DAgent().execute({
       "path": "/data/WILLY_grid",
       "output_dir": "/out/willy_inv3d",
   })

.. _agent-ensemble:

EnsembleAgent
-------------

``EnsembleAgent`` runs ensemble inversion workflows and summarizes uncertainty
from multiple predictions.  Use it when uncertainty bands or coverage metrics
are required alongside a predicted resistivity model.

.. code-block:: python
   :linenos:

   result = EnsembleAgent(n_estimators=5, epochs=50).execute({
       "path": "/data/WILLY_EDIs",
       "output_dir": "/out/willy_ensemble",
   })

   print(result.get("coverage"))

.. _agent-joint-inversion:

JointInversionAgent
-------------------

``JointInversionAgent`` runs multi-modal inversion, for example MT with TEM or
CSAMT with supporting data.  Use it when multiple geophysical modalities
constrain the same subsurface target.

.. code-block:: python
   :linenos:

   result = JointInversionAgent(modalities=["mt", "tem"]).execute({
       "path": "/data/WILLY_MT",
       "secondary_path": "/data/WILLY_TEM",
       "output_dir": "/out/willy_joint",
   })

.. _agent-anomaly-detection:

AnomalyDetectionAgent
---------------------

``AnomalyDetectionAgent`` flags anomalous data patterns.  Use it for
station-frequency anomaly screening, survey triage, or finding regions that
need manual review before inversion.

.. code-block:: python
   :linenos:

   anomalies = AnomalyDetectionAgent().execute({
       "path": "/data/WILLY_EDIs",
       "output_dir": "/out/willy_anomalies",
   })

.. _agent-model-zoo:

ModelZooAgent
-------------

``ModelZooAgent`` lists available pre-trained models, downloads checkpoints,
and runs predictions where supported.

.. code-block:: python
   :linenos:

   from pycsamt.agents import ModelZooAgent

   zoo = ModelZooAgent()

   models = zoo.execute({"action": "list"})
   print(models["models"])

   checkpoint = zoo.execute({
       "action": "download",
       "model_name": "mt1d-resnet-5layer-v1",
   })
   print(checkpoint.get("checkpoint_path"))

AI Workflow Pattern
-------------------

.. code-block:: text

   MTLoaderAgent -> DataQCAgent -> DenoisingAgent
   -> AIInversionAgent or Inv2DAgent or Inv3DAgent
   -> InterpretationAgent -> ReportAgent

Add ``EnsembleAgent`` when uncertainty is part of the objective.  Add
``JointInversionAgent`` when a secondary modality is available.
:html_theme.sidebar_secondary.remove:
