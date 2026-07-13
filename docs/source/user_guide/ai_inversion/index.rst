.. _user_guide_ai_inversion:

AI inversion
============

The AI inversion guide documents the complete learned-inversion lifecycle in
pyCSAMT: scientific assumptions, training-data generation, architecture
selection, training, inference, validation, uncertainty, PINN workflows,
agent-assisted execution, and reporting.

AI inversion is a peer of :doc:`../forward/index` and
:doc:`../models/index`. Forward modeling generates responses from known earth
models; classical inversion estimates models through numerical optimization;
AI inversion learns an inverse mapping from representative examples. These
approaches can support one another, but their assumptions and validation
requirements remain distinct.

.. admonition:: Validation is part of the model
   :class: important

   A network prediction is not trustworthy merely because inference is fast or
   a training loss is small. Field use requires representative training data,
   strict data separation, out-of-distribution checks, physical diagnostics,
   uncertainty assessment, and comparison with independent evidence.

Build the workflow
------------------

.. grid:: 1 1 2 3
   :gutter: 3
   :class-container: cta-tiles

   .. grid-item-card:: Concepts
      :link: concepts
      :link-type: doc
      :img-top: ../../_static/icons/ai.svg
      :class-card: pycsamt-card sd-text-center

      Understand surrogate inversion, supported dimensions, assumptions,
      limitations, and the relationship to forward and classical inversion.

   .. grid-item-card:: Data preparation
      :link: data_preparation
      :link-type: doc
      :img-top: ../../_static/icons/user-guide-data-loading.svg
      :class-card: pycsamt-card sd-text-center

      Build synthetic and field datasets, define features and targets, normalize
      consistently, split safely, and preserve provenance.

   .. grid-item-card:: Model selection
      :link: model_selection
      :link-type: doc
      :img-top: ../../_static/icons/model-builder.svg
      :class-card: pycsamt-card sd-text-center

      Choose 1-D, profile-based 2-D, graph-based 3-D, or physics-informed
      approaches according to geometry, data volume, and scientific goals.

   .. grid-item-card:: Training
      :link: training
      :link-type: doc
      :img-top: ../../_static/icons/cpu.svg
      :class-card: pycsamt-card sd-text-center

      Configure optimization, monitor learning, resume safely, diagnose
      overfitting, and retain reproducible checkpoints and histories.

   .. grid-item-card:: Inference
      :link: inference
      :link-type: doc
      :img-top: ../../_static/icons/launch.svg
      :class-card: pycsamt-card sd-text-center

      Apply a trained model to compatible observations, verify preprocessing,
      detect unsupported inputs, and export predictions with metadata.

   .. grid-item-card:: Validation
      :link: validation
      :link-type: doc
      :img-top: ../../_static/icons/external-validation.svg
      :class-card: pycsamt-card sd-text-center

      Evaluate held-out synthetic tests, field responses, physical consistency,
      baselines, failure modes, and independent geological evidence.

Review and extend
-----------------

.. grid:: 1 1 2 2
   :gutter: 3
   :class-container: cta-tiles

   .. grid-item-card:: Uncertainty
      :link: uncertainty
      :link-type: doc
      :img-top: ../../_static/icons/diagnostic.svg
      :class-card: pycsamt-card sd-text-center

      Quantify predictive spread, calibration, ensemble behavior, domain shift,
      and uncertainty sources not captured by the network.

   .. grid-item-card:: PINN 2-D
      :link: pinn_2d
      :link-type: doc
      :img-top: ../../_static/icons/inversion-concepts.svg
      :class-card: pycsamt-card sd-text-center

      Add physical residuals and constraints to learned 2-D inversion while
      keeping numerical assumptions and validation explicit.

   .. grid-item-card:: AI inversion agents
      :link: agents
      :link-type: doc
      :img-top: ../../_static/icons/user-guide-agents.svg
      :class-card: pycsamt-card sd-text-center

      Use pyCSAMT agents to configure and coordinate AI workflows without
      hiding the underlying science API or review requirements.

   .. grid-item-card:: Reporting
      :link: reporting
      :link-type: doc
      :img-top: ../../_static/icons/export.svg
      :class-card: pycsamt-card sd-text-center

      Package datasets, configurations, checkpoints, metrics, predictions,
      uncertainty, limitations, and approvals into an auditable model record.

Recommended reading order
-------------------------

For a new AI inversion project:

#. Read :doc:`concepts` and define why AI inversion is appropriate.
#. Prepare representative data with :doc:`data_preparation`.
#. Choose an architecture using :doc:`model_selection`.
#. Follow :doc:`training` and preserve every resolved setting.
#. Complete :doc:`validation` before field inference.
#. Apply the model with :doc:`inference` and evaluate :doc:`uncertainty`.
#. Use :doc:`reporting` before releasing predictions.

Use :doc:`pinn_2d` when the workflow explicitly includes physics-informed
constraints. Use :doc:`agents` only after understanding the corresponding
manual science workflow.

Contents
--------

.. toctree::
   :maxdepth: 1

   concepts
   data_preparation
   model_selection
   training
   inference
   validation
   uncertainty
   pinn_2d
   agents
   reporting
