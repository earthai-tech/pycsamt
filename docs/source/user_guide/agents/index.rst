.. _agents:

AI Agents
=========

pyCSAMT agents provide AI-assisted and rule-based workflow automation for
survey loading, quality control, static-shift correction, phase analysis,
inversion preparation, interpretation, reporting, and orchestration. Start
with the orchestrating chat agent for guided end-to-end runs, browse the
catalogue when you need one specialised agent, and use the family pages
for the full reference of every agent group.

Start Here
----------

.. grid:: 1 1 2 2
   :gutter: 3
   :class-container: cta-tiles

   .. grid-item-card:: Agent Master
      :link: orchestrator
      :link-type: doc
      :img-top: ../../_static/icons/user-guide-agents.svg
      :class-card: pycsamt-card sd-text-center

      The orchestrating chat agent — plans and runs multi-step workflows
      end to end, grounded by the RAG assistant.

   .. grid-item-card:: Agent overview
      :link: overview
      :link-type: doc
      :img-top: ../../_static/icons/user-guide-docs.svg
      :class-card: pycsamt-card sd-text-center

      How the agent layer is organised: roles, inputs and outputs, and
      where agents plug into the processing workflow.

   .. grid-item-card:: Agent catalogue
      :link: agent_catalogue
      :link-type: doc
      :img-top: ../../_static/icons/glossary-icon-book.svg
      :class-card: pycsamt-card sd-text-center

      Browse every specialised agent — foundation, processing, inversion,
      and more — with capabilities and entry points.

   .. grid-item-card:: Agent and LLM configuration
      :link: llm_configuration
      :link-type: doc
      :img-top: ../../_static/icons/developer-notes-icon-tools.svg
      :class-card: pycsamt-card sd-text-center

      Point pyCSAMT at your LLM provider and tune agent behaviour,
      budgets, and safety rails.

Agent Families
--------------

.. grid:: 1 1 2 3
   :gutter: 3
   :class-container: cta-tiles

   .. grid-item-card:: Foundation and survey intake
      :link: foundation_agents
      :link-type: doc
      :img-top: ../../_static/icons/user-guide-data-loading.svg
      :class-card: pycsamt-card sd-text-center

      Loader, metadata, and survey-intake agents that bring field data
      into the canonical containers.

   .. grid-item-card:: Processing and diagnostics
      :link: processing_agents
      :link-type: doc
      :img-top: ../../_static/icons/user-guide-processing.svg
      :class-card: pycsamt-card sd-text-center

      QC, denoising, static-shift, frequency, and phase-analysis agents
      for data conditioning.

   .. grid-item-card:: Forward and inversion
      :link: inversion_agents
      :link-type: doc
      :img-top: ../../_static/icons/user-guide-inversion.svg
      :class-card: pycsamt-card sd-text-center

      Forward modelling, inversion preparation, 2-D/3-D engines, and
      evaluation agents.

   .. grid-item-card:: AI and model zoo
      :link: ai_model_zoo_agents
      :link-type: doc
      :img-top: ../../_static/icons/user-guide-ai.svg
      :class-card: pycsamt-card sd-text-center

      Deep-learning inversion, ensembles, anomaly detection, and the
      pretrained model-zoo agents.

   .. grid-item-card:: Orchestration and output
      :link: orchestration_output_agents
      :link-type: doc
      :img-top: ../../_static/icons/user-guide-pipeline.svg
      :class-card: pycsamt-card sd-text-center

      Pipeline, reporting, export, and map/plot output agents that close
      out a workflow.

   .. grid-item-card:: Agent coordinator
      :link: coordinator
      :link-type: doc
      :img-top: ../../_static/icons/applications.svg
      :class-card: pycsamt-card sd-text-center

      The low-level coordinator that sequences agents, shares state, and
      arbitrates between plans.

Section Contents
----------------

.. toctree::
   :maxdepth: 1
   :hidden:

   overview
   agent_catalogue
   foundation_agents
   processing_agents
   inversion_agents
   ai_model_zoo_agents
   orchestration_output_agents
   coordinator
   orchestrator
   llm_configuration
