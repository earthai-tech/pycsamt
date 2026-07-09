.. _modeling:

Modeling
========

Everything model-side in one place: classical inversion engines
(Occam2D, ModEM, MARE2DEM), synthetic forward modelling for survey
design and intuition, and the deep-learning inverters. Start with
forward modelling when you are designing a survey or building
intuition, the model integrations when you are preparing a real
inversion run, and AI inversion when you want a fast learned surrogate.

.. grid:: 1 1 2 3
   :gutter: 3
   :class-container: cta-tiles

   .. grid-item-card:: Model integrations
      :link: models/index
      :link-type: doc
      :img-top: ../_static/icons/user-guide-models.svg
      :class-card: pycsamt-card sd-text-center

      Prepare and inspect Occam2D, ModEM, MARE2DEM, and related external
      modelling or inversion engine workflows.

   .. grid-item-card:: Forward modelling
      :link: forward/index
      :link-type: doc
      :img-top: ../_static/icons/user-guide-forward.svg
      :class-card: pycsamt-card sd-text-center

      Build controlled synthetic responses, solver configurations, model
      grids, noise models, and forward-to-inversion experiments.

   .. grid-item-card:: AI inversion
      :link: ai_inversion
      :link-type: doc
      :img-top: ../_static/icons/user-guide-ai.svg
      :class-card: pycsamt-card sd-text-center

      Use AI-assisted inversion tools, training-data concepts, prediction
      outputs, and safeguards for model-assisted interpretation.

.. toctree::
   :maxdepth: 1
   :hidden:

   models/index
   forward/index
   ai_inversion
