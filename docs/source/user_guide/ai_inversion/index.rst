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

.. toctree::
   :maxdepth: 3
   :class: pycsamt-guide-toc

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
