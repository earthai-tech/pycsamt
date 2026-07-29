pycsamt.ai
==========

AI inversion, neural network architectures, training helpers, AI processing,
model-zoo utilities, and diagnostic plotting.

.. automodule:: pycsamt.ai
   :members:
   :show-inheritance:

AI Modules
----------

The newer data contracts, geological priors, objective functions, domain-gap
tools, experiment records, and validation gates have dedicated references:

.. toctree::
   :maxdepth: 1

   ai_data
   ai_geology
   ai_losses
   ai_domain_gap
   ai_experiments
   ai_validation

.. autosummary::
   :toctree: generated

   pycsamt.ai.inversion.inv1d
   pycsamt.ai.inversion.inv2d
   pycsamt.ai.inversion.inv3d
   pycsamt.ai.inversion.ensemble
   pycsamt.ai.inversion.joint
   pycsamt.ai.inversion.calibration
   pycsamt.ai.nets.cnn1d
   pycsamt.ai.nets.drcnn
   pycsamt.ai.nets.fcn
   pycsamt.ai.nets.gcn
   pycsamt.ai.nets.resnet
   pycsamt.ai.nets.unet
   pycsamt.ai.processing.anomaly
   pycsamt.ai.processing.classify
   pycsamt.ai.processing.denoise
   pycsamt.ai.processing.qc
   pycsamt.ai.training.augment
   pycsamt.ai.training.dataset
   pycsamt.ai.training.metrics
   pycsamt.ai.training.trainer
   pycsamt.ai.plot.compare
   pycsamt.ai.plot.convergence
   pycsamt.ai.plot.diagnostics
   pycsamt.ai.plot.section
