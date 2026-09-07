pycsamt.geology
===============

General-purpose geological domain knowledge with no electromagnetic
dependency: resistivity-to-lithology classification, the literature-compiled
rock/fluid property table behind it, pluggable rock-property sources,
ground-truth borehole logs, field structural measurements (strike/dip,
trend/plunge, fault traces), and geologic swatch/pattern packs.
:doc:`interp` builds on this package (see
:class:`pycsamt.interp.ModelCalibrator`) rather than the other way round.

.. automodule:: pycsamt.geology
   :members:
   :show-inheritance:

Geology Modules
---------------

.. autosummary::
   :toctree: generated

   pycsamt.geology.borehole
   pycsamt.geology.lithology
   pycsamt.geology.rock_library
   pycsamt.geology.rock_providers
   pycsamt.geology.structural
   pycsamt.geology.patterns
