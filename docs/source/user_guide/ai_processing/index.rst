.. _user_guide_ai_processing:

AI-Assisted Processing
=======================

:mod:`pycsamt.ai.processing` applies trained models to eight
data-preparation tasks. Four already have a deterministic, rule-based
counterpart in :doc:`../emtools/index` -- quality control, noise
suppression, anomaly detection, and dimensionality classification --
and each tool here is a *complement*, not a replacement, for it. Gap
filling and uncertainty calibration have no rule-based equivalent at
all: nothing else in pyCSAMT reconstructs a genuinely missing
impedance cell or re-estimates an inversion error floor from signal
quality. Distortion triage is a third kind: it routes each station
toward whichever *existing* correction tool
(:mod:`pycsamt.emtools.ss` or :mod:`pycsamt.emtools.gb`) its regime
actually calls for, rather than complementing either directly.
:doc:`ts_denoise` is the odd one out in a different way: every other
tool here works on the frequency-domain impedance tensor, while it
works one stage further upstream, on the raw field time series before
it is even Fourier-transformed. Start with :doc:`overview` for how the
tools relate to each other, the shared estimator interface, and when
to reach for a learned model over a rule-based one.

.. toctree::
   :maxdepth: 3
   :class: pycsamt-guide-toc

   overview
   qc
   denoise
   ts_denoise
   anomaly
   classify
   imputer
   uncertainty
   distortion
