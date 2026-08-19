pycsamt.emtools
===============

Electromagnetic processing, diagnostics, tensor analysis, quality control,
static-shift correction, source-effect tools, and plotting helpers.

.. seealso::

   :doc:`../user_guide/emtools/index` for narrative, runnable examples built module
   by module (currently: :doc:`../user_guide/emtools/tf`).

.. automodule:: pycsamt.emtools
   :members:
   :show-inheritance:
   :exclude-members: RoseStyle, overlay_noise_cone, overlay_spectral_holes,
                      plot_consistency_fan, plot_coverage_psection,
                      plot_qc_quicklook, plot_snr_hist,
                      plot_xyyx_crossover_map

Processing and QC
-----------------

.. autosummary::
   :toctree: generated

   pycsamt.emtools.afmag
   pycsamt.emtools.inspect
   pycsamt.emtools.qc
   pycsamt.emtools.frequency
   pycsamt.emtools.ss
   pycsamt.emtools.remove_noise
   pycsamt.emtools.diag
   pycsamt.emtools.spectra
   pycsamt.emtools.plot
   pycsamt.emtools.overview

Tensor and Dimensionality
-------------------------

.. autosummary::
   :toctree: generated

   pycsamt.emtools.impedance
   pycsamt.emtools.tensor
   pycsamt.emtools.tf
   pycsamt.emtools.gb
   pycsamt.emtools.skew
   pycsamt.emtools.strike
   pycsamt.emtools.dimensionality
   pycsamt.emtools.anisotropy
   pycsamt.emtools.gradient_imaging

Source and Field-Zone Effects
-----------------------------

.. autosummary::
   :toctree: generated

   pycsamt.emtools.fieldzone
   pycsamt.emtools.source_array
   pycsamt.emtools.source_effects
   pycsamt.emtools.csumt
   pycsamt.emtools.lcurve
   pycsamt.emtools.legacy

Rose Diagram Styling
---------------------

.. autoclass:: pycsamt.emtools.RoseStyle
   :members:

Public QC Plot Functions
------------------------

.. autofunction:: pycsamt.emtools.overlay_noise_cone
.. autofunction:: pycsamt.emtools.overlay_spectral_holes
.. autofunction:: pycsamt.emtools.plot_consistency_fan
.. autofunction:: pycsamt.emtools.plot_coverage_psection
.. autofunction:: pycsamt.emtools.plot_qc_quicklook
.. autofunction:: pycsamt.emtools.plot_snr_hist
.. autofunction:: pycsamt.emtools.plot_xyyx_crossover_map
