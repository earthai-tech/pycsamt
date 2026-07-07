.. _emtools_gb:

pycsamt.emtools.gb — Groom-Bailey Galvanic Distortion
======================================================

:mod:`pycsamt.emtools.gb` estimates and removes frequency-independent
galvanic distortion from MT/AMT/CSAMT impedance tensors using a
Groom-Bailey-style model.  For each station, pyCSAMT fits a real 2x2
distortion matrix :math:`D` over a selected period band:

.. math::

   Z_\mathrm{obs}(f) \approx D\,Z_{2D}(f),

where :math:`Z_{2D}` is anti-diagonal at each frequency.  The fitted
matrix is decomposed into gain, twist, shear, and anisotropy-style
parameters, and can then be removed as :math:`Z_0 = D^{-1} Z_\mathrm{obs}`.

Functions
---------

- :func:`~pycsamt.emtools.gb.groom_bailey_table`
- :func:`~pycsamt.emtools.gb.apply_groom_bailey`
- :func:`~pycsamt.emtools.gb.groom_bailey_decomposition`
- :class:`~pycsamt.emtools.gb.GroomBaileyResult`

Usage
-----

Estimate distortion parameters without changing the data:

.. code-block:: python

   from pycsamt.emtools import groom_bailey_table

   gb = groom_bailey_table(
       sites,
       band=(1e-3, 10.0),
       robust=True,
   )
   print(gb[[
       "station", "twist_deg", "shear", "anisotropy",
       "rms_fit", "diagonal_ratio_before", "diagonal_ratio_after",
   ]])

Estimate and apply the correction:

.. code-block:: python

   from pycsamt.emtools import groom_bailey_decomposition

   result = groom_bailey_decomposition(
       sites,
       band=(1e-3, 10.0),
       apply=True,
   )
   corrected_sites = result.sites
   gb_table = result.table

Notes
-----

The fitted distortion matrix is real and frequency-independent, which is
the galvanic assumption behind the Groom-Bailey family of decompositions.
The implementation uses robust alternating least squares and is intended
as a reproducible preprocessing and audit step before 2-D inversion.  It
does not resolve the scalar static-shift gain ambiguity uniquely; the
reported gain is therefore a fitted matrix scale, while interpretation
should focus mainly on twist, shear, anisotropy, residual fit, and the
effect of correction on diagonal leakage.
