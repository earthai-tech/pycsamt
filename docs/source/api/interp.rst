pycsamt.interp
==============

Hydrogeological interpretation helpers, borehole calibration, uncertainty,
fusion, and report/export routines, built on top of the geological domain
knowledge in :doc:`geology`.

.. automodule:: pycsamt.interp
   :members:
   :show-inheritance:

.. note::

   ``RockDatabase``, ``RockEntry``, ``Layer``, ``StratigraphicLog``,
   ``Borehole``, and ``Interval`` remain importable from
   ``pycsamt.interp`` for backward compatibility, but their canonical
   home — and their generated reference pages — is now :doc:`geology`.

Interpretation Modules
----------------------

.. autosummary::
   :toctree: generated

   pycsamt.interp.calibrate
   pycsamt.interp.constraints
   pycsamt.interp.export
   pycsamt.interp.fusion
   pycsamt.interp.hydro
   pycsamt.interp.hydromodel
   pycsamt.interp.petrophysics
   pycsamt.interp.plot
   pycsamt.interp.timelapse
   pycsamt.interp.uncertainty
