pycsamt.api
===========

Public API configuration objects, style controls, CLI parameter helpers,
table/result view objects, and high-level reader functions.

.. automodule:: pycsamt.api
   :members:
   :show-inheritance:
   :exclude-members: SiteOrderingConfig, configure_ordering, reset_ordering, configure_topo, reset_topo

Site Ordering Configuration
---------------------------

.. autoclass:: pycsamt.api.SiteOrderingConfig
   :members:

.. autofunction:: pycsamt.api.configure_ordering
.. autofunction:: pycsamt.api.reset_ordering

The live process-wide configuration is available as
:data:`pycsamt.api.PYCSAMT_ORDERING`. See
:ref:`api-configuration` for modes, thresholds, overrides, and examples.

Topography Configuration
------------------------

.. autofunction:: pycsamt.api.configure_topo
.. autofunction:: pycsamt.api.reset_topo

API Modules
-----------

.. autosummary::
   :toctree: generated

   pycsamt.api.agents
   pycsamt.api.bunch
   pycsamt.api.cli
   pycsamt.api.cli.config
   pycsamt.api.cli.options
   pycsamt.api.cli.params
   pycsamt.api.control
   pycsamt.api.docs
   pycsamt.api.interp
   pycsamt.api.ordering
   pycsamt.api.pipe
   pycsamt.api.pipe.config
   pycsamt.api.plot
   pycsamt.api.property
   pycsamt.api.section
   pycsamt.api.station
   pycsamt.api.style
   pycsamt.api.typing
   pycsamt.api.util
   pycsamt.api.view
   pycsamt.api.view.config
   pycsamt.api.view.frame
   pycsamt.api.view.io
   pycsamt.api.view.progress
   pycsamt.api.view.result
   pycsamt.api.view.survey
   pycsamt.api.view.tables
