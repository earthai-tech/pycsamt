pycsamt.airborne
================

Format-neutral airborne electromagnetic survey containers: flight
lines and datasets built from ``EMTF`` documents, a technology/format
registry, native-I/O extension points, structural QC, and the
``AirborneSite``/``AirborneSites`` ``Sites``-shaped read path. AFMAG,
ZTEM, and MobileMT map decoded scientific arrays onto this common
model rather than each inventing its own containers.

.. seealso::

   :doc:`../user_guide/airborne/index` for narrative, runnable
   examples built page by page (data model, site view, technology/
   format registry, and structural QC).

.. automodule:: pycsamt.airborne
   :members:
   :show-inheritance:

Core Containers
---------------

.. autosummary::
   :toctree: generated

   pycsamt.airborne.base
   pycsamt.airborne.navigation

Site View
---------

.. autosummary::
   :toctree: generated

   pycsamt.airborne.site

Technology and Format Registry
------------------------------

.. autosummary::
   :toctree: generated

   pycsamt.airborne.registry
   pycsamt.airborne.io

Quality Control
---------------

.. autosummary::
   :toctree: generated

   pycsamt.airborne.qc

Shared Validation Helpers
-------------------------

.. autosummary::
   :toctree: generated

   pycsamt.airborne.validation

AFMAG Adapter
-------------

.. autosummary::
   :toctree: generated

   pycsamt.airborne.afmag
   pycsamt.airborne.afmag.adapter
   pycsamt.airborne.afmag.base
   pycsamt.airborne.afmag.datatypes
   pycsamt.airborne.afmag.constants

ZTEM Adapter
------------

.. autosummary::
   :toctree: generated

   pycsamt.airborne.ztem
   pycsamt.airborne.ztem.adapter
   pycsamt.airborne.ztem.base
   pycsamt.airborne.ztem.constants

MobileMT Adapter
----------------

.. autosummary::
   :toctree: generated

   pycsamt.airborne.mobilemt
   pycsamt.airborne.mobilemt.adapter
   pycsamt.airborne.mobilemt.base
   pycsamt.airborne.mobilemt.datatypes
   pycsamt.airborne.mobilemt.constants
