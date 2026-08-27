pycsamt.format
==============

Backend-neutral PCSF and PCSM representations for electromagnetic inversion
results. The package provides the shared model schema, binary and text I/O,
conversion adapters, multiline construction, topography integration,
point-cloud extraction, regridding, and model provenance helpers.

For format concepts, complete workflows, and browsable file examples, see
:doc:`../user_guide/models/pcsf_format`.

Public facade
-------------

The most commonly used classes and functions are re-exported directly from
:mod:`pycsamt.format`.

.. automodule:: pycsamt.format
   :members:
   :show-inheritance:

Schema and serialization
------------------------

.. autosummary::
   :toctree: generated

   pycsamt.format.schema
   pycsamt.format.io
   pycsamt.format.text

Backend adapters
----------------

The adapter namespace contains converters for supported inversion codes and
generic constructors for classical, machine-learning, and deep-learning
results that already expose their geometry and resistivity arrays.

.. autosummary::
   :toctree: generated

   pycsamt.format.adapters
   pycsamt.format.adapters.generic
   pycsamt.format.adapters.occam2d
   pycsamt.format.adapters.modem3d
   pycsamt.format.adapters.mare2dem

Geometry and visualization support
----------------------------------

.. autosummary::
   :toctree: generated

   pycsamt.format.multiline
   pycsamt.format.pointcloud
   pycsamt.format.regrid
   pycsamt.format.topography
   pycsamt.format.topo_source

Provenance
----------

.. autosummary::
   :toctree: generated

   pycsamt.format.provenance
