pycsamt.format
==============

Backend-neutral PCSF and PCSM representations for electromagnetic inversion
results, plus the human-readable PCBH borehole exchange contract. The package
provides shared schemas, serialization, conversion adapters, model/borehole
association, rendering contracts, and provenance helpers.

For format concepts, complete workflows, and browsable file examples, see
:doc:`../user_guide/models/pcsf_format`. For spatial multi-hole projects, see
:doc:`../user_guide/geology/pcbh_format`.

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

Points and interpretation exchange
----------------------------------

Small, versioned JSON documents (each with a plain-CSV sibling) that
travel alongside a PCSF model: PCPT points/targets, the PCGL
resistivity-to-geology legend, and PCGS field structural evidence.

.. autosummary::
   :toctree: generated

   pycsamt.format.pointset
   pycsamt.format.geology
   pycsamt.format.structure

Borehole exchange
-----------------

The :mod:`pycsamt.format.borehole` namespace contains PCBH schema objects,
canonical JSON I/O, CSV/LAS adapters, trajectory derivation, PCSF association,
and application-neutral 3-D render/export contracts.

.. autosummary::
   :toctree: generated

   pycsamt.format.borehole
   pycsamt.format.borehole.schema
   pycsamt.format.borehole.jsonio
   pycsamt.format.borehole.csvio
   pycsamt.format.borehole.xlsxio
   pycsamt.format.borehole.lasio
   pycsamt.format.borehole.relational
   pycsamt.format.borehole.trajectory
   pycsamt.format.borehole.pcsf
   pycsamt.format.borehole.render
   pycsamt.format.borehole.exports
