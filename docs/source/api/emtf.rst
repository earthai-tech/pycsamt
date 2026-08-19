pycsamt.emtf
============

Format-neutral electromagnetic transfer-function scientific core: the
``EMTF`` document, matrix-oriented ``TransferFunction``/``StatisticalEstimate``
objects, the EMTF datatype registry, EDI/EDI-SPECTRA/EMTF-XML
interoperability, and rotation/covariance transformation.

.. automodule:: pycsamt.emtf
   :members:
   :show-inheritance:

Core Scientific Model
---------------------

.. autosummary::
   :toctree: generated

   pycsamt.emtf.document
   pycsamt.emtf.transfer
   pycsamt.emtf.estimates
   pycsamt.emtf.datatypes
   pycsamt.emtf.validation
   pycsamt.emtf.base

Rotation and Covariance
-----------------------

.. autosummary::
   :toctree: generated

   pycsamt.emtf.orientation

EDI and TFBundle Interoperability
---------------------------------

.. autosummary::
   :toctree: generated

   pycsamt.emtf.converters.edi
   pycsamt.emtf.converters.spectra
   pycsamt.emtf.converters.bundle

EMTF XML
--------

.. autosummary::
   :toctree: generated

   pycsamt.emtf.xml.reader
   pycsamt.emtf.xml.writer
   pycsamt.emtf.xml.serializer
   pycsamt.emtf.xml.parser
   pycsamt.emtf.xml.constants
