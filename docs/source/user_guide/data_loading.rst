.. _user_guide_data:

Data loading
============

.. note:: Placeholder — content coming in v2 stable.

pyCSAMT v2 accepts EDI files as the primary data format.

.. code-block:: python

   from pycsamt.io import EDICollection

   coll = EDICollection.from_dir("data/edi/")
   print(coll)
