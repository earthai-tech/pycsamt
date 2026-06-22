.. _installation:

Installation
============

Requirements
------------

* Python **3.9** or later
* NumPy ≥ 1.22, SciPy ≥ 1.8, Matplotlib ≥ 3.5

Install from PyPI
-----------------

.. code-block:: bash

   pip install pycsamt

With AI/ML back-ends
~~~~~~~~~~~~~~~~~~~~

Install with `PyTorch <https://pytorch.org>`_ (recommended):

.. code-block:: bash

   pip install pycsamt[torch]

Install with `TensorFlow <https://tensorflow.org>`_:

.. code-block:: bash

   pip install pycsamt[tensorflow]

Install everything (both back-ends + optional scientific extras):

.. code-block:: bash

   pip install pycsamt[full]

Install from source
-------------------

.. code-block:: bash

   git clone https://github.com/earthai-tech/pycsamt.git
   cd pycsamt
   pip install -e ".[full]"

.. note::

   The inversion solvers (Occam2D, ModEM) are distributed as Fortran
   source code under ``pycsamt/models/*/``_source/``.  They are compiled
   on first use via ``f2py``.  A Fortran compiler (``gfortran``) must be
   available on the system.

Verify installation
-------------------

.. code-block:: python

   import pycsamt
   print(pycsamt.__version__)
