.. _contributing:

Contributing
============

pyCSAMT v2 is developed on `GitHub <https://github.com/earthai-tech/pycsamt>`_.

Development setup
-----------------

.. code-block:: bash

   git clone https://github.com/earthai-tech/pycsamt.git
   cd pycsamt
   pip install -e ".[full,docs,dev]"

Running the tests
-----------------

.. code-block:: bash

   pytest pycsamt/ -v --tb=short

Building the docs
-----------------

.. code-block:: bash

   cd docs
   pip install -r requirements-docs.txt
   make html
   # open build/html/index.html

Commit style
------------

* Descriptive imperative subject line (≤ 72 chars).
* Co-author line: ``Co-Authored-By: earthai-tech <earthai-tech@users.noreply.github.com>``.

Reporting issues
----------------

Use the `GitHub issue tracker <https://github.com/earthai-tech/pycsamt/issues>`_.
