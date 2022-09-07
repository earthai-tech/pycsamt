Installation
============

Installing via pip 
------------------
To install the `pycsamt` package from the Python package index into an established
Python environment, use the pip command:

.. code-block:: bash
   
   pip install pycsamt 
   pip install --user pycsamt (for window users) 
   
However, it is recommended the installation from the repository to get the latest development code.

Intalling from source 
----------------------
To install from source, clone the project with git: 

.. code-block:: bash 

   git clone https://github.com/WEgeophysics/pycsamt.git 
  
Or download the latest version from the project webpage: https://github.com/WEgeophysics/pycsamt

In the source directory use the command

.. code-block:: bash

   python setup.py install
   
   
Required Dependencies
---------------------
`pycsamt` was originally built on python 3.6. However, the last version requires at least **python 3.8**.

It calls on the core Python data analytics stack, and a third party parsing library:

* Numpy
* Scipy
* Pandas
* MTpy
* xlrd
* regex
* tqdm
* pytest
* flake8
* pyyaml
* qtpy
* netcdf4

These modules should build automatically if you are installing via `pip`. If you are building from
the source code, or if pip fails to load them, they can be loaded with the same `pip` syntax as
above.   

Optional Dependencies
---------------------
In order to plot results from the model as shown in :doc:`quick start <../quick_start>`:

* Matplotlib
* Numba 

Additional Resources
--------------------
The `pycsamt user guide <https://github.com/WEgeophysics/pyCSAMT/blob/master/docs/user_guide.pdf>`_ contains recipes that can help you get start up with pycsamt.
Moreover, to go through step by step installation, one may refer to `Step by Step Installation Guide <https://github.com/WEgeophysics/pyCSAMT/wiki/pyCSAMT-installation-guide-for-Windows--and-Linux>`_.



