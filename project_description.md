# pyCSAMT (Python  for Controlled Source Audio-frequency Magnetotellurics )
[![Documentation Status](https://readthedocs.org/projects/pycsamt/badge/?version=latest)](https://pycsamt.readthedocs.io/en/latest/?badge=latest) [![Build Status](https://travis-ci.com/WEgeophysics/pyCSAMT.svg?branch=master)](https://travis-ci.com/WEgeophysics/pyCSAMT)
     
An open source far field controlled source audio-frequency magnetotellurics 
for standard data processing , modeling and geophysical interpretation enhancement.  

## Overview 


* [Definition](#https://ui.adsabs.harvard.edu/abs/2018EGUGA..2013744L/abstract)

Controlled Source Audio-frequency Magnetotellurics or CSAMT is an enhanced frequency 
domain EM program using synchronous stacking and averaging, and Fourier integral methods to 
improve the signal to noise ratio. Later, CSAMT is well appreciated in geophysical commumity and
consequently was used as a suitable exploration method well-established in deep geological structure detection.
Today the method is broadly applied in  diverse of exploration problems such as mineral , hydrocarbon,  groundwater resources, 
as well as mapping the fault-zones etc. 

* [Purpose](#Purpose)

pyCSAMT contains bacics steps and improve CSAMT standard data processing as well as the modeling using [OCCAM2D](https://marineemlab.ucsd.edu/Projects/Occam/index.html).
The software deals with more than 150  geological structures and electrical properties of rocks to generate a pseudo-stratigraphy 2D map to enhance geophysical interpretation especially in more geological complex area ( with various tectonic accidents). 


 * [Note](#https://iopscience.iop.org/article/10.1088/1742-6596/1127/1/012021)
 
pyCSAMT, actually only works  in [far field](https://electronics.stackexchange.com/questions/487691/why-are-e-and-b-field-in-phase-in-far-field-electromagnetic-wave-propagation). Furhermore , it uses [OCCAM2D](https://marineemlab.ucsd.edu/Projects/Occam/index.html) open source sofware as modeling software. Nevertheless,
several  outputs are provided for other external modeling softwares like [MTpy](https://github.com/MTgeophysics/mtpy), [OasisMontaj](http://updates.geosoft.com/downloads/files/how-to-guides/Oasis_montaj_Gridding.pdf)
and [GoldenSoftware](https://www.goldensoftware.com/products/surfer).

## Installation 

This theme is distributed on [PyPI](https://pypi.org/project/pycsamt/) and can be installed with `pip`:
```
$ pip install pycsamt

or 

$ pip install user pycsamt

```
For more information read the full documentation on [installing pycsamt](#pycsamt) 

## Contributing

If you would like to help modify or enhance the project, you’ll find more information on [issue template](https://github.com/WEgeophysics/pyCSAMT/blob/master/ISSUE_TEMPLATE.md) file located in [git repository](https://github.com/WEgeophysics/pycsamt).