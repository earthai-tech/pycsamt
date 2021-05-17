# pyCSAMT : A Python open-source toolkit for Controlled Source Audio-frequency Magnetotellurics (CSAMT)
[![Documentation Status](https://readthedocs.org/projects/pycsamt/badge/?version=latest)](https://pycsamt.readthedocs.io/en/latest/?badge=latest) [![Build Status](https://travis-ci.com/WEgeophysics/pyCSAMT.svg?branch=master)](https://travis-ci.com/WEgeophysics/pyCSAMT)
[![Requirements Status](https://requires.io/github/WEgeophysics/pyCSAMT/requirements.svg?branch=master)](https://requires.io/github/WEgeophysics/pyCSAMT/requirements/?branch=master)
  ![GitHub](https://img.shields.io/github/license/WEgeophysics/pyCSAMT?color=blue&logo=GNU&logoColor=red) ![GitHub release (latest by date)](https://img.shields.io/github/v/release/WEgeophysics/pyCSAMT?color=orange)    

_For far-field CSAMT standard data processing, modeling and geophysical interpretation enhancement._

## Overview 


* [Definition](#https://ui.adsabs.harvard.edu/abs/2018EGUGA..2013744L/abstract)

CSAMT is an enhanced frequency domain EM program using synchronous stacking and averaging, and Fourier integral methods to 
improve the signal to noise ratio. Later, CSAMT is well appreciated in geophysical commumity and
consequently was used as a suitable exploration method well-established in deep geological structure detection.
Today the method is broadly applied in  diverse of exploration problems such as mineral , hydrocarbon,  groundwater resources, 
as well as mapping the fault-zones etc. 

* [Purpose](#Purpose)

 `pycsamt` contains basic steps and improves CSAMT standard data processing and uses  [OCCAM2D](https://marineemlab.ucsd.edu/Projects/Occam/index.html) modeling software for the fisrt release.
It deals with more than 150  geological structures and electrical properties of rocks to generate a pseudo-stratigraphy 2D map to enhance geophysical interpretation especially in more geological complex area ( with various tectonic accidents). 


 * [Note](#https://iopscience.iop.org/article/10.1088/1742-6596/1127/1/012021)
 
The software actually only works  in [far field](https://electronics.stackexchange.com/questions/487691/why-are-e-and-b-field-in-phase-in-far-field-electromagnetic-wave-propagation)
and provided several for other external modeling softwares like [MTpy](https://github.com/MTgeophysics/mtpy), [OasisMontaj](http://updates.geosoft.com/downloads/files/how-to-guides/Oasis_montaj_Gridding.pdf)
and [GoldenSoftware](https://www.goldensoftware.com/products/surfer).

## Installation 

Distributed on [PyPI](https://pypi.org/project/pycsamt/1.0.3/) and can be installed with `pip`:
`$ pip install pycsamt` or `$ pip install user pycsamt`. For more details, please refer to the [installation guide](https://github.com/WEgeophysics/pyCSAMT/wiki/pyCSAMT-installation-guide-for-Windows--and-Linux) on our [wiki page](https://github.com/WEgeophysics/pyCSAMT/wiki). 

## Contributing

If you would like to help modify or enhance the project, you'are welcome and you'll find more information on [issue template](https://github.com/WEgeophysics/pyCSAMT/blob/master/ISSUE_TEMPLATE.md) file located in [git repository](https://github.com/WEgeophysics/pycsamt).