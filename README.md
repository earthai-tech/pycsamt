# _pycsamt_: A package for audio-frequency magnetotellurics

[![Documentation Status](https://readthedocs.org/projects/pycsamt/badge/?version=latest)](https://pycsamt.readthedocs.io/en/latest/?badge=latest) [![Build Status](https://travis-ci.com/WEgeophysics/pyCSAMT.svg?branch=master)](https://travis-ci.com/WEgeophysics/pyCSAMT) 
  ![GitHub](https://img.shields.io/github/license/WEgeophysics/pycsamt?color=blue&label=licence&logo=GNU&logoColor=red) ![GitHub release (latest by date)](https://img.shields.io/github/v/release/WEgeophysics/pyCSAMT?color=orange) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5674430.svg)](https://doi.org/10.5281/zenodo.5674430)
  [![PyPI version](https://badge.fury.io/py/pycsamt.svg)](https://badge.fury.io/py/pycsamt)

## Overview 

  * **Purpose**
    
    Previously, the software was intended for controlled source audio-frequency magnetotelluric (CSAMT) data processing (hereinafter the suffix CSAMT) and mostly related
    to the groundwater exploration. The recent development is focused on the audio-magnetotelluric(AMT) methods. Indeed, the AMT methods are used broadly in diverse of exploration problems such as mineral, hydrocarbon,  groundwater resources, as well as the fault-zone mapping above the 1km depth. 
    _pycsamt_ is designed to bring a piece of solution to the problems encountered by using AMT methods. It contains steps of AMT data processing and deals with [OCCAM2D](https://marineemlab.ucsd.edu/Projects/Occam/index.html) of [DeGroot-Hedlin and Constable, 1990](https://doi.org/10.1190/1.1442303) , 
    the [MT2DInvMatlab](https://doi.org/10.1016/j.cageo.2008.10.010)  of [Lee et al., 2009](https://doi.org/10.1016/j.cageo.2008.10.010) and [ModEM](https://sites.google.com/site/modularem/download) of [Kelbert et al., 2014](https://doi.org/10.1016/j.cageo.2014.01.010)
    for the modeling purpose.
    
    It also provides processing tools for filtering and processing data( the trimming moving-average (MA), the fixed dipole-length MA ([Zonge International Engineering (Zonge, 2000)]( http://www.zonge.com/legacy/PDF_DatPro/Astatic.pdf )), the adaptative MA ( [Torres-verdìn and Bostick, 1992](https://doi.org/10.1190/1.2400625)). These filters are mostly used for fast removing the static effect especially in electromagnetic-array profiling survey. Some others filters such as "simple" for outliers removal and "PCA" can also be applied upstream for a particular data where the interferences are very strong (e.g. intenses humman activities, power lines, ...). Moreover, the  "Savitzky-Golay" filter is also added to remove high-frequency noise from data since it has the advantage of preserving the original shape and features of the signal better than other types of filtering approaches such as MA techniques (simple, exponential, cumulative, weight). 
    
 * **Note**
 
    For long periods or MT methods (below 1Hz), it is recommended to visit other suitable softwares such as  [MTpy](https://github.com/MTgeophysics/mtpy.git), [FEMT2D](https://github.com/ruboerner/FEMT2D), [razorback](https://github.com/BRGM/razorback) or consult the [MTNet](https://www.mtnet.info/main/source.html) website. Nonetheless, the sofware has a feature to generate outputs/objects for other external modeling softwares such as [MTpy](https://github.com/MTgeophysics/mtpy), [OasisMontaj](http://updates.geosoft.com/downloads/files/how-to-guides/Oasis_montaj_Gridding.pdf) and [GoldenSoftware](https://www.goldensoftware.com/products/surfer).


## Documentation 

* [Installation Guide](https://pycsamt.readthedocs.io/en/latest/installation.html?highlight=installation)
* [API Documentation](https://pycsamt.readthedocs.io/en/latest/)
* [Home Page](https://github.com/WEgeophysics/pyCSAMT/wiki)
* [User Guide](https://github.com/WEgeophysics/pyCSAMT/blob/develop/docs/pyCSAMT%20User%20Guide.pdf)

## Credits

We use or link some third-party software (beside the usual tool stack: [Numba](https://numba.pydata.org/), [Numpy](https://numpy.org/), [Scipy](https://scipy.org/), [SumPy](https://www.sympy.org/en/index.html), [Matplotlib](https://matplotlib.org/)) and are grateful for all the work made by the authors of these awesome open-source tools:
* [MTpy](https://github.com/MTgeophysics/mtpy.git)
* [Occam2D](https://marineemlab.ucsd.edu/Projects/Occam/index.html)
* [ModEM](https://sites.google.com/site/modularem/)
* Zonge Engineering softwares:
    - [AMTAVG](http://www.zonge.com/legacy/DatPro.html/)
    - [ASTATIC](http://www.zonge.com/legacy/PDF_DatPro/Astatic.pdf)

## System requirements 

* Python 3.7+ 

## Citations 

 We'd much appreciate if you consider citing the toolbox as a contribution in a published work:

> *Kouadio, K.L., Liu, R., Mi, B., Liu, C., 2022. pyCSAMT: An alternative Python toolbox for groundwater exploration using controlled source audio-frequency magnetotelluric. J. Appl. Geophys. 201, 104647. https://doi.org/10.1016/j.jappgeo.2022.104647.*
> 

A case study using the software can also be found  in the paper below: 

> *Kouadio, L. K., Liu, R., Malory, A. O., Liu, W., Liu, C., A novel approach for water reservoir mapping using controlled source audio - frequency magnetotelluric in Xingning area , Hunan Province, China. Geophys. Prospect., https://doi.org/10.1111/1365-2478.13385*
>

## Contributors

1. Department of Geophysics, School of  Info-physics and Geomatics Engineering, [Central South University](https://en.csu.edu.cn/), China. 
2. Equipe de Recherche Géophysique Appliquée, Laboratoire de Geologie Ressources Minerales et Energetiques, UFR des Sciences de la Terre et des Ressources Minières, [Université Félix Houphouët-Boigny]( https://www.univ-fhb.edu.ci/index.php/ufr-strm/), Cote d'Ivoire.

* Developer: 1, 2- Kouadio Laurent,  <etanoyau@gmail.com> / <lkouao@csu.edu.cn>.
