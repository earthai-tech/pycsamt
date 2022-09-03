# pycsamt: An open-source Python toolkit for audio-frequency magnetotellurics

[![Documentation Status](https://readthedocs.org/projects/pycsamt/badge/?version=latest)](https://pycsamt.readthedocs.io/en/latest/?badge=latest) [![Build Status](https://travis-ci.com/WEgeophysics/pyCSAMT.svg?branch=master)](https://travis-ci.com/WEgeophysics/pyCSAMT) [![Requirements Status](https://requires.io/github/WEgeophysics/pyCSAMT/requirements.svg?branch=master)](https://requires.io/github/WEgeophysics/pyCSAMT/requirements/?branch=master)
  ![GitHub](https://img.shields.io/github/license/WEgeophysics/pycsamt?color=blue&label=licence&logo=GNU&logoColor=red) ![GitHub release (latest by date)](https://img.shields.io/github/v/release/WEgeophysics/pyCSAMT?color=orange) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5674430.svg)](https://doi.org/10.5281/zenodo.5674430)
  [![PyPI version](https://badge.fury.io/py/pycsamt.svg)](https://badge.fury.io/py/pycsamt)
  ![GitHub repo size](https://img.shields.io/github/repo-size/WEgeophysics/pycsamt?color=0A4CEE&style=flat-square)
  
## Overview 

    * **Purpose**
    
    Originally, the software was intended for controlled source audio-frequency magnetotelluric (CSAMT) data processing (hence the suffix CSAMT) and was mostly related
    to the groundwater exploration. However, the developement has been advanced to encompass the audio-magnetotelluric(AMT) methods i.e the frequency above 1Hz ( such as the CSAMT, the natural 
    source AMT (NSAMT)). The software is used to broadly applied in diverse of exploration problems such as mineral, hydrocarbon,  groundwater resources, as well as the fault-zone mapping . 
    It contains steps of AMT data processing and deals with [OCCAM2D](https://marineemlab.ucsd.edu/Projects/Occam/index.html) of [DeGroot-Hedlin, C. and Constable, S., 1990](https://doi.org/10.1190/1.1442303) , 
    the [MT2DInvMatlab](https://doi.org/10.1016/j.cageo.2008.10.010)  of [Lee et al., 2009](https://doi.org/10.1016/j.cageo.2008.10.010) and [ModEM](https://sites.google.com/site/modularem/download) of [Kelbert et al., 2014](https://doi.org/10.1016/j.cageo.2014.01.010)
    for the modeling part.
    
    It also provides processing tools for filtering data and restoring amplitudes at weak frequency and the "dead-dand" tensor or scalar values especially 
    for the NSAMT where the natural signals are not under our control and suffer from frequency ranges with little or no signal. Some implemented filters are  the trimming moving-average (MA)(TMA),
    the fixed-length MA (FMA) (both proposed by [Zonge International Engineering (Zonge, 2000)]( http://www.zonge.com/legacy/PDF_DatPro/Astatic.pdf )) and the adaptative MA based on idea of [Torres-verdìn and Bostick, 1992](https://doi.org/10.1190/1.2400625) usefull for fast
    removing the static effect especially in electromagnetic-array profiling survey. Some others filters such as "simple_filter" for outliers removal and "PCA_removal" filters can also be applied 
    upstream for a particular data where the interferences are very strong (e.g. intenses humman activities , power lines, ...). Moreover, the  "Savitzky-Golay" filter is also added to remove high-frequency noise from data since it 
     has the advantage of preserving the original shape and features of the signal better than other types of filtering approaches such as MA techniques (simple, exponential, cumulative, weight) . 
     
    Besides the AMT processing, the software implements a special sub-package (geology + drilling ) referred  as `geodrill`, entirely dedicated to improve the groundwater exploration techniques. The aim is to reduce the numerous unsucessful 
    drillings mostly occured due to their wrong locations after AMT geophysical surveys. It consists to minimize the use of supplement methods to AMT which commonly increases the operating budgets 
    to right locate  the drilling thereby reducing the misinterpretation of modeling results(e.g., demarcating well the fracture zones). More details can be found in the **Citations** section. 
    
    
 * **Note**
 
    For long periods or MT methods (below 1Hz), it is recommended to visit other suitable softwares such as  [MTpy](https://github.com/MTgeophysics/mtpy.git), [FEMT2D](https://github.com/ruboerner/FEMT2D) or [razorback](https://github.com/BRGM/razorback)   
    or visit the [MTnet](https://www.mtnet.info/main/source.html) website. Nonetheless, the sofware has a feature to generate outputs/objects for other external modeling softwares such as [MTpy](https://github.com/MTgeophysics/mtpy), [OasisMontaj](http://updates.geosoft.com/downloads/files/how-to-guides/Oasis_montaj_Gridding.pdf)
    and [GoldenSoftware](https://www.goldensoftware.com/products/surfer).


## Documentation 

* [Installation Guide](https://pycsamt.readthedocs.io/en/latest/installation.html?highlight=installation)
* [Demo](https://pycsamt.readthedocs.io/en/latest/demo.html?highlight=demo) 
* [API Documentation](https://pycsamt.readthedocs.io/en/latest/)
* [Home Page](https://github.com/WEgeophysics/pyCSAMT/wiki)
* [User Guide](https://github.com/WEgeophysics/pyCSAMT/blob/develop/docs/pyCSAMT%20User%20Guide.pdf)

## Credits

We use or link some third-party software (beside the usual tool stack: numba, numpy, scipy, matplotlib) and are grateful for all the work made by the authors of these awesome open-source tools:
* MTpy: https://github.com/MTgeophysics/mtpy.git
* Occam2D: https://marineemlab.ucsd.edu/Projects/Occam/index.html
* ModEM: https://sites.google.com/site/modularem/
* Zonge Engineering softwares:
    - AMTAVG: http://www.zonge.com/legacy/DatPro.html/
    - ASTATIC: http://www.zonge.com/legacy/PDF_DatPro/Astatic.pdf

## System requirements 

* Python 3.7+ 

## Citations 

If you use pycsamt in any published work, you may consider as a contribution, citing [zenodo DOI](https://zenodo.org/record/5674430/export/hx#.YmpoXo9BxGo) or the paper below:

> *Kouadio, K.L., Liu, R., Mi, B., Liu, C., 2022. pyCSAMT: An alternative Python toolbox for groundwater exploration using controlled source audio-frequency magnetotelluric. J. Appl. Geophys. 201, 104647. https://doi.org/10.1016/j.jappgeo.2022.104647.*
> 
> *Kouadio, K.L., 2021. pyCSAMT: A Python open-source toolkit for controlled source audio-frequency magnetotelluric. https://doi.org/10.5281/zenodo.5674430.*


## Contributor

1. Department of Geophysics, School of Geosciences and Info-physics, [Central South University](http://www.zju.edu.cn/english/), China.  
2. Key Laboratory of Geosciences Big Data and Deep Resource of Zhejiang Province, School of Earth Sciences, [Zhejiang University](http://www.zju.edu.cn/english/), China.
3. Equipe de Recherche Géophysique Appliquée, Laboratoire de Geologie Ressources Minerales et Energetiques, UFR des Sciences de la Terre et des Ressources Minières, [Université Félix Houphouët-Boigny]( https://www.univ-fhb.edu.ci/index.php/ufr-strm/), Cote d'Ivoire.

* Developer: 1, 2, 3 - Kouadio Laurent a.k.a Daniel; <etanoyau@gmail.com>, <kkouao@zju.edu.cn>,



