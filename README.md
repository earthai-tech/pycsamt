# pyCSAMT : A Python open-source toolkit for Controlled Source Audio-frequency Magnetotellurics (CSAMT)

[![Documentation Status](https://readthedocs.org/projects/pycsamt/badge/?version=latest)](https://pycsamt.readthedocs.io/en/latest/?badge=latest) [![Build Status](https://travis-ci.com/WEgeophysics/pyCSAMT.svg?branch=develop)](https://travis-ci.com/WEgeophysics/pyCSAMT) [![Requirements Status](https://requires.io/github/WEgeophysics/pyCSAMT/requirements.svg?branch=develop)](https://requires.io/github/WEgeophysics/pyCSAMT/requirements/?branch=develop)
 [![Coverage Status](https://coveralls.io/repos/github/WEgeophysics/pyCSAMT/badge.svg?branch=devlop)](https://coveralls.io/github/WEgeophysics/pyCSAMT?branch=develop) ![GitHub](https://img.shields.io/github/license/WEgeophysics/pyCSAMT?color=blue&logo=GNU&logoColor=red) ![GitHub release (latest by date)](https://img.shields.io/github/v/release/WEgeophysics/pyCSAMT?color=orange) ![GitHub all releases](https://img.shields.io/github/downloads/WEgeophysics/pyCSAMT/total?color=green)

## Overview 

* **Definition**

CSAMT is geophysical method well-established  as resistivity exploration 
tool in deep geological structure detection. The method is broadly applied in  diverse of exploration problems such as mineral , hydrocarbon,  groundwater resources, 
as well as mapping the fault-zones etc. 

* **Purpose**

The software, unified open source software, contains bacics steps and improve CSAMT standard data processing as well as the modeling using [OCCAM2D](https://marineemlab.ucsd.edu/Projects/Occam/index.html).
The software also contains its inner database composed of geological structures and electrical properties of rocks, to generate  a pseudo-stratigraphy 2D map to enhance geophysical interpretation especially in more geological complex area ( with various tectonic accidents). 

* **Targets**

We hope this toolbox will help  the scientific geophysics community and those who work and grounwater exploration. In addition,  the toolbox aims  to the community of 
developers and users of that software.

 * **Note**
 
Actually pyCSAMT only works  in far field. Furhermore , it uses [OCCAM2D](https://marineemlab.ucsd.edu/Projects/Occam/index.html) open source sofware as modeling software. Nevertheless,
several  outputs are provided for other external modeling softwares like [MTpy](https://github.com/MTgeophysics/mtpy), [OasisMontaj](http://updates.geosoft.com/downloads/files/how-to-guides/Oasis_montaj_Gridding.pdf)
and [GoldenSoftware](https://www.goldensoftware.com/products/surfer).

## Documentation 
* API Documentation  : https://pycsamt.readthedocs.io/en/latest/
* Home Page : https://github.com/WEgeophysics/pyCSAMT/wiki
* Some Codes Implementation: https://github.com/WEgeophysics/pyCSAMT/wiki/How-pyCSAMT-works-%3F
* Installation Guide : https://github.com/WEgeophysics/pyCSAMT/wiki/pyCSAMT-installation-guide-for-Windows--and-Linux
* User Guide : https://github.com/WEgeophysics/pyCSAMT/blob/develop/docs/pyCSAMT%20User%20Guide.pdf


## Licence 
pyCSAMT is under GNU Lesser GPL version3 [LGPLv3](https://github.com/03-Daniel/pyCSAMT/blob/master/LICENSE.md).

## Some Features 
1. convert _*.avg_ , _*.j(.*dat)_ file  to SEG Electrical Data Interchange(EDI) file.
2. analysis and correct CSAMT static shift effects 
3. plot penetration depth 1D & 2D
4. plot pseudo-cross-section of resistivity and phase
4. plot resistivity model, plot forward response and residual model 
5. plot "pseudostratigraphy log" with true resistivity values 
6. generate true resistivity model from true resistivity values and true layer names 
7. rescale profile coordinates using  `natural|distorded`, `classic` or `equidistant` keywords 
8. read and write  EMAP or MT EDI data 
9. plot RMS, topography, azimuth and station-separation 
10. write corrected EDI,
11. write occam2D inputfiles calling [MTpy](https://github.com/MTgeophysics/mtpy.git) software. 

* **Available filters**
1. *Trimming moving average* (TMA) mostly used by [Zonge International Engineering](http://zonge.com/) .
2. *Fixed-length-dipole moving average* (FLMA) also used by [Zonge International Engineering](https://zonge.com.au/).
3. *Adaptative moving-average* (AMA) based on idea of [Torres-Verdin](https://sci-hub.se/http://dx.doi.org/10.1190/1.1443273).
4. *Removal distorsion* (`dist`)  and  *Static shift removal* (`ss`) usefull  filters to correct magnetotellurics (MT) data. 

## Units used    

* Frequency : [F] in Hz 
* Skin depth (sigma):  sigma  = 503 *sqrt([Rho]/[F]) in meters(m). 
* Apparent resistivy(Rho) : in Ω.m 
* E-field magnitude : [E]=  microvolt/meter (muv/m)
* H-field magnitude : [H] =  gamma /amp 
* Impedance Tensor [Z] in 2*2 matrices : [Z] = [E]/[H]:  km/s
* Angles : Theta in degrees clockwise 
* Location coordinates ( X =N-S , Y = E-W) in m. 
* Coordinates scaled in (UTM- Easting, Northing ) m. 
* Geomagnetic North : 0 degree azimuth 
* Step descent in m.
* Input true resistivities in Ω.m 

## A sample test using AMA and TMA  filters to correct raw *.edi files

AMA  and TMA filters can be used  to estimate average apparent resistivities at a single static-correction-reference frequency.
The following line of codes is an example to get new _*.edi_ corrected files from both filters application at each station,
refering to the EDI directory `data/edi/`.

![](https://github.com/WEgeophysics/pyCSAMT/blob/develop/quick_examples/wiki-images_quick_works/codes/demo_filter_ama_tma.PNG) 


The script below can be used to compare pseudo-cross-section of resistivity and phase of _corrected_edi outputs_ after `ama` & `tma` application  with 
_uncorrected edi_ . 

![](https://github.com/WEgeophysics/pyCSAMT/blob/develop/quick_examples/wiki-images_quick_works/codes/demo_edi_corrected.PNG) 


click [here](https://github.com/WEgeophysics/pyCSAMT/blob/develop/quick_examples/wiki-images_quick_works/codes/demo_filter_ama_tma.PNG) to see the output.

## Geophysical interpretation enhancement

If additional informations of survey area are available such as _true resistivity values_ as well as the _true layer names_, 
It's possible to used them to enhance your geophysical interpretation. There is feasibility to plot stratigraphy log 
under each station using `pycsamt.viewer.plot.Plot2d.plot_Pseudolog` or to write new resistivity model of entire survey line
using `pycsamt.geodrill.geoCore.geodrill.to_golden_software ` or `pycsamt.geodrill.geoCore.geodrill.to_asis_montaj` members from `Geodrill` module.
For instance :
 
* to plot stratigraphy log under each station,  we need to implement the command line below : 
 
![](https://github.com/WEgeophysics/pyCSAMT/blob/master/quick_examples/wiki-images_quick_works/codes/demo_plot-pseulog.PNG)

* to write new model of resistivity of survey line (here area is named  `some-where-place`), we merely need to import `.Geodrill` module  as:

![](https://github.com/WEgeophysics/pyCSAMT/blob/master/quick_examples/wiki-images_quick_works/codes/demo_geodrill.PNG) 

and after sucessfullly running, we will get the report below :

![](https://github.com/WEgeophysics/pyCSAMT/blob/master/quick_examples/wiki-images_quick_works/codes/demo_reports_geodrill.PNG)

                                                                      
* **Note** : Inversion input-files can be generated from _*.edi_ files using `pycsamt.modeling.occam2d.occam2d_write.buildingInputfiles` from `modeling.occam2d` module . 
            After applying the FDGC( Digital cartographic Standard for Geological Map Symbolization), click [here](https://github.com/WEgeophysics/pyCSAMT/blob/master/quick_examples/wiki-images_quick_works/interpretation.PNG)  to see your expected interpretation map.



## System requirements 
* Python 3.6+ 

## Contributors
  
1. Key Laboratory of Geoscience Big Data and Deep Resource of Zhejiang Province , School of Earth Sciences, Zhejiang University, China

2. Department of Geophysics, School of Geosciences and Info-physics, Central South University, China

3. Laboratoire de Géophysique Appliquée, UFR des Sciences de la Terre et des Ressources Minières, Université Félix Houphouët-Boigny, Cote d'Ivoire

* Developer's name:  [_Kouadio K. Laurent_](kkouao@zju.edu.cn), _etanoyau@gmail.com_: [1](http://www.zju.edu.cn/english/), [3](https://www.univ-fhb.edu.ci/index.php/ufr-strm/)
* Contibutors' names:
    *  [_Rong LIU_](liurongkaoyang@126.com) : [2](http://en.csu.edu.cn/)
    *  [_Binbin MI_](mibinbin@zju.edu.cn) : [1](http://www.zju.edu.cn/english/)
    *  [_Chun-ming LIU_](lifuming001@163.com): [2](http://en.csu.edu.cn/)
    *  [_Albert O. MALORY_](amalory@zju.edu.cn) :[1](http://www.zju.edu.cn/english/)
    
Any suggestion to improve the software is welcome ...

