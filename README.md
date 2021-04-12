# pyCSAMT (Python  for Controlled Source Audio-frequency Magnetotellurics )
[![Documentation Status](https://readthedocs.org/projects/pycsamt/badge/?version=latest)](https://pycsamt.readthedocs.io/en/latest/?badge=latest) [![Build Status](https://travis-ci.com/WEgeophysics/pyCSAMT.svg?branch=master)](https://travis-ci.com/WEgeophysics/pyCSAMT) [![Requirements Status](https://requires.io/github/WEgeophysics/pyCSAMT/requirements.svg?branch=master)](https://requires.io/github/WEgeophysics/pyCSAMT/requirements/?branch=master)
 [![Coverage Status](https://coveralls.io/repos/github/WEgeophysics/pyCSAMT/badge.svg?branch=master)](https://coveralls.io/github/WEgeophysics/pyCSAMT?branch=master)    


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
* Some codes implementation: https://github.com/WEgeophysics/pyCSAMT/wiki/How-pyCSAMT-works-%3F
* Installation guide : https://github.com/WEgeophysics/pyCSAMT/wiki/pyCSAMT-installation-guide-for-Windows--and-Linux


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

```
>>> from from pycsamt.ff.processing.corr import shifting
>>> for _filter in ['tma', 'ama']:
        shifting().write_corrected_edi(
                        data_fn ='data/edi', 
                        number_of_points =7.,
                        reference_frequency=8192,      # in Hz
                        number_of_skin_depth=7,  
                        dipole_length =50.,            # in meter 
                        FILTER=_filter, 
                                        )
```
The script below can be used to compare pseudo-cross-section of resistivity and phase of _corrected_edi outputs_ after `ama` & `tma` application  with 
_uncorrected edi_ . 

```
>>> from pycsamt.viewer.plot import Plot2d
>>> contourRes= 1000.                       # resistivity contour delineation in Ω.m  
>>> for path2edi_obj in ['data/edi','data/_outputEDIFiltered_AMA','data/_outputFilteredEDI_TMA']:
>>>       Plot2d().pseudocrossResPhase(fn=path2edi_obj, 
                                delineate_resistivity=contourRes,
                                )

```
click [here](https://github.com/WEgeophysics/pyCSAMT/blob/master/quick_examples/filterstests.png) to see the output.

## Geophysical interpretation enhancement

If additional informations of survey area are available such as _true resistivity values_ as well as the _true layer names_, 
It's possible to used them to enhance your geophysical interpretation. There is feasibility to plot stratigraphy log 
under each station using `pycsamt.viewer.plot.plot_Pseudolog` or to write new resistivity model of entire survey line
using `pycsamt.geodrill.geoCore.geodrill.to_golden_software ` or `pycsamt.geodrill.geoCore.geodrill.to_asis_montaj` members from `Geodrill` module.
For instance :
 
```
>>> from pycsamt.geodrill.geoCore.geodrill import Geodrill 
>>> INPUT_LAYERS = ['river water', 'fracture zone', 'augen gneiss', 'altered rocks', 'granite']  
>>> INPUT_RESISTIVITIES =[66.,70., 180., 1235. , 2202., 7000.]      # in ohm.meters 
>>> STEP_DESCENT = 200                                              # in meters. see code implementation to get more info about this parameters. 
>>> inversion_kwargs =[                                             # occam2D inversion files of survey line
                    mesh_fn: 'data/occam2D/Occam2DMesh',
                    iter_fn : 'data/occam2D/ITER17.iter',
                    model_fn : 'data/occam2D/Occam2DModel',
                    data_fn : 'data/occam2D/OccamDataFile.dat'
                    ]                                 
>>> Geodrill( **inversion_kwargs , 
                 input_resistivities=INPUT_RESISTIVITIES, 
                 input_layers =INPUT_LAYERS ,
                 step_descent =STEP_DESCENT,
                 doi ='1km'                                  # depth of investigation 
                        ).to_golden_software(filename ='some-where-place',  # survey area name
                                            to_negative_depth='True')       # export depth to negative value

```
* **Note** : Inversion input-files can be generated from _*.edi_ files using `pycsamt.modeling.occam2d.occam2d_write.buildingInputfiles` from `modeling.occam2d` module . 
            After applying the FDGC( Digital cartographic Standard for Geological Map Symbolization), click [here](https://github.com/WEgeophysics/pyCSAMT/blob/master/quick_examples/wiki-images_quick_works/interpretation.PNG)  to see your expected interpretation map.


## System requirements 
* Python 3.6+ 

## Contributors
  
a. Key Laboratory of Geoscience Big Data and Deep Resource of Zhejiang Province , School of Earth Sciences, Zhejiang University, China

b. Department of Geophysics, School of Geosciences and Info-physics, Central South University, China

c. Laboratoire de Géophysique Appliquée, UFR des Sciences de la Terre et des Ressources Minières, Université Félix Houphouët-Boigny, Cote d'Ivoire

* Developer's name: a,c. [_**Kouadio K. Laurent**_](kkouao@zju.edu.cn), _etanoyau@gmail.com_
* Contibutors' names:
    * 2. [_**Liu RONG**_](liurongkaoyang@126.com) 
    * 1. [_**BinBin MI**_](mibinbin@zju.edu.cn)
    * 2. [_**Chun-Ming LIU**_](lifuming001@163.com)
    * 1. [_**Albert O. Malory**_](amalory@zju.edu.cn) 
    
Any suggestion to improve the software is welcome ...

