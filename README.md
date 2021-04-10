# pyCSAMT (Python  for Controlled Source Audio-frequency Magnetotellurics )
[![Documentation Status](https://readthedocs.org/projects/pycsamt/badge/?version=latest)](https://pycsamt.readthedocs.io/en/latest/?badge=latest) [![Build Status](https://travis-ci.com/WEgeophysics/pyCSAMT.svg?branch=master)](https://travis-ci.com/WEgeophysics/pyCSAMT) [![Requirements Status](https://requires.io/github/pytest-dev/pytest-cov/requirements.svg?branch=master)](https://requires.io/github/pytest-dev/pytest-cov/requirements/?branch=master)
     


## Overview 

* [Definition](#Definition)

CSAMT is geophysical method well-established  as resistivity exploration 
tool in deep geological structure detection. The method is broadly applied in  diverse of exploration problems such as mineral , hydrocarbon,  groundwater resources, 
as well as mapping the fault-zones etc. 

* [Purpose](#Purpose)

The software, unified open source software, contains bacics steps and improve CSAMT standard data processing as well as the modeling using [OCCAM2D](https://marineemlab.ucsd.edu/Projects/Occam/index.html).
The software also contains its inner database composed of geological structures and electrical properties of rocks, to generate  a pseudo-stratigraphy 2D map to enhance geophysical interpretation especially in more geological complex area ( with various tectonic accidents). 

* [Targets](#Targets)

We hope this toolbox will help  the scientific geophysics community and those who work and grounwater exploration. In addition,  the toolbox aims  to the community of 
developers and users of that software.

 * [Note](#Note)
 
Actually pyCSAMT only works  in far field. Furhermore , it uses [OCCAM2D](https://marineemlab.ucsd.edu/Projects/Occam/index.html) open source sofware as modeling software. Nevertheless,
several  outputs are provided for other external modeling softwares like [MTpy](https://github.com/MTgeophysics/mtpy), [OasisMontaj](http://updates.geosoft.com/downloads/files/how-to-guides/Oasis_montaj_Gridding.pdf)
and [GoldenSoftware](https://www.goldensoftware.com/products/surfer).

## Documentation 
* API Documentation  : https://pycsamt.readthedocs.io/en/latest/
* Home Page : https://github.com/WEgeophysics/pyCSAMT/wiki
* Exercices: https://github.com/WEgeophysics/pyCSAMT/wiki/How-pyCSAMT-works-%3F
* Installation : https://github.com/WEgeophysics/pyCSAMT/wiki/pyCSAMT-installation-guide-for-Windows--and-Linux


## Licence 
pyCSAMT is under GNU Lesser GPL version3 [LGPLv3](https://github.com/03-Daniel/pyCSAMT/blob/master/LICENSE.md).

## Some Features 
1. convert *.avg, *.j(.*dat) file  to SEG Electrical Data Interchange(EDI)
2. analysis and correct CSAMT static shift effects 
3. plot penetration depth(1D, 2D ) as function of depth 
4. plot pseudo-cross section of phase and  resistivy 
4. plot resistivity model /plot roughness model and residual model 
5. plot "pseudostratigraphy log " with true resistivity values 
6. generate true resistivity model from true resistivity values and true layer names 
7. rescale profile coordinates using  `natural|distorded`, `classic` or `equidistant` keywords 
8. read /Write  EMAP or MT EDI data format 
9. plot RMS, topography , azimuth and station-separation 
10. write corrected EDI,
11. write occam2D inputfiles by straightforwardly calling [MTpy](https://github.com/MTgeophysics/mtpy.git) software. 

* **Available filters**
1. *Trimming moving average* (TMA) mostly used by [Zonge International Engineering](http://zonge.com/) .
2. *Fixed-length-dipole moving average* (FLMA) also used by [Zonge International Engineering](https://zonge.com.au/).
3. *Adaptative moving-average* (AMA) based on idea of [Torres-Verdin](https://sci-hub.se/http://dx.doi.org/10.1190/1.1443273).
4. *Removal distorsion* (`dist`) usefull to correct magnetotellurics (MT) data. 
5. *Static shift removal* (`ss`) 

## Units used    
* Apparent resistivy(Rho) : in ohm.meter 
* Frequency : [F] in Hertz 
* Skin depth (sigma):  sigma  = 503 *sqrt([Rho]/[F]) in meter  
* E-field magnitude : [E]=  microvolt/meter (muv/m)
* H-field magnitude : [H] =  gamma /amp 
* Impedance Tensor [Z] in 2*2 matrices : [Z] = [E]/[H]:  km/s
* Angles : Theta in degrees clockwise 
* Location Coordinates in meters ( X =North to south , Y = East -West)
* Coordinates scaled in meters (UTM- Easting , northing )
* Geomagnetic North : 0 degree azimuth 
* Step Descent : SD param  in meter 
* Input true resistivity in ohm-meter

## A sample test using AMA and TMA  filter to correct raw *.edi files
AMA  and TMA filter can be used  to estimate average apparent resistivities at a single static-correction-reference frequency.
We firstly corrected all _*.edi_ files at each station using the following line of codes:

```
>>> from from pycsamt.ff.processing.corr import shifting
>>> for _filter in ['tma', 'ama']:
        corr_obj= shifting().write_corrected_edi(data_fn ='edi/data', 
                             number_of_points =7.,
                             reference_frequency=8192,      # in Hz
                             number_of_skin_depth=7,  
                             dipole_length =50.,            # in meter 
                             FILTER=_filter, 
                             )
```
Then use the corrected_edi outputs for pseudo-cross-section of resistivity and phase  plots and compare  the raw outputs with 
the corrected outputs (`ama` & `tma`) 

```
>>> from pycsamt.viewer.plot import Plot2d
>>> edipath ='data/edi'                     # Current work directory assume to be /pyCSAMT
>>> contourRes= 1000.                       # resistivity contour delineation value in ohm.meters 
>>> for path2edi_obj in ['data/edi', 'data/_outputFilteredEDI_TMA','data/_outputEDIFiltered_AMA']:
>>>       Plot2d().pseudocrossResPhase(fn=path2edi_obj, 
                                delineate_resistivity=contourRes,
                                plot_style =plotStyle)

```
click [here](https://github.com/WEgeophysics/pyCSAMT/blob/master/quick_examples/filterstests.png) to see the following output : [-->ref-output](https://github.com/WEgeophysics/pyCSAMT/blob/master/quick_examples/filterstests.png)

## Geophysical interpretation enhancement 
After inversion, if some geological informations are availbale such _true resistivities_ as well as the _true layer names_, used them to enhance your geophysical interpretation. 
This is possible using `plot_Pseudolog` function.
 
```
>>> from pycsamt.viewer.plot  import Plot2d
>>> INPUT_LAYERS = ['river water', 'fracture zone', 'augen gneiss', 'altered rocks', 'granite'] # layer' names collected 
>>> INPUT_RESISTIVITIES =[66.,70., 180., 1235. , 2202., 7000.]        # in ohm.meters 
>>> STEP_DESCENT = 200                                                # in meters
>>> Plot2d().plot_Pseudolog( 
        mesh_fn= 'data/occam2D/Occam2DMesh',
        iter_fn = 'data/occam2D/ITER17.iter',
        model_fn = 'data/occam2D/Occam2DModel',
        data_fn = 'data/occam2D/OccamDataFile.dat' , 
                            station_id= 'S43', 
                            input_resistivities=INPUT_RESISTIVITIES, 
                            input_layers =INPUT_LAYERS,
                            step_descent =STEP_DESCENT,
                            doi ='1km')

```
After applying the FDGC( Digital cartographic Standard for Geological Map Symbolization), click [here](https://github.com/WEgeophysics/pyCSAMT/blob/master/quick_examples/wiki-images_quick_works/interpretation.PNG)  to see your interpretation map.


## System requirements 
* Python 3.6+ 

## Contributing 
Any suggestion to improve the software is welcome ...


*_Developer's name: **_@Daniel03_**, _etanoyau@gmail.com_
