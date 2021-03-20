# pyCSAMT (Python  for Controlled Source Audio-frequency Magnetotellurics )
[![Documentation Status](https://readthedocs.org/projects/pycsamt/badge/?version=latest)](https://pycsamt.readthedocs.io/en/latest/?badge=latest) [![Build Status](https://travis-ci.com/WEgeophysics/pyCSAMT.svg?branch=master)](https://travis-ci.com/WEgeophysics/pyCSAMT)
     


## Overview 

* [Definition](#Definition)

CSAMT is geophysical method well-established  as resistivity exploration 
tool  for the detection in deep geological structure. The method is broadly applied in a diversity of exploration problems such as mineral , hydrocarbon,  groundwater resources, 
as well as mapping the fault-zones etc. 

* [Purpose](#Purpose)

The software, unified open source software, contains bacics steps to CSAMT standard data processing  and  and try to improve it as well as the modeling using OCCAM2D.
The software also contains its inner database composed of geological structures, electrical properties of rocks, to generate  a pseudo-stratigraphy 2D map to enhance geophysical interpretation especially in more geological complex area( with various tectonic accidents). 

* [Targets](#Targets)

The first target is  the scientific geophysics community and those who work and grounwater exploration and secondly, the toolbox points out to the community of 
developers and users of that software.

 * [Note](#Note)
 
The Python toolbox gives a basic tools and  works only in far field. Furhermore , it uses  [OCCAM2D](https://marineemlab.ucsd.edu/Projects/Occam/index.html) open source sofware as modeling software , nevertheless several  outputs are provided for other external modeling softwares like [OasisMontaj](http://updates.geosoft.com/downloads/files/how-to-guides/Oasis_montaj_Gridding.pdf) and [GoldenSoftware](https://www.goldensoftware.com/products/surfer).

## Documentation 
* API Documentation  : https://pycsamt.readthedocs.io/en/latest/

## Licence 
pyCSAMT is under GNU Lesser GPL version3 [LGPLv3](https://github.com/03-Daniel/pyCSAMT/blob/master/LICENSE.md)

## Some Features 
1. convert *.avg, *.j(.data)  to SEG Electrical Data Interchange(EDI)
2. analysis and CSAMT static shift correction with TMA filter
3. plot penetration depth(1D, 2D ) as function of depth 
4. plot pseudocross section of phase and apparent resistivy 
4. plot resistivity model /plot roughness model and residual model 
5. plot "pseudostratigraphy log " with true resistivity values 
6. generate true resistivity model from true resistivity values and layer names. 
7. rescale profile coordinates 
8. read /Write EMAP EDI data format 

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

## A sample illustration using FLMA filter application 
A fixed-length-moving-average (FLMA)  filter can be used  to estimate average apparent resistivities at a single static-correction-reference frequency.
A  few line of codes yields the following output: 
```
>>> from viewer.plot import Plot2d
>>> edipath =data/                 # Current work directory assume to be os.path.join(os.environ["pyCSAMT"],'data')
>>> contourRes= 1000.              #contour delinetation value in ohm.meters 
>>> plotStyle =None                # default `imshow` can be `pcolormesh`.
>>> for path2edi_obj in [data/edi, data/correctedEdi]:
>>>       plot2d().pseudocrossResPhase(fn=path2edi_obj, 
                                delineate_resistivity=contourRes,
                                plot_style =plotStyle)
...
```
![image](https://user-images.githubusercontent.com/59920007/111862592-33f6af00-8991-11eb-994d-43039d2345bb.png)

Another filter like TMA (trimmed-moving-average) can be applied to correct apparent resistivities. MT data can also be read and corrected using  `ss` for static shift removal and `dist` for distorsion removal).

## System requirements 
* Python 3.6+ 

## Contributing 
Your suggestions are really welcome...


*_Developer's name:_ ***@Daniel03*** , _etanoyau@gmail.com_
