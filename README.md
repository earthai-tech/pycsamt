# pyCSAMT (Python  for Controlled Source Audio-frequency Magnetotelluric )
[![Documentation Status](https://readthedocs.org/projects/pycsamt/badge/?version=latest)](https://pycsamt.readthedocs.io/en/latest/?badge=latest)
     


## Overview 

* [Definition](#Definition)

CSAMT is geophysical method well-established  as resistivity exploration 
tool  for the detection in deep geological structure  within a geophysic community.
The method is broadly applied in a diversity of exploration problems such as mineral , hydrocarbon,  groundwater resources, 
as well as mapping the fault-zones etc. 

* [Purpose](#Purpose)

The software try to overcome the lack of  unified software as open source software 
for CSAMT data processing in the academic community and to enhance geophysical interpretation. 

* [Targets](#Targets)

Firstly it's the scientific geophysics community.Secondly our regards point out to the community of 
developers and users of that software.

 * [Note](#Note)
 
 For the first release, The python toolbox gives a basic tools and  works only in far field. Furhermore , it uses  [OCCAM2D](https://marineemlab.ucsd.edu/Projects/Occam/index.html) open source sofware as modeling software , nevertheless several  outputs are provided for other modeling softwares.  
## Documentation 
* API Documentation  :  

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


## A quick illustration of "pseudo-stratigraphy log" 
An inner database allows to generate a pseudo-stratigraphy log based on 2D model 
resistivity  by including  the true resistivity values,  well data and/or borehole data from field or collected from firms.
Few steps of Python shell to visualize a "pseudo-stratigraphy log" of site 01  at 1km depth (DOI).
"step descent param " is the distance in deep where the program breaks and forces inverted model values to be close to input true resistivity values fit by boreholes/well data.

```
>>> DOI = 1km 
>>> filename ='ybro'
>>> station_to_visualize =1 
>>> STEP_DESCENT = 200 
>>> INPUT_RESISTIVITIES = [312, 525, 1235, 2202., 4000., 7000.] 
>>> INPUT_LAYERS = ['alluvium', 'amphibolite','altered rock','augen gneiss', 'granite']
>>>> plot2d().plot_Pseudolog( station_id= station_to_visualize, 
                            input_resistivities=INPUT_RESISTIVITIES, 
                            input_layers =INPUT_LAYERS ,
                            step_descent =STEP_DESCENT,
                            doi =DOI)
...output :
````

![fig3](https://user-images.githubusercontent.com/59920007/109377936-923bef00-7909-11eb-97bb-ad800176d94b.png)

## System requirements 
* Python 3.6 /3.7 

## Contributing 
Your suggestions are really welcome...


*_Contributor's name_: ***@Daniel03*** , _etanoyau@gmail.com_
