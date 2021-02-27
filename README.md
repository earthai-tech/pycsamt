# pyCSAMT
Python  for Controlled Source Audio-frequency Magnetotelluric (CSAMT)

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
 
For the first release , pyCSAMT is a basic toolbox and work only  in far field. It 
does not  contain all data processing and analysis of CSAMT method.
We hope that the overview of this toolbox will bring helps to many people 
working with CSAMT method. 

## Licence 
pyCSAMT is under GNU Lesser GPL version3 [LGPLv3](https://github.com/03-Daniel/pyCSAMT/blob/master/LICENSE.md)

## Features 
1. convert *.avg, *.j(.data)  to SEG Electrical Data Interchange(EDI)
2. analysis and CSAMT static shift correction with TMA filter
3. plot penetration depth(1D, 2D ) as function of depth 
4. plot pseudocross section of phase and apparent resistivy 
4. plot resistivity model /plot roughness model and residual model 
5. plot "pseudostratigraphy log " with true resistivity values 
6. generate true resistivity model from true resistivity values and layer names. 
7. rescaling profile coordinates 
7. read /Write EMAP EDI data format 

## Units used    
* apparent resistivy(Rho) : in ohm.meter 
* Frequency : F in Hertz 
* skin depth (sigma):  sigma  = 503 *sqrt([Rho]/[F]) in meter  
* E-field magnitude : [E]=  microvolt/meter (muv/m)
* H field magnitude : [H] =  gamma /amp 
* Impedance Tensor [Z] in 2*2 matrices : [Z] = [E]/[H]:  km/s
* angles : Theta in degrees and clockwise 
* Location Coordinates in meters ( X =North to south , Y = East -West)
* coordinates scaled in meters (UTM- Easting , northing )
* Geomagnetic North : 0 degree azimuth 


## Contact
* Kouadio Laurent
* etanoyau@gmail.com
* kkouao@zju.edu.cn

## System requirements 
* Python 3.6  
* Python 3.7


## Contributing 
Contributions are welcome from anyone. Your suggestions and 
comments  are welcome to improve the software. 

## Four quick examples 
These examples are really quicks examples, we let users discover full of functionalities.
* Get pseudocross resistivity 2D map with contours 
```
>>> from viewer.plot import Plot2d 
>>> Plot2d().pseudocrossResPhase(fn = data/edi ) 
...output :
```
 ![image](https://user-images.githubusercontent.com/59920007/109303602-02993080-7876-11eb-8e84-927d6efa1184.png)

* Apply Trimming Moving Average (TMA) Filter  to correct apparent resistivity 
```
>>> from viewer.plot import Plot1d 
>>> reference_frequency, TMApoints  = 8000., 7
>>> Plot1d().plot_static_correction (data_fn = data/edi, 
                                frequency_id =reference_frequency,
                                number_of_TMApoints = TMApoints,
                                fill_between =True )
...output :
```
![image](https://user-images.githubusercontent.com/59920007/109303645-12b11000-7876-11eb-811b-8c5c47f2376f.png)
* 2D Inversion Model 
Visualize a 2D resistivity model at 1km depth   requires a few steps of Python shell:
```
>>> plot2d().plot_occam2dModel(
                     mesh_fn= data/occam2D/Occam2DMesh, 
                    	iter_fn = data/occam2D/'ITER17.iter, 
                    	model_fn = data/occam2D/Occam2DModel, 
                   	 data_fn = data/occam2D/OccamDataFile.dat, doi =’1km’)
...output :
```
![image](https://user-images.githubusercontent.com/59920007/109303739-2d838480-7876-11eb-8f38-ade2c052ddb3.png)
* Build ‘pseudo stratigraphy log’, outputs files
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

![image](https://user-images.githubusercontent.com/59920007/109303758-35432900-7876-11eb-9cfe-9b6563f31ce5.png)

