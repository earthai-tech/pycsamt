# pyCSAMT : A Python open-source toolkit for Controlled Source Audio-frequency Magnetotellurics (CSAMT)

[![Documentation Status](https://readthedocs.org/projects/pycsamt/badge/?version=latest)](https://pycsamt.readthedocs.io/en/latest/?badge=latest) [![Build Status](https://travis-ci.com/WEgeophysics/pyCSAMT.svg?branch=master)](https://travis-ci.com/WEgeophysics/pyCSAMT) [![Requirements Status](https://requires.io/github/WEgeophysics/pyCSAMT/requirements.svg?branch=master)](https://requires.io/github/WEgeophysics/pyCSAMT/requirements/?branch=master)
  ![GitHub](https://img.shields.io/github/license/WEgeophysics/pyCSAMT?color=blue&logo=GNU&logoColor=red) ![GitHub release (latest by date)](https://img.shields.io/github/v/release/WEgeophysics/pyCSAMT?color=orange) 

## Overview 

* **Definition**

    CSAMT is geophysical method well-established  as resistivity exploration 
    tool in deep geological structure detection. The method is broadly applied in  diverse of exploration problems such as mineral , hydrocarbon,  groundwater resources, 
    as well as mapping the fault-zones etc. 

* **Purpose**

    The software contains bacics steps and improve CSAMT standard data processing and deals with [OCCAM2D](https://marineemlab.ucsd.edu/Projects/Occam/index.html) for modeling part.
    It also contains its inner database composed of geological structures and electrical properties of rocks,
    based on representative chart of  Palacky (1988) and the rock and mineral property classification of Slichter and Telkes (1942)
    to generate  a pseudo-stratigraphy log for drilling operations.

 * **Note**
 
    Actually pyCSAMT only works  in far field and several  outputs are provided for other external modeling softwares such as  [MTpy](https://github.com/MTgeophysics/mtpy), [OasisMontaj](http://updates.geosoft.com/downloads/files/how-to-guides/Oasis_montaj_Gridding.pdf)
    and [GoldenSoftware](https://www.goldensoftware.com/products/surfer).

## Documentation 
* API Documentation  : https://pycsamt.readthedocs.io/en/latest/
* Home Page : https://github.com/WEgeophysics/pyCSAMT/wiki
* Some examples: https://github.com/WEgeophysics/pyCSAMT/wiki/How-pyCSAMT-works-%3F
* Installation Guide : https://github.com/WEgeophysics/pyCSAMT/wiki/pyCSAMT-installation-guide-for-Windows--and-Linux
* User Guide : https://github.com/WEgeophysics/pyCSAMT/blob/develop/docs/pyCSAMT%20User%20Guide.pdf


## Licence 
pyCSAMT is under GNU Lesser GPL version3 [LGPLv3](https://github.com/03-Daniel/pyCSAMT/blob/master/LICENSE.md).


## Units used    

* Frequency : [F] in Hz 
* Skin depth (sigma):  sigma  = 503 *sqrt([Rho]/[F]) in meters(m). 
* Apparent resistivy(Rho) : in Ω.m 
* E-field magnitude : [E]=  microvolt/meter (muv/m)
* H-field magnitude : [H] =  gamma /amp 
* Impedance Tensor [Z] in 2*2 matrices : [Z] = [E]/[H]:  km/s
* Angle : Theta in degrees clockwise 
* Location coordinates ( X =N-S , Y = E-W) in m. 
* Coordinates in (UTM- Easting, Northing ) m. 
* Geomagnetic North : 0 degree azimuth 
* Step descent in m.
* Input true resistivities in Ω.m 

## Available filters 

1. Trimming moving average (TMA) mostly used by [Zonge International Engineering](http://zonge.com/) .
2. Fixed-length-dipole moving average (FLMA) also used by [Zonge International Engineering](https://zonge.com.au/).
3. Adaptative moving-average (AMA) based on idea of [Torres-Verdin](https://sci-hub.se/http://dx.doi.org/10.1190/1.1443273).
4. MT Removal distorsion (`dist`)  and  static shift removal (`ss`) filters basically used to correct magnetotellurics (MT) data. 
                                                               
## Plot inversion misfit and geo-stratigraphy misfit (misfit G)

To plot the `misfit` from measured data and the calculated inversion data, bring the _occam response file_ (_*.rep_) and  _Occamlogfile_ (optional _*.logfile_) then 
run the script below:
 
1. Plot some fitting curves of resistivity and phase inversion after applying on observed data
the static shift correction. 
```
>>> from pysamt.modeling.occam2d import plotResponse 
>>> resPath =r'data/inversionPath'                  # path to inversion files for each line
>>> plotResponse(data_fn =resPath,
...                 stations = ['S00', 'S04', 's08', 'S12']  # sites to visualize 
...                  rms =['1.013', '1.451', '1.00', '1.069'], # rms of each line
...                  error_type ='resi' )
``` 
Click [here](https://github.com/WEgeophysics/pyCSAMT/blob/develop/quick_examples/examplefitcurves.png) to see the reference output. 

2. To plot the `misfit`of the model response from the FE algorithms: 
```
>>> from pycsamt.modeling.occam2d import getMisfit 
>>> path_data ='data/occam2D'
>>> getMisfit(response_fn = os.path.join(path_data,'RESP17.resp' ),
...         logfile=os.path.join(path_data, 'LogFile.logfile'), 
...          data_fn = path_data)
```
To see the output, click [here](https://github.com/WEgeophysics/pyCSAMT/blob/develop/quick_examples/misfit.png).

2. To evaluate the model errors `misfit G` between the the new resistivity model or stratigraphy models(NMs) from inversion models(CRMs), 
set `plot_misfit` argument to `True` . `Misfit G` computation is the best way to see whether different layers with their corresponding resistivity values
are misclassified or not. With few step of codes we can check the process:
```
>>> from pycsamt.geodrill.geoCore.geodrill import Geostratigraphy
>>> inversion_files = {'model_fn':'data/Occam2DModel', 
                       'mesh_fn': 'data/Occam2DMesh',
                        "iter_fn":'data/ITER27.iter',
                       'data_fn':'data/OccamDataFile.dat'}
>>> resistivity_values =[10,  70, 180, 1000,  3000]   # resistivity values of layers to map
>>> layer_names =['river water','sedimentary rocks', 'fracture zone',  'gravel','igneous rocks']
>>> geosObj = GeoStratigraphy(**inversion_files,
...                      input_resistivities=resistivity_values, 
...                      input_layers=layer_names)
>>> geosObj.strataModel(kind='nm', plot_misfit =True)           # 'nm':New Model
```
click [here](https://github.com/WEgeophysics/pyCSAMT/blob/develop/quick_examples/geofit.png) for reference output. 


* **Note** : 
    For CSAMT data processing and some codes implementation,
    please refer to our [wiki pages](https://github.com/WEgeophysics/pyCSAMT/wiki/How-pyCSAMT-works-%3F).

## Credits

We use or link some third-party software (beside the usual tool stack: numpy, scipy, matplotlib) and are grateful for all the work made by the authors of these awesome open-source tools:
* mtpy: https://github.com/MTgeophysics/mtpy.git
* occam2d: https://marineemlab.ucsd.edu/Projects/Occam/index.html
* zonge softwares:
    - AMTAVG: http://www.zonge.com/legacy/DatPro.html/
    - ASTATIC: http://www.zonge.com/legacy/PDF_DatPro/Astatic.pdf

## System requirements 
* Python 3.6+ 

## Contributors
  
1. Key Laboratory of Geoscience Big Data and Deep Resource of Zhejiang Province , School of Earth Sciences, Zhejiang University, China, http://www.zju.edu.cn/english/

2. Department of Geophysics, School of Geosciences and Info-physics, Central South University, China,(http://www.zju.edu.cn/english/)

3. Equipe de Recherche Géophysique Appliquée, Laboratoire de Geologie Ressources Minerales et Energetiques, UFR des Sciences de la Terre et des Ressources Minières, Université Félix Houphouët-Boigny, Cote d'Ivoire.(https://www.univ-fhb.edu.ci/index.php/ufr-strm/)

* Developer: 1-3 [_Kouadio K. Laurent ~@Daniel03_](kkouao@zju.edu.cn), <_etanoyau@gmail.com_>
* Contibutors:
    *  2- [_Rong LIU_](liurongkaoyang@126.com) 
    *  1- [_Albert O. MALORY_](amalory@zju.edu.cn)   
    *  1- [_Chun-ming LIU_](lifuming001@163.com) 


