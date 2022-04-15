# pyCSAMT: A Python open-source toolkit for Controlled Source Audio-frequency Magnetotellurics (CSAMT)

[![Documentation Status](https://readthedocs.org/projects/pycsamt/badge/?version=latest)](https://pycsamt.readthedocs.io/en/latest/?badge=latest) [![Build Status](https://travis-ci.com/WEgeophysics/pyCSAMT.svg?branch=master)](https://travis-ci.com/WEgeophysics/pyCSAMT) [![Requirements Status](https://requires.io/github/WEgeophysics/pyCSAMT/requirements.svg?branch=master)](https://requires.io/github/WEgeophysics/pyCSAMT/requirements/?branch=master)
  ![GitHub](https://img.shields.io/github/license/WEgeophysics/pyCSAMT?color=blue&logo=GNU&logoColor=red) ![GitHub release (latest by date)](https://img.shields.io/github/v/release/WEgeophysics/pyCSAMT?color=orange) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5533467.svg)](https://doi.org/10.5281/zenodo.5533467)

## Overview 

* **Definition**

    CSAMT is a geophysical method well-established  as a resistivity exploration 
    tool in deep geological structure detection. The method is broadly applied in  diverse of exploration problems such as mineral, hydrocarbon,  groundwater resources, 
    as well as mapping the fault-zones etc. 

* **Purpose**

    The software contains basic steps, uses the CSAMT standard data processing and deals with [OCCAM2D](https://marineemlab.ucsd.edu/Projects/Occam/index.html) for the modeling part.
    The idea behind the development of this toolbox is to improve the groundwater exploration techniques and fight against the numerous unsucessful drillings mostly due to their wrong
    locations after geophysical surveys. The main goal, first, consists to minimize the use of supplement methods to CSAMT which commonly increases the operating budgets.
    Secondly, it entails to right locating the drilling thereby reducing the misinterpretation of modeling results(e.g., demarcating well the fracture zones). 
    Moreover, the capability of the software to estimate the layer thicknesses with a few margin of errors,
    could indirectly help the geophysical and drilling companies to reduce their loss when purcharsing the material for borehole
    equipment (e.g., PVC pipes). Thus, to meet the global objective, the toolboox  uses the previous informations of the survey area such as the boreholes/wells and 
    geological data combined with the inversion results to generate a predicted log called a pseudostratigraphic log for drilling operations.

 * **Note**
 
    At this time, pyCSAMT only works in far field and generates several outputs for other external modeling
     softwares such as  [MTpy](https://github.com/MTgeophysics/mtpy), [OasisMontaj](http://updates.geosoft.com/downloads/files/how-to-guides/Oasis_montaj_Gridding.pdf)
    and [GoldenSoftware](https://www.goldensoftware.com/products/surfer). 
    In addition, it is important to note the software, curently, does not to solve all the CSAMT problems. The other functionalities 
    (e.g., use of Tipper data in the desert environment, also the CSEM,  etc.) will be added as
     the development of the software continues to progress.

## Documentation 

* API Documentation  : https://pycsamt.readthedocs.io/en/latest/
* Home Page : https://github.com/WEgeophysics/pyCSAMT/wiki
* Some examples: https://github.com/WEgeophysics/pyCSAMT/wiki/How-pyCSAMT-works-%3F
* Installation Guide : https://github.com/WEgeophysics/pyCSAMT/wiki/pyCSAMT-installation-guide-for-Windows--and-Linux
* User Guide : https://github.com/WEgeophysics/pyCSAMT/blob/develop/docs/pyCSAMT%20User%20Guide.pdf

## Licence 

pyCSAMT is under GNU Lesser GPL version3 [LGPLv3](https://github.com/03-Daniel/pyCSAMT/blob/master/LICENSE.md).

## Installation 

Use [PyPI release](https://pypi.org/project/pycsamt/) for quick installation 
* `$ pip install pycsamt` or
* `$ pip install --user pycsamt` (Window users)

However, it is recommended the installation from the repository to get the latest development code. 

## Available filters 

1. Trimming moving average (TMA) mostly used by [Zonge International Engineering](http://zonge.com/) .
2. Fixed-length-dipole moving average (FLMA) also used by [Zonge International Engineering](https://zonge.com.au/).
3. Adaptative moving-average (AMA) based on idea of [Torres-Verdin](https://sci-hub.se/http://dx.doi.org/10.1190/1.1443273).
4. MT Removal distorsion (`dist`)  and  static shift removal (`ss`) filters basically used to correct magnetotellurics (MT) data.

For example, to remove the static shift effect into a corrupted [SEG](https://seg.org/) EDI data, user can apply one of above filters (e.g., FLMA)
to correct the EDI data located in `data/edi` directory as the below command line (CLI) indicates. 
``` 
 $ staticshift data/edi -ft flma --ndipole 5 --dipole-length=50  
``` 
                                                               
## Plot inversion and geostratigraphy(G) misfits

For this quickstart, we assume the user has already the forward modeling files (`*.resp`, `*.dat`, `*.mesh`, `*.iter`
and `*.logfile`(optional)), and the supplement data (i.e. boreholes/wells and geological data) collected in the survey area.
 
1. Plotting some fitting curves of resistivity and phase inversion after applying on observed data
the static shift correction. For instance, we visualize the fitting curves of four survey lines with their corresponding 
RMS. The stations (e.g. `S00`,`S04`, `s08`, and `S12`) are randomly chosen of each survey line: 
```
$ fitforward  data/inversionFiles --stations S00 S04 s08 S12 --rms 1.013 1.451 1. 1.069 
```

Click [here](https://github.com/WEgeophysics/pyCSAMT/blob/develop/examples/examplefitcurves.png) to see the reference output. 

2. Plotting the `misfit` of the model response from the forward modeling (FM) of line`K1`. Set `kind` argument to `phase` for the error in phase.
```
$ misfit2d data/inversionFiles/K1.dat data/inversionFiles/K1.resp --kind=rho 
```

To see the output, click [here](https://github.com/WEgeophysics/pyCSAMT/blob/develop/examples/misfit.png).

3. Construction of the new resistivity model (NM) also called the geostratigraphy model. The latter is built from 
the additional data (e.g. borehole/well data and/or geological data) collected in the exploration area combined with the FM.
The best approach to fastly build the NM model is to gather all the data collected in the exploration area into a 
single configuration file in `*.json` or `*.yml` format. For illustrating, 
we will use `mysurveyarea_data.yml` file like: 
```
# mysurveyarea_data.yml
input_layers:
    - river water
    - sedimentary rocks
    - fracture zone
    - gravel
    - granite
    - igneous rocks
    - basement rocks
input_resistivities:
    - 10
    - 66
    - 70
    - 180
    - 1000
    - 3000
    - 7000
data_fn: data/occam2D\K1.dat
iter_fn: data/occam2D\K1.iter
mesh_fn: data/occam2D\Occam2DMesh
model_fn: data/occam2D\Occam2DModel
ptol: 0.1
beta: 5
n_epochs: 100
build: true
```
where `(data_fn, iter_fn, mesh_fn ,model_fn)` and `(ptol, beta, n_epochs, build)` are OCCAM2D inversion files resulting from FM (CRM) and the constructor parameters 
respectively. From the CLI below, NM is created. 
```
$ nm --config mysurveyarea_data.yml --show
```

4. Error Evaluation. To evaluate the geostratigraphy model error called `Misfit G` between the the NM  and CRM,
 one needs to add the following argument `--misfit` to the previous command. For instance 
 in the case where the configure file is written in `*.json` format, the CLI becomes:
```
$ nm -c mysurveyarea_data.json --show --misfit
```
Indeed, the `Misfit G` helps to determine whether the different layers with their corresponding resistivity values
are misclassified or not. An example of misfit G map can be visualized [here](https://github.com/WEgeophysics/pyCSAMT/blob/develop/examples/geofit.PNG). 

* **Note** : 
    For CSAMT data processing and an example of a deep implementation, please refer to our [wiki page](https://github.com/WEgeophysics/pyCSAMT/wiki/How-pyCSAMT-works-%3F).

## Plot the pseudostratigraphic log 

 Once the NM is built, user just needs to query the software parameters memory (SPM) to fech at each station 
 the corresponding pseudostratigraphic log. An example of the
  command to do this task (e.g. extracting the log of station `S00`) is:
```
$ pseudostratigraphic -s S00 
```
The output below with layer thicknesses estimation are displayed.

``` 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~[ PseudoStratigraphic Details: Station = S00 ]~~~~~~~~~~~~~~~~~~~~~~~~~~~~
------------------------------------------------------------------------------------------------------
|      Rank |            Stratum             |         Thick-range(m)         |     Thickness(m)     |
------------------------------------------------------------------------------------------------------
|        1. |         fracture zone          |         0.0 ----- 6.0          |         6.0          |
|        2. |             gravel             |         6.0 ----- 13.0         |         7.0          |
|        3. |            granite             |        13.0 ----- 29.0         |         16.0         |
|        4. |         igneous rocks          |        29.0 ----- 49.0         |         20.0         |
|        5. |         basement rocks         |        49.0 ----- 249.0        |        200.0         |
|        6. |         igneous rocks          |       249.0 ----- 289.0        |         40.0         |
|        7. |            granite             |       289.0 ----- 529.0        |        240.0         |
|        8. |         igneous rocks          |       529.0 ----- 699.0        |        170.0         |
|        9. |         basement rocks         |       699.0 ----- 999.0        |        300.0         |
------------------------------------------------------------------------------------------------------
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Survey Line: Occam2D files properties ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
|model = Occam2DModel     |iter  = ITER17.iter      |mesh  = Occam2DMesh      |data  = OccamDataFile.dat|
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
``` 
Click [here](https://github.com/WEgeophysics/pyCSAMT/blob/develop/examples/pseudostratigraphic_log.PNG) to see the predicted log. 

Obviously, it does not make sense to expect to drill until to reach `1km` depth. Therefore, another feature is implemented 
to help the user to only fetch from the SPM the most interesting part of the predicted log for a specific purpose. To do such task, one needs to fiddle with the `zoom` parameter. 
For instance, the CLI below with `zoom=25%` only displays the first `250m` assuming that the investigation depth  is `1000m` maximum. 
```
$ pseudostratigraphic --station=S00 --zoom 25%
```
Check the following [ouput](https://github.com/WEgeophysics/pyCSAMT/blob/develop/examples/zoom25.PNG) to see the new log. 
Futhermore,  it's also possible to provide the top (e.g. `10m`) and the bottom(e.g. `120m`) of the log for visualization as: 
```
$ pseudostratigraphic --station S00 --zoom 10 120 --fontsize 12
```

## Credits

We use or link some third-party software (beside the usual tool stack: numpy, scipy, matplotlib) and are grateful for all the work made by the authors of these awesome open-source tools:
* MTpy: https://github.com/MTgeophysics/mtpy.git
* Occam2D: https://marineemlab.ucsd.edu/Projects/Occam/index.html
* Zonge Engineering softwares:
    - AMTAVG: http://www.zonge.com/legacy/DatPro.html/
    - ASTATIC: http://www.zonge.com/legacy/PDF_DatPro/Astatic.pdf

## System requirements 
* Python 3.7+ 

## Citations 

If the software seems useful to you, you may cite the authors as a contribution.

* Kouao Laurent Kouadio, Rong Liu, Binbin Mi, Chun-ming Liu. pyCSAMT: An alternative Python toolbox for groundwater 
  exploration using controlled source audio-frequency magnetotelluric; Journal of Applied Geophysics;
  104647(2022); https://doi.org/10.1016/j.jappgeo.2022.104647.


## Contributors
  
1. Key Laboratory of Geoscience Big Data and Deep Resource of Zhejiang Province, School of Earth Sciences, [Zhejiang University](http://www.zju.edu.cn/english/), China.

2. Department of Geophysics, School of Geosciences and Info-physics, [Central South University](http://www.zju.edu.cn/english/), China.

3. Equipe de Recherche Géophysique Appliquée, Laboratoire de Geologie Ressources Minerales et Energetiques, UFR des Sciences de la Terre et des Ressources Minières, [Université Félix Houphouët-Boigny]( https://www.univ-fhb.edu.ci/index.php/ufr-strm/), Cote d'Ivoire.

* Developer: 1, 3- Kouadio K. Laurent; <etanoyau@gmail.com>, <kkouao@zju.edu.cn>,
* Contributors:
    *  2- Rong LIU; <liurongkaoyang@126.com>
    *  1- Binbin MI; <mibinbin1991@126.com>  


