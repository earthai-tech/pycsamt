# pyCSAMT: A Python watex-exploration toolbox using controlled Source Audio-frequency Magnetotellurics (CSAMT)
[![Documentation Status](https://readthedocs.org/projects/pycsamt/badge/?version=latest)](https://pycsamt.readthedocs.io/en/latest/?badge=latest) [![Build Status](https://travis-ci.com/WEgeophysics/pyCSAMT.svg?branch=master)](https://travis-ci.com/WEgeophysics/pyCSAMT)
[![Requirements Status](https://requires.io/github/WEgeophysics/pyCSAMT/requirements.svg?branch=master)](https://requires.io/github/WEgeophysics/pyCSAMT/requirements/?branch=master)
  ![GitHub](https://img.shields.io/github/license/WEgeophysics/pyCSAMT?color=blue&logo=GNU&logoColor=red) ![GitHub release (latest by date)](https://img.shields.io/github/v/release/WEgeophysics/pyCSAMT?color=orange)  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5533467.svg)](https://doi.org/10.5281/zenodo.5533467)  

_For CSAMT standard data processing, modeling and groundwater exploration (GWE) enhancement techniques._

## Overview 


* [Definition](#https://ui.adsabs.harvard.edu/abs/2018EGUGA..2013744L/abstract)

CSAMT is an enhanced frequency domain of EM program using synchronous stacking and averaging, and Fourier integral methods to 
improve the signal to noise ratio. Later, CSAMT is well appreciated in geophysical commumity and
was used as a suitable exploration method well-established in deep geological structure detection.
Today the method is broadly applied in  diverse of exploration problems such as mineral , hydrocarbon,  groundwater resources, 
as well as mapping the fault-zones etc. 

* [Purpose](#Purpose)

 `pycsamt` contains basic steps and filters for  CSAMT standard data processing and deals 
 with  [OCCAM2D](https://marineemlab.ucsd.edu/Projects/Occam/index.html) modeling software.
 The package also includes  a database geological structures and electrical properties of rocks,
 based on representative chart of  Palacky (1988) and the rock and mineral property classification of Slichter and Telkes (1942)
to generate  a pseudo-stratigraphy log for drilling operations.


 * [Note](#https://iopscience.iop.org/article/10.1088/1742-6596/1127/1/012021)
 
The software works  in [far field](https://electronics.stackexchange.com/questions/487691/why-are-e-and-b-field-in-phase-in-far-field-electromagnetic-wave-propagation)
and provided several for other external modeling softwares like [MTpy](https://github.com/MTgeophysics/mtpy), [OasisMontaj](http://updates.geosoft.com/downloads/files/how-to-guides/Oasis_montaj_Gridding.pdf)
and [GoldenSoftware](https://www.goldensoftware.com/products/surfer).

## Installation 

Use `pip` for quick installation:
```
$ pip install pycsamt| $ pip install user pycsamt

``` 
One can follow the following [installation guide](https://github.com/WEgeophysics/pyCSAMT/wiki/pyCSAMT-installation-guide-for-Windows--and-Linux).

## Quickstart 

Apply Adaptative moving-average (AMA) of [Torres-Verdin](https://sci-hub.se/http://dx.doi.org/10.1190/1.1443273) to correct  [SEG](https://seg.org/) 
Electrical Data Interchange(EDI) (e.g., EDI-files=`data/edi/*.edi`) polluted by the static shift effect by ruunning: 
``` 
$ staticshift data/edi -ft ama --nskin 3 

```
Build your [OCCAM2D](https://marineemlab.ucsd.edu/Projects/Occam/index.html) input files from EDI using the command lines (CLI) below:
```
$ occambuildinputs data/edi -mode=6 -niter 112 -cw=7 --nlayers=32 -z=1000 -zb=5000  --ifreq

``` 
Use your forward modeling (`*.resp`) and data (`*.data`) files  to plot misfit2D map  with a few step of command:  
```

$ misfit2d data/inversionFiles/K1.dat data/inversionFiles/K1.resp  

```
The most interesting part is the use of the collection of borehole/wells and geology data combined with the forward modeling to 
build the stratigraphic model of the exploration to right locate the drilling after survey. This will minimize the rate of uncessufull drilling 
and better depict the fracture zone known as the target during the GWE. Before, we prepare a `modelconfig.json` to gather all the informations collected in the exploration area like: 
```
# modelconfig.json
{
  "input_layers": [
    "river water",
    "fracture zone",
    "MWG",
    "LWG",
    "granite",
    "igneous rocks",
    "basement rocks"
  ],
  "input_resistivities": [
    10,
    66,
    70,
    100,
    1000,
    3000
  ],
  "data_fn": "data/occam2D\\OccamDataFile.dat",
  "iter_fn": "data/occam2D\\ITER17.iter",
  "mesh_fn": "data/occam2D\\Occam2DMesh",
  "model_fn": "data/occam2D\\Occam2DModel",
  "ptol": 0.1,
  "beta": 5,
  "n_epochs": 100,
  "build": true
}

```
Now, with `modelconfig.json` we can now build our stratigraphic model via a few step of command below: 
```
$ nm -c modelconfig.json --show

```
To see the error between the stratigraphic model (model predicted) and the forward modeling (occam2d model), we merly need to add `--misfit` as an argument to the previous command. 
Finally to fetch from each station the predicted log (for instance the station `S10`), we just need to run the command:
```
$ pseudostratigraphic --station=S10 --zoom=25%

```
where the `zoom` parameter indicates the most interesting part of the log. For instance `zoom=25%` shows the first 1/4 of investigation depth (DOI)
i.e. if `DOI=1000m`, only the `250 m` should be displayed.

For a deep implementation,  please refer to our [wiki page](https://github.com/WEgeophysics/pyCSAMT/wiki).
 
## Contributing

Anyone who want to enhance the project is welcome and he/she will find more informations on [issue template](https://github.com/WEgeophysics/pyCSAMT/blob/master/ISSUE_TEMPLATE.md) file located in [git repository](https://github.com/WEgeophysics/pycsamt).