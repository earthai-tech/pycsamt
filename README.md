# pyCSAMT: A Python open-source toolkit for controlled source audio-frequency magnetotelluric (CSAMT)

[![Documentation Status](https://readthedocs.org/projects/pycsamt/badge/?version=latest)](https://pycsamt.readthedocs.io/en/latest/?badge=latest) [![Build Status](https://travis-ci.com/WEgeophysics/pyCSAMT.svg?branch=master)](https://travis-ci.com/WEgeophysics/pyCSAMT) [![Requirements Status](https://requires.io/github/WEgeophysics/pyCSAMT/requirements.svg?branch=master)](https://requires.io/github/WEgeophysics/pyCSAMT/requirements/?branch=master)
  ![GitHub](https://img.shields.io/github/license/WEgeophysics/pyCSAMT?color=blue&logo=GNU&logoColor=red) ![GitHub release (latest by date)](https://img.shields.io/github/v/release/WEgeophysics/pyCSAMT?color=orange) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5674430.svg)](https://doi.org/10.5281/zenodo.5674430)
  [![PyPI version](https://badge.fury.io/py/pycsamt.svg)](https://badge.fury.io/py/pycsamt)
  ![GitHub repo size](https://img.shields.io/github/repo-size/WEgeophysics/pycsamt?color=0A4CEE&style=flat-square)

## Overview 

* **Definition**

    CSAMT is a frequency-domain electromagnetic method well-established  as a resistivity exploration tool in deep geological structure detection. The method is broadly applied in  diverse of exploration problems such as mineral, hydrocarbon,  groundwater resources, as well as mapping the fault-zones etc. 

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

* [Installation Guide](https://pycsamt.readthedocs.io/en/latest/installation.html?highlight=installation)
* [Demo](https://pycsamt.readthedocs.io/en/latest/demo.html?highlight=demo) 
* [API Documentation](https://pycsamt.readthedocs.io/en/latest/)
* [Home Page](https://github.com/WEgeophysics/pyCSAMT/wiki)
* [User Guide](https://github.com/WEgeophysics/pyCSAMT/blob/develop/docs/pyCSAMT%20User%20Guide.pdf)

## Licence 

pyCSAMT is under GNU Lesser GPL version3 [LGPLv3](https://github.com/03-Daniel/pyCSAMT/blob/master/LICENSE.md).


## Credits

We use or link some third-party software (besides the usual tool stack: numpy, scipy, matplotlib) and are grateful for all the work made by the authors of these awesome open-source tools:
* [MTpy](https://github.com/MTgeophysics/mtpy.git)
* [Occam2D](https://marineemlab.ucsd.edu/Projects/Occam/index.html)
* Zonge Engineering softwares:
    - [AMTAVG](http://www.zonge.com/legacy/DatPro.html/)
    - [ASTATIC](http://www.zonge.com/legacy/PDF_DatPro/Astatic.pdf)

## System requirements 
* Python 3.7+ 

## Citations 

If you use pyCSAMT in any published work, consider citing the paper below or the [zenodo DOI](https://zenodo.org/record/5674430/export/hx#.YmpoXo9BxGo) as a contribution.

> *Kouadio, K.L., Liu, R., Mi, B., Liu, C., 2022. pyCSAMT: An alternative Python toolbox for groundwater exploration using controlled source audio-frequency magnetotelluric. J. Appl. Geophys. 201, 104647. https://doi.org/10.1016/j.jappgeo.2022.104647.*
> 
> *Kouadio, K.L., 2021. pyCSAMT: A Python open-source toolkit for controlled source audio-frequency magnetotelluric. https://doi.org/10.5281/zenodo.5674430.*

## Contributors
  
1. Key Laboratory of Geoscience Big Data and Deep Resource of Zhejiang Province, School of Earth Sciences, [Zhejiang University](http://www.zju.edu.cn/english/), China.

2. Department of Geophysics, School of Geosciences and Info-physics, [Central South University](http://www.zju.edu.cn/english/), China.

3. Equipe de Recherche Géophysique Appliquée, Laboratoire de Geologie Ressources Minerales et Energetiques, UFR des Sciences de la Terre et des Ressources Minières, [Université Félix Houphouët-Boigny]( https://www.univ-fhb.edu.ci/index.php/ufr-strm/), Cote d'Ivoire.

* Developer: 1, 3- Kouadio K. Laurent; <etanoyau@gmail.com>, <kkouao@zju.edu.cn>,
* Contributors:
    *  2- Rong LIU; <liurongkaoyang@126.com>
    *  1- Binbin MI; <mibinbin1991@126.com>  
