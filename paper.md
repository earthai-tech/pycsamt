---
title: "The pyCSAMT package for groundwater exploration using controlled source audio-frequency magnetotelluric data "
tags:
	- Python
	- pyCSAMT
	- groundwater 
	- controlled source audio-frequency magnetotelluric 
authors:
	- name: Kouao Laurent, Kouadio
	  orcid: 0000-0001-7259-7254
	  affiliation: 1,3
	- name: Rong  Liu
	  affiliation: 2
	- name: Albert Okrah Malory
	  affiliation: 1
	- name: Bingin Mi
	  affiliation: 1
	- name: Chun-ming Liu
	  affiliation: 2
affiliations:
	- name: Key Laboratory of Geoscience Big Data and Deep Resource of Zhejiang Province, School of Earth Sciences, Zhejiang University, China.
	  index: 1
	- name: Department of Geophysics, School of Geosciences and Info-physics, Central South University, China.
	  index: 2
date: 14 March 2022
bibliography: paper.bib
---


# Introduction 

Controlled source audio-frequency magnetotelluric (CSAMT) is a frequency-domain electromagnetic method established as a 
good resistivity exploration tool for mapping the fault zones for groundwater exploration [@Asch2011; @Bernardetal1997; @Bernard1990; @Chouteau2008; @Kouadioetal2020; @Liuetal2020]. 
However, the detection of fracture zone requires additional geophysical methods to supplement the CSAMT [@Guoetal2019; @Wadi2017; @ZongeandHughes1991].
 This is expensive and despite this combination, the misinterpretation of inversion results leads to unsuccessful drillings due to the wrong location of the borehole [@Kouadioetal2020].
 We, therefore introduce pyCSAMT software to solve this problem. First, the software is used to predict the strata log at each station (pseudostratigraphic log (PS))
 by demarcating the fracture zone, and secondly estimating the layer thicknesses with less margin error useful before the drilling operations using the geological and borehole/well data collected in the survey area.


#  Key Functionalities 

pyCSAMT follows the modular approach of existing software like MTpy [@Krieger2014] and GMT [@Wessel1998], and contains an inner handler to calibrate and to scale the raw data
from different hardware into the appropriate units (SI) which recomputes the deviation errors before analysis and processing[@Mykle1996]. The software includes some electromagnetic array
profiling filters such as the trimming moving average, the fixed dipole-length moving average, and the adaptive-moving-average filter
based on the idea of [@Torres-verdìn1992] to correct the CSAMT data corrupted by the static shift effect [@Raymond1993; @Sandersetal2006; @Sandersetal1996].

The toolbox reads different CSAMT raw data formats ( e.g., \*.AVG format [@Mykle1996; @Sandersetal2006] from Zonge Engineering), \**.*DAT format proposed by [@Chave1994] and
the standard Electronic Data Interchange (EDI) file format). It is composed of three main packages with different roles: *ff*, *geodrill,* and *viewer*. The *ff* package encompasses the *core* and the
*processing* (a set of *analysis* and *processing* modules) sub*-*packages. Figure 1 shows an overview of pyCSAMT packages and sub-packages with their roles.

![pyCSAMT packages structures and the keys modules. The colors in the workflow diagram represent which parts of the software are used in each step. For example, the modules in the geodrill packages are usedfor NM construction and PS prediction ](paper_figures/pycsamt_workflow_and_packages.png){width="7.260416666666667in" height="5.303517060367454in"}


The core sub-package contains functionality to read and write CSAMT data from industry-standard formats such as AVG, DAT, and EDI including metadata from the header of the EDI file, the location, and also the impedance tensor (Z).

The processing sub-package is designed to facilitate working with \*.DAT, \*AVG and EDI data and generating inputs for existing third-party processing codes, for example, the Z codes of [@Kirkbyetal2019; @Krieger2014].

The viewer package is essentially dedicated to data and log sequences visualization (1D and 2D plots).

The modeling package of the toolbox uses the finite-element (FE) structured grid and deals with the OCCAM2D [@DeGroot-Hedlin1990] to invert the processed data. It uses the FE algorithm
developed by [@Wannamakeretal1987] to generate the OCCAM2D input and output data for the model visualization (e.g., Figure 2a). Moreover, it also provides some output files for other external modeling software like
oasis montaj of Geosoft [@GeosoftCorporation2021], and surfer of Golder Software corporations [@GoldenSoftware2021].

The geodrill package mainly deals with geological, borehole and/or well data collected from the survey area.  It also includes a geological database composed of rock properties such as the electrical properties and the minerals classification 
of [@Slichter1942] and [@Palacky1988] for new model construction (NM) from which the PS under each station is retrieved.  Thus from the forward modeling results, NM is constructed (e.g., Figure 2b) and estimated model error (e.g., Figure 2c). The  PS technique developed in  [@kouadioetal2022] 
 is used for log forecasting (e.g., Figure 2d) and thickness estimation display retrieved from the software memory (e.g., Figure 2e). Overall, the packages provide features coded in Python classes, methods, and functions.

![An example of software functionalities. a) 2D forward modeling results (error floors set at 10% apparent resistivities and 20% phase with starting model set at 150 Ω.m) with the broad delimitation of four different geological structures (S1 to S4). b) NM model generated
including the geological and boreholes data. c) The error map between the 2D inversion model and the NM. It could indicate the layer misclassification with larger errors. d) An example of PS built at station S00 at 125 m depth. e) The display layer thicknesses estimated at station S00 fetched from the memory storage thin 1km depth.](paper_figures/NMpseudostratigraphic_log_at_125_m_depth.png){width="7.490972222222222in" height="6.0in"}


# Illustrative Examples

To illustrate the use of the software, the real CSAMT data collected from the survey carried out in the Xingning area, Hunan province, China, are available on [CSAMT data in Xingning area, Hunan Province, China \|Zenodo](https://zenodo.org/record/5533467#.YVK6mnzithF). Originally,
data was acquired using the GDP-32 II hardware of Zonge International Engineering. It was recorded in AVG format (*\*.avg*) with the station location file (*\*.stn*). Moreover, a practical example has been
validated using the data of tested boreholes (Bo) and is implemented in [@kouadioetal2022] with examples scripts and workshop material. The error thickness evaluated between the predicted log and the mechanical boreholes (Bo) was less than 06 meters.

Furthermore, another published full case history paper demonstrates the key functionality contained in pyCSAMT and is also available in [@kouadioetal2022]. It can ascertain the significance scope of the software
and how the outputs of software for 3D visualization allow estimating the water reservoir thickness.

# Availability

The software is distributed under a Lesser General Public License v3.0 and is available from [WEgeophysics/pyCSAMT: Python for Controlled Source Audio-frequency Magnetotellurics (CSAMT)(github.com)](https://github.com/WEgeophysics/pyCSAMT).

# Acknowledgments 

We thank the Subsurface Imaging and Sensing Laboratory of the School of Earth Sciences of Zhejiang University and Central South University for supervising the project. The development of this Python toolbox has been
supported by the Key Laboratory of Geoscience Big Data and Deep Resource of Zhejiang Province of School of Earth Sciences at the National Zhejiang University, China.

# References 


