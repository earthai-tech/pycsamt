# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0


# This file is pycsamt/zonge/components.py ( the last configuration the structure)

"""
One-stop public façade that re-exports all high-level containers
(metadata, measurements, QC helpers, tensor stacks…).
"""

from .heads import Head # inherit from AVGComponentBase
from .properties import ( # all this module , already implement get and set 
    Hardware,
    SurveyAnnotation,
    SurveyConfiguration,
    Receiver,
    Transmitter,
    SkipFlag,
)
from .meas import ( #  they inherit from AVGComponentBase 
    CompMeas,
    Amps,
    Frequency,
)
from .resphase import Resistivity, Phase # inherit from VariationBase which also inherit from AVGComponentBase
from .survey import Station # inherit from AVGComponentBase
from .var import (
    PcEmag, # inherit from from VariationBase which also inherit from AVGComponentBas
    PcHmag,
    PcRho,
    SPhz,
    SHphz,
    SEphz,
)
from .z import Z # inherit from AVGComponentBase
from .info import DataInfo  # do not inherit from nothing . 

__all__ = [
    # header & metadata
    "Head",
    "Hardware",
    "SurveyAnnotation",
    "SurveyConfiguration",
    "Receiver",
    "Transmitter",
    "SkipFlag",

    # measurements / site-level info
    "CompMeas",
    "Amps",
    "Frequency",
    "Station",

    # data-quality / variation metrics
    "PcEmag",
    "PcHmag",
    "PcRho",
    "SPhz",
    "SHphz",
    "SEphz",

    # apparent-field estimates
    "Resistivity",
    "Phase",

    # full impedance tensor stack
    "Z",
    
    # Info 
    "DataInfo"
]
# Dont write all the components for now. I give you to have an Idea How I will structure 
# the package pycsamt/zonge/ directory. 

# Now we want to newly reogarganize, how you typycally know , each part of 
# avg file is a class object . 
# why we want in this version  , is first write the pycsamt/zonge/base.py 
# all object should inherit from AVGComponentBase rewrite to be flexible to 
# all object . 

# each object should implement a class method from_avg ( which accept the avg file,) 
# and return the specific  component informations get from the avg file. 
# the transfer class can be used to exchange between legacy avg and modern avg 

# the __init__ of each object accept the data ( which is information of the avg file 
# supposed already passed or mannually given ..)
# each class object , except should implement read ( which read the information
# collected from avg file and store as attributes of this components class )
# the method write , get the attributes stored and reconstruct the write string format 
# of the avg. This way willhelp to reconstruct AVG file later easily using each 
# object components . 

#   To achieve it , let smart revise the previous class I passed , which is not 
# intuitive and not respond to the better development guide. 
# The AVGProperties , The AVGComponentBase. 
# revise them , revise these base class more professionnaly to better fit the 
# the new structure we want to implement . 

# currently focus on the new base module after. You can also propose the best architecture to  
# achieve this professionnaly . 
#
 
