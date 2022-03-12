# -*- coding: utf-8 -*-
import os 
import sys 

if __package__ is None or __name__ == '__main__': 
    sys.path.append( os.path.dirname(os.path.dirname(__file__)))
    sys.path.insert(0, os.path.dirname(__file__))
    __package__= 'pycsamt'
    
from pycsamt.cli import (
    nm,
    pseudostratigraphic,
    cedi,
    j2edi,
    avg2edi,
    rewriteedi,
    penetration1d, 
    penetration2d, 
    rms, 
    misfit2d, 
    
)

from pycsamt.gui import (
    oc2d_bdin,
    wa2edi,
    wa2j,
    corrected_edi,
    p_moc2d,
    p2log,
    oas2f,
    gs2f,
    wi2d,
    cmake_dh ,

)



