# -*- coding: utf-8 -*-
"""
Tip to intepolate 
Created on Mon Aug  1 19:13:45 2022
Intrepolate EDI and write a new files. 

Raw Data are sometimes little messy with messing frequency values 
Best idea is to interpolate the frequency to putt all EDIs on to 
the same frequency range. 

@author: Daniel

"""

import numpy as np
from pycsamt.ff.core.edi import  Edi_collection# , Edi 
from pycsamt.utils.func_utils import get_interpolate_freqs 

edi_fn = r"C:\Users\Daniel\Desktop\Data\AMT\E1\edi_r"

savepath = r'C:\Users\Daniel\Desktop\Data\AMT\E1\edi_i'
prefix_new_edi = 'i'
# multiple interpolation 
# you can get the interpolate frequencies in the collection of edifiles 
cObjs = Edi_collection(list_of_edifiles= edi_fn)
ifreqs, nfreq = get_interpolate_freqs(cObjs.ediObjs, to_log10=True) 
new_freq = np.logspace (*ifreqs, nfreq) 
for obj in cObjs.ediObjs: 
    new_Z = obj.interpolateZ(new_freq_array=new_freq)
    obj.write_new_edifile( 
        new_Z=new_Z, savepath = savepath, 
        new_edi_fn = prefix_new_edi
        
        ) 
