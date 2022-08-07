# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 10:23:01 2022

@author: Daniel


The factors ss_x and ss_y are in resistivity scale, so the
entries of  the matrix "S" need to be given by their
square-roots!
       
>>> from pycsamt.ff.processing import Processing 
>>> edipath = '/Users/Desktop/ediout/'
>>> Processing.noiseRemoval(edi_fn =edipath)
            
"""

from pycsamt.ff.processing import Processing 

# edifile = r'C:\Users\Daniel\Desktop\Data\AMT\E1\edi_test\new_csa000.edi' #new_E1_1.edi'
# outputedi = 'rmss_E1_1.edi'

ss_x = .5 #correction factor for x component
ss_y = 1.2  # correction factor for y component
num_freq  = None #  for distorsion removal 
# ediobj = Edi(edifile )
# ediobj.remove_static_shift( output_edi=True,
#                            new_edi_fn=None)
kind ='ss' # can be 'ss' for staticshift effect or 'dt' for distortion removal 

edipath = r'C:\Users\Daniel\Desktop\Data\AMT\E1\edi_i'
savepath = r'C:\Users\Daniel\Desktop\Data\AMT\E1\edi_ss'

prefix_new_edi = 'ss'
Processing.noiseRemoval(edipath,
                        kind =kind ,
                        ss_x = ss_x ,
                        ss_y =ss_y ,
                        num_freq = num_freq, 
                        savepath = savepath,
                        new_edi_fn= prefix_new_edi 
                        )
