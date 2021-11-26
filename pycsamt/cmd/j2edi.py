# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 16:21:30 2021

    < python pycsamt j2edi --data=data/j  -s=data/j_edi>

@author: @Daniel03
"""

import os 
import sys 
import argparse 
from pycsamt.ff.core import CSAMT

PROG = os.path.basename (__file__).replace('.py', '')

def main(): 
    parser =argparse.ArgumentParser(
        prog= PROG, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter ,
        description= "Output J-EDI from Alan G. Jones (1994)  J/DAT file",
        ) 
    parser.add_argument('-d', '--data-fn', '--data',
                        type =str,
                        help ='Path of *.j/dat files collected from JEDI hardware.', 
                        dest='data_fn', 
                        required=True,
                        )
    parser.add_argument('-s', '--savepath', 
                        dest ='savepath',
                        type =str , 
                        help = 'Save output edifiles'
                        ) 

    return parser.parse_args()

if __name__== '__main__':
    args = main()
    sys.stdout.write( CSAMT().j2edi(
                               data_fn =args.data_fn,
                               savepath =args.savepath)
        )