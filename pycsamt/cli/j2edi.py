# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 16:21:30 2021

    < j2edi --data=data/j  -s=data/j_edi>
    < $ python pycsamt/cli/j2edi.py --data=data/j > 

"""

import os 
import sys 
import argparse 

# import pycsamt
from pycsamt.core import CSAMT
PROG = os.path.basename (__file__).replace('.py', '')

def main(): 
    parser =argparse.ArgumentParser(
        prog= PROG, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter ,
        description= "Output J-EDI from Alan G. Jones (1994)  J/DAT file",
        #usage =pycsamt.poof_cli_usage (pycsamt.j2edi.__doc__)
        epilog =cmd 
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
                        help = 'Save J collection destination files'
                        ) 
    
    args =parser.parse_args()
    CSAMT().j2edi(
                data_fn =args.data_fn,
                savepath =args.savepath
        )
    
cmd= """ < j2edi --data=data/j  -s=data/j_edi> | < $ python pycsamt/cli/j2edi.py --data=data/j > """

if __name__== '__main__':
    main()