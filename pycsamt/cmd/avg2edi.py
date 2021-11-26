# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 16:08:01 2021
 
    < python pycsamt avg2edi -d='data/avg/K1.AVG -p=data/avg/K1.stn -s=data/K1_edi>

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
        description= "Output SEG-EDI from Zonge International Engineering AVG file",
        ) 
    parser.add_argument('-d', '--data-fn', '--data',
                        type =argparse.FileType('r'),
                        help ='*.avg format collected from GDP-II hardware.', 
                        dest='data_fn', 
                        required=True,
                        )
    
    parser.add_argument('-p', '--profile-fn', '--profile', 
                        type =argparse.FileType('r'),
                        help ='Station location file', 
                        dest = 'profile_fn',
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
    sys.stdout.write( CSAMT().avg2edi(
                               data_fn =args.data_fn,
                               profile_fn =args.profile_fn, 
                               savepath =args.savepath)
        )