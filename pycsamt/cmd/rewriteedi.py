# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 17:21:37 2021


    < python pycsamt rewriteedi -d data/edi/new_csa000.edi -n=testedi --dtype=emap --verbose=4  >

@author: @Daniel03
"""

import os 
import sys 
import argparse 
from pycsamt.ff.core.edi import Edi 

PROG = os.path.basename (__file__).replace('.py', '')

def main(): 
    parser =argparse.ArgumentParser(
        prog= PROG, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter ,
        description= "Rewrite EDI files. Call `~.edi.Edi_collection` "\
            "to rewrite multiples files. "
        ) 

    parser.add_argument('-d','-edi', '--edi-fn', '--data-fn','--edi-filename',
                        type =argparse.FileType('r'),
                        dest='edi_fn', 
                        required=True,
                        help ='Path to EDI file to rewrite' 
                            )
    parser.add_argument('-n', '--name', '--output-filename','new-edi-name', 
                            dest ='new_edifilename', 
                            type =str , 
                            help ='Output filename used as EDI-prefix.'
                            )
    parser.add_argument('-dtype', '--data-type', '--datatype',
                        dest ='datatype', 
                        choices =('mt', 'emap'),
                        help ='Type of EDI file "mt" or "emap"'
                        )
    
    parser.add_argument('-s', '--savepath', 
                        dest ='savepath',
                        type =str , 
                        help = 'Save rewritten EDI file',
                        ) 
    parser.add_argument('-v', '--verbose', 
                        dest ='verbose',
                        default =0, 
                        action='count',
                        help = 'Control the level of verbosity. Higher more messages',
                        )

    return parser.parse_args()

if __name__== '__main__':
    args = main()
    sys.stdout.write(Edi().write_edifile(edi_fn= args.edi_fn, 
                                 savepath =args.savepath,
                          datatype=args.datatype,
                          new_edifilename =args.new_edi_name,
                          verbose = args.verbose
                          )
        )
