# -*- coding: utf-8 -*-
#       Author: Kouadio K.Laurent<etanoyau@gmail.com>
#       Licence: LGPL
"""
Created on Thu Nov 25 17:21:37 2021
    Rewrite Electromagnetic  Array Profiling (EMAP) and Magnetotelluric (MT)

    < pycsamt rewriteedi -d data/edi/new_csa000.edi -n=testedi --dtype=emap --verbose=4  >

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
            "to rewrite multiples files in SEG EMAP or MT  formats.", 
        allow_abbrev=False, 
        epilog =' | '.join(cmd)
        ) 

    parser.add_argument('data',
                        type =str,
                        help ='Path to EDI files to rewrite' 
                            )
    parser.add_argument('-n', '--name', '--output-filename','--new-edis', 
                        dest ='new_edis', 
                        type =str , 
                        help ='Output filename used as EDI-prefix.'
                        )
    parser.add_argument('-dtype', '--data-type', '--datatype',
                        dest ='datatype', 
                        choices =('mt', 'emap'),
                        help ='Type of EDI file "emap" or "mt". `emap` for'
                        ' Electromagnetic Array Profiling and `mt` for Magnetotelluric '
                        )
    
    parser.add_argument('-s', '--savepath', 
                        dest ='savepath',
                        type =str , 
                        help = 'Save rewritten EDI file',
                        ) 
    parser.add_argument('-v', '--verbose', 
                        dest ='verbose',
                        type =int, 
                        default =0, 
                        help = 'Control the level of verbosity.',
                        )

    args= parser.parse_args()
    if os.path.isfile (args.data): # if edi is a single file
        args.data = os.path.dirname (args.data) # get the directory
    if os.path.isdir(args.data): 
        args.data =[os.path.join(args.data, f) 
                    for f in os.listdir(args.data) if f.find('.edi')>=0]
        if len(args.data) ==0: 
            sys.stdout.write(
                'No edi-data found in the given path. '
                'Please provide the right path containing at least one edi file.')
        else: 
            for f in args.data : 
                Edi().write_edifile(edi_fn= f, 
                                    savepath =args.savepath,
                                      datatype=args.datatype.lower(),
                                      new_edifilename =args.new_edis,
                                      verbose = args.verbose
                                      )

cmd = ['EXAMPLE COMMANDS: < $ pycsamt rewriteedi  data/edi --verbose 2 >', 
       '< $ python pycsamt/cli/rewriteedi.py data/edi/new_csa000.edi -n=testedi -dtype=emap --verbose=4 >', 
    '< $ python pycsamt/cli/rewriteedi.py data/edi -dtype=emap  >', 
]

if __name__== '__main__':
    main()

