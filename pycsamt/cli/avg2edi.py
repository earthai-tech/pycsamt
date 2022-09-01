# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 16:08:01 2021
 
    <$ avg2edi -d='data/avg/K1.AVG -p=data/avg/K1.stn -s=data/K1_edi>
    <$ python pycsamt/cli/avg2edi.py  -d='data/avg/K1.AVG -p=data/avg/K1.stn >
     
"""

import os  
import argparse 
from pycsamt.core import CSAMT

PROG = os.path.basename (__file__).replace('.py', '')

def main(): 
    parser =argparse.ArgumentParser(
        prog= PROG, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter ,
        description= "Output SEG-EDI from Zonge International Engineering file",
        epilog= cmd 
        ) 
    parser.add_argument('-d', '--data-fn', '--data',
                        # type =str #argparse.FileType('r'),
                        help ='*.avg format collected from GDP-II hardware.', 
                        dest='data_fn', 
                        required=True,
                        )
    
    parser.add_argument('-p', '--profile-fn', '--profile', 
                        # type =str #argparse.FileType('r'),
                        help ='Station location file', 
                        dest = 'profile_fn',
                        required=True,
                        )
    parser.add_argument('-s', '--savepath', 
                        dest ='savepath',
                        # type =str , 
                        help = 'Save destination (output) edifiles'
                        ) 

    args = parser.parse_args()
    CSAMT().avg2edi(
                    data_fn =args.data_fn,
                    profile_fn =args.profile_fn, 
                    savepath =args.savepath
                    
    )
    
cmd =""" EXAMPLE COMMANDS: <$ avg2edi -d='data/avg/K1.AVG -p=data/avg/K1.stn -s=data/K1_edi>
 | <$ python pycsamt/cli/avg2edi.py  -d='data/avg/K1.AVG -p=data/avg/K1.stn > """
 
if __name__== '__main__':
    main()
