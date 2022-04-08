# -*- coding: utf-8 -*-
"""
Created on Sat Nov 27 17:03:19 2021

    <$ pycsamt misfit2d data/inversionFiles/K1.dat data/inversionFiles/K1.resp >
    
    <$ python pycsamt/cli/misfit2d.py data/inversionFiles/K1.dat 
    | data/inversionFiles/K1.resp --kind=phase >
"""
import os 
import sys 
import argparse 
from pycsamt.modeling.occam2d import getMisfit 

PROG = os.path.basename (__file__).replace('.py', '')

def main (): 
    parser =argparse.ArgumentParser(
        prog= PROG, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter ,
        # fromfile_prefix_chars='@',
        description= " Visualize the Occam2D Response Misfit",
        ) 

    parser.add_argument('data',
                        help ='Add the path to occam2d data file (e.g. *.dat)', 
                  
                        )
    parser.add_argument('response', 
                        help ='Add the path to occam2d response file (e.g. *.resp)', 
                 
                        )
    
    parser.add_argument('-k','--kind','--visualization-type',
                        dest='kind',
                        default ='rho', 
                        choices =('res', 'rho', 'phase'),
                        help ='Specificy the kind of data to plot [RHO|PHASE]', 
                        )

    args= parser.parse_args()
    getMisfit(data_fn =  args.data ,
                  resp_fn = args.response ,
                   kind=args.kind
                  )
                        

if __name__== '__main__':
    main()
    
    