# -*- coding: utf-8 -*-
"""
Created on Sat Nov 27 17:03:19 2021

    < pycsamt misfit2d -p=data/occm2d --resp= K1.resp --data=K1.dat >
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
        description= " Visualize the Occam2D Response Misfit",
        ) 
    parser.add_argument('-p', '--path', '--occam-path',
                        type =str,
                        help ='Path to Occam2D inversion files', 
                        dest='path', 
                        required=True,
                        )
    
    parser.add_argument('-d', '--data', '--data-fn',
                        type =str,
                        help ='Specify only the data filename e.g. *.dat ', 
                        dest='data_fn', 
                        required=True,
                        )
    parser.add_argument('-r', '--resp', '--response_fn',
                        type =float,
                        help ='Specify only the response filename e.g. *.resp ', 
                        dest='resp_fn', 
                        default =1., 
                        required=True,
                        )
    parser.add_argument('-k','--kind',
                        dest='kind',
                        default ='rho', 
                        choices =('res', 'rho', 'phase'),
                        help ='Specificy the kind of data to plot[RHO|PHASE]', 
                        )
    parser.add_argument('--tl','-target-line',
                        help ='Show the larget line', 
                        dest='show_target_line', 
                        action ='store_true',
                        )
    
    return parser.parse_args()

if __name__== '__main__':
    args = main()
    sys.stdout.write(getMisfit(data_fn = os.path.join(args.path,args.data_fn ),
                  resp_fn = os.path.join(args.path,args.resp_fn ),
                  kind=args.kind)
                        )
    