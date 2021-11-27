# -*- coding: utf-8 -*-
"""
Created on Sat Nov 27 16:40:47 2021
    < pycsamt rms -d data/occam2d/LogFile.logfile >
    
"""
import os 
import sys 
import argparse 
from pycsamt.viewer.plot import Plot1d
PROG = os.path.basename (__file__).replace('.py', '')

def main (): 
    parser =argparse.ArgumentParser(
        prog= PROG, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter ,
        description= " Visualize the R.M.S of Occam2D FE forward modeling",
        ) 
    parser.add_argument('-d', '--data-fn', '--data',
                        type =argparse.FileType('r'),
                        help ='Path to Occam2D log files', 
                        dest='data_fn', 
                        required=True,
                        )
    parser.add_argument('-t', '--target', '--rms-target',
                        type =float,
                        help ='Target value of root-mean-squared', 
                        dest='target', 
                        default =1., 
                        )
    parser.add_argument('--grid','-show-grid',
                        help ='Show the grid (Plot)', 
                        dest='show_grid', 
                        action ='store_false',
                        )
    parser.add_argument('--tl','-target-line',
                        help ='Show the larget line', 
                        dest='show_target_line', 
                        action ='store_true',
                        )
    
    parser.add_argument('-s', '--savefig', 
                        dest ='savefigure',
                        type =str , 
                        help = 'Save figure'
                        ) 
    
    return parser.parse_args()

if __name__== '__main__':
    args = main()
    sys.stdout.write(Plot1d().plotRMS(fn =args.data_fn, 
                        target=args.target, 
                        show_grid =args.show_grid,
                        show_target_line = args.show_target_line, 
                        savefig = args.savefigure 
                        )
    )