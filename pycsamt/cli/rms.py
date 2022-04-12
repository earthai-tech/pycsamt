# -*- coding: utf-8 -*-
"""
Created on Sat Nov 27 16:40:47 2021
    < rms  data/occam2d/LogFile.logfile >
    
"""
import os 
 
import argparse 
from pycsamt.viewer.plot import Plot1d
PROG = os.path.basename (__file__).replace('.py', '')

def main (): 
    parser =argparse.ArgumentParser(
        prog= PROG, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter ,
        description= " Visualize the R.M.S of Occam2D FE forward modeling",
        epilog= ' | '.join(cmd)
        ) 
    parser.add_argument('data',
                        type =str,
                        help ='Path to Occam2D log files.', 
                        )
    parser.add_argument('-t', '--target', '--rms-target',
                        type =float,
                        help ='Target value of root-mean-squared (RMSE)', 
                        dest='target', 
                        default =1., 
                        )
    parser.add_argument('-s', '--savefig', 
                        dest ='savefigure',
                        type =str , 
                        help = 'Save figure'
                        ) 
    parser.add_argument('--grid','--show-grid',
                        help ='Show the grid (Plot)', 
                        dest='show_grid', 
                        action ='store_true',
                        )
    parser.add_argument('--tline','--show-target-line',
                        help ='Show the larget line', 
                        dest='show_target_line', 
                        action ='store_true',
                        )
    

    
    args = parser.parse_args()
    Plot1d().plotRMS(fn =args.data, 
                        target=args.target, 
                        show_grid =args.show_grid,
                        show_target_line = args.show_target_line, 
                        savefig = args.savefigure 
                        )
    
cmd = ['EXAMPLE COMMANDS: < $ rms  data/occam2d/logFile.logfile --show-target >', 
       '< $ python pycsamt/cli/rms.py  data/occam2d/logFile.logfile --show-target  >', 
       '< $ python pycsamt/cli/rms.py  data/occam2d/logFile.logfile --grid --show-target >', 
]
if __name__== '__main__':
    main()
