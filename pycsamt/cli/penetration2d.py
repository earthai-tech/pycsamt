# -*- coding: utf-8 -*-
"""
Created on Sat Nov 27 16:25:59 2021
    < python pycsamt -d data/j --style pcolormesh>
"""

import os 
import sys 
import argparse 
from pycsamt.viewer.plot import Plot2d 

PROG = os.path.basename (__file__).replace('.py', '')

def main(): 
    parser =argparse.ArgumentParser(
        prog= PROG, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter ,
        description= "Visualize the 2D map of the frequency skin depths "
        )
    parser.add_argument('-d', '--data-fn', '--data',
                        type =str,
                        help ='Path to [EDI|AVG|J] files', 
                        dest='fn', 
                        required=True,
                        )

    parser.add_argument('-doi','--depth', '--image-depth', 
                        dest ='doi', 
                        type =float, 
                        default =1000,
                        help ='Max depth to visualize'
                        )
    parser.add_argument('-p', '--profile-fn', '--profile',
                      type =argparse.FileType('r'),
                      dest='profile_fn', 
                      required=True,
                      help ='Path to Zonge Engineering station location.'\
                          ' Set the profile only the data is avg format.'
                          )
    parser.add_argument('--snames', '--station-names', 
                        dest ='rename_station', 
                        type =list, 
                        help ='List of new stations names'
                        )

    parser.add_argument('--style', '--plot-style', 
                    dest ='plot_style',
                    default ='imshow',
                    choices =('pcolormesh', 'imshow'),
                    help = 'style of plot'
                    ) 
    parser.add_argument('-s', '--savefig', 
                        dest ='savefigure',
                        type =str , 
                        help = 'Save figure'
                        ) 


    return parser.parse_args()



if __name__== '__main__':
    args = main()
    sys.stdout.write(Plot2d().penetration2D(fn = args.data_fn, 
                              profile_fn= args.profile_stn ,
                              plot_style=args.plot_style,  
                            doi=args.doi, 
                            savefig=args.savefigure 
                            )
    )
