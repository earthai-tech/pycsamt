# -*- coding: utf-8 -*-
"""
Created on Sat Nov 27 16:04:29 2021
    < python penetration1d -d=data/edi -f 1024 3000 8000 -o=portrait>

"""

import os 
import sys 
import argparse 
from pycsamt.viewer.plot import Plot1d 

PROG = os.path.basename (__file__).replace('.py', '')

def main(): 
    parser =argparse.ArgumentParser(
        prog= PROG, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter ,
        description= "Visualize  the 1D skin depth at selected frequencies "
        ) 
    parser.add_argument('-d', '--data-fn', '--data',
                        type =str,
                        help ='Path to [EDI|AVG|J] files', 
                        dest='fn', 
                        required=True,
                        )
    
    parser.add_argument('-f', '--frequencies', '--select-frequency',
                        nargs='+',
                        type =float,
                        help ='List of frequencies to visualize', 
                        dest='selected_frequency', 
                        required=True,
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
    
    parser.add_argument('-o', '--orientation', '--orient-fig',
                        dest ='orientation', 
                        type =str , 
                        default ='landscape',
                        help ='Orientation of the savefigure'
                        )
    parser.add_argument('-s', '--savefig', 
                        dest ='savefigure',
                        type =str , 
                        help = 'Save figure'
                        ) 
    parser.add_argument('--rot', '--rotate-stations',
                        type =int , 
                        dest ='rotate_stations',
                        default = 90,
                        help = 'Rotate the station labels'
                        )

    return parser.parse_args()

if __name__== '__main__':
    args = main()

    sys.stdout.write(Plot1d().penetration1D (
        fn =args.fn,
        profile_fn= args.profile_stn, 
        selected_frequency =args.selected_frequency, 
        rename_station =args.rename_station, 
        rotate_station_names = args.rotate_stations , 
        orientation=args.orientation, 
        savefigure =args.savefigure
        )
                    
    )
