# -*- coding: utf-8 -*-
"""
Created on Sat Nov 27 16:04:29 2021
    < python penetration1d -d=data/edi -f 1024 3000 8000 -o=portrait>

"""

import os 
import sys 
import argparse 
from pycsamt.view.plot import Plot1d 

PROG = os.path.basename (__file__).replace('.py', '')

cmd = ['EXAMPLE COMMANDS: < $ penetration1d data/edi -f 1024 300 8000 -o=portrait --verbose >', 
    '< $ python pycsamt/cli/penetration1d.py data/avg/K1.avg -f=1024 -p=data/avg/K1.stn >', 
    '< $ python pycsamt/cli/penetration1d.py data/edi -f 1024 3000 8000 -o=portrait >',
]

msg ='AVG file is detected. Please provide the Zonge station file `*.stn`'

def main(): 
    parser =argparse.ArgumentParser(
        prog= PROG, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter ,
        description= "Visualize  the 1D skin depth at selected frequencies ", 
        epilog = ' | '.join(cmd), 
        allow_abbrev=False, 
        ) 
    
    parser.add_argument('data',
                        type =str,
                        help ='Path to [EDI|AVG|J] files', 
                        )
    
    parser.add_argument('-f', '--frequencies', '--select-frequency',
                        nargs='+',
                        type =float,
                        help ='List of frequencies to visualize', 
                        dest='selected_frequency', 
                        required=True,
                        )
    
    parser.add_argument('-p', '--profile-fn', '--profile',
                      type =str, #argparse.FileType('r'),
                      dest='profile_fn', 
                      help ='Path to Zonge Engineering station location.'\
                          ' Set the profile if the *.avg data format is given.'
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
    parser.add_argument('--rot', '--rotate-station',
                        type =int , 
                        dest ='rotate_station',
                        default = 90,
                        help = 'Rotate the station labels'
                        )
    
    parser.add_argument('-v', '--verbose',
                        dest ='verbose',
                        action='store_true',
                        help = 'Control the level of verbosity',
                        )
    
    args = parser.parse_args()
    
    

    if args.data.lower().endswith('.avg') and args.profile_fn is None: 
        sys.stdout.write ('penetration1d: error: AVG file is detected.'
                          'Please provide the Zonge station file `*.stn`')
    else : 
        
        Plot1d().penetration1D (
        fn =args.data,
        profile_fn= args.profile_fn, 
        selected_frequency =args.selected_frequency, 
        rename_station =args.rename_station, 
        rotate_station_names = args.rotate_station , 
        orientation=args.orientation, 
        savefigure =args.savefigure, 
        verbose = args.verbose 
        )
        

if __name__== '__main__':
    main()


                    
    
