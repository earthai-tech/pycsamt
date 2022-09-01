# -*- coding: utf-8 -*-
"""
Created on Sat Nov 27 16:25:59 2021
    < $ penetration2d data/j --style pcolormesh>
"""

import os 
import sys 
import argparse 
from pycsamt.view.plot import Plot2d 

PROG = os.path.basename (__file__).replace('.py', '')

def main(): 
    parser =argparse.ArgumentParser(
        prog= PROG, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter ,
        description= "Visualize the 2D map of the frequency skin depths ",
        allow_abbrev=False,
        epilog = ' | '.join(cmd)
        )
    parser.add_argument('data',
                        type =str,
                        help ='Path to [EDI|AVG|J] files', 
                        )
    
    parser.add_argument('-d','--depth', '--doi', '--image-depth', 
                        dest ='doi', 
                        type =float, 
                        default =1000,
                        help ='Max depth to visualize'
                        )
    parser.add_argument('-p', '--profile-fn', '--profile',
                      type =str,
                      dest='profile_fn', 
                      help ='Path to Zonge Engineering station location.'\
                          ' Set the profile if the data format (*.avg) is given.'
                          )
        
    parser.add_argument('-c', '--cmap', 
                    dest ='cm',
                    default ='twilight',
                    help = 'Matplotlib Color map' 
                    ) 
    
    parser.add_argument('-s', '--savefig', 
                        dest ='savefigure',
                        type =str , 
                        help = 'Save figure'
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



    args= parser.parse_args()
    
    if args.data.lower().endswith('.avg') and args.profile_fn is None: 
        sys.stdout.write ('penetration2d: error: AVG file is detected.'
                          'Please provide the Zonge station file `*.stn`')
    else : 
        Plot2d().penetration2D(fn = args.data, 
                                  profile_fn= args.profile_fn ,
                                  plot_style=args.plot_style,  
                                doi=args.doi, 
                                savefig=args.savefigure , 
                                cm = args.cm, 
                                )
cmd = ['EXAMPLE COMMANDS: < $ penetration2d data/edi --doi 1000 >', 
    '< $ python pycsamt/cli/penetration2d.py data/avg/K1.avg  -p=data/avg/K1.stn >', 
    '< $ python pycsamt/cli/penetration2d.py data/edi -d=1000 --cmap=jet_r >',
]

if __name__== '__main__':
    main()

