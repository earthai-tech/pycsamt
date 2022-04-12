# -*- coding: utf-8 -*-
#   Created on Thu Nov 25 15:26:19 2021

"""
        Plot the pseudostratigraphic from  NM
        
     <$ pycsamt pseudostratigraphic -s=S00 -z=25% > 
     <$ python pycsamt/cli/pseudostratigraphic.py --station=S00 --zoom=25%  > 
     

 `zoom` parameter can be a list or a float number. Please refer to
    pycsamt.geodrill.geocore.GeoStratigraphy.plotPseudostratigraphic.__doc__
    
"""
import os 
import sys 
import argparse 
# import pycsamt 
from pycsamt.geodrill.geocore import GeoStratigraphy 

PROG = os.path.basename (__file__).replace('.py', '')

def main(): 
    parser =argparse.ArgumentParser(
        prog= PROG, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter ,
        description=  'Plot the pseudostratigraphic from  NM.',
        epilog = ' | '.join(cmd),
        fromfile_prefix_chars='@',
        #usage = pycsamt.poof_cli_usage (pycsamt.pseudostratigraphic.__doc__),
        ) 
    #mygroup = parser.add_mutually_exclusive_group(required =True)
    
    parser.add_argument('-z', '--zoom', 
                    nargs ='+', 
                    dest ='zoom',
                    help ='Visualize a specific part of the pseudostratigraphic log'\
                        'Can be a list of `[top, bottom]` i.e it accepts only two values for fitting.'\
                            'Given values more than 1 are gathered into a list.')
    parser.add_argument('-s', '--station',
                        help = 'Station|Site ID. `S00` is known as the first station name.'+\
                            ' See the documentation at `GeoStratigraphy.plotPseudostratigraphic.__doc__`'\
                                ' for a deep implementation as the use of Matplotlib hatches customizing.',
                        default ='S00', 
                        ) 
    
    parser.add_argument( '--lf', '--fontsize', '--label-fontsize', 
                        dest ='fontsize',
                        type =float, 
                        default =12.,
                        help =' Label fontsize.'
                        )
    

  
    args = parser.parse_args()

    GeoStratigraphy.plotPseudostratigraphic(
        station =args.station, zoom =args.zoom, 
        annotate_kws ={'fontsize': args.fontsize} )
    
    #return parser.parse_args()
cmd = ['EXAMPLE COMMANDS: < $ pycsamt pseudostratigraphic --station=S10 --zoom=25% >', 
    '< $ python pycsamt/cli/pseudostratigraphic.py -s=s17 -z=10 120 >', 
]

if __name__=='__main__': 
    main()
    
# sys.stdout.write( GeoStratigraphy.plotPseudostratigraphic(
#     station =args.station, zoom =args.zoom, 
#     annotate_kws =args.annotate_kws )
#                   )