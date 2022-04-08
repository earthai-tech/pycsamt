# -*- coding: utf-8 -*-
#   Created on Thu Nov 25 15:26:19 2021

"""
        Plot the pseudostratigraphic from  NM
        
     <python pycsamt pseudostratigraphic -s=S00 -z=25% >

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
        epilog = GeoStratigraphy.plotPseudostratigraphic.__doc__,
        fromfile_prefix_chars='@',
        #usage = pycsamt.poof_cli_usage (pycsamt.pseudostratigraphic.__doc__),
        ) 
    
    parser.add_argument('-s', '--station', 
                        dest ='station',
                        help = 'Station|Site id to fetch the log. e.g.S00'
                        ) 
    parser.add_argument('-a', '--annotate', '--akws', 
                        dest ='annotate_kws',
                        type =dict , 
                        default ={'fontsize':12},
                        help =' Matplotlib log layout properties'
                        )
    
    parser.add_argument('-z', '--zoom', 
                        nargs ='*',
                        #action ='append', 
                        dest ='zoom',
                        help ='Visualize a specific part of the'\
                            ' pseudostratigraphic log. ')
        
        
    args = parser.parse_args()

    GeoStratigraphy.plotPseudostratigraphic(
        station =args.station, zoom =args.zoom, 
        annotate_kws =args.annotate_kws)
    
    #return parser.parse_args()

if __name__=='__main__': 
    main()
    
# sys.stdout.write( GeoStratigraphy.plotPseudostratigraphic(
#     station =args.station, zoom =args.zoom, 
#     annotate_kws =args.annotate_kws )
#                   )