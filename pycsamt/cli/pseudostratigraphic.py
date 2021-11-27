# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 15:26:19 2021

     <python pycsamt ps -s=S00 -z=25% --annotate={'fontsize':12}>

 `zoom` parameter can be a list of single float number. Please refer to
    pycsamt.geodrill.geocore.GeoStratigraphy.plotPseudostratigraphic.__doc__
"""
import os 
import sys 
import argparse 
from pycsamt.geodrill.geocore import GeoStratigraphy 

PROG = os.path.basename (__file__).replace('.py', '')

def main(): 
    parser =argparse.ArgumentParser(
        prog= PROG, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter ,
        description= "Plot the pseudostratigraphic from  NM.",
        epilog = GeoStratigraphy.plotPseudostratigraphic.__doc__,
        ) 
    
    parser.add_argument('-s', '--station', 
                        dest ='station',
                        type =str , 
                        help = 'Station!Site id to fetch the log. e.g.S00'
                        ) 
    parser.add_argument('-a', '--annotate', '--akws', 
                        dest ='annotate_kws',
                        type =dict , 
                        default ={'fontsize':12},
                        help =' Matplotlib log layout properties'
                        )
    
    parser.add_argument('-z', '--zoom', 
                        nargs ='?',
                        type ='+',
                        dest ='zoom',
                        help ='Visualize a specific part of the'\
                            ' pseudostratigraphic log. ')

    return parser.parse_args()

if __name__== '__main__':
    args = main()

    sys.stdout.write( GeoStratigraphy.plotPseudostratigraphic(
        station =args.station, zoom =args.zoom, 
        annotate_kws =args.annotate_kws )
                      )