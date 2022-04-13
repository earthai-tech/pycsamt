# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 10:17:58 2022
       
     CLI to plot pseudo-cross-section of apparent
        resistivity and phase. Deals with[AVG|J|EDI] file.
        
        <$ pseudocrossresistivityandphase data/edi --help >
        
        
@author: @Daniel03
"""

import os 
import sys 
import argparse 
from pycsamt.viewer.plot import Plot2d 

PROG = os.path.basename (__file__).replace('.py', '')
mpl_ls ='https://matplotlib.org/stable/gallery/lines_bars_and_markers/linestyles.html'
def main(): 
    parser =argparse.ArgumentParser(
        prog= PROG, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter ,
        description= "Plot Pseudo-cross resistivity and phase",
        epilog =' | '.join(cmd)
        ) 
    #my_group = parser.add_mutually_exclusive_group(required=True)
        
    parser.add_argument('data',
                        type =str,
                        help ='Path to [EDI|AVG|J] files', 
                        )

    parser.add_argument('-p', '--profile-fn',
                    type =str,
                    dest='profile_fn', 
                    help ='Path to Zonge Engineering station location.'\
                        ' Set the profile if the given data is (AVG) format.'
                        )

    parser.add_argument('-s', '--savefigure', 
                    dest ='savefig',
                    type =str , 
                    help = 'Path to save the corrected edifiles.'
                    ) 


    parser.add_argument('--rhoc', '--res-contour', '--resistivity-contour',
                        dest ='contourrho',
                        nargs='*', 
                        type=float,
                        #action ='append', 
                        help =''.join([
                            'Resistivity contour for visualization'
                            'Note that the resitivity values are not on' 
                            ' logarithm scale (ohm m ).Thus, for multiple contour'
                            ' delineation  put value are gathering onto list'
                            '(e.g: [500, 7000] ohm.m)']
                            )
                        )
    parser.add_argument('--phic', '--phase-contour',
                        dest ='contourphase', 
                        nargs='*',
                        #action='append', 
                        type=float,
                        help =''.join([
                            'Likewise the resistivity contour for visualization'
                            'phase can be can be 45 degree or else and gathering'
                            ' onto a list']
                            )
                        )
    parser.add_argument('--cline', '--contour-line',
                        dest ='contourline', 
                        default ='-',
                        help =f'Matplotlib linestyle for contour. Refer to {mpl_ls!r}'\
                            'for further details.'
                        )
    parser.add_argument('--style', '--plot_style', '--nflma',
                        dest ='ps', 
                        default ='imshow',
                        choices =('pcolormesh', 'imshow'),
                        help ='plot style. Should be [pcolormesh |imshow]', 

                        )

    args= parser.parse_args()
    
    if args.data.lower().endswith('.avg') and args.profile_fn is None: 
        sys.stdout.write ('correctedi: error: AVG file is detected.'
                          'Please provide the Zonge station file `*.stn`')
    else:
        Plot2d().pseudocrossResPhase(
            fn=args.data, 
            profile_fn=args.profile_fn, 
            delineate_resistivity=args.contourrho,
            delineate_phase=args.contourphase,
            plot_style =args.ps, 
            savefig = args.savefig, 
            contour_lines_style=args.contourline 
                                        )
cmd=['EXAMPLE CMDs: <$ pseudocrossresistivityandphase data/edi --help >', 
     '<$ python pycsamt/cli/pseudocrossresistivityandphase.py data/edi >', 
     '<$ python pycsamt/cli/pseudocrossresistivityandphase.py data/edi --rhoc 500 7000>', 
     '<$ python pycsamt/cli/pseudocrossresistivityandphase.py data/edi --phic 45 --style=pcolormesh >']

if __name__== '__main__':
    main()


