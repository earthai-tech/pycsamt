# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 09:21:55 2022
    CLI  to plot static correction . Deal with [EDI|J|AVG]files.

@author: @Daniel03
"""


import os 
import sys 
import argparse 
from pycsamt.utils.func_utils import smart_format 
from pycsamt.view.plot import Plot1d

PROG = os.path.basename (__file__).replace('.py', '')

def main(): 
    parser =argparse.ArgumentParser(
        prog= PROG, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter ,
        description= "Plot the corrected resistivity after removing the"\
            " shift effect from [EDI|AVG|J] files",
        epilog =' | '.join(cmd)
        ) 
    #my_group = parser.add_mutually_exclusive_group(required=True)
        
    parser.add_argument('data',
                        type =str,
                        help ='Path to [EDI|AVG|J] files', 
                        )
    parser.add_argument('-ft', '--filter', 
                        type =str,
                        dest='FILTER',
                        default ='ama',
                        choices=('ama', 'tma', 'flma'),
                        help ='Filter type to correct EDI files It could be one of' 
                        f" the following {smart_format(('ama', 'tma', 'flma'))} filters."  
                        )
    parser.add_argument('-p', '--profile-fn',
                    type =str,
                    dest='profile_fn', 

                    help ='Path to Zonge Engineering station location.'\
                        ' Set the profile if the given data is (AVG) format.'
                        )

    parser.add_argument('-f','-rf', '--refreq', '--reference-frequency',
                        type =float,
                        help ='Reference frequency in Hertz.', 
                        dest='reference_frequency', 
                        
                        )
    
    parser.add_argument('-s', '--savefigure', 
                    dest ='savefig',
                    type =str , 
                    help = 'Path to save figure.'
                    ) 
    
    parser.add_argument('-n', '--name', '--output-filename', 
                        dest ='filename', 
                        type =str , 
                        help ='Output filename used as EDI-prefix.'
                        )

    parser.add_argument('-dl', '--dipole-length',
                        dest ='dipole_length', 
                        type =int , 
                        default =50,
                        help ='length of fixed dipole for Hanning window.'
                        )
    #-------------------------------------------------------------------------
    parser.add_argument('--npoint', '--number-of-points', '--ntma',
                        dest ='number_of_points', 
                        type =int , 
                        default =1,
                        help ='Number of TMA filter points', 
                        metavar='NUMBER_OF_POINTS'
                        )
    
    parser.add_argument('--nskin', '--number-of-skin-depth', '--nama',
                        dest ='number_of_skin_depth', 
                        type =int, 
                        default =3,
                        help ='Number of AMA skin depths', 
                        metavar='NUMBER_OF_SKIN_DEPTHS'
                        )
    parser.add_argument('--ndipole', '--number-of-dipoles', '--nflma',
                        dest ='number_of_dipoles', 
                        type =int , 
                        default =5,
                        help ='Number of FLMA dipole length', 
                        metavar='NUMBER_OF_DIPOLES'
                        )
    parser.add_argument('--fbtw', '--fill-between', 
                        dest ='fill_between', 
                        action='store_true',
                        help ='Fill between the uncorrected and corrected resistivity.', 
                        )
    parser.add_argument('--cfbtw', '--color-fill-between', 
                        dest ='colorfillbetween', 
                        default ='thistle', 
                        help ='Color to customize the fill between the uncorrected and corrected resistivity.', 
                        )
    parser.add_argument('--ctma', '--color-tma-filter', 
                        dest ='colortma', 
                        default ='blue', 
                        help ='Customize the color for corrected TMA resistivity.', 
                        )
    parser.add_argument('--cflma', '--color-flma-filter', 
                        dest ='colorflma', 
                        default ='aqua', 
                        help ='Customize the color for corrected FLMA resistivity.', 
                        )
    
    parser.add_argument('--mk', '--marker', 
                        dest ='marker', 
                        default ='x', 
                        help ='Matplotlib marker. See https://matplotlib.org/stable/api/markers_api.html', 
                        )
    parser.add_argument('--ms', '--markersize', 
                        dest ='ms',
                        type =float, 
                        default =3., 
                        help ='Matplotlib markersize. See https://matplotlib.org/stable/api/markers_api.html', 
                        )

    parser.add_argument('--cmkedge', '--markeredgecolor', 
                        dest ='markeredgecolor', 
                        default ='k', 
                        help ='Matplotlib markersize. See https://matplotlib.org/stable/api/markers_api.html', 
                        )
    
    parser.add_argument('--ls', '--line-style', 
                        dest ='ls', 
                        default ='-', 
                        help ='Matplotlib markersize. See https://matplotlib.org/stable/gallery/lines_bars_and_markers/linestyles.html', 
                        )

    parser.add_argument('--lw', '--line-width', 
                        dest ='lw',
                        type=float,
                        default =1.5, 
                        help ='Matplotlib line style. See https://matplotlib.org/stable/gallery/lines_bars_and_markers/linestyles.html', 
                        )
    parser.add_argument('--fs', '--font-size', 
                        dest ='fs', 
                        type =float, 
                        default =2., 
                        help ='Matplotlib fontsize. See https://www.geeksforgeeks.org/change-font-size-in-matplotlib/', 
                        )
    
    #-------------------------------------------------------------------------


    args= parser.parse_args()
        
    if args.data.lower().endswith('.avg') and args.profile_fn is None: 
        sys.stdout.write ('correctedi: error: AVG file is detected.'
                          'Please provide the Zonge station file `*.stn`')
    else:
        
        if args.number_of_dipoles is not None: 
            args.number_of_points= args.number_of_dipoles 
            
        Plot1d(
            mstyle= args.marker,ms=args.ms,fs=args.fs, 
            lw=args.lw, 
            markeredgecolor= args.markeredgecolor).plot_static_correction(
            data_fn =args.data , 
            profile_fn=args.profile_fn, 
            frequency_id= args.reference_frequency,
            number_of_points=args.number_of_points, 
            fill_between =args.fill_between,
            ADD_FILTER =args.FILTER, 
            dipole_length=args.dipole_length,
            fill_between_color= args.colorfillbetween, 
            tma_color=args.colortma, 
            flma_color= args.colorflma,
            number_of_skin_depth = args.number_of_skin_depth, 
            savefig =args.savefig
            )

cmd=['EXAMPLE CMDs: <$ staticshift data/edi -ft tma --npoint 3 >', 
     '<$ python pycsamt/cli/staticshift.py data/edi --help >', 
     '<$ python pycsamt/cli/staticshift.py data/edi -ft flma --mk D --cflma=blue --ndipole=5 >']

if __name__== '__main__':
    main()

