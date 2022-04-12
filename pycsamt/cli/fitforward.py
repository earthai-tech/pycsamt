# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 11:08:23 2022

    CLI plot the Fitting resistivity curves 
    <$ pycsamt forward  --help >
    
@author: @Daniel03
"""

import os 
import argparse 
from pycsamt.modeling.occam2d import plotResponse 

PROG = os.path.basename (__file__).replace('.py', '')

def main(): 
    parser =argparse.ArgumentParser(
        prog= PROG, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter ,
        description= "Fitting forward inversion curves.",
        epilog =' | '.join(cmd)
        ) 
    #my_group = parser.add_mutually_exclusive_group(required=True)
        
    parser.add_argument('path',
                        type =str,
                        help ='Path to Occam2d forward modeling data.'\
                            ' Use responde(*.resp) and data (*.dat)files.', 
                        )


    parser.add_argument('-s', '--stations', 
                    dest ='stations',
                    nargs ='*', 
                    type =str, 
                    help = 'Station to visualize it filtting curve at a corresponding line.'\
                        " Can be ['S00', 'S04', 's08', 'S12'] for instance. Note that each station"\
                            'corespond to each individual line. Looking for the examples, '\
                                'the 04 stations fits 04 survey lines.'
                    ) 

    parser.add_argument('-e', '--error-type',
                        dest ='error',
                        choices =range(1,5), 
                        default =4, 
                        help =''.join([
                            'Errors type from raw and forward response. For instance '
                            '- `1` : visualize the misfit compute manually '
                            '- `2`:  visualize the raw error from raw occam data ' 
                            '- `3`: visualzie the error data and phase  defined as ' 
                                    '* error = (input data - forward data)/ RESI'
                            '- `4` or `residual`: Visualize only  the residual data.'
                            ' Typically default is `residual`. Get more info in '
                            ':doc:`~pycsamt.modeling.occam2d.plotResponse`']
                            )
                        )
    parser.add_argument('--rms', '--root-mean-squared', 
                    dest ='rms',
                    nargs ='*', 
                    help = 'RMS of each survey line. If RMS is given, their values must fits'\
                        ' the number of given stations to label the legend.'
                        " Can be ['1.013', '1.451', '1.00', '1.069'] for instance."
                    ) 


    args= parser.parse_args()

    plotResponse(data_fn =args.path,
                    stations = args.stations, # ['S00', 'S04'],# 's08', 'S12'],  # sites to visualize 
                     rms =args.rms, #['1.013', '1.451'],# '1.00', '1.069'], # rms of each line
                      error_type =args.error )
    
    
cmd=['EXAMPLE CMDs:<$ pycsamt fitforward --stations S01 S02 --rms 1.12 1.23 >', 
     '<$ python pycsamt/cli/fitforward.py data/inversionFiles -s s05 S07 --rms=1.14 >', 
     '<$ python pycsamt/cli/pseudocrossresistivityandphase.py data/edi --phic 45 --style=pcolormesh >']

if __name__== '__main__':
    main()
