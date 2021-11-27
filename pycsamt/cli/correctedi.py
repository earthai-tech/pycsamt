# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 16:27:24 2021

    < python pycsamt cedi -d data/edi -f flma --ndipole 150 --rf 1024  >

"""

import os 
import sys 
import argparse 
from pycsamt.ff.processing import Processing, corr

PROG = os.path.basename (__file__).replace('.py', '')

def main(): 
    parser =argparse.ArgumentParser(
        prog= PROG, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter ,
        description= "Correction and remove static"\
            " shift effect from [EDI|AVG|J] files",
        ) 
    parser.add_argument('-d', '--data-fn', '--data',
                        type =str,
                        help ='Path to [EDI|AVG|J] files', 
                        dest='data_fn', 
                        required=True,
                        )
    parser.add_argument('-f', '--filter', '--FILTER',
                        type =str,
                        dest='FILTER',
                        default ='ama',
                        choices=list(corr.TAGS.keys()),
                        help ='Filter to correct EDI files' 
                        )
    parser.add_argument('--rf', '--refreq', '--reference-frequency',
                        type =float,
                        help ='Reference frequency in Hertz.', 
                        dest='reference_frequency', 
                        required=True,
                        )
    parser.add_argument('--name', '--output-filename', 
                        dest ='filename', 
                        type =str , 
                        help ='Output filename used as EDI-prefix.'
                        )
    parser.add_argument('-p', '--profile-fn', '--profile',
                        type =argparse.FileType('r'),
                        dest='profile_fn', 
                        required=True,
                        help ='Path to Zonge Engineering station location.'\
                            ' Set the profile only the data is avg format.'
                            )

    parser.add_argument('--npoint', '--number-of-points', '--ntma',
                        dest ='number_of_points', 
                        type =int , 
                        default =1,
                        help ='Number of TMA filter points'
                        )
    
    parser.add_argument('--dl', '--dipole-length', '--ndipole',
                        dest ='dipole_length', 
                        type =int , 
                        default =50,
                        help ='length of fixed dipole for Hanning window'
                        )
    parser.add_argument('--nskin', '--number-of-skin-depth', '--nama',
                        dest ='number_of_skin_depth', 
                        type =int , 
                        default =3,
                        help ='Number of AMA skin depths'
                        )
    
    parser.add_argument('-s', '--savepath', 
                        dest ='savepath',
                        type =str , 
                        help = 'Save output edifiles corrected'
                        ) 
    parser.add_argument('--kws', 
                        type =dict , 
                        dest ='ckws',
                        default = dict(),
                        help = 'Additional keywords arguments. Please refer'\
                        ' to <~.Processing.correct_edi.__doc__> for more details.'
                        )

    return parser.parse_args()

if __name__== '__main__':
    args = main()
    sys.stdout.write( Processing().correct_edi(
                                data_fn =args.data_fn,
                                FILTER=args.FILTER, 
                                filename =args.filename,
                                reference_frequency= args.reference_frequency,
                                savepath =args.savepath, 
                                number_of_points= args.number_of_points, 
                                dipole_length= args.dipole_length, 
                                number_of_skin_depth= args.number_of_skin_depth,
                                **args.ckws
                                )
        )

