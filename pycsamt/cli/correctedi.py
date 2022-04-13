# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 16:27:24 2021

    < $ correctedi -d data/edi -f flma --ndipole 150 --rf 1024  >
    
                        or using the script : 
                            
    <$ python pycsamt/cli/correctedi.py --help > 

"""

import os 
import sys 
import argparse 
from pycsamt.utils.func_utils import smart_format 
from pycsamt.ff.processing import Processing, corr

PROG = os.path.basename (__file__).replace('.py', '')

def main(): 
    parser =argparse.ArgumentParser(
        prog= PROG, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter ,
        description= "Correction and remove static"\
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
                        choices=list(corr.TAGS.keys()),
                        help ='Filter type to correct EDI files It could be one of' 
                        f' the following {smart_format(list(corr.TAGS.keys()))} filters.'  
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
    
    parser.add_argument('-s', '--savepath', 
                    dest ='savepath',
                    type =str , 
                    help = 'Path to save the corrected edifiles.'
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
                        dest ='number_of_points', 
                        type =int , 
                        default =5,
                        help ='Number of FLMA dipole length', 
                        metavar='NUMBER_OF_DIPOLES'
                        )
    #-------------------------------------------------------------------------

    args= parser.parse_args()
        
    if args.data.lower().endswith('.avg') and args.profile_fn is None: 
        sys.stdout.write ('correctedi: error: AVG file is detected.'
                          'Please provide the Zonge station file `*.stn`')
    else:
        Processing().correct_edi(
                                data_fn =args.data,
                                FILTER=args.FILTER, 
                                filename =args.filename,
                                reference_frequency= args.reference_frequency,
                                savepath =args.savepath, 
                                number_of_points= args.number_of_points, 
                                dipole_length= args.dipole_length, 
                                number_of_skin_depth= args.number_of_skin_depth,
                                )
cmd=['EXAMPLE CMDs: <$ correctedi data/edi -ft tma --npoint 3 >', 
     '<$ python pycsamt/cli/correctedi.py data/edi -ft flma --ndipole 7 >', 
     '<$ python pycsamt/cli/correctedi.py data/edi --filter=ama --nskin=3 >']

if __name__== '__main__':
    main()

