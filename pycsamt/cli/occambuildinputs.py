# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 11:59:05 2022
      Special class to build occam2d imput files with.
      see `pycsamt.modeling.occam2d.occam2d_write.buildingInputfiles.__doc__`
      
      <$ occambuildinputs --help > 
      
@author: @Daniel03
"""

import os 
import argparse 
from  pycsamt.modeling.occam2d import occam2d_write

PROG = os.path.basename (__file__).replace('.py', '')

def main(): 
    parser =argparse.ArgumentParser(
        prog= PROG, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter ,
        description= """Build occam2d model input files. 
        [*.dat, *.model, *.mesh, *.iter, and *.startup]. Refer to 
        `:doc: `pycsamt.modeling.occam2d.occam2d_write` for more details.""",
        epilog =cmd, 
        allow_abbrev=False, 
        ) 
  
    parser.add_argument('path',
                        type =str,
                        help ='Path to where EDI files are located.'
                        )
    
    parser.add_argument( '-nl', '--nlayers',  
                    dest ='n_layers',
                    type=int, 
                    default=31., 
                    help = 'Number of model layers.'
                    )
  
    parser.add_argument( '-z', '-edoi', '--z-target',
                    dest ='investigation_depth',
                    type=float, 
                    default=1100., 
                    help = """Expected investigation depth or the target
                    of the depth expected to image in meters."""
                    )

    parser.add_argument( '-zb', '--z-bottom', '--emax-depth', 
                    dest ='z_bottom',
                    type=float, 
                    default=5000., 
                    help = """
                    Exaggerate depth in meters. Note that the exaggeration
                    depth , must be enough as possible.
                    Around 5* time  the investigation depth (doi) """
                    )
     
    parser.add_argument( '-z1', '--thick-layer1',  
                    dest ='z1_layer',
                    type=float, 
                    default=5., 
                    help = """
                    The thickness of the first layer in meters. """
                    )
    
    parser.add_argument( '-r0', '-rho0', '--starting-res', 
                    dest ='resistivity_start',
                    type=float, 
                    default=1., 
                    help = """Starting the resistivity model in log10 ohm.meters """
                    )
    

    parser.add_argument( '-niter', '--iteration-numbers', 
                    dest ='iteration_to_run',
                    type=float, 
                    default=100., 
                    help = """Number of iterations to run for building model blocks"""
                    )
    

    parser.add_argument( '-mode', '--occam-mode', 
                    dest ='occam_mode',
                    choices= ('1', '2', '5', '6'), 
                    default='6', 
                    help = """Occam model modes.
                    1 or log_all       --->   Log resistivity of TE and TM plus Tipper(MT),
                    2 or log_te_tip    --->   Log resistivity of TE plus Tipper(MT),
                    5 or log_te        --->   Log resistivity of TE,
                    6 or log_tm        --->   Log resistivity of TM, 
                    Get more info in :doc:`~pycsamt.modeling.occam2d.Data` or 
                    refer to < http://marineemlab.ucsd.edu/Projects/Occam/sharp/index.html> 
                    If you want to implement for other modes for MT data, see <https://github.com/MTgeophysics/mtpy.git`
                    """
                    ) 
    
    parser.add_argument( '-eptm', '--error-phase-tm', 
                    dest ='TM_phase_error',
                    type=float, 
                    default=20., 
                    help = """TM mode error phase in percentages."""
                    )
    
    parser.add_argument( '-ertm', '--error-rho-tm', 
                    dest ='TM_res_error',
                    type=float, 
                    default=10., 
                    help = """TM mode error resistivity in percentages."""
                    )
    
    parser.add_argument( '-epte', '--error-phase-te', 
                    dest ='TE_phase_error',
                    type=float, 
                    default=20., 
                    help = """TE mode error phase in percentages."""
                    )
    
    parser.add_argument( '-erte', '--error-rho-te', 
                    dest ='TE_res_error',
                    type=float, 
                    default=10., 
                    help = """TE mode error resistivity in percentage"""
                    )

    parser.add_argument( '-cw', '--cell-width', 
                    dest ='cell_width',
                    type=float, 
                    default=5, 
                    help = """width of cells with in station area in meters"""
                    )
    parser.add_argument( '-xpad', '--padding-x', 
                    dest ='x_pad_multiplier',
                    type=float, 
                    default=1.7, 
                    help = """Horizontal node X pad multiplier. 
                    It controls size of padding"""
                    )
    parser.add_argument( '-bt','--brick-trigger',
                    dest ='trigger',
                    type=float, 
                    default=1.7, 
                    help = """Brick trigger for block construction. 
                    [ float ] multiplier to merge model blocks at 
                    depth.  A higher number increases the number of
                    model blocks at depth."""
                    )
    
    parser.add_argument( '--strike', '--geoelectric-strike', 
                    dest ='geoelectric_strike',
                    type=float, 
                    default=0., 
                    help = """Geolectrical strike.If not given, will set to 0."""
                    )

    parser.add_argument( '--name-outdata', 
                    dest ='data_basename',
                    default ='OccamDataFile',
                    help = 'Rename the occam2d output data file.'
                    ) 
    parser.add_argument( '--name-startup', 
                    dest ='startup_basename',
                    default='Startup', 
                    help = 'Rename the occam2d output startup file.'
                    ) 
    parser.add_argument( '--ifreq', '--interpolate-frequency', 
                    dest ='interpolate_frequency',
                    action='store_true', 
                    help = 'Interpolate the frequency from EDI-files. If '\
                        'interpolate frequency is set to true, be sure to bring'\
                            ' the limit of interpolatation in logspace frequency.'
                    ) 

    parser.add_argument( '--nfreq', '--frequency-number', 
                    dest ='number_of_frequency',
                    type=int, 
                    default=17, 
                    help = 'Number of frequency for interpolations.'
                    )

    parser.add_argument( '--ifreqrange', '--frequency-range', 
                    dest ='intp_freq_logspace',
                    nargs='+', 
                    type=float, 
                    help = """Interpolate the frequency range. If 
                        'interpolate frequency is set to true, user compolsorily needs to
                        ' provide the range of interpolation in logspace log10.
                        ' e.g. (-1,4,number_of_frequency )=(0.1, 10000, 17) where the last item '
                         'is number of frequency. """
                    ) 
    
    parser.add_argument( '--savepath', 
                    dest ='savepath',
                    help = 'Path to save the occam2d inputfiles'
                    ) 

    
    args= parser.parse_args()
    occam2d_write.buildingInputfiles(
        edi_fn =args.path,
        geoelectric_strike= args.geoelectric_strike,
        interpolate_freq= args.interpolate_frequency, 
        intp_freq_logspace =args.intp_freq_logspace, 
         iteration_to_run= args.iteration_to_run, 
         resistivity_start = args.resistivity_start, 
         startup_basename=args.startup_basename,
         res_tm_err= args.TM_res_error, 
         phase_tm_err= args.TM_phase_error , 
         occam_mode= args.occam_mode, 
         n_layers = args.n_layers, 
         cell_width = args.cell_width , 
         x_pad_multiplier = args.x_pad_multiplier, 
         trigger =args.trigger, 
         z_bottom =args.z_bottom , 
         z1_layer =args.z1_layer, 
         z_target = args.investigation_depth, 
         occamDataFile = args.data_basename, 
                            )
    
cmd="""
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Note: pyCSAMT uses `MTpy` for building model inputfiles. 
    So when you build your occam2d  building files on WINDOW, you will meet a 
    little bug into the `MTpy module` such ``FileNotFoundError::
        
        [Errno 2] No such file or directory: '/tmp/profile_angle0.png` 
        relate to `line 1392: plt.savefig('/tmp/profile_angle0.png')`

    Please comment this line of code in the MTpy module and run it again or 
    create a temp directory to hold the profile angle image.

<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

EXAMPLE COMMANDS:<$ occambuildinputs --help> | <$ occambuildinputs data/edi > |
      <$ python pycsamt/cli/occambuildinputs.py data/edi -mode=6 -niter 112 -cw=7 --nlayers=32 -z=1000 -zb=5000  --ifreq >', 
"""
    
    
if __name__== '__main__':
    main()
