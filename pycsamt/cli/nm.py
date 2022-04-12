# -*- coding: utf-8 -*-
#       Author: Kouadio K.Laurent<etanoyau@gmail.com>
#       Licence: LGPL
#       Created on Thu Nov 25 13:38:29 2021
"""
    Compute the stratigraphic model (NM) from the forward modeling CRM
    
    <$ nm --data =*.file --model=*.model --iter=*.iter --data=*.dat> 
  
  where *. is the filename. In addition, rather than passing all the occam2D 
  files individually as arguments, use --crmf option instead. for instance:
      `crm_files` is a dict valuewhich contain all occam2d model files like:
          
          crm_files =dict(
                    mesh_fn = 'data/occam2d/Occam2DMesh',
                    iter_fn ='data/occam2d/ITER17.iter', 
                    model_fn ='data/occam2d/Occam2DModel', 
                    data_fn ='data/occam2d/OccamDataFile.dat')
          with TRES= [10, 70, 100, 1000, 3000]  and 
          LN= ['river water', 'fracture zone', 'MWG', 'LWG', igneous rocks] 
               
   The CLI should be written as e.g.:
       
  <$ pycsamt nm --crmf=crm_files --tres=TRES --ln=LN --beta=4 --ptol=0.2 > 
  
  <$ python pycsamt/cli/nm.py -c pycsamt/metadata/e.g.data.json --show --misfit > 
  
"""

import os 
# import sys 
import argparse 
from pycsamt.geodrill.geocore import GeoStratigraphy 

cmd = [
    'EXAMPLE COMMANDS: <$ nm -c pycsamt/metadata/e.g.data.yml --show >', ' | ',
    '<$ python pycsamt/cli/nm.py -c pycsamt/metadata/e.g.data.json --show --misfit >'
]

prog = os.path.basename (__file__).replace('.py', '')
hptol = ''. join(['likelihood error. Probabilistic error ', 
                  '  between the TRES collected from boreholes or wells', 
                  ' and the CRM from inversion.']) 
def main(): 
    
    parser =argparse.ArgumentParser(
        prog= prog, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter ,
        description ='Generate strata model (NM) from the forward modeling CRM.', 
        #usage =pycsamt.poof_cli_usage (pycsamt.nm.__doc__), 
        epilog = ''.join(cmd),
        allow_abbrev=False, 
        ) 
    my_group = parser.add_mutually_exclusive_group(required=False)
    parser.add_argument('-c', '--config', 
                    dest ='config', 
                    help =''.join([
                    'CONFIG file is the suitable approach used to wrap all the ',
                    'arguments into a single file rather than providing them individualy.',
                    'It could be *.YML or *.JSON files.', 
                    'For futher details about the arrangement of CONFIG file, type', 
                    ' the following command < $ pycsamt mconfig --show -yml > or',
                    ' < $ python pycsamt/cli/mconfig.py --show -yml >'])
                    )
    parser.add_argument('--show', 
                   dest ='show',
                   action ='store_true', 
                   help ='Visualize the strata model'
                   ) 
    
    parser.add_argument('--kind', 
                        dest ='kind',
                        default ='nm', 
                        choices= ('nm', 'crm'),
                        help=''.join([
                            'Type of model to visualize. It could be the ', 
                            'forward model `crm` or the strata model `nm`.']
                            )
                        )
   
    parser.add_argument('--misfit', 
                        dest ='misfit',
                        action='store_true',
                        help=''.join([
                           ' visualise the error (misfit G) between the `nm` and `crm`',
                           ]
                        ))
    
    my_group.add_argument('-p', '--ptol', 
                        dest ='ptol',
                        type =float,
                        default =0.1,
                        help =hptol ) 
    my_group.add_argument('-b', '--beta', 
                        dest ='beta', 
                        type =int , 
                        default=5, 
                        help ='Number of the CRM block constructor.'
                        ) 
    my_group.add_argument('-n','--n-epochs', '--n-iter',
                        dest ='n_epochs',
                        default =100, 
                        type =int, 
                        help ='Number of iteration for gradient'\
                            ' descent to converge.'
                        ) 
    my_group.add_argument('--crm',
                          dest ='crm',
                          help =''.join([
                            'Calculated resistivity model(CRM) from forward computation.',
                            'Load *.NPY/NPZ data from Numpy (https://numpy.org/doc/stable/reference/generated/numpy.load.html).', 
                            'As a reminder, CRM is composed of numpy 2D array ', 
                            '(i.e. X -horizontal nodes and Z-vertical nodes). ',
                            'If given, user does not need to supply the Occam2d ', 
                            'files (*.dat, *.iter., *.model, *.mesh).']
                            )
                                        ) 
    my_group.add_argument('--df', '--data-fn','--data', 
                        type =argparse.FileType('r'),
                        dest ='data_fn',
                        help ='Occam2d data file')
    my_group.add_argument('--iter', '--iter-fn', '-if',
                        type =argparse.FileType('r'),
                        dest='iter_fn', 
                        help ='Occam2d iteration file')
    my_group.add_argument('--mf', '--model-fn', '--model',
                        type =argparse.FileType('r'),
                        help ='Occam2d model file', 
                        dest='model_fn')
    my_group.add_argument('--mef', '--mesh-fn','--mesh', 
                        type =argparse.FileType('r'),
                        help ='Occam2d mesh file', 
                        dest = 'mesh_fn',
                        )

    my_group.add_argument('--ln','--layer-names','--LN',
                        dest ='input_layers',
                        action ='append', 
                        nargs ='+', 
                        help ='Geological rocks or layers names '\
                            'collected in the exploration area'
                        ) 
    my_group.add_argument('--tres','--input-resistivities','--TRES', 
                        dest ='input_layers',
                        action ='append', 
                        nargs ='+', 
                        help ='Electrical rock properties collected'\
                            ' in the exploration area'
                        ) 
    my_group.add_argument('--build', 
                    dest ='build', 
                    help ="Option to trigger the NM construction", 
                    action ='store_false',
                    )
   

    
    parser.add_argument('-v', '--verbosity', 
                       dest ='verbose',
                       default =0, 
                       action='count',
                       help= 'Verbosity: control the level of output messages',
                       )
    
    args = parser.parse_args()
    
    if args.config is not None:
        gObj = GeoStratigraphy.geoArgumentsParser(
            config_file = args.config)
    else:
        gObj = GeoStratigraphy(
                                crm = args.crm ,
                                data_fn = args.data_fn , 
                                mesh_fn = args.mesh_fn , 
                                iter_fn = args.iter_fn , 
                                input_layers = args.input_layers, 
                                input_resistivities = args.input_resistivities,
                                beta= args.beta , 
                                ptol =args.ptol , 
                                n_epochs = args.n_epochs , 
                                verbose = args.verbose , 
                                build =args.build, 

                                )
    if args.show:
        _= gObj.strataModel(kind=args.kind,
                                misfit_G =args.misfit)
    
if __name__== '__main__':
    main()


