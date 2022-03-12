# -*- coding: utf-8 -*-
#Created on Thu Nov 25 13:38:29 2021
"""

  <python pycsamt nm --data =*.file --model=*.model --iter=*.iter --data=*.dat> 
  
  where *. is the filename. In addition, rather than passes all the occam model
  files as arguments, use --crmf option instead. for instance:
      `crm_files` is a dict valuewhich contain all occam2d model files like:
          
          crm_files =dict(
                    mesh_fn = 'data/occam2d/Occam2DMesh',
                    iter_fn ='data/occam2d/ITER17.iter', 
                    model_fn ='data/occam2d/Occam2DModel', 
                    data_fn ='data/occam2d/OccamDataFile.dat')
          with TRES= [10, 70, 100, 1000, 3000]  and 
          LN= ['river water’, 'fracture zone’, 'MWG', ‘LWG'] 
               
   The CLI should be write as:
       
  <python pycsamt nm --crmf=crm_files --tres=TRES --ln=LN --beta=4 --ptol=0.2> 
  
"""
import os 
import sys 
import argparse 
from pycsamt.geodrill.geocore import GeoStratigraphy 

prog = os.path.basename (__file__).replace('.py', '')
hptol = ''. join(['likelihood error parameter. Probabilistic error ', 
                  '  between the TRES collected from boreholes or wells', 
                  ' and the CRM from inversion.']) 
def main(): 
    parser =argparse.ArgumentParser(
        prog= prog, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter ,
        description= "Compute the stratigraphic model NM from CRM of %(prog)s"
        ) 
    parser.add_argument('--crm', dest ='crm',
                        help = 'Calculated resistivity model  from Occam2d') 
    parser.add_argument('--odict', '--crmf', '--dict', '--oc2d',
                        dest ='oc2df',
                        type =dict , 
                        default =dict(),
                        help =''.join(['Dictionnary of occam model files.',
                                       ' Usefull rather than passing ',
                                       'individual occam2d model files.'])
                        )
    
    parser.add_argument('--df', '--data-fn','--data', 
                        type =argparse.FileType('r'),
                        dest ='data_fn',
                        help ='Occam data file')
    parser.add_argument('--iter', '--iter-fn', '-if',
                        type =argparse.FileType('r'),
                        dest='iter_fn', 
                        help ='Occam2d iteration file')
    parser.add_argument('--mf', '--model-fn', '--model',
                        type =argparse.FileType('r'),
                        help ='Occam2d model file', 
                        dest='model_fn')
    parser.add_argument('--mef', '--mesh-fn','--mesh', 
                        type =argparse.FileType('r'),
                        help ='Occam2d mesh file', 
                        dest = 'mesh_fn',
                        )
    
    parser.add_argument('-p', '--ptol', 
                        dest ='ptol',
                        type =float,
                        default =0.1,
                        help =hptol ) 
    parser.add_argument('-b', '--beta', 
                        dest ='beta', 
                        type =int , 
                        default=5, 
                        help ='Divided coefficient of the CRM block.'
                        ) 
    parser.add_argument('-n','--n-epochs', '--n-iter',
                        dest ='n_epochs',
                        default =100, 
                        type =int, 
                        help ='Number of iteration for gradient'\
                            ' descent to converge.'
                        ) 
    parser.add_argument('--ln','--layer-names','--LN',
                        dest ='input_layers',
                        required =True,
                        nargs ='+', 
                        help ='Geological rocks or layers names '\
                            'collected in the exploration area'
                        ) 
    parser.add_argument('--tres','--input-resistivities','--TRES', 
                        dest ='input_layers',
                        type =float,
                        required =True,
                        nargs ='+', 
                        help ='Electrical rock properties collected'\
                            ' in the exploration area'
                        ) 
    parser.add_argument('-v', '--verbose', 
                        dest ='verbose',
                        default =0, 
                        action='count',
                        help= 'verbosity-Control the level of out-messages')
    
    return parser.parse_args()
    
    
    
if __name__== '__main__':
    args = main()

    sys.stdout.write(GeoStratigraphy(crm = args.crm ,
                                      data_fn = args.data_fn , 
                                      mesh_fn = args.mesh_fn , 
                                      iter_fn = args.iter_fn , 
                                      input_layers = args.input_layers, 
                                      input_resistivities = args.input_resistivities,
                                      beta= args.beta , 
                                      ptol =args.ptol , 
                                      n_epochs = args.n_epochs , 
                                      verbose = args.verbose , 
                                      **args.oc2df
                                      )
                      )

