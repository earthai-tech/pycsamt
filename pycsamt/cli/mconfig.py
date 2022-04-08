# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 19:31:39 2022

@author: @Daniel03
"""
import os 
import sys 
import json 
import yaml
import argparse 

prog = os.path.basename (__file__).replace('.py', '')

cmd =['<$ pycsamt mconfig --show --json >', '\n', '|',
      '< $ python pycsamt/cli/mconfig.py --show --yml >'
] 


def display_help_config(path_to_config_file:str ,
                        ctype:str ='yml',
                        cprefix:str='e.g.data',
                        show:bool =True): 
    """ Display the content as the help for the model config file.
    :param confif_file: Path-Like object, 
        can be *.yml of *.json file 
    
    """
    cfile = f'{cprefix}.json' if ctype =='json' else f'{cprefix}.yml'
    # If the config file is not installed, ignore the error
    data , isfile = None , False

    if os.path.isdir (path_to_config_file): 
        cfile = os.path.join(path_to_config_file, cfile)
        isfile = os.path.isfile(cfile) 
    if isfile: 
        if show: 
            with open (cfile , 'r', encoding='utf8') as f: 
                data = ''.join( f.readlines())
        if not show:
            if cfile.endswith('.yml'): 
                with open(cfile) as fy: 
                    data =  yaml.load(fy, Loader=yaml.SafeLoader) 
            else : 
                with open(cfile) as fj: 
                    data =  json.load(fj)  
            data =str(data)
            
    if data is None: 
        data ='Config files (*e.g.data.YMl/*e.g.data.JSON) not found!'
        
    return data 


def main ():
    
    parser =argparse.ArgumentParser(
        prog= prog, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter ,
        description ='Show example of configuration file content for NM construction.', 
        #usage =pycsamt.poof_cli_usage (pycsamt.nm.__doc__), 
        epilog = ''.join(cmd),
        allow_abbrev=False, 
        ) 
    
    parser.add_argument('--show', 
                   dest ='show',
                   action ='store_true', 
                   help =''.join([
                       'Visualize the content of configure files. Indeed,', 
                   ' the config file grabs all the `pycsamt.geodrill.geocore.GeoStratigraphy` arguments into a single file.'
                   ]
            )
             )
    mygroup = parser.add_mutually_exclusive_group()
    
    parser.add_argument('--yml', 
                   action ='store_true', 
                   help ='stdout of the content of YAML file.'
                   )

    mygroup.add_argument('--json', 
                   action ='store_true', 
                   help ='stdout of the content of JSON file.'
                   )
    args = parser.parse_args()
    
    ftype= 'json' if args.json else 'yml' 
    
    sys.stdout.write( display_help_config(path_to_config_file='pycsamt/_mdata',
                        ctype = ftype, show =args.show)
                      )
    
if __name__=='__main__': 
    main() 

    
    
    
    
    
    
    
    
    
    
    
        