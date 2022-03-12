# -*- coding: utf-8 -*-
#       created on Fri Oct 29 12:31:01 2021
#       Author: Kouadio K.Laurent<etanoyau@gmail.com>
#       Licence: LGPL
import os
import json 
import yaml
import csv 
import warnings
import copy
import joblib 
import pickle
import datetime
import shutil 
from six.moves import urllib 
from typing import (
    Union,
    TypeVar,
    Any,
    Optional,
    Type,
    Tuple,
    Dict ,
    Sequence, 
    Iterable,
    List, 
)

import pycsamt.utils.func_utils as FU
from pycsamt.utils._csamtpylog import csamtpylog


T=TypeVar('T')
A=TypeVar('A')

V= List[Union[float, int]]

_logger=csamtpylog.get_csamtpy_logger(__name__)


def serialize_data(data:Union [V, T], 
                   filename:Optional[A]=None, 
                   force:bool =True, 
                    savepath:Optional[Union[A,str]] =None,
                    verbose:int =0): 
    """ Store a data into a binary file 
    
    :param data: Object
        Object to store into a binary file. 
    :param filename: str
        Name of file to serialize. If 'None', should create automatically. 
    :param savepath: str, PathLike object
         Directory to save file. If not exists should automaticallycreate.
    :param force: bool
        If ``True``, remove the old file if it exists, otherwise will 
        create a new incremenmted file.
    :param verbose: int, get more message.
    :return: dumped or serialized filename.
        
    :Example:
        
        >>> import numpy as np
        >>> import pycsamt.bases import serialize_data
        >>> data = np.arange(15)
        >>> file = serialize_data(data, filename=None,  force=True, 
        ...                          savepath =None, verbose =3)
        >>> file
    """
    def _cif(filename, force): 
        """ Control the file. If `force` is ``True`` then remove the old file, 
        Otherwise create a new file with datetime infos."""
        f = copy.deepcopy(filename)
        if force : 
            os.remove(filename)
            if verbose >2: print(f" File {os.path.basename(filename)!r} "
                      "has been removed. ")
            return None   
        else :
            # that change the name in the realpath 
            f= os.path.basename(f).replace('.pkl','') + \
                f'{datetime.datetime.now()}'.replace(':', '_')+'.pkl' 
            return f

    if filename is not None: 
        file_exist =  os.path.isfile(filename)
        if file_exist: 
            filename = _cif (filename, force)
    if filename is None: 
        filename ='__mymemoryfile.{}__'.format(datetime.datetime.now())
        filename =filename.replace(' ', '_').replace(':', '-')
    if not isinstance(filename, str): 
        raise TypeError(f"Filename needs to be a string not {type(filename)}")
    if filename.endswith('.pkl'): 
        filename = filename.replace('.pkl', '')
 
    _logger.info (
        f"Save data to {'memory' if filename.find('memo')>=0 else filename}.")    
    try : 
        joblib.dump(data, f'{filename}.pkl')
        filename +='.pkl'
        if verbose > 2:
            print(f'Data dumped in `{filename} using to `~.externals.joblib`!')
    except : 
        # Now try to pickle data Serializing data 
        with open(filename, 'wb') as wfile: 
            pickle.dump( data, wfile)
        if verbose >2:
            print( 'Data are well serialized using Python pickle module.`')
    # take the real path of the filename
    filename = os.path.realpath(filename)

    if savepath is  None:
        dirname ='_memory_'
        try : savepath = FU.sPath(dirname)
        except :
            # for consistency
            savepath = os.getcwd() 
    if savepath is not None: 
        try:
            shutil.move(filename, savepath)
        except :
            file = _cif (os.path.join(savepath,
                                      os.path.basename(filename)), force)
            if not force: 
                os.rename(filename, os.path.join(savepath, file) )
            if file is None: 
                #take the file  in current word 
                file = os.path.join(os.getcwd(), filename)
                shutil.move(filename, savepath)
            filename = os.path.join(savepath, file)
                
    if verbose > 0: 
            print(f"Data are well stored in {savepath!r} directory.")
            
    return os.path.join(savepath, filename) 
    
def load_serialized_data (filename:Union[A, str], verbose:int =0): 
    """ Load data from dumped file.
    :param filename: str or path-like object 
        Name of dumped data file.
    :return: Data reloaded from dumped file.

    :Example:
        
        >>> from pycsamt.bases import load_serialized_data
        >>> data = load_serialized_data(
        ...    filename = '_memory_/__mymemoryfile.2021-10-29_14-49-35.647295__.pkl', 
        ...    verbose =3)

    """
    if not isinstance(filename, str): 
        raise TypeError(f'filename should be a <str> not <{type(filename)}>')
        
    if not os.path.isfile(filename): 
        raise FileExistsError(f"File {filename!r} does not exist.")

    _filename = os.path.basename(filename)
    _logger.info(
        f"Loading data from {'memory' if _filename.find('memo')>=0 else _filename}.")
   
    data =None 
    try : 
        data= joblib.load(filename)
        if verbose >2:
            (f"Data from {_filename !r} are sucessfully"
             " reloaded using ~.externals.joblib`!")
    except : 
        if verbose >2:
            print(f"Nothing to reload. It's seems data from {_filename!r}" 
                      " are not dumped using ~external.joblib module!")
        
        with open(filename, 'rb') as tod: 
            data= pickle.load (tod)
            
        if verbose >2: print(f"Data from `{_filename!r} are well"
                      " deserialized using Python pickle module.`!")
        
    is_none = data is None
    if verbose > 0:
        if is_none :
            print("Unable to deserialize data. Please check your file.")
        else : print(f"Data from {_filename} have been sucessfully reloaded.")
    
    return data 

def parse_csv(csv_fn:A =None,
              data:Sequence [Dict[str, Type[Any]]] =None, 
              todo:Optional[str] ='reader', 
               fieldnames: Iterable[str]=None, 
               savepath:Optional[str] =None,
               header: bool=False, 
               verbose:int=0,
               **csvkws) -> T: 
    """ Parse comma separated file or collect data from CSV.
    
    :param csv_fn: csv filename,or output CSV name if `data` is 
        given and `todo` is set to ``write|dictwriter``.Otherwise the CSV 
        output filename should be the `c.data` or the given variable name.
    :param data: Sequence Data in Python obj to write. 
    :param todo: Action to perform with JSON: 
        - reader|DictReader: Load data from the JSON file 
        - writer|DictWriter: Write data from the Python object 
        and create a CSV file
    :param savepath: If ``default``  should save the `csv_fn` 
        If path does not exist, should save to the <'_savecsv_'>
        default path.
    :param fieldnames: is a sequence of keys that identify the order
        in which values in the dictionary passed to the `writerow()`
            method are written `csv_fn` file.
    :param savepath: If ``default``  should save the `csv_fn` 
        If path does not exist, should save to the <'_savecsv_'>
        default path .
    :param verbose: int, control the verbosity. Output messages
    :param csvkws: additional keywords csv class arguments 
    
    .. see also:: Read more about CSV module in:
        https://docs.python.org/3/library/csv.html or find some examples
        here https://www.programcreek.com/python/example/3190/csv.DictWriter 
        or find some FAQS here: 
    https://stackoverflow.com/questions/10373247/how-do-i-write-a-python-dictionary-to-a-csv-file
        ...
    :Example: 
        >>> import pycsamt.utils.geo_utils as GU 
        >>> import pycsamt.bases as BS
        >>> geo_kws ={'oc2d': GU.INVERS_KWS, 
                      'TRES':GU.TRES, 'LN':GU.LNS}
        # write data and save to  'csvtest.csv' file 
        # here the `data` is a sequence of dictionary geo_kws
        >>> BS.parse_csv(csv_fn = 'csvtest.csv',data = [geo_kws], 
                         fieldnames = geo_kws.keys(),todo= 'dictwriter',
                         savepath = 'data/saveCSV')
        # collect csv data from the 'csvtest.csv' file 
        >>> BS.parse_csv(csv_fn ='data/saveCSV/csvtest.csv',
                         todo='dictreader',fieldnames = geo_kws.keys()
                         )
    
    """
    todo, domsg =return_ctask(todo) 
    
    if todo.find('write')>=0:
        csv_fn = get_config_fname_from_varname(
            data, config_fname= csv_fn, config='.csv')
    try : 
        if todo =='reader': 
            with open (csv_fn, 'r') as csv_f : 
                csv_reader = csv.reader(csv_f) # iterator 
                data =[ row for row in csv_reader]
                
        elif todo=='writer': 
            # write without a blank line, --> new_line =''
            with open(f'{csv_fn}.csv', 'w', newline ='',
                      encoding ='utf8') as new_csvf:
                csv_writer = csv.writer(new_csvf, **csvkws)
                csv_writer.writerows(data) if len(
                    data ) > 1 else csv_writer.writerow(data)  
                # for row in data:
                #     csv_writer.writerow(row) 
        elif todo=='dictreader':
            with open (csv_fn, 'r', encoding ='utf8') as csv_f : 
                # generate an iterator obj 
                csv_reader= csv.DictReader (csv_f, fieldnames= fieldnames) 
                # return csvobj as a list of dicts
                data = list(csv_reader) 
        
        elif todo=='dictwriter':
            with open(f'{csv_fn}.csv', 'w') as new_csvf:
                csv_writer = csv.DictWriter(new_csvf, **csvkws)
                if header:
                    csv_writer.writeheader()
                # DictWriter.writerows()expect a list of dicts,
                # while DictWriter.writerow() expect a single row of dict.
                csv_writer.writerow(data) if isinstance(
                    data , dict) else csv_writer.writerows(data)  
                
    except csv.Error: 
        raise csv.Error(f"Unable {domsg} CSV {csv_fn!r} file. "
                      "Please check your file.")
    except: 

        msg =''.join([
        f"{'Unrecognizable file' if todo.find('read')>=0 else'Unable to write'}"
        ])
        
        raise TypeError(f'{msg} {csv_fn!r}. Please check your'
                        f" {'file' if todo.find('read')>=0 else 'data'}.")
    cparser_manager(f'{csv_fn}.csv',savepath, todo=todo, dpath='_savecsv_', 
                    verbose=verbose , config='CSV' )
    
    return data 
   
def parse_json(json_fn: Union[A, str] =None,
               data:Dict[str, Type[Any]] =None, 
               todo:Optional[str] ='load',
               savepath:Optional[str] =None,
               verbose:int =0,
               **jsonkws) -> T:
    """ Parse Java Script Object Notation file and collect data from JSON
    config file. 
    
    :param json_fn: Json filename, URL or output JSON name if `data` is 
        given and `todo` is set to ``dump``.Otherwise the JSON output filename 
        should be the `data` or the given variable name.
    :param data: Data in Python obj to serialize. 
    :param todo: Action to perform with JSON: 
        - load: Load data from the JSON file 
        - dump: serialize data from the Python object and create a JSON file
    :param savepath: If ``default``  should save the `json_fn` 
        If path does not exist, should save to the <'_savejson_'>
        default path .
    :param verbose: int, control the verbosity. Output messages
    
    .. see also:: Read more about JSON doc
            https://docs.python.org/3/library/json.html
         or https://www.w3schools.com/python/python_json.asp 
         or https://www.geeksforgeeks.org/json-load-in-python/
         ...
 
    :Example: 
        >>> import pycsamt.utils.geo_utils as GU 
        >>> import pycsamt.bases as BS
        >>> geo_kws ={'oc2d': GU.INVERS_KWS, 
                      'TRES':GU.TRES, 'LN':GU.LNS}
        # serialize json data and save to  'jsontest.json' file
        >>> BS.parse_json(json_fn = 'jsontest.json', 
                          data=geo_kws, todo='dump', indent=3,
                          savepath ='data/saveJSON', sort_keys=True)
        # Load data from 'jsontest.json' file.
        >>> BS.parse_json(json_fn='data/saveJSON/jsontest.json', todo ='load')
    
    """
    todo, domsg =return_ctask(todo)
    # read urls by default json_fn can hold a url 
    try :
        if json_fn.find('http') >=0 : 
            todo, json_fn, data = fetch_json_data_from_url(json_fn, todo)
    except:
        #'NoneType' object has no attribute 'find' if data is not given
        pass 

    if todo.find('dump')>=0:
        json_fn = get_config_fname_from_varname(
            data, config_fname= json_fn, config='.json')
        
    JSON = dict(load=json.load,# use loads rather than load  
                loads=json.loads, 
                dump= json.dump, 
                dumps= json.dumps)
    try :
        if todo=='load': # read JSON files 
            with open(json_fn) as fj: 
                data =  JSON[todo](fj)  
        elif todo=='loads': # can be JSON string format 
            data = JSON[todo](json_fn) 
        elif todo =='dump': # store data in JSON file.
            with open(f'{json_fn}.json', 'w') as fw: 
                data = JSON[todo](data, fw, **jsonkws)
        elif todo=='dumps': # store data in JSON format not output file.
            data = JSON[todo](data, **jsonkws)

    except json.JSONDecodeError: 
        raise json.JSONDecodeError(f"Unable {domsg} JSON {json_fn!r} file. "
                              "Please check your file.", f'{json_fn!r}', 1)
    except: 
        msg =''.join([
        f"{'Unrecognizable file' if todo.find('load')>=0 else'Unable to serialize'}"
        ])
        
        raise TypeError(f'{msg} {json_fn!r}. Please check your'
                        f" {'file' if todo.find('load')>=0 else 'data'}.")
        
    cparser_manager(f'{json_fn}.json',savepath, todo=todo, dpath='_savejson_', 
                    verbose=verbose , config='JSON' )

    return data 
 
def fetch_json_data_from_url (url:str , todo:str ='load'): 
    """ Retrieve JSON data from url 
    :param url: Universal Resource Locator .
    :param todo:  Action to perform with JSON:
        - load: Load data from the JSON file 
        - dump: serialize data from the Python object and create a JSON file
    """
    with urllib.request.urlopen(url) as jresponse :
        source = jresponse.read()
    data = json.loads(source)
    if todo .find('load')>=0:
        todo , json_fn  ='loads', source 
        
    if todo.find('dump')>=0:  # then collect the data and dump it
        # set json default filename 
        todo, json_fn = 'dumps',  '_urlsourcejsonf.json'  
        
    return todo, json_fn, data 
    
    
def return_ctask (todo:Optional[str]=None) -> Tuple [str, str]: 
    """ Get the convenient task to do if users misinput the `todo` action.
    
    :param todo: Action to perform: 
        - load: Load data from the config [YAML|CSV|JSON] file
        - dump: serialize data from the Python object and 
            create a config [YAML|CSV|JSON] file."""
            
    def p_csv(v, cond='dict', base='reader'):
        """ Read csv instead. 
        :param v: str, value to do 
        :param cond: str, condition if  found in the value `v`. 
        :param base: str, base task to do if condition `cond` is not met. 
        
        :Example: 
            
        >>> todo = 'readingbook' 
        >>> p_csv(todo) <=> 'dictreader' if todo.find('dict')>=0 else 'reader' 
        """
        return  f'{cond}{base}' if v.find(cond) >=0 else base   
    
    ltags = ('load', 'recover', True, 'fetch')
    dtags = ('serialized', 'dump', 'save', 'write','serialize')
    if todo is None: 
        raise ValueError('NoneType action can not be perform. Please '
                         'specify your action: `load` or `dump`?' )
    
    todo =str(todo).lower() 
    ltags = list(ltags) + [todo] if  todo=='loads' else ltags
    dtags= list(dtags) +[todo] if  todo=='dumps' else dtags 

    if todo in ltags: 
        todo = 'loads' if todo=='loads' else 'load'
        domsg= 'to parse'
    elif todo in dtags: 
        todo = 'dumps' if todo=='dumps' else 'dump'
        domsg  ='to serialize'
    elif todo.find('read')>=0:
        todo = p_csv(todo)
        domsg= 'to read'
    elif todo.find('write')>=0: 
        todo = p_csv(todo, base ='writer')
        domsg =' to write'
        
    else :
        raise ValueError(f'Wrong action {todo!r}. Please select'
                         f' the right action to perform: `load` or `dump`?'
                        ' for [JSON|YAML] and `read` or `write`? '
                        'for [CSV].')
    return todo, domsg  

def parse_yaml (yml_fn:str =None, data:Union[T , list, tuple]=None,
                todo:str ='load', savepath:Optional[str] =None,
                verbose:int =0, **ymlkws) -> T: 
    """ Parse yml file and collect data from YAML config file. 
    
    :param yml_fn: yaml filename and can be the output YAML name if `data` is 
        given and `todo` is set to ``dump``.Otherwise the YAML output filename 
        should be the `data` or the given variable name.
    :param data: Data in Python obj to serialize. 
    :param todo: Action to perform with YAML: 
        - load: Load data from the YAML file 
        - dump: serialize data from the Python object and create a YAML file
    :param savepath: If ``default``  should save the `yml_fn` 
        to the default path otherwise should store to the convenient path.
        If path does not exist, should set to the default path.
    :param verbose: int, control the verbosity. Output messages
    
    .. see also:: Read more about YAML file https://pynative.com/python-yaml/
         or https://python.land/data-processing/python-yaml and download YAML 
         at https://pypi.org/project/PyYAML/
         ...
    :Example: 
        >>> import pycsamt.utils.geo_utils as GU 
        >>> geo_kws ={'oc2d': GU.INVERS_KWS, 
                      'TRES':GU.TRES, 'LN':GU.LNS}
        # serialize yaml data 
        >>> GU.parse_yaml(data =geo_kws , todo='dump',
                          savepath =None, sort_keys=True, 
                          default_flow_style=True, verbose =1)
        ... '_saveYAML_/testgeo.yml'
        # load the data from data `testgeo.yml`
        >>> GU.parse_yaml('_saveYAML_/testgeo.yml')
        ... {'LN': ['river water','fracture zone','MWG','LWG','granite',
        ...   'igneous rocks','basement rocks'],
        ...  'TRES': [10, 66, 70, 100, 1000, 3000],
        ...  'oc2d': {'data_fn': 'data/occam2D\\OccamDataFile.dat',
        ...           'iter_fn': 'data/occam2D\\ITER17.iter',
        ...           'mesh_fn': 'data/occam2D\\Occam2DMesh',
        ...           'model_fn': 'data/occam2D\\Occam2DModel'}}
    """ 
    
    todo, domsg =return_ctask(todo)
    #in the case user use dumps or loads with 's'at the end 
    if todo.find('dump')>= 0: 
        todo='dump'
    if todo.find('load')>=0:
        todo='load'
    if todo=='dump':
        yml_fn = get_config_fname_from_varname(data, yml_fn)
    try :
        if todo=='load':
            with open(yml_fn) as fy: 
                data =  yaml.load(fy, Loader=yaml.SafeLoader)  
                # args =yaml.safe_load(fy)
        elif todo =='dump':
        
            with open(f'{yml_fn}.yml', 'w') as fw: 
                data = yaml.dump(data, fw, **ymlkws)
    except yaml.YAMLError: 
        raise yaml.YAMLError(f"Unable {domsg} YAML {yml_fn!r} file. "
                             'Please check your file.')
    except: 
        msg =''.join([
        f"{'Unrecognizable file' if todo=='load'else'Unable to serialize'}"
        ])
        
        raise TypeError(f'{msg} {yml_fn!r}. Please check your'
                        f" {'file' if todo=='load' else 'data'}.")
        
    cparser_manager(f'{yml_fn}.yml',savepath, todo=todo, dpath='_saveyaml_', 
                    verbose=verbose , config='YAML' )

    return data 
 
def cparser_manager (cfile:Union[A, str],
                     savepath:Union[Type[Any], str] =None, 
                     todo:str ='load', dpath:Optional[str] =None,
                     verbose:int =0, **pkws): 
    """ Save and output message according to the action. 
    :param cfile: name of the configuration file
    :param savepath: Path-like object 
    :param dpath: default path 
    :param todo: Action to perform with config file. Can ve 
        ``load`` or ``dump``
    :param config: Type of configuration file. Can be ['YAML|CSV|JSON]
    :param verbose: int, control the verbosity. Output messages
    
    """
    if savepath is not None:
        if savepath =='default': 
            savepath = None 
        yml_fn,_= move_cfile(cfile,savepath, dpath=dpath)
    if verbose > 0: 
        print_cmsg(yml_fn, todo, **pkws)
        
    
def get_config_fname_from_varname(data: Union[T, list, tuple],
                                  config_fname:Optional[str] =None,
                                  config:str ='.yml') -> str: 
    """ use the variable name given to data as the config file name.
    :param data: Given data to retrieve the variable name 
    :param config_fname: Configurate variable filename. If ``None`` , use 
        the name of the given varibale data 
    :param config: Type of file for configuration. Can be ``json``, ``yml`` 
        or ``csv`` file. default is ``yml``.
    :return: str, the configuration data.
    
    """
    try:
        if '.' in config: 
            config =config.replace('.','')
    except:pass # in the case None is given
    
    if config_fname is None: # get the varname 
        # try : 
        #     from varname.helpers import Wrapper 
        # except ImportError: 
        #     import_varname=False 
        #     import_varname = FU.subprocess_module_installation('varname')
        #     if import_varname: 
        #         from varname.helpers import Wrapper 
        # else : import_varname=True 
        try : 
            for c, n in zip(['yml', 'yaml', 'json', 'csv'],
                            ['cy.data', 'cy.data', 'cj.data',
                             'c.data']):
                if config ==c:
                    config_fname= n
                    break 
            if config_fname is None:
                raise # and go to except  
        except: 
            #using fstring 
            config_fname= f'{data=}'.split('=')[0]
            
    elif config_fname is not None: 
        config_fname= config_fname.replace(
            f'.{config}', '').replace(f'.{config}', '').replace('.yaml', '')
    
    return config_fname
    

def move_cfile (cfile:str , savepath:Optional[str]=None, **ckws):
    """ Move file to its savepath and output message. 
    If path does not exist, should create one to save data.
    :param cfile: name of the configuration file
    :param savepath: Path-like object 
    :param dpath: default path 
    
    :returns: 
        - configuration file 
        - out message 
    """
    savepath = FU.cpath(savepath, **ckws)
    try :shutil.move(cfile, savepath)
    except: warnings.warn("It seems the path already exists!")
    
    cfile = os.path.join(savepath, cfile)
    
    msg = ''.join([
    f'--> Data was successfully stored to {os.path.basename(cfile)!r}', 
        f' and saved to {os.path.realpath(cfile)!r}.']
        )
        
    return cfile, msg

def print_cmsg(cfile:str, todo:str='load', config:str='YAML') -> str: 
    """ Output configuration message. 
    
    :param cfile: name of the configuration file
    :param todo: Action to perform with config file. Can ve 
        ``load`` or ``dump``
    :param config: Type of configuration file. Can be [YAML|CSV|JSON]
    """
    if todo=='load': 
        msg = ''.join([
        f'--> Data was successfully stored to {os.path.basename(cfile)!r}', 
            f' and saved to {os.path.realpath(cfile)!r}.']
            )
    elif todo=='dump': 
        msg =''.join([ f"--> {config.upper()} {os.path.basename(cfile)!r}", 
                      " data was sucessfully loaded."])
    return msg 

def sPath (name_of_path:str):
    """ Savepath func. Create a path  with `name_of_path` if path not exists."""
    try :
        savepath = os.path.join(os.getcwd(), name_of_path)
        if not os.path.isdir(savepath):
            os.mkdir(name_of_path)#  mode =0o666)
    except :
        warnings.warn("The path seems to be existed !")
        return
    return savepath 




