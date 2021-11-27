# -*- coding: utf-8 -*-
#       created on Fri Oct 29 12:31:01 2021
import os
import warnings
import copy
import joblib 
import pickle
import datetime
import shutil 

from pycsamt.utils._csamtpylog import csamtpylog
_logger=csamtpylog.get_csamtpy_logger(__name__)

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

def serialize_data(data , filename:str=None,  force:bool=True, 
                         savepath:str =None, verbose:int =0): 
    """ Store a data into a binary binary file 
    
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
        >>> file = GU.serialize_data(data, filename=None,  force=True, 
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
        try : savepath = sPath(dirname)
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
    
def load_serialized_data (filename:str, verbose:int =0): 
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





