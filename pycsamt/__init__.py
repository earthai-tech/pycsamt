# -*- coding: utf-8 -*-
#  Author: Kouadio Laurent alias Daniel <etanoyau@gmail.com> 
#  GNU Licence: GPL version 3
  
"""
`pycsamt`_ : A Python Toolbox for audio-magnetotellurics 
=========================================================

The toolbox is designed to bring a piece of solution by providing fast & robust 
data processing and signal recovering of the audio-magnetotelluric(AMT) methods 
i.e the frequency above 1Hz. It also provides processing tools for filtering 
data and restoring amplitudes at weak frequency (tensor or or scalar values) 
in the "dead-dand" especially for the natural source AMT where the signals are  
not under our control and suffer from frequency ranges with little or no signal.
Originally, the software was focused on `CSAMT`_ data processing, hereinafter 
the suffix CSAMT. The development, currently encompasses all the AMT methods.
Futhermore, it also implements a special subpackage named `geodrill` for reducing
the numerous unsucessful drillings mostly occurred due to their wrong locations
after AMT geophysical surveys. 

.. |CSAMT| replace :: controlled source audio-frequency magnetotelluric
.. _documentation: https://pycsamt.readthedocs.io/en/latest/
.. _pycsamt: https://github.com/WEgeophysics/pycsamt

"""
# +---------------------------------------------------------------------------+
# | Module will take a while to install the dependencies the first time it is |
# | launched. This won't occur if the package  is installed via PyPI since the|
# | latter will take in charge the dependencies installation during the       |
# | setting up progress.                                                      |
# +---------------------------------------------------------------------------+
# define pycamt release version below
# see https://packaging.python.org/guides/single-sourcing-package-version/ 

__version__ = "1.2.0"

# load pycsamt default logging config
import sys 
import subprocess 
import logging 
import os
import tempfile

# setup the package 
if __package__ is None or __name__ == '__main__': 
    sys.path.append( os.path.dirname(os.path.dirname(__file__)))
    sys.path.insert(0, os.path.dirname(__file__))
    __package__= 'pycsamt' 
    
# configure the logger 
from pycsamt.utils._csamtpylog import csamtpylog
try: 
    csamtpylog.load_configure(os.path.join(
        os.path.abspath('.'),'pycsamt', 'utils', "p.configlog.yml"))
except: 
    csamtpylog.load_configure(os.path.join(
        os.path.abspath('.'),'utils', "p.configlog.yml"))

# set loging Level
logging.getLogger('matplotlib').setLevel(logging.WARNING)

epsg_dict = {
    28350: ['+proj=utm +zone=50 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs', 50],
    28351: ['+proj=utm +zone=51 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs', 51],
    28352: ['+proj=utm +zone=52 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs', 52],
    28353: ['+proj=utm +zone=53 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs', 53],
    28354: ['+proj=utm +zone=54 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs', 54],
    28355: ['+proj=utm +zone=55 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs', 55],
    28356: ['+proj=utm +zone=56 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs', 56],
    3112:  [
        '+proj=lcc +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=134 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs',
        0 ],
    4326: ['+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs', 0]
}

def is_installing (
    module:str ,
    upgrade:bool =True ,
    DEVNULL:bool =False, 
    action:bool =True,
    verbose:int =0,
    **subpkws
 )-> bool : 
    """ Install or uninstall a module using the subprocess under the hood.
    
    install the dependencies if failed using PyPI. This might delay the software
    to run for the fisrt time use.
    
    :param module: str, module name 
    :param upgrade:bool, install the lastest version.
    :param verbose:output a message 
    :param DEVNULL: decline the stdoutput the message in the console 
    :param action: str, install or uninstall a module 
    :param subpkws: additional subprocess keywords arguments.
    
    :return:  module successfull installed or not. 
    :rtype: bool 
    
    :Example --> install and uninstall the `tqdm` packages: :: 
        
        >>> from pycsamt import is_installing
        >>> is_installing(
            'tqdm', action ='install', DEVNULL=True, verbose =1)
        >>> is_installing(
            'tqdm', action ='uninstall', verbose =1)
    """
    _logger= csamtpylog.get_csamtpy_logger(__name__)
    #implement pip as subprocess 
    # refer to https://pythongeeks.org/subprocess-in-python/
    if not action: 
        if verbose > 0 :
            print("---> No action `install`or `uninstall`"
                  f" of the module {module!r} performed.")
        return action  # DO NOTHING 
    
    success=False 

    action_msg ='uninstallation' if action =='uninstall' else 'installation' 

    if action in ('install', 'uninstall', True) and verbose > 0:
        print(f'---> Module {module!r} {action_msg} will take a while,'
              ' please be patient...')
        
    cmdg =f'<pip install {module}> | <python -m pip install {module}>'\
        if action in (True, 'install') else ''.join([
            f'<pip uninstall {module} -y> or <pip3 uninstall {module} -y ',
            f'or <python -m pip uninstall {module} -y>.'])
        
    upgrade ='--upgrade' if upgrade else '' 
    
    if action == 'uninstall':
        upgrade= '-y' # Don't ask for confirmation of uninstall deletions.
    elif action in ('install', True):
        action = 'install'

    cmd = ['-m', 'pip', f'{action}', f'{module}', f'{upgrade}']

    try: 
        STDOUT = subprocess.DEVNULL if DEVNULL else None 
        STDERR= subprocess.STDOUT if DEVNULL else None 
    
        subprocess.check_call(
            [sys.executable] + cmd, stdout= STDOUT, stderr=STDERR,
                              **subpkws)
        if action in (True, 'install'):
            # freeze the dependancies
            reqs = subprocess.check_output(
                [sys.executable,'-m', 'pip','freeze'])
            [r.decode().split('==')[0] for r in reqs.split()]
            _logger.info( f"{action_msg.capitalize()} of `{module}` "
                         "and dependancies was successfully done!") 
        success=True
        
    except: 
        _logger.error(f"Failed to {action} the module =`{module}`.")
        
        if verbose > 0 : 
            print(f'---> Module {module!r} {action_msg} failed. Please use'
                f' the following command: {cmdg} to manually do it.')
    else : 
        if verbose > 0: 
            print(f"{action_msg.capitalize()} of `{module}` "
                      "and dependancies was successfully done!") 
        
    return success
# install basic dependencies 
itqdm =False 
try : 
    import tqdm 
except ImportError: 
    itqdm =is_installing (
        'tqdm',
        DEVNULL =True)
else :
    itqdm =True 
    
# auto intall numba and the suitable version.
# make sure that numby work at the present date (09/22)
# with numpy version 1.12 or less   
inumba = False 
 
try : 
    import numba   
except ImportError: 
    inumba =is_installing (
        'numba',
        DEVNULL =True)
else :
    inumba  =True 
    
#Install MTpy if not install yet. 
imtpy =False 
try : 
    import mtpy  
except ImportError: 
    imtpy =is_installing  (
        'mtpy',
        DEVNULL =True)
else :
    imtpy =True 
    
# set path for demo and testing module
# Get the repository dir like https://github.com/WEgeophysics/pycsamt /
#C:github.com/WEgeophysics/pyCSAMT  
PYCSAMT_ROOT = os.path.normpath(
    os.path.abspath(
        os.path.dirname(
            os.path.dirname(__file__)
        )
    )
)

EDI_DATA_DIR = os.path.normpath(
    os.path.join(PYCSAMT_ROOT , 'data/edi'))
AVG_DATA_DIR = os.path.normpath(
    os.path.join(PYCSAMT_ROOT , 'data/avg'))
J_DATA_DIR = os.path.normpath(
    os.path.join(PYCSAMT_ROOT , 'data/j'))
DRILL_DATA_DIR = os.path.normpath(
    os.path.join(PYCSAMT_ROOT , 'data/drill_examples_files'))
OCCAM2D_DATA_DIR = os.path.normpath(
    os.path.join(PYCSAMT_ROOT , 'data/occam2d'))
STN_DATA_DIR = os.path.normpath(
    os.path.join(PYCSAMT_ROOT , 'data/stn_profiles'))
CONFIG_DATA_DIR = os.path.normpath(
    os.path.join(PYCSAMT_ROOT , 'data/_conffiles'))

SYSTEM_TEMP_DIR = tempfile.gettempdir()
NEW_TEMP_DIR=tempfile.mkdtemp(prefix="pycsamt_tmpdir_")








