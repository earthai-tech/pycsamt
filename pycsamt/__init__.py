# define pycamt release version below
# see https://packaging.python.org/guides/single-sourcing-package-version/ 
__version__ = "1.1.2"

# load pycsamt default logging config
import sys 
import logging 
import os
import tempfile

from pycsamt.utils._csamtpylog import csamtpylog
from pycsamt.utils.func_utils import  subprocess_module_installation

if __package__ is None or __name__ == '__main__': 
    sys.path.append( os.path.dirname(os.path.dirname(__file__)))
    sys.path.insert(0, os.path.dirname(__file__))
    __package__= 'pycsamt'
    

try: 
    csamtpylog.load_configure(os.path.join(
        os.path.abspath('.'),'pycsamt', 'utils', "p.configlog.yml"))
except: 
    csamtpylog.load_configure(os.path.join(
        os.path.abspath('.'),'utils', "p.configlog.yml"))

# set loging Level
logging.getLogger('matplotlib').setLevel(logging.WARNING)

# intall basic module tqdm 
itqdm =False 
try : 
    import tqdm 
except ImportError: 
    itqdm =subprocess_module_installation (
        'tqdm',
        DEVNULL =True)
else :
    itqdm =True 
    
# auto intall numba  if not installed yet  
inumba = False 
 
try : 
    import numba   
except ImportError: 
    inumba =subprocess_module_installation (
        'numba', 
        DEVNULL =True)
else :
    inumba  =True 
    
# auto install mtpy 
imtpy =False 
try : 
    import mtpy  
except ImportError: 
    imtpy =subprocess_module_installation (
        'mtpy',
        DEVNULL =True)
else :
    imtpy =True 
    
epsg_dict = {
    28350: ['+proj=utm +zone=50 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs', 50],
    28351: ['+proj=utm +zone=51 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs', 51],
    28352: ['+proj=utm +zone=52 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs', 52],
    28353: ['+proj=utm +zone=53 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs', 53],
    28354: ['+proj=utm +zone=54 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs', 54],
    28355: ['+proj=utm +zone=55 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs', 55],
    28356: ['+proj=utm +zone=56 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs', 56],
    3112: [
        '+proj=lcc +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=134 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs',
        0],
    4326: ['+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs', 0]
}

# Get the repository dir like https://github.com/WEgeophysics/pycsamt /C:github.com/WEgeophysics/pyCSAMT  
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

