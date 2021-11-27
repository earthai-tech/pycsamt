# define mtpy release version below
# see https://packaging.python.org/guides/single-sourcing-package-version/ 
__version__ = "1.1.0"

# load pycsamt default logging config
import sys 
import os 
import logging 

from pycsamt.utils._csamtpylog import csamtpylog

if __package__ is None or __name__ == '__main__': 
    sys.path.append( os.path.dirname(os.path.dirname(__file__)))
    sys.path.insert(0, os.path.dirname(__file__))
    __package__= 'pycsamt'

csamtpylog.load_configure(os.path.join(
    os.path.abspath('.'),'pycsamt', 'utils', "p.configlog.yml"))

# set loging Level
logging.getLogger('matplotlib').setLevel(logging.WARNING)

# intall basic module tqdm 
itqdm =False 
try : 
    import tqdm 
except ImportError: 
    from pycsamt.utils.func_utils import  subprocess_module_installation
    itqdm =subprocess_module_installation ('tqdm', DEVNULL =True)
else :
    itqdm =True 