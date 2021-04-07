# define mtpy release version below
# see https://packaging.python.org/guides/single-sourcing-package-version/ 
__version__ = "1.0.03"

# load pycsamt default logging config
import os 
import logging 
import sys
from pycsamt.utils._csamtpylog import csamtpylog

csamtpylog.load_configure(os.path.join(os.path.abspath('.'),
                                       'pycsamt', '_logfile',
                                       "main_logging_configfile.yml"))

# set loging Level
logging.getLogger('matplotlib').setLevel(logging.WARNING)

# let set the systeme path find memory dataBase
sys.path.insert(0, os.path.abspath('.'))  
sys.path.insert(0, os.path.abspath('./viewer')) 
sys.path.insert(0, os.path.abspath('./geodrill'))

# if __name__=='__main__': 

    # print(sys.path)