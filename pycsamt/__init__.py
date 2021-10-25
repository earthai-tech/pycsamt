# define mtpy release version below
# see https://packaging.python.org/guides/single-sourcing-package-version/ 
__version__ = "1.0.9"

# load pycsamt default logging config
import sys 
import os 
import logging 

from pycsamt.utils._csamtpylog import csamtpylog

if __package__ is None or __name__ == '__main__': 
    sys.path.append( os.path.dirname(os.path.dirname(__file__)))
    sys.path.insert(0, os.path.dirname(__file__))
    __package__= 'pycsamt'

csamtpylog.load_configure(os.path.join(os.path.abspath('.'),
                                        'pycsamt', '_logfile',
                                        "main_logging_configfile.yml"))

# set loging Level
logging.getLogger('matplotlib').setLevel(logging.WARNING)

