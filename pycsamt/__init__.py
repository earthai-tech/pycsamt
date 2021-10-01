# define mtpy release version below
# see https://packaging.python.org/guides/single-sourcing-package-version/ 
__version__ = "1.0.7"

# load pycsamt default logging config
import os 
import logging 

from pycsamt.utils._csamtpylog import csamtpylog

csamtpylog.load_configure(os.path.join(os.path.abspath('.'),
                                        'pycsamt', '_logfile',
                                        "main_logging_configfile.yml"))

# set loging Level
logging.getLogger('matplotlib').setLevel(logging.WARNING)

