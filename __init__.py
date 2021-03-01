
# version 
__version__ ='1.0.01'
# import logging and load pyCSAMT default logging config

import logging 

from csamtpy.utils._csamtpylog import csamtpylog

csamtpylog.load_configure()

# set loging Level
logging.getLogger('matplotlib').setLevel(logging.WARNING)

