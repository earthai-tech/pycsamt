
# version 
__version__ ='1.0.01'
# import logging and load pyCSAMT default logging config

import logging 
import sys
import os 

from csamtpy.utils._csamtpylog import csamtpylog

csamtpylog.load_configure()

# set loging Level
logging.getLogger('matplotlib').setLevel(logging.WARNING)

# let set the systeme path find memory dataBase
 
sys.path.insert(0, os.path.abspath('.'))  
sys.path.insert(0, os.path.abspath('..')) 
sys.path.insert(0, os.path.abspath('../..'))
# sys.path.insert(0, os.path.abspath('.//csamtpy'))
# sys.path.insert(0, os.path.abspath('.//viewer'))
# sys.path.insert(0, os.path.abspath('.//geodrill'))

if __name__=='__main__': 

    print(sys.path)