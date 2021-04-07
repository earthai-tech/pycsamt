# import mtpy , if not exist try to install it
# process  to install mtpy using subprocess

import sys, warnings
import subprocess
from pycsamt.utils._csamtpylog import csamtpylog

SUCCESS_IMPORT_MTPY = False

try :
    import mtpy 
    
except :
    # implement pip as a subprocess:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install',
    'mtpy'])
    # process output with an API in the subprocess module:
    reqs = subprocess.check_output([sys.executable, '-m', 'pip',
    'freeze'])
    installed_packages = [r.decode().split('==')[0] for r in reqs.split()]
        # get list of intalled dependancies 
    csamtpylog().get_csamtpy_logger().info(
        'MTpy module was successfully intalled with  dependancies.')
    
    SUCCESS_IMPORT_MTPY =True  
 
else : 
    mtpy_version = [int(ss) for ss in mtpy.__version__.split('.')]
    if mtpy_version[0] == 1:
        if mtpy_version[1] < 1:
            warnings.warn('Note: need mtpy version 1.1.0 or higher to write imput occam2D files '
                          'might not work propertly.', ImportWarning)
            csamtpylog().get_csamtpy_logger().warning(
                'Note: need mtpy version 1.14.0 to write imput occam2D files'
                ' or might not work properly.')
            SUCCESS_IMPORT_MTPY  =False 
            
        else : 
            msg = 'MTpy is already installed , trying import `occam2d` module !'
            csamtpylog().get_csamtpy_logger().info( ''.join([msg, 
                                '`occam2d from MTpy was successfully imported !']))
            SUCCESS_IMPORT_MTPY  =True 
    
