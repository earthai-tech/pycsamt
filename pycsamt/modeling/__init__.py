# import mtpy , if not exist try to install it
# using subprocess module

import warnings

from pycsamt.utils._csamtpylog import csamtpylog

SUCCESS_IMPORT_MTPY = False

try :
    import mtpy 
    
except :
    from pycsamt.utils import func_utils 
    
    SUCCESS_IMPORT_MTPY= func_utils.subprocess_module_installation('mtpy')
    msg = 'was successfully intalled.' \
        if SUCCESS_IMPORT_MTPY else 'installation was failed!'
    csamtpylog().get_csamtpy_logger().info(f"MTpy {msg}")

else : 
    mtpy_version = [int(ss) for ss in mtpy.__version__.split('.')]
    if mtpy_version[0] == 1:
        if mtpy_version[1] >= 1:
            SUCCESS_IMPORT_MTPY  =True 
        elif mtpy_version[1] < 1:
            warnings.warn(
                'Note: need mtpy version 1.1.0 or higher to write imput  '
                 'occam2D files might not work propertly.', ImportWarning)
            csamtpylog().get_csamtpy_logger().warning(
                'Note: need mtpy version 1.1.0 or more to write imput'
                '  occam2D files or might not work properly.')
            SUCCESS_IMPORT_MTPY  =False 
            

    
