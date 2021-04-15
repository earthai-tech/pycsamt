
import os 
from tests.__init__ import csamtpylog

def remove_control(rm_file, extensions=['.xlsx', '.csv', '.dat', '.bln'],
                   **kws) :
    """
    small tip to remove file in the export directory if user doest
    know it specific extension.
    rm_file: is only the name of fmile. It doesnt include the extension 
    of files. 
    
    """
    type_of_file =kws.pop('type_of_file', 'File')
    
    extensions=[fex.lower() for fex in extensions]
    
    for file_extension in extensions :
        if file_extension.find('.') <0 : 
            file_extension ='.'+file_extension
        #build file 
        rm_file = ''.join([rm_file, file_extension])
        
        try : 
            if os.path.exists(rm_file):
                os.remove(rm_file) 
        except : 
            csamtpylog().get_csamtpy_logger().error(
                'File <{0}> does not exist!'.format(rm_file))
            f=0
        else :
            csamtpylog().get_csamtpy_logger().info(
                '{0}: <{1}> is sucessfully removed !'.
                format(type_of_file, os.path.basename(rm_file)))
            f=1
            
            return rm_file 
        
        if f==0 : return 