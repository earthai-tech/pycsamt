
import os
import unittest 
from tests.modeling.__init__ import diff_files


def compare_diff_files(refout, refexp):
    """
    Compare diff files like expected files and output files generated after 
    runnning scripts. 
    :param refout: list of reference output files generated after runing scripts
    :type refout: list 
    
    :param refexp: recreated expected files for comparison 
    :param refexp: list 

    """
    for outfile , expfile in zip(sorted(refout), 
                                   sorted(refexp)):
            unittest.TestCase.assertTrue(os.path.isfile(outfile),
                                "Ref output data file does not exist,"
                                "nothing to compare with"
                                )
            
            print(("Comparing", outfile, "and", expfile))
            
            is_identical, msg = diff_files(outfile, expfile, ignores=['Date/Time:'])
            print(msg)
            unittest.TestCase.assertTrue(is_identical, 
                            "The output file is not the same with the baseline file.")