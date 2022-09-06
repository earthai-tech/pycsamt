# -*- coding: utf-8 -*-
"""
Description:
    Test to write corrected edi. Collected uncorrected edi , read edi and find 
    the data type either : MT or EMAP  and apply filters like `ama` , `flma` 
    `tma`, `ss` or `dist`. 
    {`ama` : adaptative moving average, 
    `flma`: fixed-length-dipole moving average 
    `tma`: trimming moving average 
    `ss` : static shift removal  and 
    `dist`: distorsion removal}

References:
    scripts/write_corrected_edi.py
Created on Wed Apr 14 13:15:19 2021

@author: @Daniel03
"""
import os
import  unittest 

from pycsamt.processing import corr

from tests import  (make_temp_dir,
                    TEST_TEMP_DIR,
                    EDI_DATA_DIR)
from tests import survey_testname
from tests.modeling.__init__ import (csamtpylog,
                                     reset_matplotlib)
from tests.processing.__init__ import compare_diff_files

class TestEDICOR(unittest.TestCase):
    """
    Writing corrected edi files from raw edifiles 
    * call edifiles from '/data/edi/' directory 
    * outputs 'temp/TESTEDICOR' from `'nibykro_survey' area
    
    """
    _filters =['ama', 
               'tma', 
               'flma',
               'ss', 
        ]
    params_filters = {'number_of_points':7, 
                      'dipole_length':50., 
                      "number_of_skin_depth": 7, 
                      'reduce_res_factor_x':1, 
                      'reduce_res_factor_y': 1 , 
                      }
    @classmethod 
    def setUpClass(cls):
        """
        Reset building matplotlib plot and generate tempdir inputfiles 
        
        """

        cls._temp_dir = make_temp_dir(cls.__name__)

    def setUp(self): 
        reset_matplotlib()
        if not os.path.isdir(TEST_TEMP_DIR):
            print('--> outdir not exist , set to None !')
            csamtpylog().get_csamtpy_logger().error('Outdir does not exist !')
        
    def wrap_single_processing_filter(self, FILTER,edi_path, 
                                      filename =None, save_path=None, **kws):
        """
        Single method used to wrap a simple processing filter and apply it for
        all filters 
        return name of filter successfully run. 
        .. note:: reference frequency should be comoputed automatically
        
        """
        # new edi output filenames 
  
        csamtpylog().get_csamtpy_logger().info(
            'Ready to enter %s filter and parameters '
            'in processing module'% FILTER)
        try :
            corr.Processing().write_corrected_edi(data_fn = edi_path, 
                                            FILTER=FILTER,
                                            filename=filename, 
                                            savepath =save_path ,
                                         number_of_points =kws['number_of_points'],
                                         number_of_skin_depth=kws['number_of_skin_depth'], 
                                         dipole_length =kws['dipole_length'], 
                                         reduce_res_factor_x=kws['reduce_res_factor_x'], 
                                         reduce_res_factor_y = kws['reduce_res_factor_y'], 
                                        )
  
        except : 
            csamtpylog().get_csamtpy_logger().error(
                'something wrong when enter %s filter'%FILTER)
            
        else : 
            return FILTER 
        
    def test_write_corrected_edi(self):
        """
        Write corrected edi with differents filters 
        
        """
        
        # let collected the number of edifiles 
        input_edi_files = os.listdir(EDI_DATA_DIR)
        expected_edifiles= ['{0}.S{1:02}.edi'.format(survey_testname, nn ) for nn
                            in range(len(input_edi_files))]
    
        summary_outputs={}
        save_edi_corrected = os.path.join(TEST_TEMP_DIR, self.__class__.__name__)
        
        for ff, filter_used in enumerate(self._filters) :
            csamtpylog().get_csamtpy_logger().info(
                ' {0}-{1} applying to correct edi.'.format(ff+1, filter_used))
            try :
                successf=self.wrap_single_processing_filter(filter_used,
                                                   filename=survey_testname ,
                                                   edi_path=EDI_DATA_DIR, 
                                                   save_path= save_edi_corrected, 
                                                   **self.params_filters)
            
            except: 
                csamtpylog().get_csamtpy_logger().error(
                    'Problem occurs while setting {0} filter ' 
                    'to correct edifiles'.format(filter_used))
            else : 
            # compare outputfiles and expected files 
                expected_edifiles=[os.path.join(save_edi_corrected, expfile) 
                                   for expfile in expected_edifiles]
                self.assertEqual(len(os.listdir(save_edi_corrected)), len(expected_edifiles), 
                                 'reference ouputs  ={0} are different from'
                                 ' expected output = {1}.'.format(
                                     len(os.listdir(save_edi_corrected)), 
                                     len(expected_edifiles)))
                # if ok , compare the files 
                compare_diff_files(refout=[os.path.join(save_edi_corrected, outfile)
                                           for outfile in os.listdir(save_edi_corrected)] ,
                                   refexp = expected_edifiles)
                # if ok , then keep values 
                summary_outputs[successf]=(len(os.listdir(save_edi_corrected)),
                                              sorted([os.path.join(save_edi_corrected, file)
                                               for file in os.listdir(save_edi_corrected)]) )
                
        # compare all outputs if there are the same 
        self.assertEqual(len(self._filters), len(summary_outputs), 
                         'all filters are not applied corrected !Filters' 
                         'successfully run are {0}.'.format(list(summary_outputs.keys())))
        # compared outputs and numbers 
        # letcollected the number of edifiles 
        input_edi_files = os.listdir(EDI_DATA_DIR)
        expected_edifiles= ['{0}.S{1:02}.edi'.format(survey_testname, nn ) for nn
                            in range(len(input_edi_files))]
        expectedDICT= {ff:(len(expected_edifiles), [os.path.join(
            save_edi_corrected, expfile) for expfile in expected_edifiles])
                           for ff in self._filters}
     
        self.assertDictContainsSubset(
            expectedDICT, summary_outputs, 
            'Expected dict containers are seriously different from'  
            'reference output dictionnary')
        


if __name__=='__main__': 

    unittest.main()
        
        
        
        
        
        
        
        