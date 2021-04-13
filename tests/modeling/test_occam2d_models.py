# -*- coding: utf-8 -*-
# ! /usr/bin/env python
"""
Description:
    Test create input files required to read  occam2d inversion model files, which includes:
    ('Occam2DMesh',    'Occam2DModel',    'Startup', 'OccamDataFile.dat')
    These input files are created from standard edi data file set.

References:
    scripts/read_occam2d_res_model.py

CreationDate:   25/03/2021
Developer:      etanoyau@gmail.com/kkouao@zju.edu.cn

    
"""
from tests.modeling.__init__ import reset_matplotlib, csamtpylog, diff_files

import os , datetime

import  unittest 

from pycsamt.modeling.occam2d import occam2d_write,Iter2Dat 


from tests import EDI_DATA_DIR, OC2D_DIR, make_temp_dir,TEST_TEMP_DIR 

class TestWriteOccam2DInputfiles(unittest.TestCase):
    """
    Test to generate Occam2D inputfiles such as 
    `Startup` , 'Occam2DMesh', 'Occam2DModel', 'OccamDataFile.dat' files.
    pyCSAMT try to intall `MTpy` packages if not exists'. If an errors occurs, 
    may try to install `MTpy` mannually. 
    
    """

    _expected_edi_dir = EDI_DATA_DIR 
    _expected_save_dir = TEST_TEMP_DIR 
    geo_electric_strike = 34.
    
    
    @classmethod 
    def setUpClass(cls):
        """
        Reset building matplotlib plot and generate tempdir inputfiles 
        
        """
        reset_matplotlib()
        cls._temp_dir = make_temp_dir(cls.__name__)

    def setUp(self): 
        
        if not os.path.isdir(self._expected_edi_dir):
            print('--> outdir not exist , set to None !')
            csamtpylog().get_csamtpy_logger().error('Outdir does not exist !')
        

    def test_occam2d_buildingInputfiles(self): 
        """
        Test occam2d inputbuilding files . 
        First check if code succeffully run and secundly check whether the 
        expected files much the number of files generated as  well as 
        the ouputfiles names .
        
        """
        for file  in [self._expected_edi_dir , self._expected_save_dir]: 
            self.assertIsInstance(file, str,
                                  'Expected paths must be PathLike str !')
            
        try : 
            stations_list = [os.path.join(self._expected_edi_dir, file) 
             for file in os.listdir(self._expected_edi_dir) if 
                             self._expected_edi_dir.endswith('.edi')]
        except : 
            csamtpylog.get_csamtpy_logger().error( 
                'something wrong when call edifiles in %s' % self._expected_edi_dir )
        else : 
            print('...> {0} stations successfully  read'.format(len (stations_list)))
            csamtpylog.get_csamtpy_logger().info( 
                'Now write occam2d inputfuiles in %s' % self._expected_edi_dir )
            occam2d_write.buildingInputfiles(edi_fn= self._expected_edi_dir,
                                   savepath =os.path.join(self._expected_save_dir,
                                                          self.__class__.__name__), 
                                   geoelectric_strike= self.geo_electric_strike)
            
        # now assert expected files and output files 
        # check whether the tem dir files if create
        self.assertTrue(os.path.exists(os.path.join(TEST_TEMP_DIR , self.__class__.__name__ )),
                                       'Temp dir path does not exists')
        
        # now control existing files 
        out_tem_dir = os.path.join(TEST_TEMP_DIR ,self.__class__.__name__ )

        output_data_file =  [os.path.join(out_tem_dir, file) for file in os.listdir(out_tem_dir)]
    
        output_data_file = sorted( output_data_file) 
        
        self.assertEqual(len(output_data_file), 4,
                         'Number of Occam2D inputfiles  built is %s ,'
                         ' far from expected 4 files.'%len(output_data_file))
        
        # now check the name of generate files and compared to expected files 
        for ii, afile in enumerate( sorted(('Occam2DMesh', 'Occam2DModel',
                                            'OccamDataFile.dat', 'Startup'))):
            expected_data_file = os.path.join(out_tem_dir, afile) 
            
            self.assertTrue(os.path.isfile(expected_data_file),
                                "Ref output data file does not exist, nothing to compare with"
                                )
        
            print(("Comparing", output_data_file[ii], "and", expected_data_file))
            
            is_identical, msg = diff_files(output_data_file[ii], 
                                           expected_data_file, ignores=['Date/Time:'])
            print(msg)
            self.assertTrue(is_identical, "The output file is not the same with the baseline file.")
          
class OtherModelingTest(unittest.TestCase) :
    """
    Others modeling tests. Collect specifiles occam model files names. 
    inputfiles are (ITER17.iter,LogFile.logfile, 'Occam2DMesh',
                    'Occam2DModel','Startup', 'OccamDataFile.dat', 'RESP17.resp' )

    """

    _expected_save_dir = TEST_TEMP_DIR 
    
    oc2testfiles =[os.path.join(OC2D_DIR, file) 
                   for file in os.listdir(OC2D_DIR)]
    
    @classmethod 
    def setUpClass(cls):
        """
        Reset building matplotlib plot and generate tempdir inputfiles 
        
        """
        reset_matplotlib()
        cls._temp_dir = make_temp_dir(cls.__name__)

    def setUp(self): 
        
        if not os.path.isdir(self._expected_save_dir):
            print('--> outdir not exist , set to None !')
            self._expected_save_dir =None 
            
    def ascertain_oc2d_inputfiles(self):
        """
        Assert existing oc2d files in OC2D_DIR . 
        
        """
        # check number of existing occam2d files 
        self.assertEqual(len(self.oc2testfiles),7,
                         'Existing files number is ={0}, Expected 7 files.')
        # checker whether the files joined are effective 
        for file in self.oc2testfiles : 
            self.assertTrue(os.path.isfile(file), 'Symlink for occam2dfile is wrong ')
        
        # collect occamfiles and populate attraibutes 
        for ifile in ('ITER17.iter', 'Occam2DMesh','Occam2DModel','LogFile.logfile',
                      'OccamDataFile.dat', 'Startup', 'RESP17.resp' ):
            for  tfile in self.oc2testfiles : 
                if tfile.lower().find(ifile.lower())>=0 :
                    self.assertTrue(os.path.isfile(tfile),
                                    'Symlink for occam2dfile is wrong ')
                    if os.path.basename(tfile).lower().find('mesh')>=0:
                        self.oc2d_mesh = tfile 
                    elif os.path.basename(tfile).lower().find('dataf')>=0:
                        self.oc2d_data = tfile
                    elif  os.path.basename(tfile).endswith('iter'): 
                        self.oc2d_iter = tfile 
                    elif os.path.basename(tfile).lower().find('model')>=0:
                        self.oc2d_model = tfile
                    elif 'startup' in os.path.basename(tfile).lower(): 
                        self.oc2d_startup = tfile 
                    elif os.path.basename(tfile).endswith('.resp')>=0:
                        self.oc2d_resp = tfile
                    elif os.path.basename(tfile).lower().find('logf')>=0:
                        self.oc2d_logfile = tfile
                        
                    break
                    
    def test_write_xyz_modelfile(self):
        """
        Write x,y,z  model file called Bo Yang file. 
        Expected 2 files iterxx.yy.bln & iterxx.yy.dat 

        """
        
        self.ascertain_oc2d_inputfiles()
        
        try : 
            i2d_obj = Iter2Dat(mesh_fn=self.oc2d_mesh, 
                iter_fn = self.oc2d_iter, 
                model_fn =self.oc2d_model, 
                data_fn =self.oc2d_data)
        except :
            csamtpylog.get_csamtpy_logger().error( 
            'Something wrong when initializing iter2dat object.' )
        else :
            # collect iteration number  and iterroughness to generate filename 
            # if file is not provided
            self.iternum = i2d_obj.iter_num 
            self.iter_roughness= i2d_obj.iter_roughness
            try :
                i2d_obj.write_iter2dat_file(filename =None,
                                   doi=1000., # or 1km
                                   savepath=os.path.join(self._expected_save_dir,
                                                         self.__class__.__name__ )
                                   )
            except :
                  csamtpylog.get_csamtpy_logger().error( 
                    'something wrong when writing model'
                    ' `x`, `y`,`z` files. May check either arguments params.')
            else : 
                #check new outdir , if exists.
                self.assertTrue(os.path.isdir(
                    os.path.join(self._expected_save_dir, self.__class__.__name__ )), 
                    'outdir doesnt not exist, should be created automatically.')
            
        output_xyz_dir = os.path.join(self._expected_save_dir, 
                                      self.__class__.__name__ )

        output_xyz_files = sorted([os.path.join(output_xyz_dir, file)
                                   for file in os.listdir(output_xyz_dir)])

        # normally expected file should be 2
        self.assertEqual(len(output_xyz_files ), 2,
                     'Number of "xyz" model files  built is %s ,'
                     ' far from expected 2 files.'%len(output_xyz_files))
    
 
        for iterf in [self.iternum , self.iter_roughness]:
            self.assertIsNotNone(iterf,
                              'Iteration number or iteration roughness is None'
                              'Unable to use to generate xyz filename.')
            
        filename ='iter{0}.{1}{2}'.format(int(self.iternum), self.iter_roughness,
                                  datetime.datetime.now().month)
        
        expected_files =[os.path.join(output_xyz_dir , file) 
                         for file in (filename +'.dat', filename +'.bln')]
         # ow looping files generated and tech it with expected extension files 
        for ii, (afile, efile) in enumerate(zip(sorted(output_xyz_files), 
                                       sorted(expected_files))):#('bln', 'dat'):
            self.assertTrue(os.path.isfile(afile),
                                "Ref output data file does not exist, nothing to compare with"
                                )
            
            print(("Comparing", afile, "and", efile))
            
            is_identical, msg = diff_files(afile, 
                                           efile, ignores=['Date/Time:'])
            print(msg)
            self.assertTrue(is_identical, "The output file is not the same with the baseline file.")
            
    
            
            
            
if __name__=='__main__':
    # ttt= OtherModelingTest()
    # ttt.test_write_xyz_modelfile()
    unittest.main()
    
   