# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 13:27:27 2021

    Description:
        Test to core module . Containers of modules '`avg`, `cs` , `j`, `edi`'
        Test ouputfiles from rewriting and generating files , which includes:
        Reference input data are from AVG_DATA_DIR, EDI_DATA_DIR.
        ('K1.AVG',    'K1.stn'),    'edi')
    
    References:
        scripts $ /read_occam2d_res_model.py
        .. _module-core::`pycsamt.ff.core`
        
@author: @Daniel03

"""
from tests.modeling.__init__ import reset_matplotlib, csamtpylog, diff_files

import os 

import  unittest 

from pycsamt.ff.core import ( avg, edi, cs)

from tests import EDI_DATA_DIR, AVG_DATA_DIR, make_temp_dir,TEST_TEMP_DIR , STN_DIR 
from tests import survey_testname


class TestAVG(unittest.TestCase):
    """
    Test avg module for Reading and rewriting files 
        *. write zonge astatic file type2 to zonge avg plainty file (type1) with 
            eddor propagations
        * write avg to j files 
        * write avg to edifiles .
    """
    

    @classmethod 
    def setUpClass(cls):
        """
        Reset building matplotlib plot and generate tempdir inputfiles 
        
        """
        reset_matplotlib()
        cls._temp_dir = make_temp_dir(cls.__name__)

    def setUp(self): 
        
        if not os.path.isdir(TEST_TEMP_DIR):
            print('--> outdir not exist , set to None !')
            csamtpylog().get_csamtpy_logger().error('Outdir does not exist !')
        
    
    def test_write_avgf2tof1(self):
        """
        Test write avg file 2 (Astatic file) to avg file 1 plainty file with 
        error propagations. 
        will check in avg outidr ,differents AVG files and will write only
        astatic files.
        
        """
 
        avg_files =[os.path.join(AVG_DATA_DIR,file) for file 
                    in  os.listdir(AVG_DATA_DIR) if file.endswith('.AVG')]
        
        # control  files if exists 
        cc=0
        for afile in avg_files : 
            self.assertTrue(os.path.isfile(afile), 
                            "It seems file  ={0} does not exist !".
                            format(os.path.basename(afile)))
            cc +=1 
        if cc >1 : verb ='are','s' 
        else :verb ='is',''
        
        print('---> Number of `avg` file%s'
              ' collected %s = %s.'% (verb[1], verb[0], cc))
        
        tem_asfiles =[]     # collected astatic files 
        for file in avg_files :
            try:
                avg_obj = avg.Avg().avg_write_2_to_1(data_fn= file , 
                                   savepath = os.path.join(TEST_TEMP_DIR,
                                                           self.__class__.__name__))
            except :
                # assume that file is not zonge 
                try: 
                    if avg_obj.checker ==1 : 
                        csamtpylog().get_csamtpy_logger().info(
                            'File supposed to read and write is not ASTATIC file ! ')
                except:
                    csamtpylog().get_csamtpy_logger().error(
                        'It seems errors occurs while reading `avg`files !')
            else :
                # readinf fine 
                tem_asfiles.append(file)
                
        print('--> Finally number of `astatic` file found  and successfully'
              ' read is {0}/{1} files at all.'.format(len(tem_asfiles), len(avg_files)))
              
        # get expected filenames from astatic files
        expected_files = [''.join([os.path.basename(temfile)[:-4],'_2_to_1.avg'])
                          for temfile in tem_asfiles]
        
        #now check the outfiles 
        save_outdir = os.path.join(TEST_TEMP_DIR, self.__class__.__name__)
        expected_files =[os.path.join(save_outdir, file) for file in expected_files]
        # check the number of output files and the number of expected files 
        
        self.assertEqual(len(expected_files), len(os.listdir(save_outdir)), 
                         'Different size ascertained with expected files and output files.'
                         'Expected files and output files size are = {0} & {1} respectively.'.
                         format(len(expected_files), len(os.listdir(save_outdir))))
    
        
        for outfile, exfile in zip(sorted(os.listdir(save_outdir)), sorted(expected_files)): 
            
            self.assertTrue(os.path.isfile(os.path.join(save_outdir,outfile)),
                                "Ref output data file does not exist, nothing to compare with."
                                )
        
            print(("Comparing", os.path.join(save_outdir,outfile), "and", exfile))
            
            is_identical, msg = diff_files(os.path.join(save_outdir,outfile),
                                           exfile, ignores=['Date/Time:'])
            print(msg)
            self.assertTrue(is_identical, "The output file is not the same with the baseline file.")
       
    def test_write_edi(self):
        """
        Test to write edi files `emap` or `mt` data file,  Will test write EDI using
        specific `K1.AVG` and `K1.stn` files. with no filter applying
        will compare number of edifiles write to number of station embedded in 
        `zonge avg`file. 
        
        """
        # defile parameter 
        try :
            avg_obj = avg.Avg()
            avg_obj.avg_to_edifile(data_fn= os.path.join(AVG_DATA_DIR,'K1.AVG'),
                                   profile_fn = os.path.join(AVG_DATA_DIR, 'K1.stn'), 
                                   savepath =os.path.join(
                                       TEST_TEMP_DIR, self.__class__.__name__), 
                                   reference_frequency= None, 
                                   apply_filter='tma', 
                                   number_of_points=7, 
                                   number_of_skin_depth=3.) 
        except : 
            
            csamtpylog().get_csamtpy_logger().error(
                'It seems something wrong happened when calling `avg_to_edifile` method. ')
        else :
            # collect list of sations 
            expected_files = avg_obj.Data_section.Station.names
      
            outputdir = [file for file in 
                os.listdir(os.path.join(TEST_TEMP_DIR, self.__class__.__name__)) 
                if file.endswith('.edi')]
            
        self.assertEqual(len(expected_files), len(outputdir), 
                         'Different size occured in expected output edifiles.'
                         'Expected edi size is  = {0} rather than the number'
                         '  of survey sites = {1}.'.format(len(outputdir), 
                                                                   len(expected_files)))
            
        print('Number of survey sites  is ={0} ,is identical with number'
              '  of output edi files'.format(len(expected_files)))
            
    def test_write_j(self):
        """
        write j file from avg files . use default extension j as , '*.dat'
        Test will set 'write avgfile infos to True so check whether all the codes 
        qre well written .Can change the `j` output . *Default*  is`.dat`.
        """
       
        try :
            avg_obj =avg.Avg()
            avg_obj.avg_to_jfile(avg_data_fn=os.path.join(AVG_DATA_DIR,'K1.AVG'), 
                                        station_fn=os.path.join(AVG_DATA_DIR, 'K1.stn'),
                                        j_extension='.dat',
                                        savepath=os.path.join(
                                            TEST_TEMP_DIR, self.__class__.__name__),
                                        writeInfos=True)
        except : 
            csamtpylog().get_csamtpy_logger().error(
                'Unable to write `avg` file to `j(*.dat)` format.'
                'Something wrong happened !')
        else : 
            #collect number of stations and compared to output files 
            sites_names = avg_obj.Data_section.Station.names
            
        expected_j_files = [jfile for jfile in 
                os.listdir(os.path.join(TEST_TEMP_DIR, self.__class__.__name__)) 
                if jfile.endswith('.dat')]
        
        self.assertEqual(len(expected_j_files), len(sites_names), 
                         'Different size occured in expected output edifiles.'
                         'Expected `j(*.dat)` size is  = {0} rather than the number'
                         '  of survey sites = {1}.'.format(len(expected_j_files), 
                                                          len(sites_names)))
            
        print('Number of survey sites  is ={0} ,is identical with number'
              '  of output edi files.'.format(len(expected_j_files)))
  
class TestEDI(unittest.TestCase):
    """
    Test edifile  write edifiles 
    
    """  
    @classmethod 
    def setUpClass(cls):
        """
        Reset building matplotlib plot and generate tempdir inputfiles 
        
        """
        reset_matplotlib()
        cls._temp_dir = make_temp_dir(cls.__name__)

    def setUp(self): 
        
        if not os.path.isdir(TEST_TEMP_DIR):
            print('--> outdir not exist , set to None !')
            csamtpylog().get_csamtpy_logger().error('Outdir does not exist !')
        
    def test_rewrite_edi(self):
        """
        Rewrite edi test.
        `data_section` param  could be `emap` or `mt` detect automaticcally.
        set to false to let detect automatically. To force rewrittying edi into 
        a specific format , set`data_section` param to the corresponding data type.
        
        """

        list_of_edifiles = [os.path.join(EDI_DATA_DIR, edifile) 
                   for edifile in os.listdir(EDI_DATA_DIR) if edifile.endswith('.edi')]

        edi_head_id , edi_type = [[] for i in range(2)] # initialise empty list to 
        # collect edi head id as well as edi type  
        try :
            for edi_file in list_of_edifiles :
                self.assertIsInstance(edi_file, str, '`edi file` must be a string.'
                                      'Unable to recognize this file as an OS.Path_Like obj!.')
                
                edi_obj =edi.Edi(edi_filename= edi_file )
                edi_obj.write_edifile(savepath =os.path.join(
                                                TEST_TEMP_DIR, self.__class__.__name__),
                                      datatype=None,
                                      new_edifilename =survey_testname)
                edi_head_id.append(edi_obj.Head.dataid)
                edi_type.append(edi_obj.typefile)
                
        except : 
                csamtpylog().get_csamtpy_logger().error('An error occurs while reading and'
                                                        'Rewritting `edi file.')
        else : 
            # now  get the number of edi head id and create expcted
            # generate edi output filename 
            expected_edi_output_files = [os.path.join(os.path.join(
                                                TEST_TEMP_DIR, self.__class__.__name__),
                        survey_testname + '.{0}.edi'.format(edi_id))
                for edi_id in edi_head_id
                                        ]
            
        # now collect edioutfiles generated 
        save_edi_outdir =os.path.join(TEST_TEMP_DIR, self.__class__.__name__)
        save_edi_files = [os.path.join(save_edi_outdir, edi)
                          for edi in os.listdir(save_edi_outdir)]
        
        # fist check assert if number generate edifiles is the same than input files 
        self.assertEqual(len(list_of_edifiles), len(save_edi_files), 
                         'Difference found in expected edi output size = {0}'
                         ' with inpufiles size ={1}. '.
                         format(len(save_edi_files), len(list_of_edifiles)))

        try : 
            if all(edi_type) is True: 
                if edi_type [0] =='mt' : 
                    resp ='magnetotelluric <{0}> data'.format(
                        edi_type[0].upper()+'SECT')
                elif edi_type[0] =='emap':
                    resp= 'electromagnetic profiling array '+\
                    '<{0}> '.format(edi_type[0].upper()+'SECT')
                    
                mess ='from the same data type = {0}'.format(resp)
                
                print('---> {0} Survey edis are deeply tested and all  are '
                      ' from = {1} type.'.format(len(list_of_edifiles), mess))
            else :
                mess = 'Different edis type found while testing'+\
                    ' {0} edis files.'.format(len(list_of_edifiles))
                print('---> '+ mess)
        except : 
            csamtpylog().get_csamtpy_logger().info(
                'It seems different data type not  found become'
                ' something wrong happened while checking!')
 
        for iedif , outedif in zip(sorted(expected_edi_output_files), 
                                   sorted(save_edi_files)):
            self.assertTrue(os.path.isfile(outedif),
                                "Ref output data file does not exist,"
                                "nothing to compare with"
                                )
            
            print(("Comparing", iedif, "and", outedif))
            
            is_identical, msg = diff_files(iedif, outedif, ignores=['Date/Time:'])
            print(msg)
            self.assertTrue(is_identical, 
                            "The output file is not the same with the baseline file.")


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
    
    
class TestPROFILE(unittest.TestCase):
    """
    Test profile. Two specific test is proposed 
    1. rewrite station profile 
    2. readjustment or rescale station profile . 
        STN_DIR is '/data/stn_profiles'
    
    """
    # Define some optional params 
    ELEV = None  # Dont because `stn` file already contained elevation values 
        # if given must be (ndarray,1)  (optional params)   
    compute_AZ = True   # is value is` TRUE`,will add  azimuth computation. Most 
                        # likely 'stn' file doesnt not contain  azimuth value.
                        # we set to True  for testing all functionnalities. 
    USERNAME = 'BomoAG.'
    NEW_STNNAME =survey_testname [:-4] +'civ' #output new coordinates values 
    

    @classmethod 
    def setUpClass(cls):
        """
        Reset building matplotlib plot and generate tempdir inputfiles 
        
        """
        reset_matplotlib()
        cls._temp_dir = make_temp_dir(cls.__name__)

    def setUp(self): 
        
        if not os.path.isdir(TEST_TEMP_DIR):
            print('--> outdir not exist , set to None !')
            csamtpylog().get_csamtpy_logger().error('Outdir does not exist !')
        
   
    def test_create_station_profile(self):
        """
        to create stn profile
        
        """
        #create profile object
        refin, refout =[[] for i in range(2)]
        list_of_station_profiles =[os.path.join(STN_DIR, pfile)
                                   for pfile in os.listdir(STN_DIR) if 
                                   pfile.endswith('.stn')]
        
        for psfile in list_of_station_profiles : 
            
            self.assertTrue(os.path.isfile(psfile) ,
                            'Ref input data must be a station profile file !.')
            try :
                profile_obj =cs.Profile(profile_fn= psfile)
            except :
                csamtpylog().get_csamtpy_logger().error(
                    'Something wrong happening when reading profile object ! ')
            else :
                outputname = survey_testname + os.path.basename(psfile)[:int(
                    len(os.path.basename(psfile))/2)]
                try : 
                    
                    profile_obj.rewrite_station_profile ( easting = profile_obj.east, 
                                                         northing= profile_obj.north,
                                                         elevation = profile_obj.elev, 
                                                    area_name =self.NEW_STNNAME, 
                                                    username =self.USERNAME, 
                                                    add_azimuth =self.compute_AZ, 
                                                    savepath =os.path.join(
                                                        TEST_TEMP_DIR, 
                                                        self.__class__.__name__) , 
                                                    output_name =outputname)
    
                except : 
                    csamtpylog().get_csamtpy_logger().error(
                    'Could not rewrite profile_obj !')
                else :
                    #collect files read successfully
                    refin.append(psfile)
                    
                refout.append(outputname)
                
        try : 
            self.assertEqual(len(refin), len(refout),
                             'Size of `refin` and `ref out` are not the same.')
        except :
            csamtpylog().get_csamtpy_logger().error(
                'Diffennce size found betwen input stations size ={} files and'
                'Expected output files with size ={}.'.format(len(refin), len(refout)))
        else :
            print('---> ref input station file and ref output size'
                  ' match perfectly , is = {}'.format(len(refout)))
                  
        # collect expected files 
        save_outdir = os.path.join(TEST_TEMP_DIR, self.__class__.__name__)
        
        save_stn_files = [os.path.join(save_outdir, ofile) 
                          for ofile in os.listdir(save_outdir)]  
        # let create expected stn files 
        expected_stn_files = [os.path.join(save_outdir, ouname +'.stn')
                              for ouname in refout]
        self.assertListEqual(save_stn_files, expected_stn_files,
                              'Ref output and ref expected list must have the same size.'
                              'ref outsize is = %s while ref expectsize '
                              'is =%s'%(len(save_stn_files),len(expected_stn_files)))
        
        # now compared files 
        compare_diff_files(refout=save_stn_files, refexp=expected_stn_files)
        
if __name__=='__main__':

    unittest.main()
                
    
    
    
    
    
    
    
    
    
    
    