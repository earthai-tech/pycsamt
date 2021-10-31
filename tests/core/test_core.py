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


import os
import datetime

import  unittest 

from pycsamt.ff.core import ( avg, edi, cs,j)
from tests.modeling.__init__ import (reset_matplotlib, csamtpylog, diff_files)
from tests import EDI_DATA_DIR, AVG_DATA_DIR, make_temp_dir 
from tests import J_DATA_DIR, TEST_TEMP_DIR , STN_DIR
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
                                        profile_fn=os.path.join(AVG_DATA_DIR, 'K1.stn'),
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
            print('--> outdir does not exist , set to None !')
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
                for edi_id in edi_head_id ]
                                        
            
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

    def test_write_j2edi(self):
        """
        Convert AG.Jones files *j to SEG EDI files. we will not test all jfiles 
        We select one file among the existing jfiles in jpath '/data/j' 
        
        """
        import random 
        
        jfiles = os.listdir(J_DATA_DIR)
        # we only selected one file from all jpath
        selected_j = random.sample(jfiles, 1)  # return a list 
        jf =os.path.join(J_DATA_DIR, selected_j[0])
        csamtpylog().get_csamtpy_logger().info('{0} is successufully selected '
                                                ' for a jedi test.'.format(jf))
        # now check whether j as a file  
        try : 
            os.path.isfile(jf)
        except : 
            csamtpylog().get_csamtpy_logger().error(
                '{0} is not a file. Please check your input file assumed '
                'to be a AG. Jones files *.dat'.format(jf))
        else : 
            save_outdir = os.path.join(TEST_TEMP_DIR, self.__class__.__name__)
            j.J_collection().j2edi(jfn=jf, savepath =save_outdir )
            
        # Now test whether the files is the right file generated
        # get the name of the file and created expected output file name 
        # remember Jfile is MT data file 
        exp_jname = os.path.join(save_outdir, 
                                  'S00_mt.{0}.edi'.format(
                                   datetime.datetime.now().year))
        
        #compare both files 
        refout = os.path.join(save_outdir, 
                              os.listdir(save_outdir)[0])
        compare_diff_files(refout=[refout], refexp=[exp_jname])
       

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
    
    
class TestProfile(unittest.TestCase):
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
    
    # get a list of stn file in STN_DIR
    list_of_station_profiles =[os.path.join(STN_DIR, pfile)
                                  for pfile in os.listdir(STN_DIR) if 
                                  pfile.endswith('.stn')]

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
        To create stn profile.
        
        """
        #create profile object
        refin, refout =[[] for i in range(2)]
        # assert if this file is really a readeable file.
        for psfile in sorted(self.list_of_station_profiles) : 
            try : 
                self.assertTrue(os.path.isfile(psfile) ,
                                'Ref input data must be a station profile file !.')
            except : 
                csamtpylog().get_csamtpy_logger().error(
                    'Station profile {} is not zonge .stn file.')
            else : 
                refin.append(psfile)
                
        #now create for each file a specila obj 
        try : 
            p_objs = [cs.Profile(profile_fn =ps_obj) for ps_obj in refin]
        except : 
            csamtpylog().get_csamtpy_logger().error(
                    'Something wrong happening when reading'
                    ' station profile source file ! ')
        # create outputnames to create expected files 
        outp_names = [survey_testname + os.path.basename(psfile)[:int(
                    len(os.path.basename(psfile))/2)] 
                    for psfile in refin]
        # then loop for all object 
        for pp, p_obj in enumerate(p_objs ):
            try : 
                p_obj.rewrite_station_profile (
                                easting = p_obj.east, 
                                northing= p_obj.north,
                                elevation = p_obj.elev, 
                                area_name =self.NEW_STNNAME, 
                                username =self.USERNAME, 
                                add_azimuth =self.compute_AZ, 
                                savepath =os.path.join(
                                    TEST_TEMP_DIR, 
                                    self.__class__.__name__) , 
                                output_name =outp_names[pp]
                                )
        
            except : 
                csamtpylog().get_csamtpy_logger().error(
                'Could not rewrite station profile ,it seems something wrong !'
                ' happened while reading `rewrite_station_profile` method.')

            #collect files read successfully
            refout.append(outp_names[pp])
           
        # collect expected files 
        save_outdir = os.path.join(TEST_TEMP_DIR, self.__class__.__name__)
        refout =[os.path.join(save_outdir, ouname) for ouname in refout]
           
        # controle whether all .stn files found are successfully converted as 
        # station profile_obj
        if len(refout) != len(refin):
            print('---> Size of `refin` and `ref out` are not the same.'
                  'refin size = {0} & ref out size = {1}.'.
                  format(len(refin), len(refout)))
            try :
                self.assertLess(len(refout), len(refin),
                                'Impossible to get number of .stn file read succeessfully '
                                'greater than the number of files collected.')
            except :
                csamtpylog().get_csamtpy_logger().error(
                    'Difference size found between inputs stations size ={}  and'
                    'expected outputs size ={}.'.format(len(refin), len(refout)))
    
        elif len(refin) == len(refout):
            print('---> ref input station file and ref output size'
                  ' match perfectly , is = {0}'.format(len(refout)))
                  
        #---> collect outputfiles and expected files 
        save_stn_files = [os.path.join(save_outdir, ofile) 
                          for ofile in os.listdir(save_outdir)] 
        
        expected_stn_files = [os.path.join(save_outdir, ouname +'.stn')
                              for ouname in refout] # expected files creating

        self.assertEqual(len(save_stn_files), len(expected_stn_files),
                              'Ref output and ref expected  must have the same size.'
                              'ref outsize is = %s while ref expect size '
                              'is =%s'%(len(save_stn_files),len(expected_stn_files)))
        
        # now compared differences between files 
        compare_diff_files(refout=save_stn_files, refexp=expected_stn_files)
        
    def test_scale_station_profile(self):
        """
        Test to scaled coordinates values from `x`and  `y`
        UTM coordinates and rewrite new scaled stn files.
        
        """
         
        scaled_utm__x_coordinates=-300238.702 
        scaled_utm__y_coordinates= -2369.252
        
        rewrite_scaled_stn_file =True 
        
        # collected files in stn directory 
        readablefiles, success_read_objs =[],[]
        for readable_file in sorted(self.list_of_station_profiles): 
            try : 
                self.assertTrue(os.path.isfile(readable_file),
                                ' %s is not a readable file !'% readable_file)
            except :
                csamtpylog().get_csamtpy_logger().error(
                    'Error finding readable profile file.!')
            else :
                #collect readable file 
                readablefiles.append(readable_file)
                
        #collect pObjs from readables files 
        try :
            pObjs=[cs.Profile(profile_fn=readfs) for readfs in readablefiles]
        except :
            csamtpylog().get_csamtpy_logger().error(
                'Error collection profile objs from readable stn files collected.'
                "It's obvious some readables files are not zonge station profile file.")

        # read readfiles and save it to outdir 
        tem_save_files = os.path.join(TEST_TEMP_DIR, self.__class__.__name__)
        
        for ii, pr_obj in enumerate(pObjs):
            try : 
                pr_obj.reajust_coordinates_values(x=scaled_utm__x_coordinates ,
                                                  y=scaled_utm__y_coordinates , 
                                                stn_fn=readablefiles[ii] ,
                                                rewrite = rewrite_scaled_stn_file, 
                                                savepath=tem_save_files )
            
            except :csamtpylog().get_csamtpy_logger().error(
                    'Problem occurs while trying to scaled sites coordinates !')
            else :
                success_read_objs.append(readablefiles[ii])
        
        if len(success_read_objs) != len(readablefiles):
            self.assertLess(len(success_read_objs), len(readablefiles), 
                            'In fact obj read must be greater than'
                            ' readables files collected!')
        elif len(success_read_objs) == len(readablefiles):
            print('All readable files successfully read !')
            
            
        # create expected files only 
        expected_files =[os.path.join(tem_save_files, 
                                os.path.basename(stn_fn).split('.')[0] + '_reaj.stn')
                                for stn_fn in success_read_objs]
        scaled_stn_files= [os.path.join(tem_save_files, file) 
                           for file in os.listdir(tem_save_files)
                           if file.endswith('_reaj.stn')]
        
        #now compare  difference between expected files & outputfiles
        self.assertEqual(len(scaled_stn_files), len(expected_files),
                         'Difference sizes found between'
                         ' ref outs = {0} & expected files = {1}'.
                         format(len(scaled_stn_files),len(expected_files)))
        
        compare_diff_files(refout=scaled_stn_files, refexp=expected_files)

if __name__=='__main__':

    unittest.main()
    
    
                
    
    
    
    
    
    
    
    
    
    
    