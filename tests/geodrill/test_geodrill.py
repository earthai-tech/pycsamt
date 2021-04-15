# -*- coding: utf-8 -*-
"""
Description:
    Test to write new_model-resistivity  with oasis- montaj, golden software 
   and final test to write DrillHole data 
    the data type either : MT or EMAP  and apply filters like `ama` , `flma` 
    `tma`, `ss` or `dist`. 
    {`ama` : adaptative moving average, 
    `flma`: fixed-length-dipole moving average 
    `tma`: trimming moving average 
    `ss` : static shift removal  and 
    `dist`: distorsion removal}
Created on Wed Apr 14 16:14:00 2021

references modules : `pycsamt.geodrill.geoCore.geodrill.Geodrill`
                    :members:`.to_oasis_montaj`
                             `.to_golden_software`
                    `pycsamt.geodrill.geoCore.geodrill.Drill`
                    `pycsamt.geodrill.geoCore.geodrill.Drill`.Geosurface`
                    
Refencences :
    Scripts: /write_geo_model_files.py
        /write_geo_model_files_to_oasis.py
        /write_geosurface_files.py
        /write_borehole_data.py

@author: @Daniel03

"""

from tests.modeling.__init__ import csamtpylog, reset_matplotlib
from tests.processing.__init__ import compare_diff_files
from tests.geodrill.__init__ import remove_control
import os, datetime

import  unittest
import pytest 

from pycsamt.geodrill.geoCore.geodrill import (Geodrill, Drill, Geosurface)

from tests import  make_temp_dir, TEST_TEMP_DIR ,I2DAT_DIR, OC2D_DIR, AVG_DATA_DIR
from tests import OAS_DIR, DRILL_PARSER_DIR

from tests import survey_testname


class TestGEODRILL(unittest.TestCase):
    """
    Writing new resistivity model by adding new resistivity values and 
    new layers names . 
    Use both option to outputfiles and to check whether the script work as well
    """
    #set the main parameters 
    
    main_geo_params={
        'doi':'1km',
        'step_descent': 200., 
        'input_resistivities' : [66.,70., 180., 1000., 3000., 10000., 20000.], 
        'input_layers' : ['river water', 
                          'fracture zone' ,
                          'granite ',
                          'altered rock', 
                          'granite']
                    }
    
    #-> first option : use occam2d modelfiles   ---------
    oc2d_inversion_kwargs={'data_fn':'OccamDataFile.dat', 
                      'mesh_fn': 'Occam2DMesh', 
                      'model_fn':'Occam2DModel' , 
                      'iter_fn': 'ITER17.iter', 
                      }
    for key, modef in oc2d_inversion_kwargs.items(): 
        oc2d_inversion_kwargs[key]= os.path.join(OC2D_DIR, modef)
        
    # second option :use x,y,z model files 
    i2d_files_kwargs={
                    'iter2dat_fn' : os.path.join(I2DAT_DIR,'K1.iter.dat'),#'iter17.2412.dat'),
                    'bln_fn': os.path.join(I2DAT_DIR,'K1.bln')#'iter17.2412.bln')
                    }

    station_profile =os.path.join(AVG_DATA_DIR, 'K1.stn')
    
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
            
            
    def create_geo_obj(self, **kws): 
        """
        kws dict can be either `oc2d_inversion_kwargs` parameters or
        i2d_files_kwargs parameters to figure out the two option of files 
        assumed to be read by geodrill module.
        create geo obj for two members `to_oasis` and `to_golden_software`
        return model_object. 
        
        """
        # create saveoutdir 
        self.save_geo_outdir =os.path.join(TEST_TEMP_DIR, 
                                           self.__class__.__name__)

        if 'model' in list(kws.keys()):
            csamtpylog().get_csamtpy_logger().info(
                'Creating ocam2d objects from OCCAM2D model files')
            
        elif 'iter2dat' in list(kws.keys()):
            csamtpylog().get_csamtpy_logger().info(
                'Creating i2d object from iter2dat `x`, `y`, `z` model files '
                )
        try : 
            geo_obj =  Geodrill( **kws)  
        except : 
            csamtpylog().get_csamtpy_logger().error(
                'Trouble occurs while setting geo-model'
                ' object from `Geodrill` module.Try to fix it'
                )
        else : return geo_obj 
        
    
    def test_to_geolden_software(self):
        """
        Test writing oasis montaj new model with new resistivity files
        should return 4 outputs:
            * `nibykro_survey.{month}_aver.dat`-model averaged 
            * `nibykro_survey.{month}_rr.dat` :-roughness model 
            * `nibykro_survey.{month}_sd.dat` :step descent model
            * `nibykro_survey.{month}_bln.rr`: station location file.
    
        """

        for KWS, scale, option_depth in zip([self.oc2d_inversion_kwargs , self.i2d_files_kwargs], 
                                               ['m', 'km'], 
                                               [True, False]):
            try :
                geo_obj = self.create_geo_obj(**KWS, 
                                              **self.main_geo_params)

                geo_obj.to_golden_software(filename =survey_testname , 
                                           savepath = self.save_geo_outdir, 
                                           scale =scale, 
                                           to_negative_depth=option_depth) 
            except : 
                if 'model' in list(KWS.keys()): 
                    msg='Occam2d inputfiles'
                else :msg ='x,y,z model files'
                
                csamtpylog().get_csamtpy_logger().error(
                    'It seems there is something to fix in {}.'
                    ' Could not read propertly files.'.format(msg))
            else :
                #collected expected files , must be located in 
                gf_expectedfiles =[os.path.join(self.save_geo_outdir, gff) 
                                   for gff in ['{0}.{1}{2}'.format(
                        survey_testname, datetime.datetime.now().month,index) 
                    for index in ['_aver.dat',
                                  '_rr.dat',
                                  '_sd.dat', 
                                  '_yb.bln']]
                                    ]
                
                gf_outputfiles =[os.path.join(self.save_geo_outdir, outgff) for outgff in 
                                 os.listdir(self.save_geo_outdir) 
                                 if (outgff.endswith('.dat') or outgff.endswith('.bln')) ]
                #now compare files 
                self.assertEqual(len(gf_expectedfiles), len(gf_outputfiles), 
                                 'Different size found between reference output ={0} '
                                 ' & expected files = {1}.'.
                                 format(len(gf_outputfiles),len(gf_expectedfiles) ))
                
                compare_diff_files(refout=gf_outputfiles,
                                   refexp=gf_expectedfiles )
                
    def test_to_oasis_montaj(self) :
        """
        Test model outfiles files wrinte in excel sheet in xls.
        output extension files : `xlsx or `csv`. 
        coordinates scaled (X,Y):`(0., 0.), (-300238.702 , -2369.252)`
        export to log10 resistivities :Test (`True` & `False`)
        inputfiles :Test i2dfiles and occamfiles.
        return : K1.stn.2.main._cor_oas for xlsx  and 
                : 
        
        """
        for KWS, xy_corr, log10rho, file_type in zip(
                [self.oc2d_inversion_kwargs , self.i2d_files_kwargs], 
                                               [(0., 0.), (-300238.702 , -2369.252)],            
                                               [True, False], 
                                               ['xlsx', 'csv']):

            try : 
                geo_obj = self.create_geo_obj(**KWS, 
                                              **self.main_geo_params)
    
                geo_obj.to_oasis_montaj (profile_fn =self.station_profile, 
                                 to_negative_depth = log10rho, 
                                 savepath =self.save_geo_outdir, 
                                 filename =survey_testname , 
                                 scalled_east_north=(xy_corr[0],xy_corr[1]), 
                                 to_log10 =log10rho, 
                                 writeType=file_type, 
                                 output_s_XY=False) # not output scale file 
            except : 

                if 'model' in list(KWS.keys()): 
                    msg='Occam2d inputfiles'
                else :msg ='x,y,z model files'
                
                csamtpylog().get_csamtpy_logger().error(
                    'It seems there is something to fix in {}.'
                    ' Could not read propertly files.'.format(msg))
            else :
            #collected expected files when extension is xlsx
                if file_type =='xlsx':
                    oas_expectedfiles =[os.path.join(self.save_geo_outdir, 
                                            '{0}.{1}.main._cor_oas.xlsx'.format(
                                    survey_testname,datetime.datetime.now().month)
                                    )]
                    oas_outputfiles =[os.path.join(self.save_geo_outdir, ff) 
                                      for ff in os.listdir(self.save_geo_outdir)
                                      if ff.endswith('_oas.xlsx')]
                   
                    
                elif file_type =='csv':
                    
                    oas_expectedfiles =[os.path.join(self.save_geo_outdir, oasfile) 
                                        for oasfile in ['{0}.{1}{2}'.format(
                                                survey_testname,
                                                datetime.datetime.now().month, efile)
                                        for efile in ['_sd_cor_oas.csv', 
                                                     '_aver_cor_oas.csv', 
                                                     '_rr_cor_oas.csv', 
                                                     '.main._cor_oas.csv']]
                                        ]
                    oas_outputfiles = [os.path.join(self.save_geo_outdir, csvfile)
                                       for csvfile in os.listdir(self.save_geo_outdir)
                                       if csvfile.endswith('_oas.csv')]
                    
                    compare_diff_files(refout = oas_outputfiles ,
                                       refexp= oas_expectedfiles )
                    
                self.assertEqual(len(oas_expectedfiles), len(oas_outputfiles), 
                                 'Difference raises between expected files = {0} &'
                                 ' reference outputs = {1}.'.format(
                                     len(oas_expectedfiles), len(oas_outputfiles)))
                    

    def test_geosurface(self): 
        """
        Test build geosurface map at (38m, 100.) depth
        outputfiles is in `xlsx` and `csv`.
        """
                                
        # call geosurface object 
        self.save_geo_outdir =os.path.join(TEST_TEMP_DIR, 
                                           self.__class__.__name__)
        depths_to_image =[38., 100.]
        
        csamtpylog().get_csamtpy_logger().info(
            'Reading to export geosurface data at {0}m & {1}'.format(*depths_to_image) )
        
        try :
            csamtpylog().get_csamtpy_logger().info(
                'Trigging  geosurface object!')
            geo_surface_obj = Geosurface( path =OAS_DIR , 
                                         depth_values = depths_to_image , 
                                         )
        except : 
            csamtpylog().get_csamtpy_logger().error(
                'Building geosurface obj failed. Something wrong happenned!'
                'Check your oasis files')
            
        for xformat in ['.csv', '.xlsx']:
            try :
                geo_surface_obj.write_file(fileformat = xformat , 
                                            savepath =self.save_geo_outdir )
            except : 
                csamtpylog().get_csamtpy_logger().error(
                    'Export geosurface file to {0} failed !'
                    'Something wrong happend reading oasfiles.'.format(xformat))
            else:
                #collected the input geo and created expected file whithout 
                # extension format' xlsx and '.csv`
                make_filename = ''.join([file.replace(
                    '_cor_oas','').replace('.csv','').replace('.xlsx','') 
                                    for file in os.listdir(OAS_DIR)]) 
                
                if xformat=='.xlsx':
                    gs_expectedfiles = [os.path.join(self.save_geo_outdir,
                                        make_filename + '_gs{0}{1}'.format(
                                        datetime.datetime.now().month, xformat ) )]
                    
                elif xformat=='.csv':
                    
                    gs_expectedfiles = [os.path.join(self.save_geo_outdir, 
                                                     ''.join([make_filename + str(xval),
                                        '_gs{0}.'.format(datetime.datetime.now().month),
                                        'csv']))
                                 for xval in list(geo_surface_obj.geosurface_dico.keys())]
    
                
                gs_outputfiles = [os.path.join(self.save_geo_outdir, file )
                                  for file in os.listdir(self.save_geo_outdir)
                                  if (file.endswith(xformat) 
                                      and file.find('_cor_oas')<0)]
                
                self.assertEqual(len(gs_expectedfiles), len(gs_outputfiles), 
                                 'Difference found between reference output = {0}'
                                 ' & expected outputs = {1}.'.format(
                                     len(gs_outputfiles),len(gs_expectedfiles)))
                
                if xformat =='.csv': # compare files 
                    compare_diff_files(refout = gs_outputfiles ,
                                       refexp = gs_expectedfiles )
                    
    @pytest.mark.skip(reason='Test succeeded on Windox env. with Python 3.7'
                      'but required latest version of pandas library on Linux env.')            
    def test_make_drillhole (self): 
        """
        Test to generate a new DH file. 
        .. note: When parser file is provided , the praser file will read as 
        as main even `build_mannually_welldata` is set to `True`.
        Therefore , to force `buid borehle mannually` by entering data step by step
        set 'build_manually_welldata' to True and set `well_filename` to None.
        
        """
        parser_file = os.path.join(DRILL_PARSER_DIR, 'nbleDH.csv')
        savepath = os.path.join(TEST_TEMP_DIR, self.__class__.__name__)
        
        for dh_type in [ 'collar', 'Geology','Sample','Elevation',
                       'Azimuth', '*']: 
            
            filename=os.path.join(savepath, os.path.basename(parser_file).lower(
                            ).replace('.csv','').replace('.xlsx',''))
            
            # remove the file after succesfully run
            remove_control(rm_file=filename, type_of_file ='Borehole')
  
            kind_of_data2output=dh_type
            # we already test mannually, it's run well , than we test the outputs 
             # of all dh_type 
            buid_borehole_manually =False
             
            try :
                borehole_obj = Drill (well_filename= parser_file, 
                               build_manually_welldata= buid_borehole_manually)
            except : 
                csamtpylog().get_csamtpy_logger().error(
                    'Build borehole failed!  Unable to create borehole obj.')
            else : 
                # then read 
                borehole_obj.writeDHData(data2write=kind_of_data2output, 
                                         savepath = savepath, writeType='.xlsx')
                refout = os.path.join(savepath, os.path.basename(parser_file).lower(
                            ).replace('.csv','').replace('.xlsx','')+'.xlsx')
                self.assertEqual(''.join([filename, '.xlsx']), refout, 
                            'Difference found between reference output = {0} &'
                            'Expected file = {1}.'.format(refout,''.join([filename, '.xlsx'])))
                
                
              
if __name__=='__main__': 
    gt = TestGEODRILL()
    # gt.test_to_geolden_software() 
    # gt.test_to_oasis_montaj()   
    # gt.test_geosurface ()   
    gt.test_make_drillhole()
    # unittest.main()


                

    
        
        
        
        
        
        
        
        
            
            