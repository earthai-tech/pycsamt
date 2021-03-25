# -*- coding: utf-8 -*-
# ! /usr/bin/env python
"""
Description:
    Test create input files required to read  occam2d inversion model files, which include:
    ('Occam2DMesh',    'Occam2DModel',    'Startup', 'IterXX.iter')
    These input files are created from standard edi data file set.

References:
    scripts/read_occam2d_res_model.py

CreationDate:   25/03/2021
Developer:      kkouao@zju.edu.cn

inspired from script of Zhang Fei of Geosciences Australia: fei.zhang@ga.gov.au
    
"""
from . import reset_matplotlib
from . import diff_files

import os
from unittest import TestCase

# import csamtpy.modeling.occam2d as occam2d
from viewer.plot import Plot2d

from tests import EDI_DATA_DIR, OC2D_DIR, make_temp_dir,TEST_TEMP_DIR 



class TestOccam2D(TestCase):
    @classmethod
    def setUpClass(cls):
        reset_matplotlib()
        cls._temp_dir = make_temp_dir(cls.__name__)

    def setUp(self):

        # set the dir to the output from the previously correct run
        self._expected_output_dir = os.path.join(OC2D_DIR, 'Occam2d')

        if not os.path.isdir(self._expected_output_dir):
            self._expected_output_dir = None

        # directory to save readed outputfiles t files
        self._output_dir = make_temp_dir('Occam2d', self._temp_dir)

    def test_fun(self):
        """

        :return:
        """
        
        outdir = self._main_func(edipath=EDI_DATA_DIR)

        for afile in ('Occam2DMesh', 'Occam2DModel', 'Occam2DStartup'):
            output_data_file = os.path.join(outdir, afile)
            self.assertTrue(os.path.isfile(output_data_file), "output data file not found")

            expected_data_file = os.path.join(self._expected_output_dir, afile)

            self.assertTrue(os.path.isfile(expected_data_file),
                            "Ref output data file does not exist, nothing to compare with"
                            )

            print(("Comparing", output_data_file, "and", expected_data_file))

            is_identical, msg = diff_files(output_data_file, 
                                           expected_data_file, ignores=['Date/Time:'])
            print(msg)
            self.assertTrue(is_identical, "The output file is not the same with the baseline file.")

    def test_plot_model_and_responses(self):
        """
            test function
            :return: T/F
            """
        # imaging depth : Maximum depth investigation 
        doi = '1km'                 #  can be float like 1000 = 1km 
        
        #plot style 
        plotStyle ="pcolormesh"            # if None Default is 'imshow', can be 
                                    #["pcolormesh"]
        
        savefigure = os.path. join(TEST_TEMP_DIR, 'k1.png')
        #-----OCCAM 2D output data files ---------
        # path to occam Data file
        path_to_occam_data='OccamDataFile.dat'
        
        # path to occam Mesh file 
        path_to_occam_mesh = 'Occam2DMesh'
        
        # path to occam Model file 
        path_to_occam_model = 'Occam2DModel'
        # path to Occam Iteration file 
        path_to_occam_iter='ITER10.iter'
        
        # show a report 
        see_report =False      # generate groundwater report           
        
        figsize =[9,9]
        # call plot obj 
        plot2d_obj = Plot2d(fig_size =figsize)#,fig_dpi = 600. )
        plot2d_obj.plot_occam2dModel(mesh_fn=os.path.join(OC2D_DIR, path_to_occam_mesh), 
                            iter_fn = os.path.join(OC2D_DIR, path_to_occam_iter), 
                            model_fn =os.path.join(OC2D_DIR, path_to_occam_model ), 
                            data_fn =os.path.join(OC2D_DIR, path_to_occam_data), 
                            doi= doi, 
                            savefig =savefigure, 
                            plot_style =plotStyle, 
                            show_report = see_report )
