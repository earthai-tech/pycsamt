# -*- coding: utf-8 -*-
#       Created on Tue Aug 30 14:47:01 2022
#       Author: Kouadio K.Laurent<etanoyau@gmail.com>
#       Licence: GPL

"""
Module ModEM 
=============
.. |ModEM| replace:: Modular system for inversion of electromagnetic 

A |ModEM| geophysical data module . It used  `MTpy`_ module of Lars Krieger and 
Jared R. Peacock ( https://github.com/MTgeophysics/pycsamt.git) modem module as 
dependancy. You may cite the authors when you use these codes. For a deep 
implementation, it is recommended to use full `MTpy` package especially when 
the data is  an MT composed of Tipper and Spectra. 

References 
-----------
.. [1] Kirkby, A., Zhang, F., Peacock, J., Hassan, R., Duan, J., 2019. The MTPy
     software package for magnetotelluric data analysis and visualisation. 
     J. Open Source Softw. 4, 1358. https://doi.org/10.21105/joss.01358
.. [2] Kirkby, A.L., Zhang, F., Peacock, J., Hassan, R., Duan, J., 2020. 
    Development of the open-source MTPy package for magnetotelluric data analysis 
    and visualisation. https://doi.org/10.11636/132198
.. [3] Krieger, L., Peacock, J.R., 2014. MTpy: A Python toolbox for magnetotellurics.
     Comput. Geosci. 72, 167–175. https://doi.org/10.1016/j.cageo.2014.07.013


Created on Tue Aug 30 14:47:01 2022

@author: Daniel

"""

from pycsamt._csamtpylog import csamtpylog 
from pycsamt.core.edi import Edi_collection 
from pycsamt.utils.exceptions import ModEMError

_logger =csamtpylog.get_csamtpy_logger(__name__)

try:
    from mtpy.modeling.modem import ( 
        data , 
        model ,
        data_model_analysis ,
        convariance,
        residual,
        control_fwd, 
        control_inv,  
        station
    )
except  ModuleNotFoundError as e:
    msg=("This error happens because {}. You try to import some ModEM"
         " utilities while some dependencies are not installed. Prior"
         " intall 'MTpy' for dealing with all dependencies at once or"
         " install the related missing module.")
    raise ModEMError (msg.format(str(e).lower()))
except : 
    raise ImportError(
        "Loading failed. Unable to import `MTpy` package. You can get"
        " 'MTpy' at <https://github.com/MTgeophysics/pycsamt.git>`.")
    
__all__= ['DataModelAnalysis', 
          'Data', 
          'Model', 
          'Covariance',
          'ControlFwd', 
          'ControlInv',
          'Stations'
          ]

class Residual(residual.Residual):
    """
    class to contain residuals for each data point, and rms values for each
    station
    
    ====================== ====================================================
    Attributes/Key Words   Description
    ====================== ====================================================
    work_dir
    residual_fn            full path to data file
    residual_array         numpy.ndarray (num_stations) structured to store
                           data.  keys are:
                               * station --> station name
                               * lat --> latitude in decimal degrees
                               * lon --> longitude in decimal degrees
                               * elev --> elevation (m)
                               * rel_east -- > relative east location to
                                               center_position (m)
                               * rel_north --> relative north location to
                                               center_position (m)
                               * east --> UTM east (m)
                               * north --> UTM north (m)
                               * zone --> UTM zone
                               * z --> impedance tensor residual (measured - modelled)
                                       (num_freq, 2, 2)
                               * z_err --> impedance tensor error array with
                                       shape (num_freq, 2, 2)
                               * tip --> Tipper residual (measured - modelled)
                                       (num_freq, 1, 2)
                               * tipperr --> Tipper array with shape
                                       (num_freq, 1, 2)
    rms
    rms_array              numpy.ndarray structured to store station
                           location values and rms.  Keys are:
                               * station --> station name
                               * east --> UTM east (m)
                               * north --> UTM north (m)
                               * lat --> latitude in decimal degrees
                               * lon --> longitude in decimal degrees
                               * elev --> elevation (m)
                               * zone --> UTM zone
                               * rel_east -- > relative east location to
                                               center_position (m)
                               * rel_north --> relative north location to
                                               center_position (m)
                               * rms --> root-mean-square residual for each
                                         station
    rms_tip
    rms_z
    ====================== ====================================================
    """
    def __init__(self, **kws):
        super ().__init__( **kws)
        self._logging=csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        
        for key in list(kws.keys()): 
            setattr(self, key, kws[key])
            
                    
class ControlInv (control_inv.ControlInv): 
    """ Read and write control file for how the inversion starts and how it 
    is run"""
    
    def __init__(self, **kws):
        super ().__init__( **kws)
        self._logging=csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        
        for key in list(kws.keys()): 
            setattr(self, key, kws[key])
            
class ControlFwd (control_fwd.ControlFwd): 
    """ Read and write control file from Forward modeling 
        This file controls how the inversion starts and how it is run """
        
    def __init__(self, **kws):
        super ().__init__( **kws)
        self._logging=csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        
        for key in list(kws.keys()): 
            setattr(self, key, kws[key])
            
            
class Covariance (convariance.Covariance): 
    """
    Read and write covariance files
    
    Parameters 
    -----------
    grid_dimensions: ndarray ((Nx, Ny, Nz))
        Dimension of Grid is the same like :attr:`Model.res_model`` shape 
        Nx : number of grid along x-axis 
        Ny: Number of grid along y-axis 
        Nz: Number of grid along z-depth. 
        
    """
    def __init__(self, grid_dimensions=None, **kws):
        super ().__init__( grid_dimensions= grid_dimensions, **kws)
        self._logging=csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        
        for key in list(kws.keys()): 
            setattr(self, key, kws[key])
            
            
class DataModelAnalysis (data_model_analysis.DataModelAnalysis): 
    """
    Extract info from a pair of files (namely .dat and .rho) of modem 
    re-write the data into other formats such as csv.
    Inversion results Get a slice of the model data for analysis and plotting 
    visualization.

    The output CSV file include StationName, Lat, Long, X, Y, Z, Log(Resistivity)
    where (X,Y,Z) are relative distances in meters from the mesh's origin.
    Projection/Coordinate system must be known in order to associate 
    (Lat, Long) to (X, Y)
    
    
    Parameters 
    ----------
    
    filedat: str , 
        Path to |ModEM| Data file. 
    filerho: str 
        Path to |ModEM| resistivity file
    plot_orient:str 
        plot orientation. Can be ['ew','ns', 'z']. Default is *ew* 
        
    
    """
    def __init__(self, filedat, filerho, plot_orient='ew', **kws):
        super ().__init__( filedat, filerho, plot_orient=plot_orient, **kws)
        self._logging=csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        
        for key in list(kws.keys()): 
            setattr(self, key, kws[key])

class Data (data.Data): 
    """ 
    Data will read and write .dat files for ModEM and convert a WS data file
    to ModEM format.

    ..note: :: the data is interpolated onto the given periods such that all
               stations invert for the same periods.  The interpolation is
               a linear interpolation of each of the real and imaginary parts
               of the impedance tensor and induction tensor.
               See pycsamt.core.mt.MT.interpolate for more details

    Arguments
    ------------
        **edi_list** : list
                       list of full paths to .edi files you want to invert for

    ====================== ====================================================
    Attributes              Description
    ====================== ====================================================
    _dtype                 internal variable defining the data type of
                           data_array
    _logger                python logging object that put messages in logging
                           format defined in logging configure file, see MtPyLog
                           for more information
    _t_shape               internal variable defining shape of tipper array in
                           _dtype
    _z_shape               internal variable defining shape of Z array in
                           _dtype
    center_position        (east, north, evel) for center point of station
                           array.  All stations are relative to this location
                           for plotting purposes.
    comp_index_dict        dictionary for index values of component of Z and T
    station_locations      Stations object
    data_array             numpy.ndarray (num_stations) structured to store
                           data.  keys are:
                               * station --> station name
                               * lat --> latitude in decimal degrees
                               * lon --> longitude in decimal degrees
                               * elev --> elevation (m)
                               * rel_east -- > relative east location to
                                               center_position (m)
                               * rel_north --> relative north location to
                                               center_position (m)
                               * east --> UTM east (m)
                               * north --> UTM north (m)
                               * zone --> UTM zone
                               * z --> impedance tensor array with shape
                                       (num_freq, 2, 2)
                               * z_err --> impedance tensor error array with
                                       shape (num_freq, 2, 2)
                               * tip --> Tipper array with shape
                                       (num_freq, 1, 2)
                               * tipperr --> Tipper array with shape
                                       (num_freq, 1, 2)
    data_fn                full path to data file
    data_period_list       period list from all the data
    edi_list               list of full paths to edi files
    error_type_tipper      [ 'abs' | 'floor' ]
                           *default* is 'abs'
    error_type_z           [ 'egbert' | 'mean_od' | 'eigen' | 'median']
                           *default* is 'egbert_floor'
                                * add '_floor' to any of the above to set the
                                  error as an error floor, otherwise all
                                  components are give weighted the same

                                * 'egbert'  sets error to
                                            error_value_z * sqrt(abs(zxy*zyx))
                                * 'mean_od' sets error to
                                            error_value_z * mean([Zxy, Zyx])
                                            (non zeros)
                                * 'eigen'   sets error to
                                            error_value_z * eigenvalues(Z[ii])
                                * 'median'  sets error to
                                            error_value_z * median([Zxx, Zxy, Zyx, Zyy])
                                            (non zeros)
                           A 2x2 numpy array of error_type_z can be specified to
                           explicitly set the error_type_z for each component.

    error_value_z          percentage to multiply Z by to set error
                           *default* is 5 for 5% of Z as error
                           A 2x2 numpy array of values can be specified to
                           explicitly set the error_value_z for each component.

    error_value_tipper     absolute error between 0 and 1.
    fn_basename            basename of data file. *default* is 'ModEM_Data.dat'
    formatting             ['1' | '2'], format of the output data file, *default* is '1'
    header_strings         strings for header of data file following the format
                           outlined in the ModEM documentation
    inv_comp_dict          dictionary of inversion components
    inv_mode               inversion mode, options are: *default* is '1'
                               * '1' --> for 'Full_Impedance' and
                                             'Full_Vertical_Components'
                               * '2' --> 'Full_Impedance'
                               * '3' --> 'Off_Diagonal_Impedance' and
                                         'Full_Vertical_Components'
                               * '4' --> 'Off_Diagonal_Impedance'
                               * '5' --> 'Full_Vertical_Components'
                               * '6' --> 'Full_Interstation_TF'
                               * '7' --> 'Off_Diagonal_Rho_Phase'

    inv_mode_dict          dictionary for inversion modes
    max_num_periods        maximum number of periods
    model_epsg             epsg code for model projection, provide this to
                           project model to non-utm coordinates. Find the epsg
                           code for your projection on
                           http://spatialreference.org/ref/ or google search
                           epsg "your projection"
    model_utm_zone         alternative to model_epsg, choose a utm zone to
                           project all sites to (e.g. '55S')
    mt_dict                dictionary of mtpy.core.mt.MT objects with keys
                           being station names
    period_buffer          float or int
                           if specified, apply a buffer so that interpolation doesn't
                           stretch too far over periods
    period_dict            dictionary of period index for period_list
    period_list            list of periods to invert for
    period_max             maximum value of period to invert for
    period_min             minimum value of period to invert for
    period_buffer          buffer so that interpolation doesn't stretch too far
                              over periods. Provide a float or integer factor, 
                              greater than which interpolation will not stretch.
                              e.g. 1.5 means only interpolate to a maximum of
                              1.5 times each side of each frequency value
    rotate_angle           Angle to rotate data to assuming 0 is N and E is 90
    save_path              path to save data file to
    units                  [ [V/m]/[T] | [mV/km]/[nT] | Ohm ] units of Z
                           *default* is [mV/km]/[nT]
    wave_sign_impedance    [ + | - ] sign of time dependent wave.
                           *default* is '+' as positive downwards.
    wave_sign_tipper       [ + | - ] sign of time dependent wave.
                           *default* is '+' as positive downwards.
    ====================== ====================================================

    ========================== ================================================
    Methods                    Description
    ========================== ================================================
    center_stations            Center station locations to the middle of cells,
                               might be useful for topography.
    change_data_elevation      At each station in the data file rewrite the
                               elevation, so the station is on the surface,
                               not floating in air.
    compute_inv_error          compute the error from the given parameters
    convert_modem_to_ws        convert a ModEM data file to WS format.
    convert_ws3dinv_data_file  convert a ws3dinv file to ModEM fomrat,
                               **Note** this doesn't include tipper data and
                               you need a station location file like the one
                               output by pycsamt.modeling.ws3dinv
    fill_data_array            fill the data array from mt_dict
    filter_periods             Select the periods of the mt_obj that are in
                               per_array. used to do original freq inversion.
    get_header_string          reset the header sring for file
    get_mt_dict                get mt_dict from edi file list
    get_parameters             get important parameters for documentation
    get_period_list            make a period list to invert for
    get_relative_station_locations     get station locations from edi files
    project_stations_on_topography     This method is used in add_topography().
                                       It will Re-write the data file to change
                                       the elevation column. And update
                                       covariance mask according topo elevation
                                       model.
    read_data_file             read in a ModEM data file and fill attributes
                               data_array, station_locations, period_list,
                               mt_dict
    write_data_file            write a ModEM data file
    write_vtk_station_file     write a vtk file for station locations.  For now
                               this in relative coordinates.
    ========================== ================================================

    :Example 1 --> create inversion period list: ::

        >>> import os
        >>> import pycsamt.modeling.modem as modem
        >>> edi_path = r"/home/mt/edi_files"
        >>> edi_list = [os.path.join(edi_path, edi) \
                        for edi in os.listdir(edi_path)\
                        if edi.find('.edi') > 0]
        >>> md = modem.Data(edi_list, period_min=.1, period_max=300,\
                            max_num_periods=12)
        >>> md.write_data_file(save_path=r"/home/modem/inv1")

    :Example 2 --> set inverions period list from data: ::

        >>> import os
        >>> import pycsamt.core.cs
        >>> import pycsamt.modeling.modem as modem
        >>> edi_path = r"/home/mt/edi_files"
        >>> edi_list = [os.path.join(edi_path, edi) \
                        for edi in os.listdir(edi_path)\
                        if edi.find('.edi') > 0]
        >>> md = modem.Data(edi_list)
        >>> #get period list from an .edi file
        >>> mt_obj1 = cs.CSAMT(edi_list[0])
        >>> inv_period_list = 1./mt_obj1.Z.freq
        >>> #invert for every third period in inv_period_list
        >>> inv_period_list = inv_period_list[np.arange(0, len(inv_period_list, 3))]
        >>> md.period_list = inv_period_list
    >>> md.write_data_file(save_path=r"/home/modem/inv1")

    :Example 3 --> change error values: ::

        >>> import pycsamt.modeling.modem as modem
        >>> mdr = modem.Data()
        >>> mdr.read_data_file(r"/home/modem/inv1/ModEM_Data.dat")
        >>> mdr.error_type = 'floor'
        >>> mdr.error_floor = 10
        >>> mdr.error_tipper = .03
        >>> mdr.write_data_file(save_path=r"/home/modem/inv2")

    :Example 4 --> change inversion type: ::

        >>> import pycsamt.modeling.modem as modem
        >>> mdr = modem.Data()
        >>> mdr.read_data_file(r"/home/modem/inv1/ModEM_Data.dat")
        >>> mdr.inv_mode = '3'
        >>> mdr.write_data_file(save_path=r"/home/modem/inv2")

    :Example 5 --> rotate data: ::

        >>> md.rotation_angle = 60
        >>> md.write_data_file(save_path=r"/home/modem/Inv1")
        >>> # or
        >>> md.write_data_file(save_path=r"/home/modem/Inv1", \
                               rotation_angle=60)
    """
 
    def __init__(self, edi_list,  period_min = None, period_max = None, 
                 max_num_periods = None, data_period_list = None,  **kws): 
        self._logging=csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        edi_list = Edi_collection(edi_list).edifiles 
        super ().__init__(edi_list= edi_list , period_min = period_min,
                          period_max = period_max, max_num_periods = max_num_periods,
                          data_period_list = data_period_list,  **kws)
        
        
        for key in list(kws.keys()): 
            setattr(self, key, kws[key])
            
    
class Model (model.Model):
    """ 
    make and read a FE mesh grid

    The mesh assumes the coordinate system where:
        x == North
        y == East
        z == + down

    All dimensions are in meters.

    The mesh is created by first making a regular grid around the station area,
    then padding cells are added that exponentially increase to the given
    extensions.  Depth cell increase on a log10 scale to the desired depth,
    then padding cells are added that increase exponentially.

    Arguments
    -------------
        **station_object** : pycsamt.modeling.modem.Stations object
                            .. seealso:: pycsamt.modeling.modem.Stations

    Examples
    -------------

    :Example 1 --> create mesh first then data file: ::

        >>> import pycsamt.modeling.modem as modem
        >>> from  pycsamt.processing import get_full_frequency 
        >>> from pycsamt.utils.func_utils import get_ediObjs 
        >>> import os
        >>> # 1) make a list of all .edi files that will be inverted  or 
        >>> # read ediObjs 
        >>> edi_path = r"data/edis"
        >>> ediObjs = get_ediObjs (edipath) # or get the list of .edi files 
        >>> edi_list = list(map(lambda o: o.edifile , ediObjs))
        >>> # 2) Make a Stations object
        >>> stations_obj = modem.Stations()
        >>> stations_obj.get_station_locations(edi_list)
        >>> # 3) make a grid from the stations themselves with 200m cell spacing
        >>> mmesh = modem.Model(stations_obj)
        >>> # change cell sizes
        >>> mmesh.cell_size_east = 20
        >>> mmesh.cell_size_north = 20
        >>> mmesh.ns_ext = 3000 # north-south extension
        >>> mmesh.ew_ext = 2000 # east-west extension of model
        >>> mmesh.make_mesh()
        >>> # check to see if the mesh is what you think it should be
        >>> mmesh.plot_mesh()
        >>> # all is good write the mesh file
        >>> mmesh.write_model_file(save_path=r"data/Inv1")
        >>> # create data file
        >>> f = get_full_frequency(ediObjs)
        >>> # create object with filling all arguments (1) 
        >>> md = modem.Data(edi_list,period_min = 1/f.max() , period_max = 1/f.min() ,
                            max_num_periods = len(f) //2 ,model_utm_zone ='49N',
                            data_period_list = 1/f )
        >>> # create object and set attribute arguments 
        >>> md.modem.Data(edi_list) 
        >>> md.data_period_list = 1/f  
        >>> md.max_num_periods = len(f) //2
        >>> md.model_utm_zone ='49N'  
        >>> md.period_min = 1/f.max() ; md.period_max = 1/f.min()
        >>> md.station_locations = mmesh.station_locations
        >>> md.write_data_file(save_path=r"data/Inv1)

    :Example 2 --> Rotate Mesh: ::

        >>> mmesh.mesh_rotation_angle = 60
        >>> mmesh.make_mesh()

    .. note:: ModEM assumes all coordinates are relative to North and East, and
             does not accommodate mesh rotations, therefore, here the rotation
             is of the stations, which essentially does the same thing.  You
             will need to rotate you data to align with the 'new' coordinate
             system.

    ==================== ======================================================
    Attributes           Description
    ==================== ======================================================
    _logger              python logging object that put messages in logging
                         format defined in logging configure file, see MtPyLog
                         more information
    cell_number_ew       optional for user to specify the total number of sells
                         on the east-west direction. *default* is None
    cell_number_ns       optional for user to specify the total number of sells
                         on the north-south direction. *default* is None
    cell_size_east       mesh block width in east direction
                         *default* is 500
    cell_size_north      mesh block width in north direction
                         *default* is 500
    grid_center          center of the mesh grid
    grid_east            overall distance of grid nodes in east direction
    grid_north           overall distance of grid nodes in north direction
    grid_z               overall distance of grid nodes in z direction
    model_fn             full path to initial file name
    model_fn_basename    default name for the model file name
    n_air_layers         number of air layers in the model. *default* is 0
    n_layers             total number of vertical layers in model
    nodes_east           relative distance between nodes in east direction
    nodes_north          relative distance between nodes in north direction
    nodes_z              relative distance between nodes in east direction
    pad_east             number of cells for padding on E and W sides
                         *default* is 7
    pad_north            number of cells for padding on S and N sides
                         *default* is 7
    pad_num              number of cells with cell_size with outside of
                         station area.  *default* is 3
    pad_method           method to use to create padding:
                         extent1, extent2 - calculate based on ew_ext and
                         ns_ext
                         stretch - calculate based on pad_stretch factors
    pad_stretch_h        multiplicative number for padding in horizontal
                         direction.
    pad_stretch_v        padding cells N & S will be pad_root_north**(x)
    pad_z                number of cells for padding at bottom
                         *default* is 4
    ew_ext               E-W extension of model in meters
    ns_ext               N-S extension of model in meters
    res_scale            scaling method of res, supports
                           'loge' - for log e format
                           'log' or 'log10' - for log with base 10
                           'linear' - linear scale
                         *default* is 'loge'
    res_list             list of resistivity values for starting model
    res_model            starting resistivity model
    res_initial_value    resistivity initial value for the resistivity model
                         *default* is 100
    mesh_rotation_angle  Angle to rotate the grid to. Angle is measured
                         positve clockwise assuming North is 0 and east is 90.
                         *default* is None
    save_path            path to save file to
    sea_level            sea level in grid_z coordinates. *default* is 0
    station_locations    location of stations
    title                title in initial file
    z1_layer             first layer thickness
    z_bottom             absolute bottom of the model *default* is 300,000
    z_target_depth       Depth of deepest target, *default* is 50,000
    ==================== ======================================================


    ==================== ======================================================
    Methods              Description
    ==================== ======================================================
    add_topography_to_model2    if air_layers is non-zero, will add topo: read
                                in topograph file, make a surface model.
                                Call project_stations_on_topography in the end,
                                which will re-write the .dat file.
                                If n_airlayers is zero, then cannot add topo
                                data, only bathymetry is needed.
    assign_resistivity_from_surfacedata     assign resistivity value to all
                                            points above or below a surface
                                            requires the surface_dict attribute
                                            to exist and contain data for
                                            surface key (can get this
                                            information from ascii file using
                                            project_surface)
    get_parameters       get important model parameters to write to a file for
                         documentation later.
    interpolate_elevation2  project a surface to the model grid and add
                            resulting elevation data to a dictionary called
                            surface_dict. Assumes the surface is in lat/long
                            coordinates (wgs84)
    make_mesh            makes a mesh from the given specifications
    make_mesh_from_center   The mesh is built by first finding the center of
                            the station area. Then cells are added in the north
                            and east direction with width cell_size_east and
                            cell_size_north to cover all the station area.
    make_z_mesh          Create a mesh grid for vertical Earth layers.
    make_z_mesh_exp      version of make_z_mesh method in order to create
                         exp-increasing cell sizes from the top layer
    make_z_mesh_new      new version of make_z_mesh. make_z_mesh and M
    plot_mesh            plots mesh to make sure everything is good
    plot_mesh_xy         add mesh grid lines in xy plan north-east map
    plot_mesh_xz         display the mesh in North-Depth aspect
    plot_topograph       display topography elevation data together with
                         station locations on a cell-index N-E map
    print_mesh_params    print out the mesh-related paramas
    print_model_file_summary    print out the summary of the model file
    project_surface      project a surface to the model grid and add resulting
                         elevation data to a dictionary called surface_dict.
                         Assumes the surface is in lat/long coordinates (wgs84),
                         if not, need to supply the epsg of the surface xy
                         points
    read_dem_ascii       read in dem which is ascii format
    read_model_file      read an initial file and return the pertinent
                         information including grid positions in coordinates
                         relative to the center point (0,0) and starting model.
    read_ws_model_file   reads in a WS3INV3D model file
    write_model_file     writes an initial model file that includes the mesh
    write_vtk_file       write a vtk file to view in Paraview or other
    ==================== ======================================================
    
    """
    def __init__(self, stations_object = None, data_object = None, **kws ): 
        super ().__init__(stations_object= stations_object, data_object = data_object, 
                          **kws)
        self._logging=csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        
        for key in list(kws.keys()): 
            setattr(self, key, kws[key])
    
    
class Stations(station.Stations):
    """
    station locations class

    .. note:: If the survey steps across multiple UTM zones, then a
             distance will be added to the stations to place them in
             the correct location.  This distance is
             _utm_grid_size_north and _utm_grid_size_east.  You should
             these parameters to place the locations in the proper spot
             as grid distances and overlaps change over the globe.
             **This is not implemented yet**
    """    
    def __init__(self, **kws):
        self._logging=csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        super().__init__(**kws)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    