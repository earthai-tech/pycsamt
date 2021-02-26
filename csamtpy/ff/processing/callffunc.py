# -*- coding: utf-8 -*-
"""
===============================================================================
    Copyright © 2021  Kouadio K.Laurent
    
    This file is part of pyCSAMT.
    
    pyCSAMT is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    pyCSAMT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
    
    You should have received a copy of the GNU Lesser General Public License
    along with pyCSAMT.  If not, see <https://www.gnu.org/licenses/>.

===============================================================================

    .Typical call files functions  
    .call files from Avg, AGSO, else... 
    
Created on Thu Nov 26 20:55:39 2020

@author: @Daniel03

"""
import os,warnings 
# import warnings
import numpy as np
import pandas as pd
from csamtpy.utils import exceptions as CSex
from csamtpy.utils._csamtpylog import csamtpylog
from csamtpy.utils import func_utils as func
_logger=csamtpylog.get_csamtpy_logger(__name__)



#*** file from _geocodes folder :AGSO & AGSO.STCODES ******


def _code_strata():
    return ['code', 'label', 'name', 'pattern',
                  'size', 'density', 'thickness', 'color']


def agso_data (): 
    """
    Geological data codes processing
    :: will later deprecated :: 
    """
    
    agsofiles =['AGSO.csv', 'AGSO_STCODES.csv']
    
    for file in agsofiles : 
        path_to_agsofile=os.path.join(os.environ["pyCSAMT"], 'geodrill',
                                                 "_geocodes", file)
        with open(path_to_agsofile,'r', encoding='utf8') as f:
                fgeo=f.readlines()
        for ss, stratum in enumerate(fgeo):
            stratum=stratum.strip('\n').split(',')
            fgeo[ss]=stratum
        if file =='AGSO.csv':
            fagso=fgeo
        elif file =='AGSO_STCODES.csv':
            fagso_stcodes=fgeo
            fagso_stcodes[0]=fagso[0]
        else :
            _logger.warn("Geocodes files < AGS0 & AGS0_STCODES> Not found !")
            CSex.pyCSAMTError_file_handling('Geocodes files :< AGS0 & AGS0_STCODES> not found.'\
                                            'check the right path to _geocodes folder. ')
    return fagso, fagso_stcodes


def set_stratum_on_dict():
    """
    Process to put geocodes_strata and geocodes_structures into dictionnaries  
    better way to go on metaclasses merely. Thus each keys of dictionary will be 
    its own object. 

    Returns
    -------
    strata_dict : dict
        Disctionnary of geostrata .
    structures_dict : dict
        Dictionnary of geostructures.

    """
    code_strata=_code_strata()
    drop_bar=[' / ',' ', '/','-']
    attr_dict={}
    strata_dict, structures_dict={},{}
    for iiagso in range(2):
        fstrata=agso_data()[iiagso][1:]
        
        for ss , elem in enumerate(fstrata):
            for jj, attribute in enumerate(code_strata):
                attr_dict["{0}".format(attribute)]=elem[jj]
                stratum_description=elem[2]
            for db in drop_bar :
                if db in stratum_description : 
                    stratum_description=stratum_description.replace(db,'_')
            stratum_description=stratum_description.lower()
            if stratum_description.startswith('"') or stratum_description.endswith('"'):
                stratum_description=stratum_description[1:-1]
            if iiagso ==1 :
                structures_dict['{0}'.format(stratum_description)]=attr_dict
            else :
                strata_dict['{0}'.format(stratum_description)]=attr_dict
            attr_dict={}

    return strata_dict, structures_dict

#***  end call file from _geocodes folder :AGSO & AGSO.STCODES ******

#***start Avg files extract **** 
    
def _validate_avg_file(avg_dataLines):
    
    """
    Function to validate and readZonge Avg FILE :Plainty file _F1 and the Astatic file 
    _F2. Both files can be read and separatedd by the program. if Avg file doesnt match 
    any of them , the programm will generate errors . 
    
    params : 
    -------
        ** avag_dataLines: list 
                    list of avgData for all lines. 
    
    returns :
        _F1 & _F2: int , 
            type of AVG Zonge file .
  
                
    """
        
    remark =['\\ AMTAVG','$ ASPACE', '$ XMTR '] 
    count=0
    
    import inspect 
    
    func_name = inspect.getframeinfo(inspect.currentframe())[2]
    _logger.info ('Check the Avg_Datalinesfile %s',func_name )
      
    if avg_dataLines is not ['']:
        filehead="".join([ss for ss in avg_dataLines[:3]])
        #------------------checker section ------------------------
        #--->  check if astatic file is on the head_section 
        if filehead.find('astatic'.upper())>1: 
            _FT=2
        elif 'astatic'.upper() not in filehead: 
            for ss in remark :
                if filehead.find(ss) > 0 :count +=1
            #---> that mean all the remak elemt are on the head_section 
            if count >=2 :
                _FT=1
        else:
            _logger.error('Wrong Avg file.  Must respect the typical zonge file.')
            raise CSex.pyCSAMTError_avg_file('Type of AVg file provide is Wrong. No match'\
                                             'the typical Zonge engeneering file.')       
    else : 
        _logger.error('No Zonge Engeneering files exists. No data found on that curent AVG file.'\
                      'Check the better path to Avgfiles right file.')
        raise CSex.pyCSAMTError_avg_file('Check the The type of Avg file.')
            #-------------------end checker -----------------------------
    return  _FT
        



def straighten_cac2CSfile (data_array, component_column_section = None):
    """
         Sometimes head_sections of file _F2(CAC2CSAMT) provided is little 
         different colunms section name  according to different version . 
         it 's better to filter  and to check before returning the 
         correct informations we need.
    
    params : 
    -------
        data_array : ndarray 
            data from AVG astatic file 
        
         component_column_section : list 
             astactic file column comps provided . 
    returns : 
        array_we_need : ndarray
            same infos present in the plainty /1 Avg file 
        array_other_comp :  pd.Core.DataFrame 
             infos include through Astatic softwares . very usefull therefore 
             we keep it .
        
             
    """
    _comp_not_in_F1, _comp_we_need=['Z.mwgt','Z.pwgt','Z.mag',
                                     'SRes','E.wght',
                                     'H.wght','Z.%err',
                                     'Z.perr','Gdp.Blk',
                                     'Gdp.Chn', 'Gdp.Time'],['skp','Freq',
                                                            'Tx.Amp', 'E.mag',
                                                            'B.mag','B.phz',
                                                            'Z.phz','ARes.mag',
                                                            'E.%err','E.perr','B.%err',
                                                            'B.perr','Z.perr',
                                                            'ARes.%err','E.phz']           
                                            
    
    if component_column_section is None : 
        raise CSex.pyCSAMTError_inputarguments('head-section must be on list of string items.')
        
    tem_del, other_del=[],[]
    #strip comps _column to check....
    component_column_section=func._strip_item(item_to_clean=component_column_section)
    
    #---> put on dataFrame the component wee_need 
    df=pd.DataFrame(data=data_array, columns= component_column_section)
    for ss in component_column_section : 
        if ss not in _comp_we_need: 
            tem_del.append(ss)
        if ss in _comp_we_need :
            other_del.append(ss)
    df_we_need=df.drop(tem_del, axis=1)
    df_we_need.reset_index(drop=True, inplace=True)
    array_we_need= df_we_need.to_numpy()
    #--->> _keep on other dataframe the component we dont need.
    # other_comp_list =func._cross_eraser(data=component_column_section, to_del=tem_del,
    #                                     deep_cleaner=False)
    df_other_comp=df.drop(other_del, axis=1)
    df_other_comp.reset_index(drop=True, inplace=True) 

    
    #--> add skpi value if not in component section : 
    if 'skp'  not in [ii.lower() for ii in component_column_section]: # to spare a manner 'Skp' is written.
        skp_value =2   #  give a certitude value that data were good quality : 
        skp_array=np.full((array_we_need.shape[0],1),skp_value)
        array_we_need=np.concatenate((skp_array,array_we_need), axis=1)

    
    return array_we_need , df_other_comp


def truncated_data (data, number_of_reccurence, **kwargs):
    """
    Function to truncate all data according to number of frequency.
    
    Parameters
    ----------
    data : list, or nd.array 
        data must be truncate.
    number_of_freq : int
        number of frequency imaged.

    Returns
    -------
    loc_list : list
        data truncated on list.

    """
    if type(data) is list : 
        data =np.array(data) 
    value_lengh= data.shape[0]
    # index_loc = [np.int(ss) for ss in np.linspace(0, value_lengh, number_of_station)]
    index_loc = [ss for ss in np.arange(0,value_lengh, number_of_reccurence)]

    loc_list=[]
    for ss, value  in enumerate(index_loc) : 
        if ss ==len(index_loc)-1 : 
            loc_list.append(data[value:])
            break
        loc_list.append(data[value:index_loc[ss+1]])
            
    return loc_list 
    

def _numbering_station(number_of_station, number_of_freq): 
    """
    small function to numbering stations : 
    params : 
        number_of_station :int , 
            number of station found on the site. 
        number_of freq : int , 
            number of frequency found for survey at each station.
    return : 
        numbsta : station numbered 
        poly_sta : multplie station for each frequency for each stations.
    """
   
    numbsta =['S{:02}'.format(ii) for ii in range (number_of_station)]
    poly_sta=numbsta*number_of_freq 
    poly_sta.sort()
    
    return numbsta , poly_sta


def zstar_array_to_nan (zstar_array, nan_value=np.nan, keep_str =False): 
    """
    Parameters
    ----------
    zstar_array : ndarray
        array contain unfloat converter value. the unconverter value can be a '*'
    nan_value : float or np.nan type
        the nan_value could be any value either int, float or str.  i The default is np.nan.
    keep_str : bool, optional
        keep the str item on your array. f keep_str is set to 
        false and the type nan_value is str , the program will force 'keep_str_'  to True 
        to allow converter . The default is False.

    Returns
    -------
    TYPE ; ndarray 
        zstrar_array converted .

    """
    
    for kk , value in enumerate(zstar_array): 
        try : zstar_array[kk] =float(value)
        except :zstar_array[kk]=nan_value
    if type(nan_value) is str : keep_str=True
    if keep_str : return zstar_array 
    return np.array([float(kk) for kk in zstar_array])

   
def get_array_from_reffreq ( array_loc, freq_array,reffreq_value, stnNames=None):
    """ 
    Get array value at special frequency
    params : 
        array_loc : dict , dictionnary of stations :array_value :eg: S00:(ndarray,1) rho_values
        freq_array : (ndarray,1) : frequency array for CSAMT survey 
        reffreq_value : int or float : the value of frequency user want to get the value 
        stnNames : list , list of stations names . 
    
    return : 
        (ndarray,1) : an array of all station with reffreq_value . 
        eg : reffreq_value =1024. it return all value of the array at 1024Hz frequency . 
    
    """
    if stnNames is None : raise CSex.pyCSAMTError_station('You may at least specify '\
                        'the array or list of stations-Names or station id.')
    
    arrayObj,freqObj,rfObj,stnObj =array_loc ,freq_array,reffreq_value,stnNames

    def get_reffreq_index(freq_array, reffreq_value): 
        """ 
        get the index of reference index. From this index ,All array will filter data at this reffreq
        value . 
        
        """
    
        for ii, freq in enumerate(freq_array): 
            if freq == reffreq_value:
                index_rf = ii
                break
        return index_rf
    
    return np.array([values[get_reffreq_index(freq_array=freqObj, reffreq_value=rfObj)] \
                for stn in stnObj for keys, values in arrayObj.items() if stn==keys])
        
def relocate_on_dict_arrays(data_array, number_of_frequency, station_names =None): 
    """ 
    Put data arrays on dictionnary where keys is each station and value the array of that station.
    if station_names is None , program will create name of station. if station_names is given ,
    function will sorted stations names . please make sure to provide correctly station according 
    the disposal you want . 
    params : 
        number_of_frequency : array_like (ndarray,1) : array of frequency during survey 
        station_names: list of array_like : list of station . 
        
    return : dict , key :stations , value : infos at data stations . 
    
    """
        
    if station_names is not None : 
        if type(station_names) is list :
            assert len(station_names) == (data_array.size / number_of_frequency),\
                CSex.pyCSAMTError_station('Number of Station provided must be <{0}> not <{1}>'.format(np.int(data_array.size / number_of_frequency)), 
                                            len(station_names))
        else :
            assert station_names.size == data_array.size / number_of_frequency, \
            CSex.pyCSAMTError_station('Number of Station provided must be <{0}> not <{1}>'.format(np.int(data_array.size / number_of_frequency)), 
                                            station_names.size)
            station_names =station_names.tolist()
            
        NUM_STN=sorted(station_names)
    elif station_names is None  :
        station_length =np.int(data_array.size / number_of_frequency )
        
        NUM_STN=sorted((_numbering_station(number_of_station=station_length ,
                                    number_of_freq=number_of_frequency))[0])

    LIST_BROK = truncated_data(data=data_array, number_of_reccurence=number_of_frequency)   
    
    return { key: value for key, value in zip (NUM_STN, LIST_BROK)}

def dipole_center_position (dipole_position=None):
    """
    Generaly positions are taken at each electrode of dipole to that to easy correct data  for ploting and for
    noise correction , we adjust coordinate by taking the center position that means , 
    the number of points will be substract to one.
    params : dipole_postion : array_like (ndarray,1): postion array at each electrodes. 
    
    return : array_centered_value . 
    
    """
    if dipole_position is None : raise CSex.pyCSAMTError_inputarguments('NoneType can not be compute.'\
                                                                        ' Please provide the right values.')
    try : dipole_position= np.array([float(pp) for pp in dipole_position])
    except : raise CSex.pyCSAMTError_float('Cannot convert position values. '\
                                           'Position must be a number not str.Please check your data!.')
    # temp_pospp =[(dipole_position[pp]+dipole_position[pp-1])/2 for \
    #              pp, pos in enumerate( dipole_position) if ( pp>0 and pp <= dipole_position.size) ] 
    return np.array([(dipole_position[pp]+dipole_position[pp-1])/2 for \
                 pp, pos in enumerate( dipole_position) if ( pp>0 and pp <= dipole_position.size) ])
        


# if __name__=='__main__':
    
#     # print(numfreq, num_sta)
#     ch =['2', '7', '9', 'a', 'g', 'e', '*', '9','2', '7', '9', '*', 'g', 'e', 'geth', '9']
#     array =np.array(ch)
#     test=zstar_array_to_nan(zstar_array=array, keep_str=False)

    
    
    
    
    
    
    
    
    
    