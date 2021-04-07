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

.. _module-Plot_Utils:: `pycsamt.utils.plot_utils` 
    :synopsis: helpers functions for visualization 
    
Created on Tue Dec 29 19:18:44 2020

@author: @Daniel03
"""
import re , os, warnings
import numpy as np
from pycsamt.etc.infos import _sensitive as SB
import matplotlib as mpl 
import matplotlib.cm as cm 
import matplotlib.pyplot as plt
from pycsamt.utils import exceptions as CSex 
from pycsamt.utils.decorator import deprecated, redirect_cls_or_func


def share_props_for_each_plot(number_of_plot = 3 ,  **kwargs):
    """
    Function to set properties for each plot. Easy to customize line and markers.
    Function can add other properties which are not in kwargs.keys(). 
    It will set according the number of subplotsplots 
    we assume subplot are define on one columns 
    
    :param number_of_plot: number of subsplot you want show.
    :type number_of_plot: int 
    
    :returns: dictionarry of labels and properties. 
    :rtype: dict 
    
    """
    
    def set_multiple(labelprops, default=None): 
        """ 
        set multiple labels 
        
        """
        if not isinstance(labelprops, (tuple, list, str)):
            labelprops=[labelprops] 
            
        # labelprops= [";", 'y']
        # chech is the label properties is str , user can forget to put on list , function  will check 
        # at least if there is comma ", " or double point ":" then it will put it on list . 
        if  isinstance (labelprops, str ): 
            # if re.match(r',', labelprops) is not None  : labelprops= labelprops.split(',')
            if (labelprops.find(',') > 0) and (labelprops.find(':') >0): 
                labelprops=labelprops.replace(',',":")
                labelprops = labelprops.split(':')
            elif labelprops.find(',') > 0:labelprops= labelprops.split(',')
            elif labelprops.find(':') >0 : labelprops = labelprops.split(':')
            else :labelprops=[labelprops]
            
        #set the default properties if not use the first element of provide properties.
        
        if default is not None :  lbproperty = [default for ll in range(number_of_plot)]
        else : lbproperty = [labelprops[0] for ll in range(number_of_plot)]

        if len(labelprops)  < number_of_plot : 
           lbproperty[:len(labelprops)]=labelprops
        elif len(labelprops)  > number_of_plot :
            lbproperty = labelprops[:number_of_plot]
        elif len(labelprops)  ==  number_of_plot : 
            lbproperty = labelprops
        
        # else :raise CSex.pyCSAMTError_plot_tip('Number of argument of label properties'\
        #                                         ' provided doesnt not match the number of subplots. Please rectify the numve')
        return lbproperty 

    # linewidth =
    linewidth=set_multiple(labelprops=kwargs.pop('lw', 2.) , default=2.)
    
    linestyle=set_multiple(labelprops=kwargs.pop('ls',2.) , default=2.)
    color=set_multiple(labelprops=kwargs.pop('color', 'k') , default='k')
    
    marker=set_multiple(labelprops=kwargs.pop('marker', 'o'), default='o')
   
    markerfacecolor=set_multiple(labelprops=kwargs.pop('markerfacecolor', 'k'), default='k')
    markeredgecolor =set_multiple(labelprops=kwargs.pop('markeredgecolor','k') , default='k')
    
    rotate_angle_xticklabels=set_multiple(labelprops=kwargs.pop ('xtick_label_rotation',45.), default=45.)
    rotate_angle_yticklabels=set_multiple(labelprops=kwargs.pop('ytick_label_rotation',None), default=None)
    
    x_ticks_labelsize =set_multiple(labelprops=kwargs.pop('xtick_label_size', 12.) , default=12.)
    y_ticks_labelsize=set_multiple(labelprops=kwargs.pop('ytick_label_size', 12.) , default=12.)

    # dict of labelproperties.
    kws ={'lw':linewidth, 
          'ls':linestyle, 
          'marker':marker, 
          'color':color, 
          'markerfacecolor':markerfacecolor, 
          'markeredgecolor':markeredgecolor, 
          'xtick_label_rotation':rotate_angle_xticklabels, 
          'ytick_label_rotation':rotate_angle_yticklabels, 
          'xtick_label_size':x_ticks_labelsize , 
          'ytick_label_size':y_ticks_labelsize}

    #can add other properties which are not in kwargs .keys ()
    for keys in list(kwargs.keys()):
        if keys not in list(kws.keys()):
            kws.__setitem__(keys, set_multiple(labelprops=kwargs[keys], default=kwargs[keys]))
        
    
    return kws

def find_path (path =None, ptol =0.7):
    """
    check path and return filepath , edipath or jpath .
    
    :param path: full  path  to `edi`, `avag` or `j` file or directory 
    :type path: str 
    
    :param ptol:  tolerance that given by the program to judge if the number 
                    of typical file [EDI|J] 
                    to declare as path found is either "edipath" or "jpath" 
                    if none ,return None . less or equal to 1.
                    
    :type ptol: float 
    
    :returns: specific path 
    :rtype: str 
    """
    if path is None : raise CSex.pyCSAMTError_plot_tip('Can not find path to read.please provided a datapath . ')
    if path is not None : 
        if os.path.isfile (path) is True : return 'isfile'
        elif os.path.isdir(path) is True :
            if os.listdir(path) is not None : 
                ex = [file for file in os.listdir(path) if os.path.splitext(file)[1] =='.edi']
                if len(ex)/len(os.listdir(path))>= ptol :return 'edipath'
                elif len(ex)/len(os.listdir(path)) < ptol : 
                    m=[]
                    try : m= [file for file in os.listdir(path) if SB.which_file(filename = os.path.join(path, file)) =='j' ]
                    except : pass 
                    if len(m)/len(os.listdir(path)) >= ptol :return 'jpath'
                        
        return 
                
def get_stationid (stations ,  station_id) :  
    """
    Tip to get station id from user by input either integer of station name .
    
    :param stations:  list of stations known
    :type stations: list
    
    :param station_id: staion expect to plot.
    :type station_id: list, str, or int 
            
    :returns: constructed list for plotting
    :rtype: array_like 
    
    :Example:
        
        >>> from pycsamt.utils import plot_utils as punc
        >>> teslist = ['S{0:02}'.format(ii) for ii in range(23)]
        >>> ss = punc.get_stationid (stations=teslist ,  station_id=('S04',13))
        >>> print(ss)
    """                
    
    sid=[]
    if not isinstance(station_id, (list, tuple)) : station_id =[station_id]
    
    for iid in station_id : 
        if isinstance (iid, str):
           for ss , stn in enumerate(stations ): 
               if iid  ==stn : sid.append(stn)
        if isinstance(iid, (int, float)): 
            if int(iid) > len(stations): 
                raise CSex.pyCSAMTError_plot_tip('maximum station number is <{0}>.'\
                                                 'Can not plot beyond this limit. Your stations go to <{1},...,{2}>'.\
                                                     format(len(stations), stations[0], stations[-1]))
                  
            sid.append(stations[int(iid) -1]) #  the stations numberred start from 1 to len(stations)

    return sid 
    
def get_frequency_id (freq_array , frequency_id ): 
    """
    function to get id of frequency . Frequency to plot
    
    :param freq_array:  array of frequency
    :type freq_array:  nd.array,1 
    
    :param frequency_id:  frequency to plot .
    :type frequency_id: list or float  
    
    :returns: new close frequency id
    :rtype: float|list
    """
    freqID =[] 
    if isinstance(freq_array, list) : freq_array=np.array(freq_array)
    if isinstance(frequency_id, (str, float, int)): frequency_id =[frequency_id]
    
    for freqid in frequency_id:
        if isinstance(freqid ,str) :
            try : freqid =np.float(freqid)
            except : raise CSex.pyCSAMTError_plot_tip('Frequency value must be an integer of float not str ' )
        if freq_array.min() > freqid > freq_array.max() : 
            warnings.warn ('Can not find the id <{0}>. ID provided is out of the range.Frequency trend to: <{1}|{2}Hz to {3}|{4}Hz>.'.\
                                                    format(freqid, np.where(freq_array==freq_array.min()),freq_array.min(),
                                                           np.where(freq_array==freq_array.max()), freq_array.max()))
            raise CSex.pyCSAMTError_plot_tip('frequency value <{0}Hz> is out of the range'.format(freqid))
        freqID.append(freqid)
        
    return freqID

def slice_matrix (base_matrix , freq_array, doi=2000):
    """
    Function get a matrix and give new matrice slice from limit value  

    :param base_matrix: arrays (yaxis lenghth, station_length)
    :type base_matrix: ndarray
    
    :param freq_array:  frequency array range 
    :type freq_array: ndarray,1 
    
    :param doi: expect to be the depth of investigation from 
                    which data muts be selected in `m` or `km` 
    :type doi:float 
    
    :returns: matrix sliced according to doi
    :rtype: ndarray 
    """
    fm=0
    if isinstance (doi, str):
        try : doi =float(doi)
        except :
            if doi.find('km')>=0 :
                ndoi = doi.replace('km','')
                fm=1
            elif doi.find('m')>=0 :ndoi = doi.replace('m','')
            try : doi =float(ndoi)
            except :
                warnings.warn ('Can not convert <{0}> into float number. '\
                               'Depth of investigation value must be float number.'.format(doi))
                raise CSex.pyCSAMTError_plot_tip('Value provided <{0}> must be float number not <{1}>'.format(doi, type(doi)))
            else :
                if fm==1 : doi *= 1e3   # convert into m
        
    doi =float(doi) # for consistency , convert to float
    indexMatrix = np.where(base_matrix >= doi)[0][0] # take the first lien that value is greater than 

    return base_matrix[:indexMatrix, ::], freq_array[:indexMatrix], doi 
    
def delineate_curve ( dict_loc , value , atol=0.2, replace_value = np.nan): 
    """
    function to delineate value of rho and phase .
    
    :param dict_loc:  dictionnary composed of keys = stations id and values
    :type  dict_loc: dict 
    
    :param value: value to delineate curve . for single value. 
    :type  value: float |list 
    
    :param atol: tolerance parameter <=1 . Most the param is closest to 0 , 
                 most the selected data become severe.default is 0.2
    :type atol: float   
                     
    :param replace_value: could be None or np.nan , Default is np.nan
    :type replace_value: float or else
                    
    :returns: dict of delineate data
    :rtype: dict

    :Example:
        
        >>> from pycsamt.utils import plot_utils as punc
        >>> ts = np.array([0.1, 0.7, 2, 3, 8, 1000, 58,55, 85, 18])
        >>> to =np.array([51, 78, 0.25, 188, 256, 7])
        >>> tt ={'S00': ts, 'S01':to, 'S02': np.array([0, 5, 125, 789])}
        >>> print(punc.delineate_curve(dict_loc = tt, value=[50,70], atol = 0.1))
    """

    if isinstance(value, (str, int, float)): 
        try : value =float(value)
        except: 
            warnings.warn('Delineation value is float number not <{0}>'.format(type(value)))
            raise CSex.pyCSAMTError_plot_tip('Input delineation value is unacceptable. Please enter a float number.')
        else : value =[value]
        
    elif isinstance(value, (list, tuple())): 
        if type(value) ==tuple : value =list(value)
        try : 
            value =sorted([float(ss) for ss in value]) 
        except : 
            warnings.warn('Delineation value is float number not <{0}>'.format(type(value)))
            raise CSex.pyCSAMTError_plot_tip('Input delineation value is unacceptable. Please enter a float number.')
    try : 
        atol =float(atol)
        if atol > 1  or atol <0: 
            warnings.warn ('Tolerance value is <=1 (less or egal than 1). Please enter again a new value.')
            raise CSex.pyCSAMTError_plot_tip('Tolerance value is assumed to be 0<= atol <=1 .Please check your value.')
    except : 
        warnings.warn('Tolerance value <{0}> provided is wrong. Must be a float number 0<= atol <=1  not a <{1}>'.format(atol, type(value)))
        raise CSex.pyCSAMTError_plot_tip('Unacceptable type tolerance value. Please enter a float number 0<= atol <=1.')
    
    # loop dictionnary an dkey value .
    delineate_dict ={}

    for value_to_find  in value :
        delineate_data_curve = []
        for keys, value_in_dict in sorted(dict_loc.items()):
            if value_to_find in value_in_dict: 
                 delineate_data_curve.append(value_to_find)      # if find , the best one , append it .
            elif  value_to_find not in value_in_dict: 
                tem=[]
                min_value = value_to_find - atol * value_to_find # find min max value closest to valueto find
                max_value = value_to_find + atol * value_to_find
                for new_turn_value  in value_in_dict :
                    if min_value <= new_turn_value <= max_value  : 
                        tem.append(new_turn_value)       
                if tem==[] : 
                    delineate_data_curve.append(replace_value)         # append None value outside the tolerance param.
                elif tem !=[]: 
                    dx = np.abs(np.array(tem) - value_to_find)
                    minimum = dx.min()
                    valueabs = np.array(tem)[np.where (dx == minimum)] 
                    delineate_data_curve.append(float(valueabs) )
       
        delineate_dict['{0}'.format(value_to_find)] =np.array(delineate_data_curve)
  
    return delineate_dict
                    
def delineate_sparseMatrix (dict_loc, delineate_dict, replace_value =np.nan ):
    """
    Build from delineate dict a matrix according to frequency length . value doesnt 
    exist in the delineate dict will be repalce by replacevalue . Default is `nan`.
    
    :param delineate_dict: delineation value 
    :type delineate_dict: dict
    
    :param dict_loc: dictionnary composed of keys = stations id and values
    :type dict_loc: dict 
    
    :param replace_value: value to replece other value  like build a sparse matrix 
    :type replace_value: float 
    
    :returns: dit of sparse matrix 
    :rtype: dict 
    """
    dictMatrix ={}
    dict_loc =sorted(dict_loc.items())
    for dd, ( keys , dvalues) in enumerate(dict_loc): 
        
        for kk, (keyi, vv) in enumerate(sorted(delineate_dict.items())[dd]): 
        #     if vv[dd] !=np.nan : 
            os= np.where(dvalues==vv[dd], dvalues, replace_value)
            dictMatrix[keyi] = os
            if vv[dd] == replace_value : 
                dictMatrix[kk]= np.repeat(replace_value, dvalues.shape[0])
            
    return dictMatrix

def resetting_ticks ( get_xyticks,  number_of_ticks=None ): 
    """
    resetting xyticks  modulo , 100
    
    :param get_xyticks:  xyticks list  , use to ax.get_x|yticks()
    :type  get_xyticks: list 
    
    :param number_of_ticks:  maybe the number of ticks on x or y axis 
    :type number_of_ticks: int
    
    :returns: a new_list or ndarray 
    :rtype: list or array_like 
    """
    if not isinstance(get_xyticks, (list, np.ndarray) ): 
        warnings.warn ('Arguments get_xyticks must be a list not <{0}>.'.format(type(get_xyticks)))
        raise CSex.pyCSAMTError_plot_tip('<{0}> found. "get_xyticks" must be a list or (nd.array,1).'.format(type(get_xyticks)))
    
    if number_of_ticks is None :
        if len(get_xyticks) > 2 : 
            number_of_ticks = int((len(get_xyticks)-1)/2)
        else : number_of_ticks  = len(get_xyticks)
    
    if not(number_of_ticks, (float, int)): 
        try : number_of_ticks=int(number_of_ticks) 
        except : 
            warnings.warn('"Number_of_ticks" arguments is the times to see the ticks on x|y axis.'\
                          ' Must be integer not <{0}>.'.format(type(number_of_ticks)))
            raise CSex.pyCSAMTError_plot_tip('<{0}> detected. Must be integer.')
        
    number_of_ticks=int(number_of_ticks)
    
    if len(get_xyticks) > 2 :
        if get_xyticks[1] %10 != 0 : 
            get_xyticks[1] =get_xyticks[1] + (10 - get_xyticks[1] %10)
        if get_xyticks[-2]%10  !=0 : 
            get_xyticks[-2] =get_xyticks[-2] -get_xyticks[-2] %10
    
        new_array = np.linspace(get_xyticks[1], get_xyticks[-2], number_of_ticks )
    elif len(get_xyticks)< 2 : 
        new_array = np.array(get_xyticks)
 
    return  new_array
        
        
def resetting_colorbar_bound(cbmax , cbmin, number_of_ticks = 5, logscale=False): 
    """
    Function to reset colorbar ticks more easy to read 
    
    :param cbmax: value maximum of colorbar 
    :type cbmax: float 
    
    :param cbmin: minimum data value 
    :type cbmin: float  minimum data value
    
    :param number_of_ticks:  number of ticks should be 
                            located on the color bar . Default is 5.
    :type number_of_ticks: int 
    
    :param logscale: set to True if your data are lograith data . 
    :type logscale: bool 
    
    :returns: array of color bar ticks value.
    :rtype: array_like 
    """
    def round_modulo10(value): 
        """
        round to modulo 10 or logarithm scale  , 
        """
        if value %mod10  == 0 : return value 
        if value %mod10  !=0 : 
            if value %(mod10 /2) ==0 : return value 
            else : return (value - value %mod10 )
    
    if not(number_of_ticks, (float, int)): 
        try : number_of_ticks=int(number_of_ticks) 
        except : 
            warnings.warn('"Number_of_ticks" arguments is the times to see the ticks on x|y axis.'\
                          ' Must be integer not <{0}>.'.format(type(number_of_ticks)))
            raise CSex.pyCSAMTError_plot_tip('<{0}> detected. Must be integer.')
        
    number_of_ticks=int(number_of_ticks)
    
    if logscale is True :  mod10 =np.log10(10)
    else :mod10 = 10 
       
    if cbmax % cbmin == 0 : 
        return np.linspace(cbmin, cbmax , number_of_ticks)
    elif cbmax% cbmin != 0 :
        startpoint = cbmin + (mod10  - cbmin % mod10 )
        endpoint = cbmax - cbmax % mod10  
        return np.array([round_modulo10(ii) for ii in np.linspace(startpoint, endpoint, number_of_ticks)])
    
    
def slice_csamt_matrix ( block_matrix , station_offsets, depth_offsets, offset_MinMax =(0, 1000),  doi='2000m' ) : 
    """
    Using Wannamaker FE elements mesh to define rho matrix blocks , need after inversion to slice the model 
    resistivity according the offset and depth we need . This function is easy tool to slice matrix and to keep
    the part we need . station offset , depth and model resistivity 
    
    Parameters
    -----------
        * block_matrix : ndarray(station_offsets.shape[0], 
                    matrix of station depth Resistivity model
                    depth_offsets.shape[0])
        
        * depth_offset : array_like  
                   depth of investigation after generating by mesh file :>z_nodes . 
                   
        * station_offsets : array_like 
                  station _offsets : offset generate by mesh_file :>x_nodes .
                  
        * offset_MinMax : tuple  
              the interval of data to keep . eg if station location start by 0 : 
                  off[0] = min and off[-1]=max (min, max):--> index 0 :
                  minimum value of station location -->index 1 : maximum value of station location
                  *default* is (0,1000)
                        
        * doi : str , float 
                investigation depth ,  migth be [m|km]. 
                If value is provided is float number , it might take value 
                as a default unit 'meter'. i.e : 1000="1000m"
                
    Returns
    ---------
        tuple 
          new sliced station offset , new sliced depth offset  , new_matrix block , 
    """
    def depth_of_investigation(doi): 
        """
        :returns: doi in meter whether you provide value in kilometer 
        :rtype:float 
    
        """
        fm=0
        if isinstance (doi, str):
            try : doi =float(doi)
            except :
                if doi.find('km')>=0 :
                    ndoi = doi.replace('km','')
                    fm=1
                elif doi.find('m')>=0 :ndoi = doi.replace('m','')
                try : doi =float(ndoi)
                except :
                    warnings.warn ('Can not convert <{0}> into float number. '\
                                   'Depth of investigation value must be float number.'.format(doi))
                    raise CSex.pyCSAMTError_plot_tip('Value provided <{0}> must be float number not <{1}>'.format(doi, type(doi)))
                else :
                    if fm==1 : doi *= 1e3   # convert into m
            
        return float(doi) # for consistency , convert to float
    
    doi= depth_of_investigation(doi)
    
    # --> get  the depth  depth indice 
    if isinstance(station_offsets, (list, tuple)) : station_offsets=np.array(station_offsets)
    if isinstance(depth_offsets, (list, tuple)) : depth_offsets= np.array(depth_offsets)
    if isinstance(offset_MinMax, (list, tuple)) : offset_MinMax= np.array(offset_MinMax)
    
    #offsets min and maximal station offsets and depth 
    offset_min,  offset_max = offset_MinMax.min(), offset_MinMax.max()#
    offset_start, offset_end = station_offsets.min(), station_offsets.max()    
    
    #check wheren the value provide is on the offset 
    if offset_min < offset_start : 
        warnings.warn ('Povided minimal offset is out of the range !Offset minimal for plottig is ={0}'.format(offset_start))
        raise CSex.pyCSAMTError_plot_tip('Povided minimal offset out of the range !Offset minimal for plottig is ={0}'.format(offset_start))
    if offset_max > offset_end  : 
        warnings.warn ('Povided maximal offset is out of the range !Offset maximal for plottig is ={0}'.format(offset_end))
        raise CSex.pyCSAMTError_plot_tip('Povided maximal offset is out of the range !Offset maximal for plottig is ={0}'.format(offset_end))
    
    if doi > depth_offsets.max() : 
        mess ='depth <{0}> is out of depth range. Maximum depth for  is={1}.'.format(doi, depth_offsets.max())
        warnings.warn(mess)
        raise CSex.pyCSAMTError_plot_tip(mess)
    # make sure that station offset and and depth are not  in decrease oder
    sflip,dflip=0,0
    if depth_offsets[0] > depth_offsets[-1] : 
        depth_offsets =depth_offsets [::-1]
        sflip =1
    if station_offsets [0] > station_offsets [-1] : 
        station_offsets= station_offsets[::-1]
        dflip=1
        
    # check whether index values are present on depth array 
    ind_depth_s, = np.where(depth_offsets==0.)
    ind_depth_e, = np.where(depth_offsets==doi)

    #if not : 
    if len(ind_depth_s)==0 : 
        for ii, dd in enumerate(depth_offsets): 
                if dd > 0. : 
                    dmax , dmin = abs(dd-0.), abs(depth_offsets[ii-1]) 
                    if dmax < dmin : 
                        ind_depth_s =ii 
                        break
                    else :
                        ind_depth_s = ii-1
                        break
    if len(ind_depth_e) ==0 : 
        for ji, de in enumerate(depth_offsets): 
                if de > doi : 
                    dmax , dmin = abs(de-doi), abs(depth_offsets[ji-1]-doi) 
                    if dmax < dmin : 
                        ind_depth_e = ji 
                        break
                    else : 
                        ind_depth_e = ji-1
                        break

      #--> get for the station offset
    #checker whether station index are presents 
    ind_offs_s, = np.where(station_offsets==offset_min)
    ind_offs_e, = np.where(station_offsets==offset_max)
     
    # if not : 
    if len(ind_offs_s)==0 : 
        for jj, offs in enumerate(station_offsets): 
                if offs > offset_min : 
                    #claculate the difference betwen 
                    maxoff = abs(offs- offset_min) 
                    minoff = abs(station_offsets[jj-1]- offset_min)
                    if maxoff < minoff:
                        ind_offs_s = jj
                        break
                    else : 
                        ind_offs_s = jj-1
                        break
    if len(ind_offs_e) ==0 : 
         for je, offs in enumerate(station_offsets): 
                if offs> offset_max : 
                    #claculate the difference betwen 
                    maxoff = abs(offs-offset_max) 
                    minoff = abs(station_offsets[je-1]-offset_max)
                    if maxoff < minoff: 
                        ind_offs_e = je
                        break
                    else :
                        ind_offs_e = je -1 
                        break
     
    # convert to integer and leave array for indexing 
    ind_depth_s, ind_depth_e = int(ind_depth_s), int(ind_depth_e) 
    ind_offs_s, ind_offs_e= int(ind_offs_s), int(ind_offs_e)
    
    # for add +1 for the last index so to take the last value
    ind_offs_e +=+1
    ind_depth_e +=1 
 
    new_station_offsets= station_offsets[ind_offs_s: ind_offs_e]
    new_depth_offsets= depth_offsets[ind_depth_s:ind_depth_e]

    # Flip indexing if values array are sorted in highest to lower (depth even stations locations)
    if sflip ==1 :
        # flips station -index 
        ind_offs_s = - ind_offs_e
        ind_offs_e = - ind_offs_s
        new_station_offsets= new_station_offsets[::-1]
   

    if dflip ==1 : 
        # flips depth-index 
        ind_depth_s = - ind_depth_e
        ind_depth_e = - ind_depth_s
        new_depth_offsets=new_depth_offsets[::-1]
        

        
    new_block_matrix = block_matrix[ind_depth_s:ind_depth_e, ind_offs_s:ind_offs_e] 
    
    return   (new_station_offsets,  new_depth_offsets, new_block_matrix)
            
def controle_delineate_curve(res_deline =None , phase_deline =None ): 
    """
    fonction to controle delineate value given  and return value ceilling .
    
    :param  res_deline:  resistivity  value todelineate. unit of Res in `ohm.m`
    :type  res_deline: float|int|list  
    
    :param  phase_deline:   phase value to  delineate , unit of phase in degree
    :type phase_deline: float|int|list  
    
    :returns: delineate resistivity or phase values 
    :rtype: array_like 
    """
    fmt=['resistivity, phase']
 
    for ii, xx_deline in enumerate([res_deline , phase_deline]): 
        if xx_deline is  not None  : 
            if isinstance(xx_deline, (float, int, str)):
                try :xx_deline= float(xx_deline)
                except : raise CSex.pyCSAMTError_plot_tip('Value <{0}> to delineate <{1}> is unacceptable.'\
                                                      ' Please ckeck your value.'.format(xx_deline, fmt[ii]))
                else :
                    if ii ==0 : return [np.ceil(np.log10(xx_deline))]
                    if ii ==1 : return [np.ceil(xx_deline)]
  
            if isinstance(xx_deline , (list, tuple, np.ndarray)):
                xx_deline =list(xx_deline)
                try :
                    if ii == 0 : xx_deline = [np.ceil(np.log10(float(xx))) for xx in xx_deline]
                    elif  ii ==1 : xx_deline = [np.ceil(float(xx)) for xx in xx_deline]
                        
                except : raise CSex.pyCSAMTError_plot_tip('Value to delineate <{0}> is unacceptable.'\
                                                      ' Please ckeck your value.'.format(fmt[ii]))
                else : return xx_deline

def depth_of_investigation(doi): 
    """
     Depth of investigation converter 

    :param doi: depth of investigation  if value float is provided , it will considered as 
                default units in meter
    :type doi: str|float 
    
    :returns doi:value in meter
    :rtype: float           
    """
    fm=0
    if isinstance (doi, str):
        try : doi =float(doi)
        except :
            if doi.find('km')>=0 :
                ndoi = doi.replace('km','')
                fm=1
            elif doi.find('m')>=0 :ndoi = doi.replace('m','')
            try : doi =float(ndoi)
            except :
                warnings.warn ('Can not convert <{0}> into float number. '\
                               'Depth of investigation value must be float number.'.format(doi))
                raise CSex.pyCSAMTError_plot_tip('Value provided <{0}> must be float number not <{1}>'.format(doi, type(doi)))
            else :
                if fm==1 : doi *= 1e3   # convert into m
        
    return float(doi) # for consistency , convert to float


def build_new_station_id (station_id , new_station_name ): 
    """
    Fonction to build new station id including station name provided , if the length provided 
        doesnt match the length of station id 
        
    :param station id: new sites names 
    :type station id: list 
    
    :param mess:  message for debugging, Default is None 
    :type  mess: str  
             
    :Example:
        
        >>> from pycsamt.utils import plot_utils as punc
        >>> ts =[   28.,  200.,  400.,  600.,  800., 1000., 1200., 1400., 1600., 1800., 1807.]
        >>> sto = ['S{0}'.format(i) for i in range(7)]
        >>> print(sto)
        >>> stn, fu =punc.build_new_station_id(station_id = sto, new_station_name =ts)
        >>> print(stn, fu)
    """
    mess=None
    get_len_new_station_id =len(new_station_name)
    if not isinstance(new_station_name, (list, tuple, np.ndarray)) : 
                new_station_name=[str(new_station_name)] # convert value to string 
    else : 
        new_station_name=[str(stn) for stn in new_station_name] # for consistenncy rebuild sites name and keep on str value 
    if get_len_new_station_id <  len(station_id ): 
        mess ='Length of sites names provided ={0}. It doesnt match the length of  station id ={1}.'\
            ' new sites names will build including your sites names provided.'.\
                format(len(new_station_name), len(station_id ))
        warnings.warn (mess)

        
        station_id =  new_station_name + station_id[get_len_new_station_id:]
    
    elif get_len_new_station_id > len(station_id) : 
        mess ='Length of sites names provided ={0}. It doesnt match the length of  station id ={1}.'\
            ' new sites names will be reduced according the length of station id.'.\
                format(len(new_station_name), len(station_id ))
        warnings.warn (mess)
        
        station_id = new_station_name[:len(station_id) ]

    else : station_id = new_station_name 
    
    return [str(id_) for id_ in station_id] , mess

def get_conductive_and_resistive_zone (data, site_names, purpose ='groundwater',**kws): 
    
    """
    function to get the probability of conductive and resistive zone . It is not absolutely True 
    but give an overview of decison . It is not sufficient to declare that the zone is 
    favorable for any drill , but just work with probability 

    :param data: resistivity datata of survey area 
    :type data: ndarray 
    
    :param site_name: list of sites names 
    :type site_name: list 
    
    :param purpose: type of exploration , default is `groundwater` 
    :type purpose: str 
    
    :returns: report of exploration area 
    :rtype: str 
    """
    def fmt_text (data_text, fmt='~', leftspace = 3, return_to_line =77) : 
        """
        Allow to format report with data text , fm and leftspace 
        :param  data_text: a long text 
        :type  data_text: str  
            
        :param fmt:  type of underline text 
        :type fmt: str

        :param leftspace: How many space do you want before starting wrinting report .
        :type leftspae: int 
        
        :param return_to_line: number of character to return to line
        :type return_to_line: int 
        """
    
        return_to_line= int(return_to_line)
        begin_text= leftspace *' '
        text= begin_text + '~'*(return_to_line +7) + '\n'+ begin_text
    
        ss=0
        
        for  ii, num in enumerate(data_text) : # loop the text 
            if ii == len(data_text)-1 :          # if find the last character of text 
                #text = text + data_text[ss:] + ' {0}\n'.format(fmt) # take the remain and add return chariot 
                text = text+ ' {0}\n'.format(fmt) + begin_text +'~'*(return_to_line+7) +'\n'         # add a fmt strinmg to finish 
     
                break 
            if ss == return_to_line :                        # if enconter 50character then return chariot 
                if data_text[ii+1] !=' ' : 
                    text = '{0}- {1} \n {2} '.format( text,  fmt, begin_text + fmt ) #     if the word is brken then add "-" and return chariot 
                else : 
                    text ='{0} {1} \n {2} '.format(text, fmt, begin_text+fmt ) 
                ss=0
            text += num    # add charatecter  
            ss +=1
    
        return text 
    
    # other functions arguments
    
    model_offsets=kws.pop('model_offsets', None)
    site_offsets = kws.pop('site_offsets', None)
    
    
    purpose_value = np.log10(
        np.mean(
                [3.1e2, 1e0]))  # we assume the resistivity for groundater is ranged
                                        # between [3.1e3 1e0] = elecetrical props of fresh water
    
    average_rho_data =data.mean(axis = 0) # return an ndarray of average _rho  
    mins= average_rho_data.min()

    # imins,  =np.where(average_rho_data==mins) # get the position of minimum value
    imins = average_rho_data.argmin()

    maxs =average_rho_data.max()
    # imaxs, =np.where(average_rho_data==maxs) # get the position of maximum  value 
    imaxs = average_rho_data.argmax()
 
    # get the average value on the set of each line 
    mean_set_= average_rho_data.mean()         # find the average value ....
    
    # mean
    # check the mimum and compare to purpose value 
     # print('{0:=^77}'.format('Occam Response plot infos'))
     
    print('---> Average Rho on survey area is = {0} {1}'.format(np.power(10, mean_set_), 'Ω.m.'))
    
    if  imins < len(site_names) :
        print('---> Probably very conductive zone is = {0} with rho = {1} Ω.m.'.\
              format(site_names[int(imins)], np.power(10,mins)))
    if  imaxs < len(site_names) :
        print('---> Probably very resistive zone is = {0} with rho = {1} Ω.m.'.\
              format(site_names[int(imaxs)], np.power(10, maxs))) 
    
    # make a ratio minimum value and maxmim ration 
    temcon, temres=[[] for nn in range(2)] # build empty list for min sites names and 
    for jj, (index, iindex) in enumerate(zip ([imins, imaxs], [mins, maxs])):
        if index > len(site_names): # if the value is greater than sites names length , we assume that is them occam 
        # model objet then bring model file 
            if model_offsets is None  or site_offsets is None: 
                mess =''.join([' !Could not find Occam 2D model resistivity value . Could not provide a report.', 
                               'Too many minimum value find ! Need Occam model resistivity data and model station data', 
                               'to provide a report.'])
                warnings.warn(mess)
        
            icc, offs_value = find_closest_station(offset_indice= index , 
                                                     model_offsets = model_offsets,
                                                     site_offsets =site_offsets)
            if jj ==0 : 
                snn = site_names[int(icc)] # build conductive ssites names 
                if snn not in temcon:   # sometime the names of station is close , keep one , if not keep all 
                                        # else pass 
                    temcon.append((snn, iindex))  # keep station 
                    
            if jj ==1 : 
                snn = site_names[int(icc)] # build conductive ssites names
                if snn not in temres:              
                    temres.append((snn, iindex))
        
    if len(temcon) >0 : 
        for stn, rhoval  in temcon : 
            if len(temcon)>1 : # if many conductive zone has ben selected then : 
                print('--> {0} Probably  conductive zones has been selected :'.format(len(temcon)))
            
            print('---> Probably very conductive zone is = {0} with rho = {1} Ω.m.'.\
                  format(stn, np.power(10,rhoval)))
                
            imins =site_names.index(stn) # get the real minimum and the value of resistivities
            mins = rhoval
    if len(temres) > 0: 
        for stn, rhoval in temres : 
            if len(temres)>1 : # if many conductive zone has ben selected then : 
                print('--> {0} Probably  resistives zones has been selected :'.format(len(temres)))
            print('---> Probably very resistive zone is = {0} with rho = {1} Ω.m.'.\
                  format(stn, np.power(10,rhoval)))
                
            imaxs =site_names.index(stn)
            maxs = rhoval
   
            
    print('- However :')
    print('---> {0:<27} {1} {2}. '.format('Minimum ratio is',  '=',round( mins/mean_set_,7)))
    print('---> {0:<27} {1} {2}.'.format('Maximum ratio is',  '=', round(maxs/mean_set_, 7)))
    
    if purpose =='groundwater':
        p_ratio = (purpose_value / mins)*100        # presumed a prefect zone (conductive)
        # p_ratio_max : maximal trend of condictive zone 
        #p_ratio _min : minimu trend of conductive zone : 
            
        p_ratio_max, p_ratio_min  = np.log10(3.1e2)/mins, np.log10(1e0) / mins 
         
        # pseudo_trust = (mins/mean_set_ / p_ratio) * 100
        pseudo_trust = (mean_set_/ mins) *100 
        if pseudo_trust >  np.sqrt(2)* p_ratio_max *100 : 
             mess = ''.join([' resistive medium and could be ',
                             '  hydrogeologically UNfavorable zone.' ])
             fm='However,'
             # mess = 'hydrogeologically UNfavorable'
             
        elif pseudo_trust < np.sqrt(2)* p_ratio_min *100  : 
            mess ='medium too conductive and must be aware of that zone.'
            fm='However,'
            
        else : 
             mess ='hydrogeologically favorable zone and could be a right drilling point.' 
             fm=''
        
        print()
        print( '!IMPORTANT NOTE! :')
        bfreport =''.join(['This present report is',
                           '  for {0}  exploration purpose.',  
                           'If your survey purpose if far from',
                           ' {0} exploration , please ignore it.'])
        bfreport =bfreport.format(purpose)
        
        fttext = fmt_text(data_text = bfreport, fmt='-', leftspace=3)
        
        print(fttext)
        
        print( 'Report:')

        report=''.join([
            'In the case of {0} exploration, along the survey line,',
            ' site {1} seems to be the best conductive zone with reliability ',
            ' ratio estimates to  = {2} % compared to the threshold = {3} %.',
            ' {4} Site {1} is considered as {5}. '
                    ])
        report= report.format(purpose, 
                              site_names[int(imins)],
                              np.around(pseudo_trust, 2), 
                              np.around(p_ratio,2),
                              fm, 
                              mess)
        caution =''.join(['The present  report is NOT absolute. It might be',
                          ' taken with cautiousness because in groundwater ', 
                          ' exploration, many other parameters are considering',
                          ' before a possible drilling decison. We hope in  ', 
                          'the next release, introduce an intelligent drill ', 
                          'estimator `StarGain` to give better precision of ',
                          ' absolute point  with tolerance value by considering', 
                          '  multiple other parameters. '])
        
        # generate a reading report 
        textr = fmt_text(data_text = report + caution )
        print(textr)

        
def fmt_text (data_text, fmt='~', leftspace = 3, return_to_line =77) : 
    """
    Allow to format report with data text , fm and leftspace 

    :param  data_text: a long text 
    :type  data_text: str  
        
    :param fmt:  type of underline text 
    :type fmt: str

    :param leftspae: How many space do you want before starting wrinting report .
    :type leftspae: int 
    
    :param return_to_line: number of character to return to line
    :type return_to_line: int 
    """

    return_to_line= int(return_to_line)
    begin_text= leftspace *' '
    text= begin_text + fmt*(return_to_line +7) + '\n'+ begin_text

    
    ss=0
    
    for  ii, num in enumerate(data_text) : # loop the text 
        if ii == len(data_text)-1 :          # if find the last character of text 
            #text = text + data_text[ss:] + ' {0}\n'.format(fmt) # take the remain and add return chariot 
            text = text+ ' {0}\n'.format(fmt) + begin_text +fmt*(return_to_line+7) +'\n'         # add a fmt strinmg to finish 
 
            break 
        if ss == return_to_line :                        # if enconter 50character then return chariot 
            if data_text[ii+1] !=' ' : 
                text = '{0} {1}- \n {2} '.format( text,  fmt, begin_text + fmt ) #     if the word is brken then add "-" and return chariot 
            else : 
                text ='{0} {1} \n {2} '.format(text, fmt, begin_text+fmt ) 
            ss=0
        text += num    # add charatecter  
        ss +=1

    return text 

def find_closest_station( offset_indice ,  model_offsets , site_offsets )  : 
    """ 
    Get the indice of the closest offset 
    
    Parameters
    -----------
        * offset_indice: int 
             if the indix of the offset at the selected resistivity point. 
        * model_offset: array_like  
            is a large band of x_nodes resistivities  gnerated by mesh files 
        
        * sites_offsets: array_like 
             the data set offset from Occam Data file 
        
    Returns 
    --------
        indexoff : int
             index of data offset 
        get_offs : float
               value of the offset at that index 
    """
    offs_value= model_offsets[int(offset_indice)] #  keep the offset value form model 

    for  ii , iof in enumerate(site_offsets) :  # loop the data offset and find the offset in model more close to data offset
        if iof > offs_value :  # compared both offsets values 
            minabs = abs(iof - offs_value) # find the distance between offset from model and offset from ocam data 
            try : 
                maxabs = abs (site_offsets[ii-1]-offs_value) # if ii -1 not possible to find then take absolue distance from ... 
            except : 
                maxabs= abs(iof -offs_value)            # ... from data offsets 
            if minabs <= maxabs :  # if minimum offset is found at index ii of offset iof than take the index and the iofs value 
                indexoff, get_offs  = ii, iof
                break
            else :  
                indexoff, get_offs = ii-1, site_offsets[ii-1] # if ii-1 exist then get the previous offset from data offset .
                break
      # return indice of data offset so to find the sites names and the offset at that sites  index 
    
    return indexoff, get_offs
    
def find_local_maxima_minima (array): 
    """
    function to find minimum local and maximum local  on array 
    
    :param data: value of array to find minima , maxima 
    :type data: array_like 
    
    :returns: tuple of index of minima  and maxima local index and 
             array of minima maxima  value
    :rtype: array_like
    
    :Example: 
        
        >>> from pycsamt.utils import func_utils as func
        >>> from pycsamt.utils import plot_utils as punc
        >>> ts =np.array([2, 3, 4, 7, 9, 0.25, 18, 28, 86, 10, 5])
        >>> te=np.array([0.2, .3, 18, 1.6, 0.2, 0.6, 0.7, 0.8, 0.9, 1., 23.])
        >>> th =func.concat_array_from_list([ts, te], concat_axis=1)
        .. tth =th.T
        ... print(tth)
        ... print(punc.find_local_maxima_minima(tth[0]))
    """  
    localmm, indexlmm=[[] for ii in range(2)]
    
    error_mess ='Can not convert value to float number. Arry to find locals'\
                ' minmum and maximum must be an aray of float numbers.'
    if isinstance(array, (list, tuple)): # put on list 
        try :  array=np.array([float(ss) for ss in array])
        except : raise CSex.pyCSAMTError_parameter_number(error_mess)
            
    if isinstance(array, (float, str, int)):
        try :
            array =np.array([float(array)])
        except : CSex.pyCSAMTError_inputarguments(error_mess)
    
    if len(array) == 1 : 
        return (np.array([0]), array )
            
    d0, f=array[0], 0               # get the initial value , use the initial at starting points 
    localmm.append(d0)              #... can be either local maximum or local minimum 
    indexlmm.append(0)              # append index of starting point d0
    for ii, value in enumerate(array):  # let go for maximum data 
            if value >= d0 : #0.3       # try to find maximum |speudo maximum but get the minimum index when
                if f==2 :
                    localmm.append(array[ii-1])
                    indexlmm.append(ii-1)
                f=1
                d0=value                # reinitialize pseudo maximum value 
            elif value < d0 :           # # try to find minimum|speudo minimum but get the maximum index when
                if f== 1 : 
                    localmm.append(array[ii-1])
                    indexlmm.append(ii-1)
                    
                f=2
                d0=value                # # reinitialize pseudo minimum value 
            if value ==array[-1]:       # at the end take the last value to close the loop whether is minimum or maximum 
                localmm.append(value)  
                indexlmm.append(ii)
                
    return np.array(indexlmm) ,np.array(localmm)            
        
def average_rho_with_locals_minmax(array): 
    """
    How to compute mean value between local maxima and local minima 
    and keep locals minima and maxima value on the final data 
    
    :param array: data to compute the local minima and local maxima 
    :type array: array_like 
    
    :returns: array mean with data local value averaged 
    :rtype: array_like
    
    :Example: 
        
        >>> from pycsamt.utils import plot_utils as punc
        >>> mean1= punc.average_rho_with_locals_minmax(tth[0])
        >>> mean2= punc.average_rho_with_locals_minmax(tth[1])
        >>> print(mean1)
        >>> print(mean2)
    """
    tem=[]
    indexmm, _= find_local_maxima_minima(array) # find the local and maxium local index
    
    if len(indexmm) == 1 :              # once it as a singe value  
        return array[int(indexmm)]      # for consistency
    
    s_index =indexmm[1]                 #   next index as str
    
    for ii, index in enumerate(indexmm): #[ 0,  4,  5,  8, 10] #s_index = 4
        if index == indexmm[-1]:
                                    
            tem.append(np.array(array[index]))    # put single value on array for conacatenation 
            break
        
        if ii !=0 :                     # start a counter at index >0
            s_index = indexmm[ii+1]
        
        ndiff = s_index - index         # take the next index and substract the present index s_index 
                                        # ndiff = 4-0 (indexmm[1]-index (indexmm[0])) = 4 >1 
        if ndiff == 1 :
            tem.append(np.array(array[index]))   #add local minimum array[0]= 2
            
        if ndiff >1 : 
            mean_mm = array[index:s_index+1].mean() # [0:4+1]take the mean between [index0 : index 9 +1]--> [0:5]
            #moy = array[ii :indexmm[ii+1]] = array[0]-array[5]= moy (array[2:5] 2:9 /mean
            tem.append(np.array(array[index]))
            nn=np.repeat(mean_mm, ndiff-1) # 25/5 =5 ( 5, 3) 2 5 5 5 9 (ndiff =4-0-1)=3
            tem.append(nn)
               
    return np.hstack(tuple(tem))

def average_rho_in_deeper (dep_array, rho_array, step_descent ) :
    """
    function to average rho in deep according to the value provided . In fact 
    averaged rho in shorter depth distance allow us to understand the conductive zone. 
    in approximately. The most conductive zone is detected as the zone with 
    lower resistivities values . But fixing values as averaged rho , can build a 
    specific strata that could match this zone .
    
    Parameters
    ------------
        * dep_array : array_like 
                the imaged depth (doi)
            
        * rho_array: array_like
                resistivity array  
            
        * step_descent : float 
                value to step descent 
    
    Returns
    ---------
        array_like 
             rho average for each station  # dep_averaged for each station

    :Example:
        
        >>> import numpy as np
        >>> from pycsamt.utils import plot_utils as punc
        >>> pseudo_depth=np.array([  0. ,  6. , 13. , 20.  ,29. , 39.,  49. ,
        ...                       59. , 69. , 89., 109. ,129., 149. ,179.,
        ...                  209. ,249. ,289., 339., 399., 459. ,529. ,609., 
        ...                  699., 799., 899., 999.])
        >>> pseudo_depth = np.arange(0, 1220, 20)
        >>> rho = np.random.randn(len(pseudo_depth))
        >>> rho_aver, dep_aver= average_rho_in_deep(dep_array=pseudo_depth,
        ...                                        rho_array=rho, 
        ...                                        step_descent=20.)
        >>> rho_2,dep_2 =   average_rho_in_deeper (dep_array= pseudo_depth,
        ...                                         rho_array=rho,
        ...                                         step_descent=1000) 
        >>>  print(pseudo_depth)
    """
    v,r , dm, rm=[[] for i in range(4)]
    value = step_descent
    
    for index , depth in enumerate(dep_array):

        if depth <= step_descent :      # value less than step descent must be averaged 
            v.append(depth)             # keep resistivities values onto list 
            r.append(rho_array[index])

        if depth > step_descent :       # if the next value is greater than the previous one 
        
            if v !=[]:                  # ccheck if kist is not enmpty 
            
                dm.append(np.repeat(np.array(v).mean(), len(v)))  # rebuild resistivities values with rho averaged 
                rm.append(np.repeat(np.array(r).mean(), len(r)))
                step_descent += value               #increment the next descent to step of descent 
                v=[depth]       # initialise new list by adding the index value greater one 
                r=[rho_array[index]]
              
        if depth ==dep_array[-1]:
            if len(v)==1 :                  # it length last value ==1 , means is the last value of depth
                dm.append(dep_array[index])
                rm.append(rho_array[index])
            elif len(v) !=1 :               # averaged the reamin rho values  
                dm.append(np.repeat(np.array(v).mean(), len(v)))
                rm.append(np.repeat(np.array(r).mean(), len(r)))
      

    return np.hstack(tuple(rm)),np.hstack(tuple(dm))


@deprecated("function deprecated , it's does not cover the total depth when last value is outside the model depth range. ")
@redirect_cls_or_func(average_rho_in_deeper,"Function replaced , it redirects to {average_rho_in_deeper}.")
def average_rho_in_deep (dep_array, rho_array, step_descent)   :
    """
    function to average rho in deep according to the value provided . In fact 
    averaged rho in smalled depth distance allow us to understand the conductive zone. 
    in approximated. The most conductive zone is detected as the zone with 
    lower resistivities values . But fixing values as averaged rho , can build a 
    specific strata that could match the zone .
    
    Parameters
    -----------
        * dep_array : array_like
                the imaged depth (doi)
        * resistivity : array_like 
                resistivity ate the station 
        * step_descent : float 
                value to step descent 
    Returns 
    --------
    array_like 
         rho average for each station # dep_averaged for each station
  
    :Example: 
        
        >>> from pycsamt.utils import plot_utils as punc
        >>> pseudo_depth = np.arange(0, 1201.,100)
        >>> print(pseudo_depth)
        >>> dep_aver= dep(dep_array=pseudo_depth, value=100)
        >>> print(dep_aver)
        >>> print(pseudo_depth.shape)
        >>> print(dep_aver.shape)
    """
    tem, rm=[[] for i in range(2)]
    step_descent  =depth_of_investigation(step_descent)
    
    if isinstance(step_descent, (str, int)): 
        try : 
            step_descent =float(step_descent)
        except : 
            raise CSex.pyCSAMTError_plot_geoinputargument(
                'Could not convert depth value ={} to float. Please check your value.'.format(step_descent))
    
    
    if step_descent < dep_array.min(): 
        raise CSex.pyCSAMTError_inputarguments(
            'Value provided ={0} m is less than the minimum depth ={1} m.'.format(step_descent, dep_array.min()))
    
    if step_descent > dep_array.max(): 
        raise CSex.pyCSAMTError_inputarguments(
            'Value provided is = {0} m is greater than maximum depth ={1}m.'.format(step_descent, dep_array.max()))
        
    mm=0 
    _init_depth =step_descent
  
    for index, depth in enumerate(dep_array):  # index = 0

        if depth >= step_descent : # if value is greater than the maximum depth for the fist time 
            dep_averaged = dep_array[mm:index].mean() #[0 :20] 20 is ouut : 0--19 
            rho_averaged=  rho_array[mm:index].mean()
            tem.append(np.repeat(dep_averaged, index - mm)) #3-0=3
            rm.append(np.repeat(rho_averaged, index - mm)) 
            mm= index               # mm= 3 (20)
           
            step_descent += _init_depth    # value = 20+20 = 40
    
            
            
        if depth == dep_array[-1]:               # last value then check if value is modulo the init depth 
            
            if depth % _init_depth == 0 or depth < step_descent : # assume the last depth is less than the next assigned depth 

                dep_averaged = dep_array[mm:].mean()
                rho_averaged=  rho_array[mm:].mean()
                tem.append(np.repeat(dep_averaged, index -mm +1)) # add +1 because of python nanems start count from 0
                rm.append(np.repeat(rho_averaged, index -mm +1))
        

        
    return np.hstack(tuple(rm)), np.hstack(tuple(tem))
    

def get_station_id_input_resistivities(station_rho_value, number_of_layer=None): 
    """
    Get a special station input resistivities is much benefit and 
    more close to reality of plot . Indeeed , it take only the maximum value of resistivities 
    below the site and the minimum , then cut out 7  Considering all the model data and 
    choose the max resistivities and the minim resistivities to build automatic 
    resistivities whom COULD match the deth is less sure .
    Geeting a input resitivities to aplom the site , give a merly and better interpretation. 

    :param station_rho_value: value of resistivities under the 
                                site thin the maximum depth 
    :type station_rho_value: array_like
    
    :param number_of_layer: number of layer to top to bottom.
    :type number_of_layer: int 
    """
    if number_of_layer is None : number_of_layer =7.
    if isinstance(number_of_layer, (list, tuple, np.ndarray)): 
        try : 
            number_of_layer =int(number_of_layer[0])
        except : 
            raise CSex.pyCSAMTError_plot_tip(
                'Could not convert number of  layer ={0} value into integer.'.format(tuple(number_of_layer)))
            
    if isinstance(number_of_layer, (float, str)): 
        try : int(number_of_layer)
        except : CSex.pyCSAMTError_plot_tip(
                'Could not convert number of layer ={} to integer.'.format(number_of_layer))
        else : number_of_layer=int(number_of_layer)
        
    if isinstance(station_rho_value, (list, tuple)): 
        try : 
            station_rho_value= np.array(station_rho_value)
            if station_rho_value.dtype not in ['float', 'int']: 
                station_rho_value= station_rho_value.astype('float64')
                
        except : raise CSex.pyCSAMTError_plot_geoinputargument(
                'Could not convert resistivities values '
                'to float number. Please provide a right values.')
                  
                                                               
    
    # get maxrho and min rho 
    return np.linspace(station_rho_value.min(), station_rho_value.max(), number_of_layer)
        
def build_resistivity_barplot(depth_values , res_values): 
    """
    Allow to build bar plot resistivity function of  investigation depth .
    
    Parameters
    ----------
        * depth_values : aray_like 
                model investigation depth 
                
        * res_values : array_like 
                model_resistivities at each depth values 
            
    Returns 
    ---------
        d : array_like 
            resistivity barplot depth .
        r : array_like 
               specific structure resistivities 
        sumd : float
            checker number that cover in fact the total depth. 
            this numbe must absolutely match the total depth .

    """
    mess = "".join(["Depth and resistivity arrays have no the same Length .",
                    " Length Depth is ={0} and length depth is = {1}".format(len(depth_values), len(res_values))
                    ])
    if len(depth_values) !=len(res_values) : 
        warnings.warn(mess)
        raise CSex.pyCSAMTError_geodrill_inputarguments(
            'Error size. Depth and resistivity must have the same length.')
    
    # lop the resistivities array 
    v,r=[],[]
    # mm = res_values [0]
    for ii , res in enumerate(res_values):
        if res == res_values [-1] : # if res egal the last value s, check the previous 
            if res_values [ii] != res_values[-2] : # then append the last depth index to close the loop 
                v.append(depth_values[-1])
                r.append(res)
                break
            if res == res_values [-2] : # then clean the last index  of depth and replace it by the last one 
                # v[-1]= depth_values [-1]
                v.append(depth_values [-1])
                r.append(res)
                
        elif res !=res_values[ii+1] : 
            v.append(depth_values [ii])
            r.append(res)
            
    # thens substract the difference between depth values so all the sum must be egal to the maximum depth 
    d=[v[0]]        # because the minimum depth start by 0 , then append  first value depth 
    
    for dd , depth in enumerate(v): 
        if dd >0 : 
            d.append(depth - v[dd-1]) # 
            
    # check the sum if it covers the total depth .
    d= np.array(d)
    sumd = d.sum()
    # sanitize array and eliminate all negative values and zeros 
    d[d<0]=0.           # then sanitize 
    d =d[d!=0]
    
    return d, r, sumd

def annotate_tip(layer_thickness , layer_names): 
    """
    A tip to group text with the same resistivities one layer when the layer 
    are successively the same .
    
    :param layer_thickness: thickness of layer
    :type layer_thickness: array_like |list 
    
    :param layer_names: names of layers , geological structures names 
    :type layer_names: list or array_like 
    
    :returns: v ,  layer  thickness in depth 
    :rtype: array_like 
    
    :returns: ni ,  list : name of layer 
    :rtype: array_like 
  
    :Example: 
        
        >>> from pycsamt.utils import plot_utils as punc
        >>> rocks = ['Massive sulfide', 'Igneous rocks', 'Igneous rocks', 'Igneous rocks',
             'Igneous rocks', 'Igneous rocks', 'Igneous rocks', 'Igneous rocks', 
             'Massive sulfide', 'Igneous rocks', 'Igneous rocks', 'Igneous rocks'] 
        >>> resprops = [0.0, 49.0,69.0,89.0,109.0,129.0,
                        149.0,179.0,249.0,699.0,799.0,899.0]
        >>> thickness, lnames = punc.annotate_tip(layer_thickness=resprops, layer_names=rocks)
        >>>print(thickness)
        >>> print(lnames)
        >>> v, ni
        ...  [24.5, 149.0, 474.0, 799.0]
        ...  ['Massive sulfide', 'Igneous rocks', 'Massive sulfide', 'Igneous rocks']
    """
    if len(layer_thickness)==1 : 
        return layer_thickness, layer_names
    if len(layer_thickness)==0 : 
        mess = 'None value is provided, Please specify at least the thickness of one truth layer.'
        warnings.warn(mess)
        raise CSex.pyCSAMTError_plot_tip(mess)
    
    top = layer_thickness[0]
    bottom = layer_thickness[1]
    ni=[layer_names[0]]
    v=[]
    for ii, thickvalues in enumerate(layer_thickness): 
        if ii != 0 :
            if ii == len(layer_thickness)-1 :  # controle whether the layer is the same name with previous ones 
                if layer_names[ii].lower()== layer_names[ii-1].lower() : 
                    mean = top + (thickvalues - top )/2 # we reached the bottom (top +center)
                    v.append(mean)
                    
                if layer_names[ii].lower() != layer_names[ii-1].lower() :  # then do the mean with previous one 
                    mean =top + (layer_thickness[ii-1]-top)/2 # mean for the previous layers 
                    v.append(mean)
                    mean = layer_thickness[ii-1] + (thickvalues- layer_thickness[ii-1])/2  #top +center 
                    v.append(mean)
                    ni.append(layer_names[ii])                # then append the last names 
            
            elif layer_names[ii].lower() != layer_names[ii-1].lower() : # if the name is the same that the previous values 
                bottom = thickvalues # the bottom of previous layer  is the top of layer ii
                mean = (bottom - top )/2 + top         # center the annotate text (top + center)
                v.append(mean)                  # keep values and layer names 
                ni.append(layer_names[ii])                # append the name of layer 
                top =bottom                     # top becomes the bottom 


    return v, ni  
    
    
   

# if __name__=='__main__':
    
#     RES=np.array([1.200000e-01, 1.200000e-01 ,1.200000e-01, 1.200000e-01, 1.200000e-01,
#              1.200000e-01, 1.200000e-01, 1.138263e+04, 1.138263e+04, 2.276514e+04,
#              3.414765e+04, 4.553016e+04 ,3.414765e+04, 2.276514e+04, 1.138263e+04,
#              1.138263e+04, 1.200000e-01, 1.200000e-01, 1.200000e-01, 1.200000e-01,
#              1.200000e-01, 1.200000e-01, 1.200000e-01, 1.138263e+04 ,2.276514e+04,
#              6.829517e+04])
#     pseudo_depth=np.array([  0. ,  6. , 13. , 20.  ,29. , 39.,  49. , 59. , 69. , 89., 109. ,129., 149. ,179.,
#                       209. ,249. ,289., 339., 399., 459. ,529. ,609., 699., 799., 899., 999.])
    
#     barT,bar_res,  sumd = build_resistivity_barplot(depth_values=pseudo_depth ,
#             res_values=RES)
    
 
#     rocks = ['Massive sulfide', 'Igneous rocks', 'Igneous rocks', 'Igneous rocks',
#              'Igneous rocks', 'Igneous rocks', 'Igneous rocks', 'Igneous rocks', 
#              'Massive sulfide', 'Igneous rocks', 'Igneous rocks', 'Igneous rocks'] 

#     resprops = [0.0, 49.0,69.0,89.0,109.0,129.0,149.0,179.0,249.0,699.0,799.0,899.0]
    
#     thickness, lnames = annotate_tip(layer_thickness=resprops, layer_names=rocks)


    

    

    
    

    
    
    
    