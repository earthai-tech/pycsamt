# -*- coding: utf-8 -*-
# Copyright (c) 2021 Kouadio K. Laurent, on Mon Oct 25 15:05:40 2021
#       This module is a part of pycsamt utils packages
#       released under a LGL- licence.
#       @author-email:<etanoyau@gmail.com>
import os
import itertools
import warnings
import shutil 
import copy 
from six.moves import urllib 
from pprint import pprint 
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt

from pycsamt.bases import sPath 
import pycsamt.utils.func_utils as FU
import pycsamt.utils.plot_utils as PU
import pycsamt.utils.exceptions as CSex
from pycsamt.utils._csamtpylog import csamtpylog
_logger=csamtpylog.get_csamtpy_logger(__name__)

PATH = 'data/occam2D'
k_ =['model', 'iter', 'mesh', 'data']

try : 
    INVERS_KWS = {
        s +'_fn':os.path.join(PATH, file) 
        for file in os.listdir(PATH) 
                  for s in k_ if file.lower().find(s)>=0
                  }
except :
    INVERS=dict()


TRES=[10, 66,  70, 100, 1000, 3000]# 7000] #[10,  70, 100, 1000,  3000]
#[10, 66, 70, 100, 1000, 2000, 3000, 7000, 15000 ]      
LNS =['river water','fracture zone', 'MWG', 'LWG', 
      'granite', 'igneous rocks', 'basement rocks']


def lns_and_tres_split(ix,  lns, tres):
    """ Indeed lns and tres from `GeoStratigraphy` model are updated. 
    
    Then splitting the `lns` and `tres` from the topped up values is necessary.
    Kind to resetting `tres` and `ln `back to original and return each split
    of inputting layers and TRES and the automatic rocks topped up
    during the NM construction.
    
    :param ix: int 
        Number of autorocks added 
    :param lns: list 
        List of input layers
    :param tres: list 
        List of input resistivities values.
    """
    if ix ==0: return  lns, tres,[], []
    return  lns[:-ix], tres[:-ix], lns[-ix:], tres[-ix:]   
   

def get_closest_gap (value, iter_obj, status ='isin', 
                          condition_status =False, skip_value =0 ):
    """ Get the value from the minimum gap found between iterable values.
    
    :param value: float 
        Value to find its corresponding in the `iter_obj`
    :param iter_obj: iterable obj 
        Object to iterate in oder to find the index and the value that match 
        the best `value`. 
    :param condition_status:bool 
        If there is a condition to skip an existing value in the `iter_obj`, 
        it should be set to ``True`` and mention the `ship_value`. 
    :param skip_value: float or obj 
        Value to skip existing in the `iter_obj`. 
        
    :param status:str 
        If layer is in the databse, then get the electrical property and 
        from that properties find the closest value in TRES 
        If layer not in the database, then loop the database from the TRES 
        and find the auto rock name from resistivity values in the TRES
        
    :return: 
        - ix_close_res: close value with its index found in` iter_obj`
    :rtype:tuple 
    
    """
    minbuff= np.inf 
    ix_close_res =None
    in_database_args = ['isin' , 'in', 'on', 'yes', 'inside']
    out_database_args= ['outoff' , 'out', 'no', 'isoff']
    if status.lower() in in_database_args:
        status ='isin'
    elif status.lower() in out_database_args: 
        status ='isoff'
    else: 
        raise ValueError(f"Given argument `status` ={status!r} is wrong."
                         f" Use arguments {in_database_args} to specify "
                         "whether rock name exists in the database, "
                         f"otherwise use arguments {out_database_args}.")

    for i, v in enumerate(iter_obj):
        if condition_status : 
            if v==skip_value:continue # skip 
        if status=='isin': 
            try: iter(value)
            except :e_min = abs(v - value)
            else : e_min = np.abs(v - np.array(value)).min() 
        # reverse option: loop all the database 
        elif status=='isoff':
            try: iter(v)
            except :e_min = abs(value - v)
            else :e_min = np.abs(value - np.array(v)).min() 
                
        if e_min <= minbuff : 
            if status =='isoff':
                ix_close_res = (i,  value) # index and value in database  
            else:ix_close_res = (i,  v)  # index and value in TRES 
            
            minbuff = e_min 
        
    return ix_close_res

def fit_rocks(logS_array, lns_, tres_):
    """ Find the pseudo rock name at each station from the pseudovalue intres 
    
    :param logS_array: array_like of under the station resistivity value 
    :param lns_: array_like of the rocks or the pseudolayers (automatick)
    :param tres_: array_like of the TRES or the pseudo value of the TRES 
    
    :returns: list of the unik corresponding  resistivity value at each 
            station  and its fitting rock names.
            
    :Example: 
        
        >>> import pycsamt.utils.geo_utils as GU
        >>> obj= GU.quick_read_geos()
        >>> pslns , pstres,  ps_lnstres= make_strata(obj)
        >>> logS1 =obj.nmSites[0] # station S0
        >>> fit_rock(logS1, lns_= pslns, tres_= pstres)
    """
    # get the log of each stations 
    # now find the corresponding layer name from close value in 
    #pseudotres 
    unik_l= np.unique(logS_array)
    unik_fitted_rocks=list()
    for k in range(len(unik_l)): 
        ix_best,_= get_closest_gap ( value = unik_l[k], 
                                       iter_obj =tres_ )
        unik_fitted_rocks.append(lns_[ix_best])
    # now build log blocks 
    fitted_rocks =list()
    for value in logS_array : 
        ix_value, =np.where(unik_l==value) # if index found
        fitted_rocks.append(unik_fitted_rocks[int(ix_value)])
        
    return fitted_rocks 

def assert_station(id, nm =None):
    """ Assert station according to the number of stations investigated
    :param id: int or str, station number. The station counter start from 01 
        as litteral count except whn provided value in string format 
        following the letter `S`. For instance : `S00` =1
    :param nm: matrix of new stratiraphy model built. 
    :return: Index at specific station
    :Example:
        >>> import pycsamt.utils.geo_utils as GU
        >>> geoObj = GU.quick_read_geos() 
        >>> assert_station(id=47, nm=geoObj.nmSites)
        ...46
        
    """
    nstations = nm.shape[1]
    id_= PU.station_id(id)
    
    if id_> nstations : 
        msg ='Site `S{0:02}` is out of the range. Max site= `S{1:02}`.'
        msg=msg.format(id_, nstations-1)
        msg+='The last station `S{0:02}` shoud be used'\
            ' instead.'.format(nstations-1)
        warnings.warn(msg, UserWarning)
        _logger.debug(msg)
        id_= nstations 
        
    return id_
        
def find_distinct_items_and_indexes(items, cumsum =False ):
    """ Find distincts times and their indexes.
    
    :param items: list of items to get the distincts values 
    :param cumsum: bool, cummulative sum when items is a numerical values
    :returns:  
        - distinct _indexes unique indexes of distinct items 
        - distinct_items: unik items in the list 
        - cumsum: cumulative sum of numerical items
        
    :Example: 
        >>> import pycsamt.utils.geo_utils as GU
        >>> test_values = [2,2, 5, 8, 8, 8, 10, 12, 1, 1, 2, 3, 3,4, 4, 6]
        >>> ditems, dindexes, cumsum = GU.find_distinct_items_and_indexes(
            test_values, cumsum =True)
        >>> cumsum 
    """
    if isinstance(items, (tuple, np.ndarray, pd.Series)):
        items =list(items)
     
    if cumsum : 
        try : 
            np.array(items).astype(float)
        except: 
            warnings.warn('Cumulative sum is possible only with numerical '
                          f'values not {np.array(items).dtype} type.')
            cumsum_=None
        else: cumsum_= np.cumsum(items)
        
    else: cumsum_=None 

    s, init = items[0], 0
    distinct_items= list()
    ix_b=[]
    for i, value in enumerate(items):
        if value ==s : 
            if i ==len(items)-1: 
                ix_b.append(init)
                distinct_items.append(s)
            continue 
        elif value != s: 
            distinct_items.append(s)
            ix_b.append(init)
            s= value 
            init= i
           
    return  ix_b, distinct_items, cumsum_
        
def grouped_items( items, dindexes, force =True ):   
    """ Grouped items with the same value from their corresponding
    indexes.
    
    :param items: list of items for grouping.
    :param dindexes: list of distinct indexes 
    :param force: bool, force the last value to broken into two lists.
                Forcing value to be broke is usefull when the items are string.
                Otherwise, `force`  param should be ``False`` when dealing 
                numerical values.
    :return: distinct items grouped 
    
    :Example:  
        >>> import pycsamt.utils.geo_utils as GU
        >>> test_values = [2,2, 5, 8, 8, 8, 10, 12, 1, 1, 2, 3, 3,4, 4, 6]
        >>> dindexes,* _ = GU.find_distinct_items_and_indexes(
            test_values, cumsum =False)
        >>> GU.grouped_items( test_values, dindexes)
        ...  [[2, 2], [5], [8, 8, 8], [10], [12], [1, 1],
        ...      [2], [3, 3], [4, 4], [6]]
        >>> GU.grouped_items( test_values, dindexes, force =False)
        ... [[2, 2], [5], [8, 8, 8], [10], [12], [1, 1],
            [2], [3, 3], [4, 4, 6]]
    """
    gitems =list() 
    
    def split_l(list0): 
       """ split list to two when distinct values is found
       for instance list [3, 3, 4, 4]--> [3, 3], [4,4]"""
       for i, v in enumerate(list0): 
           if i ==0: continue 
           if v != list0[i-1]: 
               return [list0[:i], list0[i:]] 

    for k , val in enumerate(dindexes): 
        if k== len(dindexes)-1:
            # check the last value and compare it to the new 
            gitems.append(items[val:])
            # get the last list and check values 
            l= gitems[-1]
            if force:
                if len(set(l)) !=1: 
                    try:
                        gitems = gitems[:-1] + split_l(l)
                    except: 
                        raise TypeError(
                            'can only concatenate list (not "NoneType") to list.'
                            ' Please check your `items` argument.')
            break 
        gitems.append(items[val:dindexes[k+1]])
    # if there a empty list a the end then remove it     
    if len(gitems[-1]) ==0: 
        gitems = gitems [1:]
        
    return gitems

def fit_stratum_property (fittedrocks , z , site_tres):
    """ Separated whole blocks into different stratum and fit their
    corresponding property like depth and values of resistivities
    
    :param fittedrocks: array_like of layers fitted from the TRES 
    :param z: array like of the depth 
    :param site_tres: array like of the station TRES 
    
    :returns: 
        - s_grouped: Each stratum grouped from s_tres 
        - site_tres_grouped: The site resistivity value `site_tres` grouped 
        - z_grouped: The depth grouped (from the top to bottom )
        - z_cumsum_grouped: The cumulative sum grouped from the top to bottom
        
    :Example: 
        
        >>> import pycsamt.geodrill.geocore as GC
        >>> obj= GC.quick_read_geomodel()
        >>> logS1 = obj.nmSites[:, 0] 
        >>> sg, stg, zg, zcg= fit_stratum_property (
            obj.fitted_rocks, obj.z, obj.logS)
    """
    
    # loop the fitted rocks and find for each stratum its depth and it values
    # find indexes of distinct rocks 
    dindexes,* _ = find_distinct_items_and_indexes(fittedrocks)
    strata_grouped = grouped_items( fittedrocks , dindexes )
    # do it for tres 
    site_tres_grouped = grouped_items(site_tres, dindexes)
    # do it for depth 
    cumsumz = np.cumsum(z)
    z_grouped = grouped_items(z, dindexes, force =False)
    zcumsum_grouped = grouped_items( cumsumz, dindexes, force=False)
    
    return strata_grouped, site_tres_grouped, z_grouped , zcumsum_grouped

def get_s_thicknesses(grouped_z, grouped_s, display_s =True, station=None):
    """ Compute the thickness of each stratum from the grouped strata from 
    the top to the bottom.
    
    :param grouped_z: depth grouped according its TRES 
    :param grouped_s: strata grouped according to its TRES 
    :param s_display: bool, display infos in stdout 
    
    :returns: 
        - thick : The thickness of each layers 
        - strata: name of layers 
        - status: check whether the total thickness is equal to the 
            depth of investigation(doi). Iftrue: `coverall= 100%
            otherwise coverall is less which mean there is a missing layer 
            which was not probably taking account.
            
    :Example: 
        
        >>> import pycsamt.geodrill.geocore as GC
        >>> obj= GC.quick_read_geomodel()
        >>> sg, _, zg, _= fit_stratum_property (obj.fitted_rocks,
        ...                                    obj.z, obj.nmSites[:, 0]  )
        >>> get_s_thicknesses( zg, sg)
        ... ([13.0, 16.0, 260.0, 240.0, 470.0],
        ...     ['*i', 'igneous rocks', 'granite', 'igneous rocks', 'granite'],
        ...     'coverall =100%')
    """
    # get the distincs layers 
    
    pstrata = [stratum[0] if stratum[0] != '$(i)$' else '*i' 
               for stratum in grouped_s ]
    # pstrata =[stratum for s in pstrata else '']
    b= grouped_z[0][0] #take the first values 
    thick =[]
    thick_range =[]
    for k , zg in enumerate(grouped_z): 
        th = abs(max(zg) - b) # for consistency 
        thick.append( th) 
        str_=f"{b} ----- "
        b= max(zg)
        thick_range.append(str_ + f"{b}")
    
    doi = grouped_z[-1][-1]
    if sum(thick) == doi: 
        status = "coverall =100%"
    else : 
        status = f"coverall= {round((sum(thick)/doi)*100, 2)}%"
    
    if display_s: 
        display_s_infos(pstrata , thick_range, thick, 
                        station =station )
    # keep back the stamp 
    pstrata =['$(i)$' if s== '*i' else s for s in  pstrata]
    return thick , pstrata,  status 

def display_s_infos( s_list, s_range, s_thick, **kws):
    """ Display strata infos at the requested station.
    
    :param s_list: the pseudostratigraphic details list 
    :param s_range: the pseudostratigraphic strata range 
    :param s_thick: the pseudostratigraphic  thicknesses
    :param kws:additional keywords arguments.
    """
    linestyle =kws.pop('linestyle', '-')
    mullines = kws.pop('linelength', 102)
    station = kws.pop('station', None)
    
    if station is not None: 
        if isinstance(station, (float,int)): 
            station = 'S{0:02}'.format(station)
        ts_= '{'+':~^'+f'{mullines}'+'}'
        print(ts_.format('[ PseudoStratigraphic Details: '
                         'Station = {0} ]'.format(station)))
        
    print(linestyle *mullines)
    headtemp = "|{0:>10} | {1:^30} | {2:^30} | {3:^20} |"
    if '*i' in  s_list :
        str_i ='(*i=unknow)'
    else: str_i =''
    headinfos= headtemp.format(
        'Rank', f'Stratum {str_i}', 'Thick-range(m)', 'Thickness(m)')
    print(headinfos)
    print(linestyle *mullines)
    for i, (sn, sthr, sth) in enumerate(zip(s_list, s_range, s_thick)):
        print(headtemp.format(f"{i+1}.", sn, sthr, sth))
        
    print(linestyle *mullines)
    
def assert_len_lns_tres(lns, tres): 
    """ Assert the length of LN and Tres"""
    msg= "Input resistivity values <TRES> and input layers <LN> MUST" +\
          " have the same length. But {0} and {1} were given respectively. "
    is_the_same = len(tres) == len(lns)
    return is_the_same , msg.format(len(tres), len(lns))   

 
def _sanitize_db_items (value, force =True ): 
    """ Sanitize Database properties by removing the parenthesis and 
    convert numerical data to float. 
    
    :param value: float of list of values to sanitize.
    :param force:If `force` is `True` will return value whtout parenthesis 
        but not convert the inside values
        
    :return: A list of sanitized items 
    
    :Example:
        
        >>> import pycsamt.utils.geo_utils as GU
        >>> test=['(1.0, 0.5019607843137255, 1.0)','(+o++.)',
        ...          '(0.25, .0, 0.98)', '(0.23, .0, 1.)']
        >>> GU._sanitize_db_items (test)
        ...[(1.0, 0.5019607843137255, 1.0),
        ...    '+o++.', (0.25, 0.0, 0.98), (0.23, 0.0, 1.0)]
        >>> _sanitize_db_items (test, force =False)
        ... [(1.0, 0.5019607843137255, 1.0), 
             '(+o++.)', (0.25, 0.0, 0.98), (0.23, 0.0, 1.0)]
    """

    if isinstance(value, str): 
        value=[value]
    def sf_(v):
        """Sanitise only a single value"""
        if '(' and ')' not in  v:
            try : float(v) 
            except : return v 
            else:return float(v)
        try : 
            v = tuple([float (ss) for ss in 
                 v.replace('(', '').replace(')', '').split(',')])
        except : 
            if force:
                if '(' and ')' in v:
                    v=v.replace('(', '').replace(')', '')
        return v

    return list(map(lambda x:sf_(x), value))


def base_log( ax, thick, layers, *, ylims=None, hatch=None, color=None )  : 
    """ Plot pseudo-stratigraphy basemap and return ax 
    
    :param ax: obj, Matplotlib axis 
    :param thick: list of the thicknesses of the layers 
    :param layers: list of the name of layers
    :param hatch: list of the layer patterns
    :param color: list of the layer colors
    
    :return: ax: matplotlib axis properties 
    """
    if ylims is None: 
        ylims=[0, int(np.cumsum(thick).max())]
    ax.set_ylim(ylims)
    th_data = np.array([np.array([i]) for i in thick ]) 

    for ii, data in enumerate(th_data ): 
        next_bottom = sum(th_data [:ii]) +  ylims[0]
        ax.bar(1,
               data,
               bottom =next_bottom, 
               hatch = hatch[ii], 
               color = color[ii],
               width = .3)

    ax.set_ylabel('Depth(m)', fontsize = 16 , fontweight='bold')
    
    pg = [ylims[0]] + list (np.cumsum(thick) + ylims[0])

    ax.set_yticks(pg)
    ax.set_yticklabels([f'${int(i)}$' for i in  pg] )
    
    ax.tick_params(axis='y', 
                   labelsize= 12., 
                        )
    # inverse axes 
    plt.gca().invert_yaxis()
    return ax 

def annotate_log (ax, thick, layers,*, ylims=None, colors=None, 
                    set_nul='*unknow', bbox_kws=None, 
                    set_nul_bbox_kws=None, **an_kws): 
    """ Draw annotate stratigraphic map. 
    
    :param ax: obj, Matplotlib axis 
    :param thick: list of the thicknesses of the layers 
    :param layers: list of the name of layers
    :param set_nul: str 
        `set the Name of the unknow layers. Default is `*unknow`. Can be 
        changed with any other layer name. 
    :param bbox_kws:dict,  Additional keywords arguments of Fancy boxstyle 
        arguments
    :param set_nul_bbox_kws: dict, customize the boxstyle of the `set_nul` 
        param. If you want the bbox to be the same like `bbox_kws`, we need 
        just to write ``idem`` or `same`` or ``origin``.
        
    :return: ax: matplotlib axis properties 
    """
    
    xinf , xsup =-1, +1
    xpad =  .1* abs(xinf)/2
    ax.set_xlim([xinf, xsup])
    if ylims is None: 
        ax.set_ylim([0, int(np.cumsum(thick).max())]) #### 1 check this part 
        ylim0=0.
    else :
        ax.set_ylim(ylims)
        ylim0=ylims[0]
    # inverse axes 
    plt.gca().invert_yaxis()
    # if ylims is None:
    pg = np.cumsum([ylim0] + thick)# add 0. to thick to set the origin  #### 2
    # take values except the last y from 0 to 799 
    v_arrow_bases = [(xinf + xpad, y) for y in  pg ]
    v_xy = v_arrow_bases[:-1]
    v_xytext = v_arrow_bases[1:]
    # build the pseudo _thickness distance between axes 
    for k, (x, y) in enumerate(v_xy):
        ax.annotate('', xy=(x, y), xytext =v_xytext[k],
                    xycoords ='data',
                    # textcoords ='offset points',
                    arrowprops=dict(arrowstyle = "<|-|>", 
                                  ),  
                horizontalalignment='center',
                verticalalignment='top',                         
                )
    # ------------make horizontal arraow_[properties]
    # build the mid point where starting annotations 
    mid_loc = np.array(thick)/2 
    # if ylims is None: 
    center_positions =  pg[:-1] + mid_loc
    # else :center_positions =  pg + mid_loc
    h_xy = [ (xinf + xpad, cp) for cp in center_positions]
    h_xytext = [(0, mk ) for mk in center_positions ]
    
    # control the color 
    if colors is not None: 
        if isinstance(colors, (tuple, list, np.ndarray)):
            if len(colors) != len(thick): colors =None 
        else : colors =None 
    # build the pseudo _thickness distance between axes 
    if bbox_kws is None:  
         bbox0=  dict(boxstyle ="round", fc ="0.8", ec='k')
    if set_nul_bbox_kws in ['idem', 'same', 'origin']: 
        bbox_i = bbox0
        
    if not isinstance (set_nul_bbox_kws, dict): 
         set_nul_bbox =  None     
    if set_nul_bbox is None:
        bbox_i= dict (boxstyle='round', 
                     fc=(.9, 0, .8), ec=(1, 0.5, 1, 0.5))

    layers=[f"${set_nul}$" 
            if s.find("(i)")>=0 else s for s in layers ]
    
    for k, (x, y) in enumerate(h_xy):
        if layers[k] ==f"${set_nul}$" : 
            bbox = bbox_i
        else: bbox = bbox0
        if colors is not None:
            bbox ['fc']= colors[k]
        ax.annotate( f"{layers[k]}",  
                    xy= (x, y) ,
                    xytext = h_xytext[k],
                    xycoords='data', 
                    arrowprops= dict(arrowstyle='-|>', lw = 2.), 
                    va='center',
                    ha='center',
                    bbox = bbox, **an_kws
                 )

    return ax 


def pseudostratigraphic_log (thick, layers, station, *,
                    zoom =None, hatch=None, color=None, **annot_kws) : 
    """ Make the pseudostratigraphic log with annotate figure.
    
    :param thick: list of the thicknesses of the layers 
    :param layers: list of the name of layers
    :param hatch: list of the layer patterns
    :param color: list of the layer colors
    :parm zoom: float, list. If float value is given, it considered as a 
            zoom ratio and it should be ranged between 0 and 1. 
            For isntance: 
                - 0.25 --> 25% plot start from 0. to max depth * 0.25 m.
                
            Otherwise if values given are in the list, they should be
            composed of two items which are the `top` and `bottom` of
            the plot.  For instance: 
                - [10, 120] --> top =10m and bottom = 120 m.
                
            Note that if the length of `zoom` list is greater than 2, 
             the function will return all the plot and 
             no errors should raised. 

    :Example: 
        >>> from pycsamt.utils.geo_utils as GU   
        >>> layers= ['$(i)$', 'granite', '$(i)$', 'granite']
        >>> thicknesses= [59.0, 150.0, 590.0, 200.0]
        >>> hatch =['//.', '.--', '+++.', 'oo+.']
        >>> color =[(0.5019607843137255, 0.0, 1.0), 'b', (0.8, 0.6, 1.), 'lime']
        >>> GU.pseudostratigraphic_log (thicknesses, layers, hatch =hatch ,
        ...                   color =color, station='S00')
        >>>  GU.pseudostratigraphic_log ( thicknesses,
                                         layers,
                                         hatch =hatch, zoom =0.25,
                                         color =color, station='S00')
    """
    import matplotlib.gridspec as GridSpec
    
    is_the_same, typea_status, typeb_status= _assert_list_len_and_item_type(
        thick, layers,typea =(int, float, np.ndarray),typeb =str)
    if not is_the_same: 
        raise TypeError("Layers' thicknesses and layer names lists must have "
                        "the same. But {len(thick)} and {len(layers)} were given.")
    if not typea_status: # try to convert to float
        try : 
            thick =[float(f) for f in thick ]
        except :raise ValueError("Could not convert to float."
                                 " Please check your thickness values.")
    if not  typeb_status: layers =[str(s) for s in layers] 
    
    # get the dfault layers properties hatch and colors 
    # change the `none` values if exists to the default values
    #for hatch and colors
    # print(color)
    hatch , color = set_default_hatch_color_values(hatch, color)
    #####INSERT ZOOM TIP HERE#######
    ylims =None
    if zoom is not None: 
        ylims, thick, layers, hatch, color = zoom_processing(zoom=zoom, 
             data= thick, layers =layers, hatches =hatch,colors =color)
    #####################################
    fig = plt.figure( f"PseudoLog of Station ={station}",
                     figsize = (10,14), 
                      # dpi =300 
                     )
    plt.clf()
    gs = GridSpec.GridSpec(1,2,
                       wspace=0.05,
                       left=.08,
                       top=.85,
                       bottom=0.1,
                       right=.98,
                       hspace=.0,
                       ) 
    doi = sum(thick) 
    axis_base = fig.add_subplot(gs[0, 0],
                           ylim = [0, int(doi)] ####### check this part 
                           )
                
    axis_annot= fig.add_subplot(gs[0, 1],
                                sharey=axis_base) 
    axis_base = base_log(ax = axis_base, 
                         thick=thick, 
                         ylims=ylims, 
                         layers=layers,
                         hatch=hatch,
                         color=color)
    
    axis_annot = annotate_log(ax= axis_annot,
                              ylims=ylims, 
                             thick=thick,
                             layers=layers,
                             colors =color,
                             **annot_kws)
    
    for axis in (axis_base, axis_annot):
        for ax_ in ['top','bottom','left','right']:
            if ax_ =='left': 
                if ylims is None:
                    axis.spines[ax_].set_bounds(0, doi)
                else :  axis.spines[ax_].set_bounds(ylims[0], ylims[1])
                axis.spines[ax_].set_linewidth(3)
            else : axis.spines[ax_ ].set_visible(False)
  
        axis.set_xticks([]) 
        
    fig.suptitle( f"PseudoStratigraphic log of Station ={station}",
                ha='center',
        fontsize= 7* 2., 
        verticalalignment='center', 
        style ='italic',
        bbox =dict(boxstyle='round',facecolor ='moccasin'), 
        y=0.90)
    
    plt.show()

def _assert_list_len_and_item_type(lista, listb, typea=None, typeb=None):
    """ Assert two lists and items type composed the lists 
    
    :param lista: List A to check the length and the items type
    :param listb: list B to check the length and the items type
    :param typea: The type which all items in `lista` might be
    :param typeb: The type which all items in `listb` might be
    
    :returns: 
        - the status of the length of the two list ``True`` or ``False``
        - the status of the type of `lista` ``True`` if all items are the 
            same type otherwise ``False``
        - idem of `listb`
        
    :Example: 
        >>> thicknesses= [59.0, 150.0, 590.0, 200.0]
        >>> hatch =['//.', '.--', '+++.', 'oo+.']
        >>> _assert_list_len_and_item_type(thicknesses, hatch,
        ...                                   typea =(int, float, np.ndarray),
        ...                                    typeb =str))
        ... (True, True, True)
    """
    try: import __builtin__ as b
    except ImportError: import builtins as b
    
    def control_global_type(typ):
        """ Check the given type """
        import pandas as pd 
        builtin_types= [t for t in b.__dict__.values()
                     if isinstance(t, type)] 
        conv_type = builtin_types+ [np.ndarray, pd.Series,pd.DataFrame]
        if not isinstance( typ, (tuple, list)):
            typ =[typ]
        # Now loop the type and check wether one given type is true
        for ityp in typ:
            if ityp not in conv_type: 
                raise TypeError(f"The given type= {ityp} is unacceptable!"
                                " Need a least the following types "
                                f" {FU.smart_format([str(i) for i in conv_type])}"
                                " for checking.")
        return True
    
    is_the_same_length  = len(lista) ==len (listb)
    lista_items_are_the_same, listb_items_are_the_same =False, False
    def check_items_type(list0, type0): 
        """ Verify whether all items  in the list are the same type"""
        all_items_type = False
        is_true = control_global_type(type0)
        if is_true:
             s0= [True if isinstance(i0, type0) else False for i0 in list0 ]
             if len(set(s0)) ==1: 
                 all_items_type = s0[0]
        return  all_items_type
        
    if typea is not None :
        lista_items_are_the_same = check_items_type (lista, typea)
    if typeb is not None :
        listb_items_are_the_same = check_items_type (listb, typeb)
        
    return (is_the_same_length , lista_items_are_the_same,
            listb_items_are_the_same )
    

def set_default_hatch_color_values(hatch, color, dhatch='.--', 
                                   dcolor=(0.5019607843137255, 0.0, 1.0),
                                   force =False): 
    """ Set the none hatch or color to the default value. 
    :param hatch: str or list of layer patterns 
    :param color: str or list of layers colors
    :param dhatch: default hatch 
    :param dcolor: default color 
    :param force: Return only single tuple values otherwise put the RGB tuple
        values  in the list. For instance::
            -if ``False`` color =[(0.5019607843137255, 0.0, 1.0)]
            - if ``True`` color = (0.5019607843137255, 0.0, 1.0)
    :Example: 
        >>> from pycsamt.utils.geo_utils as  GU.
        >>> hatch =['//.', 'none', '+++.', None]
        >>> color =[(0.5019607843137255, 0.0, 1.0), None, (0.8, 0.6, 1.),'lime']
        >>> set_default_hatch_color_values(hatch, color))
    """
    fs=0 # flag to reconvert the single RGB color in tuple 
    def set_up_(hc, dhc):
        if isinstance(hc, (list, tuple, set, np.ndarray)): 
            hc =[dhc if h in (None, 'none') else h for h in hc ]
        return hc 
    
    if isinstance(hatch, str): hatch=[hatch]
    if isinstance(color, str): color=[color]
    elif len(color)==3 : 
        try:iter(color[0]) # check iterable value in tuple
        except : 
            try : color =[float(c) for c in color]
            except : raise ValueError("wrong color values.")
            else: fs=1
    hatch = set_up_(hatch, dhatch)
    color = set_up_(color, dcolor) 
    if force:
        if len(color) ==1: color = color[0]
    if len(hatch) ==1 : hatch = hatch[0]
    if fs==1: color = tuple(color)
    return hatch, color 

def print_running_line_prop(obj, inversion_software='Occam2D') :
    """ print the file  in stdout which is currently used
    " for pseudostratigraphic  plot when extracting station for the plot. """
    
    print('{:~^108}'.format(
        f' Survey Line: {inversion_software} files properties '))
    print('|' + ''.join([ "{0:<5} = {1:<17}|".format(
        i, os.path.basename( str(getattr(obj, f'{i}_fn'))))  
     for i in ['model', 'iter', 'mesh', 'data']]) )
    print('~'*108)


#*** Manage the geological rocks from files:AGSO & AGSO.STCODES ******
__agso_properties =dict(
    GIT_REPO = 'https://github.com/WEgeophysics/pyCSAMT', 
    GIT_ROOT ='https://raw.githubusercontent.com/WEgeophysics/pyCSAMT/master/',
    props_dir = 'pycsamt/geodrill/_geocodes',
    props_files = ['AGSO.csv', 'AGSO_STCODES.csv'], 
    props_codes = ['code', 'label', 'name', 'pattern','size',
            'density', 'thickness', 'color']
    )

def set_agso_properties (download_files = True ): 
    """ Set the rocks and their properties from inner files located in 
        <  'pycsamt/geodrill/_geocodes'> folder."""
        
    msg= ''.join([
        "Please don't move or delete the properties files located in", 
        f" <`{__agso_properties['props_dir']}`> directory."])
    mf =list()
    __agso= [ os.path.join(os.path.realpath(__agso_properties['props_dir']),
                           f) for f in __agso_properties['props_files']]
    for f in __agso: 
        agso_exists = os.path.isfile(f)
        if not agso_exists: 
            mf.append(f)
            continue 
    
    if len(mf)==0: download_files=False 
    if download_files: 
        for file_r in mf:
            success = fetching_data_from_pycsamt_repo(props_files = file_r, 
                      savepath = os.path.join(
                          os.path.realpath('.'), __agso_properties['props_dir'])
                      )
            if not success:
                msg_ = ''.join([ "Unable to retreive the geostructure ",
                      f"{os.path.basename(file_r)!r} property file from",
                      f" {__agso_properties['GIT_REPO']!r}."])
                warnings.warn(f"Geological structure file {file_r} "
                              f"is missing. {msg_}") 
                _logger.warn( msg_)
                raise CSex.pyCSAMTError_file_handling(
                    f"No property file {os.path.basename(file_r)!r}"
                    f" is found. {msg}.")
    for f in __agso:
        with open(f,'r' , encoding ='utf8') as fs: 
            yield([stratum.strip('\n').split(',')
                    for stratum in fs.readlines()])
            
def mapping_stratum(download_files =True): 
    """ Map the rocks properties  from  _geocodes files and fit 
    each rock to its properties. 
    
    :param download_files: bool 
        Fetching data from repository if the geostrutures files are missing.
    :return: Rocks and structures data  in two diferent dictionnaries
    """
    # get code description _index 
    ix_= __agso_properties['props_codes'].index('name')
    def mfunc_(d): 
        """ Set individual layer in dict of properties """
        _p= {c: k.lower() if c not in ('code', 'label', 'name') else k 
                 for c,  k in zip(__agso_properties['props_codes'], d) }
        id_= d[ix_].replace('/', '_').replace(
            ' ', '_').replace('"', '').replace("'", '').lower()
        return id_, _p 
    rock_and_structural_props =list()
    for agso_data in tuple(set_agso_properties(download_files)): 
        # remove the header of the property file
        rock_and_structural_props.append(
            dict(map( lambda x: mfunc_(x), agso_data[1:])))
     
    return   tuple(rock_and_structural_props)


def fetching_data_from_pycsamt_repo(repo_file, savepath =None ): 
    """ Try to retreive data from github repository.
    
    :param repo_file: str or Path-like object 
        Give the full path from the repository root to the file name.
        For instance, we want to retrieve the file 'AGSO.csv' which is located 
        in <pycsamt/geodrill/_geocodes> directory then the full path 
        is: --> 'pycsamt/geodrill/_geocodes/AGSO.csv'
        
    :return:`status`: Either ``False` for failed downloading 
            or ``True`` for successfully downloading
    """
    fmsg =['... Please wait for the second attempt...',
          '... Wait for the last attempt...']
    status=False 
    git_repo = __agso_properties['GIT_REPO']
    git_root = __agso_properties['GIT_ROOT']
    # max attempts =3 : 
    for i in range(3):
        if i ==0 :
            print("---> Please wait while fetching"
                  f" {repo_file!r} from {git_repo!r}...")
        else:print(fmsg [i-1])
        try : 
            try :
                urllib.request.urlretrieve(git_root,  repo_file )
            except: 
                with urllib.request.urlopen(git_root) as response:
                    with open( repo_file,'wb') as out_file:
                        data = response.read() # a `bytes` object
                        out_file.write(data)
        except TimeoutError: 
            if i ==2: print("---> Established connection failed "
                      " because connected host has failed to respond.")
        except:pass 
        else : 
            print(f"---> Downloading {repo_file!r} from "
                  f"{git_repo!r} was successfully done!")
            status=True
            break 
    if status: print(f"---> Downloading {repo_file!r} from {git_repo!r} "
                 "was successfully done!")
    else: print(f"---> Failed to download {repo_file!r} from {git_repo!r}!")
    # now move the file to the right place and create path if dir not exists
    if savepath is not None: 
        if not os.path.isdir(savepath): 
            sPath (savepath)
        shutil.move(os.path.realpath(repo_file), savepath )
    if not status:pprint(connect_reason )
    
    return status

def get_agso_properties(config_file =None, orient ='series'): 
    """ Get the geostructures files from <'pycsamt/geodrill/_geocodes'> and 
    set the properties according to the desire type. When `orient` is 
    ``series`` it will return a dictionnary with key equal to 
    properties name and values are the properties items.
    
    :param config_file: Path_Like or str 
        Can be any property file provided hat its obey the disposal of 
        property files found in   `__agso_properties`.
    :param orient: string value, ('dict', 'list', 'series', 'split',
        'records’, ''index') Defines which dtype to convert
        Columns(series into).For example, 'list' would return a 
        dictionary of lists with Key=Column name and Value=List 
        (Converted series). For furthers details, please refer to
        https://www.geeksforgeeks.org/python-pandas-dataframe-to_dict/
        
    :Example: 
        >>> import pycsamt.utils.geo_utils as GU
        >>> data=GU.('pycsamt/geodrill/_geocodes/AGSO_STCODES.csv')
        >>> code_descr={key:value for key , value in zip (data["CODE"],
                                                       data['__DESCRIPTION'])}
    """
    msg= ''.join(["<`{0}`> is the software property file. Please don't move "
        " or delete the properties files located in <`{1}`> directory."])
    
    pd_pos_read ={".csv":pd.read_csv, ".xlsx":pd.read_excel,
             ".html":pd.read_json,".sql" : pd.read_sql,".json":pd.read_json
             } 
    ext='none'
    if config_file is None: 
        config_file = os.path.join(os.path.realpath('.'), os.path.join(
                       __agso_properties['props_dir'],
                       __agso_properties ['props_files'][0]))
    if config_file is not None: 
        is_config = os.path.isfile(config_file)
        if not is_config : 
            if os.path.basename(config_file) in __agso_properties['props_files']:
                _logger.error(f"Unable to find  the geostructure property" 
                              f"{os.path.basename(config_file)!r} file."
                              )
                warnings.warn(msg.format(os.path.basename(config_file) , 
                                         __agso_properties['props_dir']))
            raise FileExistsError(f"File `{config_file}`does not exist")
            
        _, ext = os.path.splitext(config_file)
        if ext not in pd_pos_read.keys():
            _logger.error(f"Unable to read {config_file!r}. Acceptable formats"
                          f" are {FU.smart_format(list(pd_pos_read.keys()))}.")
            raise TypeError(
                f"Format {ext!r} cannot be read. Can only read "
                    "{FU.smart_format(list(pd_pos_read.keys()))} files."
                )
    agso_rock_props = pd_pos_read[ext](config_file).to_dict(orient)
    if ('name' or 'NAME') in agso_rock_props.keys(): 
        agso_rock_props['__DESCRIPTION'] = agso_rock_props ['name']
        del agso_rock_props['name']
        
    return  agso_rock_props

#***  end  AGSO & AGSO.STCODES management******

def map_bottom (bottom, data, origin=None): 
    """Reduce the plot map from the top assumes to start at 0. to the
    bottom value.
    
    :param bottom: float, is the bottom value from
        which the plot must be end 
    :param data: the list of layers thicknesses in meters
    :param origin: The top value for plotting.by the default 
        the `origin` takes 0. as the starting point
        
    :return: the mapping depth from the top to the maximum depth.
            - the index indicated the endpoint of number of layer 
                for visualizing.
            - the list of pairs <top-bottom>, ex: [0, bottom]>
            - the value of thicknesses ranged from  0. to the bottom 
            - the coverall, which is the cumul sum of the value of
                the thicknesses reduced compared to the total depth.
     Note that to avoid raising errors, if coverall not equal to 100%,
     will return safety default values (original values).
     
    :Example: 
        >>> ex= [ 59.0, 150.0, 590.0, 200.0]
        >>> map_bottom(60, ex)
        ... ((2, [0.0, 60], [59.0, 1.0]), 'coverall = 100.0 %')
    """
    
    cumsum_origin = list(itertools.accumulate(data)) 
    if origin is None: origin = 0.
    end = max(cumsum_origin)
    wf =False # warning flag
    coverall, index =0., 0
    wmsg = ''.join([ "Bottom value ={0} m might be less than or ",
                          "equal to the maximum depth ={1} m."])
    t_to_b = list(itertools.takewhile(lambda x: x<= bottom,
                                      cumsum_origin))
    
    if bottom > end :bottom , wf = end, True 
    elif bottom ==0 or bottom < 0: 
        bottom , wf = data[0], True 
        to_bottom=([origin , bottom], [bottom])
    elif bottom < data[0] : 
        to_bottom = ([origin ,bottom], [bottom]) 
    elif len(t_to_b) !=0 :
        # add remain extent bottom values
        if max(t_to_b) < bottom : 
            add_bottom = [abs (bottom - max(t_to_b))] 
            to_bottom = ([origin, bottom], data[:len(t_to_b)] + add_bottom )
        elif max(t_to_b) ==bottom :
            to_bottom= ([origin, sum(t_to_b)],  t_to_b)
        index =len(to_bottom[1])   # get the endpoint of view layers 
    if bottom ==end : # force to take the bottom value
        to_bottom= ([origin, bottom], data)
        index = len(data)
        
    if wf:
        warnings.warn(wmsg.format(bottom, sum(data)), UserWarning)
        wf =False # shut down the flag
    coverall=  round(sum(to_bottom[1])/ bottom, 2)
    cov = f"coverall = {coverall *100} %"
    if coverall !=1.: 
        to_bottom = (len(data), [0., sum(data)], data)
    else : to_bottom = get_index_for_mapping(index, to_bottom )
    
    return to_bottom, cov 
    
def get_index_for_mapping (ix, tp): 
    """ get the index and set the stratpoint of the top or the endpoint 
    of bottom from tuple list. The index should be used to map the o
    ther properties like `color` or `hatches`"""
    tp=list(tp)
    # insert index from which to reduce other property
    tp.insert(0, ix)
    return  tuple(tp )
    
    
def map_top (top, data, end=None): 
    """ Reduce the plot map from the top value to the bottom assumed to 
    be the maximum of investigation depth 
    
    :param top: float, is the top value from which the plot must be starts 
    :param data: the list of layers thicknesses in meters
    :param end: The maximum depth for plotting. Might be reduced, 
        but the default takes the maximum depth investigation depth 
    
    :return: 
         the mapping depth from the top to the maximum depth.
        - the index indicated the number of layer for visualizing to 
                from the top to the max depth.
        - the list of pairs <top-bottom>, ex: [top, max depth]>
        - the value of thicknesses ranged from 0. to the bottom 
        - the coverall, which is the cumul sum of the value of
            the thicknesses reduced compared to the total depth.
            It allows to track a bug issue.
            
        Note that if coverall is different 100%, will return the 
        default values giving values. 
        
    :Example: 
        >>> ex= [ 59.0, 150.0, 590.0, 200.0] # layers thicknesses 
        >>> map_top(60, ex)
        ... ((3, [60, 999.0], [149.0, 590.0, 200.0]), 'coverall = 100.0 %')
    """
    wmsg = ''.join([ "Top value ={0} m might be less than ",
                    "the bottom value ={1} m."])
    cumsum_origin = list(itertools.accumulate(data)) 
    if end is None: end = max(cumsum_origin)
    # filter list and keep value in cumsum 
    #greater  or equal to top values 
    data_= copy.deepcopy(data)
    if top < 0: top =0.
    elif top > end : 
        warnings.warn(wmsg.format(top, sum(data)), UserWarning)
        top = sum(data[:-1])
    t_to_b = list(filter(lambda x: x > top, cumsum_origin)) 
    index =0
    coverall =0.
    if len(t_to_b) !=0:
        if min (t_to_b)> top : # top = 60  --> [209.0, 799.0, 999.0]
            #get index of the min value from the cums_origin 
            # find 209 in [59.0, 209.0, 799.0, 999.0] --> index = 1
            index= cumsum_origin.index(min(t_to_b))
            #  find the  value at index =1 in data 
            #[ 59.0, 150.0, 590.0, 200.0]=> 150
             # reduce the downtop abs(59 - 60) = 1
            r_= abs(sum(data[:index]) - top )
            # reduced the data at index  1 with r_= 1
            reduce_top = abs(data [index] - r_)  # data[1]=150-1: 149 m 
            # replace the top value `150` in data with the reduced value
            data[index] = reduce_top  # 150==149
            from_top= ([top, end],data [index:] )# [ 149, 590.0, 200.0]
        elif min(t_to_b) == top: 
            index = cumsum_origin.index(min(t_to_b))
            from_top = ([top, end], data[index:])
            r_ = abs(sum(data[:index]) - top )
        
        coverall = round((sum(data[index :] +data[:index ])
                          + r_)/ sum(data_),2)
        
    cov = f"coverall = {coverall *100} %"
    if coverall !=1.:
        from_top = (index, [0., sum(data_)], data_)
    else:from_top= get_index_for_mapping(index, from_top )
        
    return from_top, cov 

def smart_zoom(v):
    """ select the ratio or value for zooming. Don't raise any error, just 
    return the original size. No shrunk need to be apply when error occurs.

    :param v: str, float or iterable for zoom
            - str: 0.25% ---> mean 25% view 
                    1/4 ---> means 25% view 
            - iterable: [0, 120]--> the top starts a 0.m  and bottom at 120.m 
            note: terable should contains only the top value and the bottom 
                value.
    :return: ratio float value of iteration list value including the 
        the value range (top and the bottom values).
    :Example: 
        >>> import pycsamt.utils.geo_utils as GU
        >>> GU.smart_zoom ('1/4')
        ... 0.25
        >>> GU.smart_zoom ([60, 20])
        ... [20, 60]
    """
    def str_c (s):
        try:
            s=float(s)
        except:
            if '/' in s: # get the ratio to zoom 
                s= list(map(float, sorted(s.split('/'))))
                s=round(min(s)/max(s),2)
            elif '%' in s: # remove % and get the ratio for zoom
                s= float(s.replace('%', ''))*1e-2
                if s >1.:s=1.
            else: s=1.
            if s ==0.:s=1.
        return s
    
    msg="Ratio value `{}` might be greater than 0 and less than 1."
    emsg = f" Wong given value. Could not convert {v} to float "
    is_iterable =False 
    try:iter(v)
    except :pass 
    else : 
        if isinstance(v, str): v= str_c(v)
        else: is_iterable = True
    if not is_iterable: 
        try:  float(v)
        except ValueError:
            s=False 
            try:v = str_c(v)
            except ValueError :warnings.warn(emsg)
            else :s=True # conversion to zoom ratio was accepted
            if not s:
                warnings.warn(emsg)
                v=1.
        else: 
            if v > 1. or v ==0.:
                warnings.warn(msg.format(v)) 
                v=1.
    elif is_iterable : 
        if len(v)>2:
            warnings.warn(f" Expect to get size =2 <top and bottom values>. "
                          f"Size of `{v}` should not  be greater than 2,"
                          f" but {len(v)} were given.", UserWarning)
            v=1.
        try : v= list(map(float, sorted(v)))
        except:  
            warnings.warn(emsg)
            v=1.
    return v

def frame_top_to_bottom (top, bottom, data ): 
    """ Compute the value range between the top and the bottom.
    
    Find the main range value for plot ranged between the top model and 
        the bottom. 
    :param top: is the top value where plot might be started 
    :param bottom: is the bottom value where plot must end 
    :param data: is the list of thicknesses of each layers 
    :param cumsum_origin: is the list of cumul sum of the `data` 
    
    :returns: 
        - the index for other properties mapping. It indicates the
            index of layer for the top and the index of layer for the 
            bottom for visualizing.
        - the list of pairs top-bottom , ex: [10, 999.0] 
                where tp -> 10 and bottom ->999. m
        - the value of thicknesses ranged from the top to the bottom 
            For instance:  [49.0, 150.0, 590.0, 200.0] where 
            - 49 is the thockness of the first layer 
            - 200 m is the thickness of the 
        - the coverall allow to track bug issues.The thickness of layer 
            for visualizing should ne the same that shrank. If not the same 
            the mapping was not successfully done. Therefore coverall 
            will be different to 100% and function will return the raw data
            instead of raising errors. 
            
    :Example: 
        >>> import pycsamt.utils.geo_utils as GU
        >>> layer_thicknesses = [ 59.0, 150.0, 590.0, 200.0]
        >>> top , bottom = 10, 120 # in meters 
        >>> GU.frame_top_to_bottom( top = top, bottom =bottom,
                                data = layer_thicknesses)
        ...(([0, 2], [10, 120], [49.0, 61.0]), 'coverall = 100.0 %')
 
    """
    if top > bottom :
        warnings.warn( f"Top value ={top} should be less than"
                      " the bottom ={bottom} ")
        top=0.
    if top ==bottom :top , bottom = 0.,  sum(data) 
    if top <0 : top =0.
    if bottom > sum(data): 
        warnings.warn( f"Bottom value {bottom} should be less than"
                      f" or equal to {sum(data)}")
        bottom =sum(data)
        
    # test data = [ 59.0, 150.0, 590.0, 200.0]
    data_safe = copy.deepcopy(data)
    data_ = copy.deepcopy(data)
    # get the value from the to to the bottom 
    tm,*_ = map_top (top, data = data )
    ixt, _, tm = tm # [49.0, 150.0, 590.0, 200.0]
    bm, *_= map_bottom(bottom, data = data_ )
    ixb, _, bm = bm  # [59.0, 150.0, 391.0])
    #remove the startpoint and the endpoint from top and bottom 
    sp = [tm[0]]  # 149.
    del tm[0]       #--> [150.0, 590.0, 200.0]
    ep = [bm[-1]]         # 391.
    del bm[-1]      # --> [59.0, 150.0]
    # computethe intersection of the two list 
    inter_set_map_tb = set(tm).intersection(set(bm))
    # set obj classification is sometimes messy, so let loop 
    # to keep the layer disposalthe same like the safe data value
    # if len(inter_set_map_tb )!=0: 
    inter_set_map_tb=[v for v in data_safe if v in inter_set_map_tb]
    top_bottom = sp + inter_set_map_tb + ep 
    # compute coverall to track bug issues 
    coverall = round(sum(top_bottom)/ abs(top-bottom ), 2)
    cov = f"coverall = {coverall *100} %"
    top_bottom = ([ixt, ixb], [top, bottom], top_bottom)
    if coverall !=1.:
        top_bottom = ([0, len(data_safe) ],
                      [0., sum(data_safe )], data_safe )
 
    return top_bottom, cov

def zoom_processing(zoom, data, layers =None, 
                    hatches=None, colors =None): 
    """ Zoom the necessary part of the plot. 
    
    If some optionals data are given such as `hatches`, `colors`, `layers`,
    they must be the same size like the data.
    
    :param zoom: float, list. If float value is given, it's cnsidered as a 
            zoom ratio than it should be ranged between 0 and 1. 
            For isntance: 
                - 0.25 --> 25% plot start from 0. to max depth * 0.25 m.
                
            Otherwise if values given are in the list, they should be
            composed of two items which are the `top` and `bottom` of
            the plot.  For instance: 
                - [10, 120] --> top =10m and bottom = 120 m.
                
            Note that if the length of list  is greater than 2, the function 
            will return all the plot and  no errors should raised.
    :param data: list of composed data. It should be the thickness from 
        the top to the bottom of the plot.
        
    :param layers: optional, list of layers that fits the `data`
    :param hatches: optional, list of hatches that correspond to the `data` 
    :param colors: optional, list of colors that fits the `data`
    
    :returns: 
        - top-botom pairs: list composed of top bottom values
        - new_thicknesses: new layers thicknesses computed from top to bottom
        - other optional arguments shrunk to match the number of layers and 
            the name of exact layers at the depth.
    
    :Example: 
        >>> import pycsamt.utils.geo_utils as GU
        >>> layers= ['$(i)$', 'granite', '$(i)$', 'granite']
        >>> thicknesses= [59.0, 150.0, 590.0, 200.0]
        >>> hatch =['//.', 'none', '+++.', None]
        >>> color =[(0.5019607843137255, 0.0, 1.0), None, (0.8, 0.6, 1.),'lime'] 
        >>> zoom_processing(zoom=0.5 , data= thicknesses, layers =layers, 
                              hatches =hatch, colors =color) 
        ... ([0.0, 499.5],
        ...     [59.0, 150.0, 290.5],
        ...     ['$(i)$', 'granite', '$(i)$'],
        ...     ['//.', 'none', '+++.'], 
        ...     [(0.5019607843137255, 0.0, 1.0), None, (0.8, 0.6, 1.0)])
        
    """
    def control_iterable_obj( l, size= len(data)): 
        """ Control obj size compared with the data size."""
        if len(l) != size : 
            warnings.mess(f"The iterable object {l}  and data "
                        f" must be same size ={size}. But {len(l)}"
                        f" {'were' if len(l)>1 else 'was'} given.")
            l=None 
        return l 
 
    zoom = smart_zoom(zoom) # return ratio or iterable obj 
    y_low, y_up = 0., sum(data)
    ix_inf=0
    try : 
        iter(zoom)
    except :
        if zoom ==1.:# straightforwardly return the raw values (zoom =1)
            return [y_low, y_up], data, layers,  hatches, colors
 
    if isinstance(zoom, (int, float)): #ratio value ex:zoom= 0.25
        # by default get the bootom value and start the top at 0.
        y_up = zoom * sum(data)  # to get the bottom value  as the max depth 
        bm, *_=map_bottom(y_up, data = data )
        ix_sup, _, maptopbottom = bm  # [59.0, 150.0, 391.0])
       
        y = [y_low, y_up ]

    if isinstance(zoom , (list, tuple, np.ndarray)): 
        top, bottom = zoom 
        tb, *_= frame_top_to_bottom (top, bottom, data )
        ixtb, y, maptopbottom= tb 
        ix_inf, ix_sup = ixtb  # get indexes range
        
    if hatches is not None:
        hatches = control_iterable_obj(hatches)
        hatches =hatches [ix_inf:ix_sup]
    if colors is not None:
        colors= control_iterable_obj(colors)
        colors = colors[ix_inf:ix_sup]
    if layers is not None:
        layers= control_iterable_obj(layers)
        layers = layers [ix_inf:ix_sup]
        
    return y, maptopbottom, layers, hatches , colors 


##############connection git error ##########################
connect_reason ="""<ConnectionRefusedError><No connection could  '
            be made because the target machine actively refused it>.
            There are some possible reasons for that:
         1. The server is not running. Hence it wont listen to that port. 
             If it's a service you may want to restart the service.
         2. The server is running but that port is blocked by Windows Firewall
             or other firewall. You can enable the program to go through 
             firewall in the Inbound list.
        3. There is a security program on your PC, i.e a Internet Security 
            or Antivirus that blocks several ports on your PC.
        """   
