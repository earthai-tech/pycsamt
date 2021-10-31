# -*- coding: utf-8 -*-
# Copyright (c) 2021 Kouadio K. Laurent, on Mon Oct 25 15:05:40 2021
#       This module is a part of pycsamt utils packages
#       released under a LGL- licence.
#       @author: K.KL alias Daniel03<etanoyau@gamil.com>
import os
import warnings 
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt

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
        
        >>> import pycsamt.utils.geo_utils as GU
        >>> obj= GU.quick_read_geos()
        >>> logS1 = obj.nmSites[:, 0] 
        >>> obj = make_strata(obj)
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
        
        >>> import pycsamt.utils.geo_utils as GU
        >>> obj= GU.quick_read_geos()
        >>> obj= make_strata(obj)
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
        """Santise only a single value"""
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


def base_log( ax, thick, layers, hatch=None, color=None )  : 
    """ Plot pseudo-stratigraphy basemap and return ax 
    
    :param ax: obj, Matplotlib axis 
    :param thick: list of the thicknesses of the layers 
    :param layers: list of the name of layers
    :param hatch: list of the layer patterns
    :param color: list of the layer colors
    
    :return: ax: matplotlib axis properties 
    """

    th_data = [np.array([i]) for i in thick ]
    # doi = sum(thick)
    for ii, data in enumerate(th_data ): 
        next_bottom = sum(th_data [:ii])
        ax.bar(1,
               data,
               bottom =next_bottom, 
               hatch = hatch[ii], 
               color = color[ii],
               width = .3)
    # inverse axes 
    plt.gca().invert_yaxis()

    ax.set_ylabel('Depth(m)', fontsize = 16 , fontweight='bold')

    pg = [0.]+ list(np.cumsum([thick]))
    ax.set_yticks(pg)
    ax.set_yticklabels([f'${int(i)}$' for i in  pg] )
    ax.tick_params(axis='y', 
                   labelsize= 12., 
                        )

    return ax 

def annotate_log (ax, thick, layers, colors=None, set_nul='*unknow',
                    bbox_kws=None, set_nul_bbox_kws=None, **an_kws): 
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
    ax.set_ylim([0, int(np.cumsum(thick).max())])
    # inverse axes 
    plt.gca().invert_yaxis()
    # add 0. to thick to set the origin  
    pg = np.cumsum([0] + thick)
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
    # build the mid point where staring annotations 
    mid_loc = np.array(thick)/2 
    center_positions =  pg[:-1] + mid_loc
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
                    hatch=None, color=None, **annot_kws) : 
    """ Make the pseudostratigraphic log with annotate figure.
    
    :param thick: list of the thicknesses of the layers 
    :param layers: list of the name of layers
    :param hatch: list of the layer patterns
    :param color: list of the layer colors
    
    :Example: 
        >>> from pycsamt.utils.geo_utils as GU   
        >>> layers= ['$(i)$', 'granite', '$(i)$', 'granite']
        >>> thicknesses= [59.0, 150.0, 590.0, 200.0]
        >>> hatch =['//.', '.--', '+++.', 'oo+.']
        >>> color =[(0.5019607843137255, 0.0, 1.0), 'b', (0.8, 0.6, 1.), 'lime']
        >>> GU.pseudostratigraphic_log (thicknesses, layers, hatch =hatch ,
        ...                   color =color, station='S00')
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
                           ylim = [0, int(doi)]
                           )
                
    axis_annot= fig.add_subplot(gs[0, 1],
                                sharey=axis_base) 
    axis_base = base_log(ax = axis_base, 
                         thick=thick, 
                         layers=layers,
                         hatch=hatch,
                         color=color)
    
    axis_annot = annotate_log(ax= axis_annot,
                             thick=thick,
                             layers=layers,
                             colors =color,
                             **annot_kws)
    
    for axis in (axis_base, axis_annot):
        for ax_ in ['top','bottom','left','right']:
            if ax_ =='left': 
                axis.spines[ax_].set_bounds(0, doi)
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
        >>> color =[(0.5019607843137255, 0.0, 1.0), None, (0.8, 0.6, 1.),  'lime']
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


def __build_ps__token(obj):
    """ Build a special token for each GeoStratigraphic model. Please don't 
    edit anything here. Force editing is your own risk."""
    import random 
    random.seed(42)
    __c =''.join([ i for i in  [''.join([str(c) for c in obj.crmSites.shape]), 
     ''.join([str(n) for n in obj.nmSites.shape]),
    ''.join([l for l in obj.input_layers]) + str(len(obj.input_layers)), 
    str(len(obj.tres))] + [''.join(
        [str(i) for i in [obj._eta, obj.beta, obj.doi,obj.n_epochs,
                          obj.ptol, str(obj.z.max())]])]])
    __c = ''.join(random.sample(__c, len(__c))).replace(' ', '')                                               
    n= ''.join([str(getattr(obj, f'{l}'+'_fn'))
                         for l in ['model', 'iter', 'mesh', 'data']])
    n = ''.join([s.lower() 
                 for s in random.sample(n, len(n))]
                ).replace('/', '').replace('\\', '')
    
    return ''.join([n, __c]).replace('.', '')


def print_running_line_prop(obj, inversion_software='Occam2D') :
    """ print the file  in stdout which is currently used
    " for pseudostratigraphic  plot when extracting station for the plot. """
    
    print('{:~^108}'.format(
        f' Survey Line: {inversion_software} files properties '))
    print('|' + ''.join([ "{0:<5} = {1:<17}|".format(
        i, os.path.basename( str(getattr(obj, f'{i}_fn'))))  
     for i in ['model', 'iter', 'mesh', 'data']]) )
    print('~'*108)



#*** file from _geocodes folder :AGSO & AGSO.STCODES ******
def _code_strata():
    return ['code', 'label', 'name', 'pattern',
                  'size', 'density', 'thickness', 'color']

def agso_data (): 
    """
    Geological data codes processing
    
    .. deprecated:: function  will later deprecated 
    
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
            raise CSex.pyCSAMTError_file_handling(
                'Geocodes files :< AGS0 & AGS0_STCODES> not found.'
                 'check the right path to _geocodes folder.')
    return fagso, fagso_stcodes


def set_stratum_on_dict():
    """
    Process to put geocodes_strata and geocodes_structures into dictionnaries  
    better way to go on metaclasses merely. Thus each keys of dictionary will be 
    its own object. 

    Returns
    -------
        strata_dict : dict
            Disctionnary of geostrata 
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
            if stratum_description.startswith('"') or\
                stratum_description.endswith('"'):
                stratum_description=stratum_description[1:-1]
            if iiagso ==1 :
                structures_dict['{0}'.format(stratum_description)]=attr_dict
            else :
                strata_dict['{0}'.format(stratum_description)]=attr_dict
            attr_dict={}

    return strata_dict, structures_dict

#***  end call file from _geocodes folder :AGSO & AGSO.STCODES ******

