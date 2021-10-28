# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 15:05:40 2021

@author: @Daniel03
"""
import os
import warnings 
import copy
import numpy as np
import pandas as pd 

import pycsamt.utils.func_utils as FU
import pycsamt.utils.plot_utils as PU
from pycsamt.geodrill.geoCore.geodrill import GeoStratigraphy
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

def quick_read_geos(lns=LNS, tres=TRES):
    """Quick read and build the geostratigraphy model (NM) 
    
    :param lns: list of input layers 
    :param tres: list of input true resistivity values 
    
    :Example: 
        >>> import pycsamt.utils.geo_utils as GU 
        >>> obj= GU.quick_read_geos()
        >>> GU.fit_tres(obj.input_layers, obj.tres, obj.auto_layer_names)
    """
    if len(INVERS_KWS) ==0: 
        _logger.error("NoneType can not be read! Need the basics Occam2D"
                         f" inversion {FU.smart_format(k_)} files.")

        raise ValueError("NoneType can not be read! Need the basics Occam2D"
                         f" inversion {FU.smart_format(k_)} files.")
        
    geosObj = GeoStratigraphy( input_resistivities=tres, 
                      input_layers=lns,**INVERS_KWS)
    geosObj._createNM()
    
    return geosObj 

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
   
def fit_tres(lns, tres, autorocks, force=False, **kws): 
    """ Read and get the resistivity values from tres that match the 
     the given layers.
     
    Find the layers and  their corresponding resistivity values from the 
    database especially when values in the TRES and LN are not the same
    length. It's therefore not possible to match each value to its
    correspinding layer name. Therefore the best approach is to read the
    TRES and find the layer name in the database based on the closest value.

    :param lns: list of input layers 
    :param tres: list of input true resistivity values 
    :param autorocks: list of the autorocks found when building the new model.
    :param force: bool, force fitting resistivity value with the rocks in 
            the database whenever the size of rocks match perfectly 
            the number of the rocks. Don't do that if your are sure that the 
            TRES provided fit the  layers in LNS.
    :param kws: is database column property. Default is
        `['electrical_props', '__description']`
        
    :returns: new pseudolist contains the values of rocks retrived from 
        database as well as it closest value in TRES.
    """

    def flip_back_to_tuple(value , substitute_value, index=1): 
        """convert to tuple to list before assign values and 
          reconvert to tuple  after assignment for consistency. 
          `flip_back_to_tuple` in line in this code is the the same like :
                newTRES[ii] = list(newTRES[ii])
                newTRES[ii][1] = val
                newTRES[ii] = tuple (newTRES[ii]) 
          """ 
        value = list(value)
        if index is not None: 
            value[index] = substitute_value
        else : value = substitute_value
        return tuple (value) 

     
    ix = len(autorocks)
    lns0, tres0, rlns, rtres= lns_and_tres_split(ix,  lns, tres)
    if len(lns0) > len(tres0): 
        msg= ''.join(['Number of given layers `{0}` should not be greater ',
                      ' than the number of given resistivity values` {1}`.'])
        msg= msg.format(len(lns0), len(tres0))
        
        n_rock2drop = len(tres0)-len(lns0) 
        msg += f" Layer{'s' if abs(n_rock2drop)>1 else ''} "\
            f"{FU.smart_format(lns0[n_rock2drop:])} should be ignored."
    
        lns0 = lns0[: n_rock2drop]
        warnings.warn(msg)
        _logger.debug(msg)
       
    if sorted([n.lower() for n in lns0]
              ) == sorted([n.lower() for n in lns]): 
        if not force: 
            return lns0, tres0, [(a, b) for a , b in zip(lns0, tres0)]
        
    r0 =copy.deepcopy(tres0)
    # for consistency, lowercase the layer name
    # get the properties [name and electrical properties]  
    # from geoDataBase try to build new list with none values 
    # loop for all layer and find their index then 
    # their elctrical values 
    #           if name exist in database then:
    #           loop DB layers names 
    #           if layer is found then get it index 
    lns0 =[ln.lower().replace('_', ' ') for ln in lns0 ]
    _gammaRES, _gammaLN = GeoStratigraphy._getProperties(**kws)

    newTRES =[None for i in tres0]
    temp=list()
    for ii, name in enumerate(lns0) : 
        if name in _gammaLN: 
            ix = _gammaLN.index (name) 
            temp.append((name,_gammaRES[ix])) 
            
    # keep the lns0 rocks that exists in the database 
    # and replace the database value by the one given 
    #in tres0 and remove the tres value with 
    # unknow layer by its corresponding value.
    if len(temp)!=0: 
        for name, value in temp: 
            ix, val = get_closest_gap (value= value, iter_obj=tres0)
            newTRES[ix]= (name, val) 
            tres0.pop(ix) 
    # try to set the values of res of layer found in 
    # the database is not set = 0 by their corresponding
    # auto -layers. if value is in TRES. We consider 
    #that the rocks does not exist and set to None
    for ii, nvalue in enumerate(newTRES):
        try: 
            iter(nvalue[1])
        except:
            if nvalue is not None and nvalue[1]==0. :
                newTRES[ii]= None 
            continue 
        else: 
            # if iterable get the index and value of layers
            # remove this values in the tres 
            ix, val = get_closest_gap (value=nvalue[1], iter_obj=tres0)
            newTRES[ii] = flip_back_to_tuple (newTRES[ii], val, 1) 
            tres0.pop(ix) 
            
    for ii, nvalue in enumerate(tres0):
        ix,_val=  get_closest_gap (value=nvalue,status ='isoff', 
                                   iter_obj=_gammaRES, 
                          condition_status =True, skip_value =0 )
        # get the index of this values in tres
        index = r0.index (_val) 
        newTRES[index] = (_gammaLN[ix], nvalue)
        
    # create for each tres its pseudorock name 
    # and pseudorock value
    pseudo_lns = [a [0] for a in newTRES] + rlns 
    pseudo_tres = [b[1] for b in newTRES] + rtres 
    newTRES += [(a, b) for a , b in zip(rlns, rtres)]
    
    return pseudo_lns , pseudo_tres , newTRES 

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


def fit_default_layer_properties(layers, dbproperties_= ['hatch', 'colorMPL']): 
    """ Get the default layers properties  implemented in database. 
     
    For instance get the hatches and colors from given layers implemented in 
    the database by given the database `dbproperties_`.
    
    :param layers: str or list of layers to retrieve it properties
        If specific property is missing , ``'none'`` will be return 
    :param db_properties_: str, list or database properties 
    :return: property items sanitized
    
    :Example: 
        
        >>> import pycsamt.utils.geo_utils as GU
        >>> GU.fit_default_layer_properties(
        ...    ['tuff', 'granite', 'evaporite', 'saprock']))
        ... (['none', 'none', 'none', 'none'],
        ...     [(1.0, 1.0, 0.0), (1.0, 0.0, 1.0), (0.0, 0.0, 1.0),
        ...     (1.0, 0.807843137254902, 1.0)])
    """
    # for consistency check again and keep the DB properties untouchable.
    dbproperties_= ['colorMPL' if g.lower().find('mpl') >=0 else 
                    'FGDC' if g.lower()=='fgdc'else g.lower() 
                    for g in dbproperties_]
    if isinstance(layers , str): layers =[layers]
    assert_gl = ['yes' if isinstance(ll, str) else 'no' for ll in layers]
    if not len(set(assert_gl)) ==1: 
        raise TypeError("Wrong given layers. Names should be a string!")
    if 'name' or'__description' not in  dbproperties_: 
        dbproperties_.insert(0, 'name')
    
    __gammaProps = GeoStratigraphy._getProperties(dbproperties_)
    
    r_props =[['none' for i in layers] for j in range(len(__gammaProps)-1)]
    for k  , l in enumerate(layers): 
        if l  in __gammaProps[0] :
            ix= __gammaProps[0].index(l)
            for kk, gg in enumerate(r_props) : 
                gg[k]= __gammaProps[1:][kk][ix]
                
    r_props = [_sanitize_db_items(r_props[k], force=True )
               for k in range(len (r_props))]
    return tuple(r_props)

 
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



