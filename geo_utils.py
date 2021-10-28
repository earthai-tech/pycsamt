# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 15:05:40 2021

@author: @Daniel03
"""
import os
import warnings 
import numpy as np
from pycsamt.geodrill.geoCore.geodrill import GeoStratigraphy

path=r'F:\ThesisImp\occam2D\invers+files\inver_res\K1'
inversion_files = {'model_fn':'Occam2DModel', 
                    'mesh_fn': 'Occam2DMesh',
                    "iter_fn":'ITER17.iter',
                    'data_fn':'OccamDataFile.dat'
                    }
input_resistivity_values =[10, 66, 70, 100, 1000, 2000, 
                                3000, 7000, 15000 ] 
input_layer_names =['river water', 'fracture zone', 'granite']
inversion_files = {key:os.path.join(path, vv) for key,
                    vv in inversion_files.items()}
  
#     with np.errstate(divide='ignore'):
#         ss= np.array(inpt2) /np.array(input_resistivity_values )
# #         print(ss )
geosObj = GeoStratigraphy(**inversion_files,
                  input_resistivities=input_resistivity_values, 
                  input_layers=input_layer_names)

geosObj._createNM()
# crm = geosObj.crmSites
# nrm = geosObj.nmSites

lns = geosObj.input_layers 
tres= geosObj.tres
autoln= geosObj.auto_layer_names



def reset_lns_and_tres(ix,  lns, tres):
    """ Resetting tres and ln back to original and return the layer names exempt 
    of automatic rocks topped up during the NM construction.
    :param ix: int 
        Number of autorocks added 
    :param lns: list 
        List of input layers
    :param tres: list 
        List of input resistivities values.
    """
    if ix ==0: return tres, lns, [], []
    return  lns[:-ix], tres[:-ix], lns[-ix:], tres[-ix:]   


def make_strata(svector, ln, tres, ptol=0.1): 
    """ Build the strata log to fit the resistivity values of the TRES"""
    # get the uniques values in 
    unik_value = np.unique(svector)
    newLN =['' for i in range(len(unik_value))]
    e_min = ptol+1
    for ii, itres in enumerate(unik_value):
        for jj, ires in enumerate(tres): 
            e_min  = abs(itres -ires)
            if e_min <= ptol: 
                newLN[ii] = ln[jj] 

    return unik_value,  newLN

    
def fit(lns, tres, autorocks,  **kws): 
    """ read and get the resistivity values from tres that match the 
     the given layers if layer already exists in the geoDataBase."""
     
    ix = len(autorocks)
    lns0, tres0, rlns, rtres= reset_lns_and_tres(ix,  lns, tres)
    import copy
    r0 =copy.deepcopy(tres0)
    # for consistency, lowercase the layer name
    # get the properties [name and electrical properties]  from geoDataBase 
    # try to build new list with none values 
    # loop for all layer and find their index then their elctrical values 
    #           if name exist in database then:
    #           loop DB layers names 
    #           if layer is found then get it index 

    lns =[ln.lower() for ln in lns ]
    _gammaRES, _gammaLN = geosObj._getProperties(**kws)
    newTRES =[None for i in tres0]
    
    for ii, name in enumerate(lns0) : 
        if name in _gammaLN: 
            ix = _gammaLN.index (name) 
            newTRES[ii] =(name,_gammaRES[ix]) 
            
    # try to set thevalues of res of layer found in the database 
    for ii, nvalue in enumerate(newTRES):
        
        try: 
            iter(nvalue[1])
        except:
            # now if value is in TRES is not set = 0.
            #consider that the rocks does not exist and set to None
            if nvalue is not None and nvalue[1]==0. :
                newTRES[ii]= None 
            continue 
        else: 
            # if iterable get the index and value of layer names 
            # remove this values in the tres 
            ix, val = get_value_from_mingap (value=nvalue[1], iter_obj=tres0)
            newTRES[ii][1] = val
            tres0.pop(ix) 
    
    for ii, nvalue in enumerate(tres0):
        ix,_val=  get_value_from_mingap (value=nvalue,status ='isoff', 
                                   iter_obj=_gammaRES, 
                          condition_status =True, skip_value =0 )
        
        # get the index of this values in tres
        index = r0.index (_val) 
        newTRES[index] = (_gammaLN[ix], nvalue)
    
    # create for each tres its pseudorocks names and pseudorock values
    pseudo_lns = [a [0] for a in newTRES] + rlns 
    pseudo_tres = [b[1] for b in newTRES] + rtres 
    newTRES += [(a, b) for a , b in zip(rlns, rtres)]
    
    return pseudo_lns , pseudo_tres , newTRES 




def get_value_from_mingap (value, iter_obj, status ='isin', 
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
        - ix_close_res: close value with index found in` iter_obj`
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


# print(tres)
# ln0, tres0, noneLN, noneTres = assert_length_LNTRES(len(autoln), ln, tres)  

# print(STRUCT)
# print(noneTres)

# stratares,  stratanames= make_strata(nrm [:, 0], ln0, tres0)
# def fit(nm, ln, tres): 
#     """ Fit for each LN it corresponding values in Tres"""
#     ...
new_tres = fit_tres(lns, tres, autorocks= autoln)
print(new_tres)
# lns= ['river water', 'fracture zone', 'granite', 'massive sulphide', 'sedimentary rocks']
# tres= [1.0, 1.8195439355418688, 1.845098040014257, 2.0, 3.0, 3.3010299956639813, 3.4771212547196626, 3.845098040014257, 4.176091259055681, -0.616462374, 0.5451962386363636]
# print(fit_tres(lns, tres))

# # geosObj.stratigraphyModel(kind ='nm',misfit_G=False)

# lns0, tres0, _= reset_lns_and_tres(2,  lns,tres)
# print(get_value_from_mingap (value=(30.2, 631.0), iter_obj=tres0, 
#                               condition_status =False, skip_value =0 ))
# _gammaRES, _gammaLN = geosObj._getProperties()
# print(find_auto_roc_res(itres= 4.176091259055681 , e_properties= _gammaRES , n_properties= _gammaLN))

# print(get_value_from_mingap (value=4.176091259055681,status ='isoff', iter_obj=_gammaRES, 
#                               condition_status =True, skip_value =0 ))

# print(find_closest_res(dbres =(30.2, 631.0), tres = tres0))



# print(lns0, tres0)
# STRUCT= findGeostructures(tres0)
# print(STRUCT)
# print(fit_tres(lns0, tres0))

# print('ln=', ln)
# print('tres=', tres)
# print(geosObj._getProperties())
# print('ln0=', ln0)
# print('tres0=', tres0)
# print('strataRES=', stratares)
# print('strataNames=',stratanames)