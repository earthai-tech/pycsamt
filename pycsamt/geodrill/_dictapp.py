# -*- coding: utf-8 -*-
#       Copyright © 2021  Kouadio K.Laurent, Licence: LGPL
#       @author: KouaoLaurent alias @Daniel03 <etanoyau@gmail.con>
#       Created on Tue Oct 13 14:52:08 2020
"""
use : 
    for sql dataBase works 
    dico to sql_database . 

"""

# from datetime import datetime, time, date, timezone 

class SqlQ(object):
    """
    Build sql _requests - Not use for others
    you may populate request for your purpose. the request in not general , 
    the user must change the request according its will . 
    
    """
    
    sql_req = ['select Formation, WellName, Depth, GR, ILD_log10, DeltaPHI, PHIND from nofacies_db '\
                   'join blind_stuart_db using(WellName)',
                   
               'select * from nofacies_db join blind_stuart_db using(WellName) order by a_mort',
               
               'select WellName from nofacies_db intersect select WellName from blind_stuart_db',
               
               'select WellName from nofacies_db except select WellName from blind_stuart_db',
               
               'select WellName from blind_stuart_db except select WellName from nofacies_db',
               
               'select distinct WellName from nofacies_db union select WellName from blind_stuart_db ',
               
               'select s2.Formation,s2.WellName, s2.Depth, s2.GR, s2.ILD_log10, s2.DeltaPHI, s2.PHIND,'\
                       ' s2.PE, s2.NM_M, s2.RELPOS, s1.Depth,s1.Facies,s1.LithLabel from nofacies_db as s2,'\
                        ' blind_stuart_db as s1 where s2.WellName=s1.WellName and s1.Depth=s2.Depth',
                       ]
class Glob(object): 
    """
    Spaces of variables and fonctions pseudo-globales .
    dictionnary can be set outside the container class Glob,
    following the dicoT datastructuration . 
    
    :Example: 
        
        >>> value_DB =[('id_new', 'i', 'new_vision'),
        ...   ('infoTab','k','no comment'), 
        ...   ('collar','d','collarDH')]
        >>> Glob.dicoT.__setitem__('TableDB_set',value_DB )
        
    """
    
    
    # struture of database, dictionnary of Table and fields # 
    
    dicoT={ "key_table_one": [("main_f","k",'id_','type_of_value','posi','num_of_value','date_time',"primary key",'comment'),
                                ("main_f","k",'id_','type_of_value','posi','num_of_value','date_time',"primary key",'comment'),
                                ("main_f","k",'id_','type_of_value','posi','num_of_value','date_time',"primary key",'comment'),
                                ("main_f","k",'id_','type_of_value','posi','num_of_value','date_time',"primary key",'comment'),
                                ("main_f","k",'id_','type_of_value','posi','num_of_value','date_time',"primary key",'comment')],
    
           "key_table_two": [("main_f","k",'id_','type_of_value','posi','num_of_value','date_time',"primary key",'comment'),
                                ("main_f","k",'id_','type_of_value','posi','num_of_value','date_time',"primary key",'comment'),
                                ("main_f","k",'id_','type_of_value','posi','num_of_value','date_time',"primary key",'comment'),
                                ("main_f","k",'id_','type_of_value','posi','num_of_value','date_time',"primary key",'comment'),
                                ("main_f","k",'id_','type_of_value','posi','num_of_value','date_time',"primary key",'comment')]}
    
