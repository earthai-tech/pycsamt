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

Created on Tue Oct 13 14:52:08 2020

@author: Daniel03

use : 
    for sql dataBase works 
    dico to sql_database . 

"""


# from datetime import datetime, time, date, timezone 



class Glob(object): 
    """
    Spaces of variables and fonctions pseudo-globales .
    dictionnary can be set outside the container class Glob,
    following the dicoT datastructuration . 
    
     eg : 
         value_DB =[('id_new', 'i', 'new_vision'),
            ('infoTab','k','no comment'), 
            ('collar','d','collarDH')]
         Glob.dicoT.__setitem__('TableDB_set',value_DB )

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
    
