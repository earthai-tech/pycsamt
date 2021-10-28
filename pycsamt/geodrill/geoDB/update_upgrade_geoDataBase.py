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

    .Script to UPDATE and UPGRADE GeoDataBase. GeoDataBase is the Core of geodrill module.
     More than 150 rocks as well as their properties have already set. User can
     add more informations to  increase the chance to plot a truth  pseudostratigraphy log
     with their right pattern , colors and resistivity values.  In addition , provide 
     the right resistivity value of  differents geological rocks,  may increase the accuracy
     to find the closest layer when only input resistivity values  are set (_see_module `Geodrill`)
    . pyCSAMT geoDataBase is flexible , can be edited at least you  are aware of what you are doing.
   
     GeoDatabse deals with all geological informations and  The electrical properties of rocks
     their pattern and the Digital Cartography map symbolization. 
     The database code colums 
          ::: ('code','label','__description','pattern', 'pat_size',	'pat_density',
            'pat_thickness','rgb', 'electrical_props', 'hatch', 'colorMPL', 'FGDC') ::: 
              
           ! `__description` can be replaced  by `name` !
     --Have fun ---
 
Created on Tue Feb 16 11:37:47 2021

@author: @Daniel03
"""

from pycsamt.geodrill.geoDB.sql_recorder import GeoDataBase 


# code and Label are automatically , dont need to add info . if provided , it will be rejected 
GET_INFO =False                                 # see the information already set in DataBase 
                                                # before nany update and upgrade , it is better to set this to true to see 
                                                # if the structure exist and what are the corresponding values in database .
UPDATE =  False                                 # if update provide the geostructure name and the columns to modify 
UPGRADE =False                                  # if upgrade provide all  informations 

name_geo =  'granite'                           # name of stucture 
pattern_geo = 258.                              # In fact pattern are linked with hatch . the size of elements 
pat_size_geo = 1.5                              # size of hatch symbol (float) 
pat_density_geo =.23                            # density of hatch symbol (float)
pat_thickness = 2.                              # thick ness of hatch symbol (float)
rgb_geo ='R128G52B26'                           # RGB colors : eg RG28 or hexadecima #000ff  or 'blue'. program will convert to "RGB" except alpha
hatch_geo ='+.//.Oo'                            # matplotlib pattern build with symbolisation map ["/", "\\", "|", '-', '+', 'x', 'o', 'O', '.', '*'] 
electrical_props = [2e3 ,  15e3]                 # range of electrical property of rocks (maximum and minium)
colorMPL = (1.0, 0.25, 0.)                      # if rgb is provided , this key will set automatically , dont need to add 
FGDC ='none'                                    # Dont need to fill this part , Digital cartographic 
                                                # Standard for Geological  Map Symbolisation, will be fill later  at least you 
                                                # know what you do.
    

    
# fill the update dictionnary  and set UPDATE to TRUE  # update informationalready exist in database 
update_geo_dict = {'name': name_geo, 
                   'electrical_props':electrical_props}
# fill the upgrade dictionnary  and set UPGRADE to TRUE  # set new information to dataBase 
upgrade_geo_dict = {
                    'name': name_geo, 
                    'pattern': pattern_geo, 
                    'pat_size':pat_size_geo, 
                    'pat_density':pat_density_geo, 
                    'pat_thickness': pat_thickness,
                    'electrical_props':electrical_props , 
                    'hatch': hatch_geo , 
                    'colorMPL': colorMPL, 
                    }

# GET INFO ALREADY EXISTS 
if GET_INFO is True : 
    GeoDataBase()._reminder_geo_recorder(name_geo) # set the name you want to check  eg : GeoDataBase()._reminder_geo_recorder('amphibolite') 
#update INFO ALREADY EXISTS
if UPDATE is True : 
    rgeoDB= GeoDataBase()._update_geo_structure(**update_geo_dict)

#UPGRADE NEW STRUCTURE
if UPGRADE is True : 
    geodatabase_obj= GeoDataBase._add_geo_structure( **upgrade_geo_dict )
        

                
