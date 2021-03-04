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
.. _module-Request::`geodrill.geoDB.sqlrequests`
    :synopsis: Deal with Boreholes and well data , requests to DataBase
    
Created on Tue Oct 13 17:36:35 2020

@author: @Daniel03

"""

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
