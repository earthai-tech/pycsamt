# -*- coding: utf-8 -*-
"""
    .script to rewrite jfile. 
    
Created on Tue Jan  5 12:11:07 2021

@author:K.L ~ @Daniel03
"""
import os 
from pycsamt.ff.core.j import J_collection as J 

#--- > set jfiles path 
path ='data/j'

#---> add your path to save files
savepath = None 
#---> add survey_name on the file 
surveyname ='Niable'
#---> add your extension file to get the format of file you want. Defaut is '.dat' 
jextension ='.dat'

#--> create J_object 
j_obj =J()
#---> call list of files <----
jfiles_list =[os.path.join(path, jfile) for jfile in os.listdir(path)] 
#--> rewrite j
j_obj.rewrite_jfiles(list_of_jfiles=jfiles_list,
                     savepath =savepath,
                     survey_name =surveyname, 
                     j_extension=jextension)



