# -*- coding: utf-8 -*-
"""
    .script to rewrite jfile. 
    
Created on Tue Jan  5 12:11:07 2021

@author: @Daniel03

"""
import os 
from pycsamt.ff.core.j import J_collection as J 

#--- > set jfiles path 
path =r'C:\Users\Administrator\Desktop\pyCS_datasets\testj'

#---> add your path to save files , if None , savepath is your current work directory 
savepath = None
#---> add survey_name on the file 
surveyname ='Niable'
#---> add your extension file to get the format of file you want . Defaut is '.dat' 
jextension ='.dat'

#---> set to True to get info about the script .Defalut is False 

see_documentation =False 
#--> set to False to only read info about rewrite obj.
REWRITE= True  

#--> create J_object 
j_obj =J()
if see_documentation : 
    help(j_obj.rewrite_jfiles())

if REWRITE : 
    #---> call list of files <----
    jfiles_list =[os.path.join(path, jfile) for jfile in os.listdir(path)] 
    #--> rewrite j
    j_obj.rewrite_jfiles(list_of_jfiles=jfiles_list,
                         savepath =savepath,
                         survey_name =surveyname, 
                         j_extension=jextension)



