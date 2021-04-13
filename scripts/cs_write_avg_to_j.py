# -*- coding: utf-8 -*-
"""
     .Script to write straightforwardly your avg file to A.G.Jones J-format. Both type of 
     AVG file (plainty file or AsTatic file can be converted into Jformat. )
     
Created on Wed Jan  6 13:49:02 2021

@author: @Daniel03
"""

import os 

from pycsamt.ff.core.avg import Avg 

#---> path to your avg file 
path = r'C:\Users\Administrator\Desktop\test\avg_toj'
#---> savepath 
savepath = r'C:\Users\Administrator\Desktop\test\avg_toj\jfilenew'

#--> Zonge avgfile name 
avgfile ='LCS01.AVG'
#--> station profile name 
station_file = 'K1.stn'
#--> j_extension _name 
j_extension ='.dat'

#---> survey_name 
surveyName =None
#---> add your avg info in your jformat output file 
write_avg_file_infos = False

#---> see documentation : set to True 
see_doc =False
#---> set to True and run cell to write files. 
WRITE =True 

#set UTM zone 
utmZone ='49N'
#----> set path to files 
path_to_avg=os.path.join(path, avgfile)
path_to_station_file = os.path.join(path, station_file)
#--->
avg_obj =Avg()
if WRITE : 
    avg_obj.avg_to_jfile(avg_data_fn=path_to_avg, 
                                station_fn=path_to_station_file,
                                j_extension=j_extension,
                                savepath=savepath  ,
                                writeInfos=write_avg_file_infos, 
                                utm_zone =utmZone , 
                                survey_name =surveyName)
if see_doc :help(avg_obj.avg_to_jfile)