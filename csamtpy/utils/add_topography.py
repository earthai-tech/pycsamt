# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 14:36:23 2020

@author: @Daniel03

sets up input files for running 2d occam inversions adding elevation profile 

"""
import os  
import numpy as np 
import glob 


class Topo: 
    """
    Build numpy array  to add elevation profile and elevation model. 
    
    Attributes :
    ===========
        * add_elevation_profile : np:ndarray (2, num_elev_points)
            x=relative coordinates (offest ) , y= elevation   
        * add_elevation_model : np.nadarry (3,num_elev_points) 
            x=easting (UTM), northing(UTM) , elevation 
        * path_elev : is the path where the profiles files is located .
            ex: os.chdir(os.path.dirname())
                os.path.join(os.path.dirname(), K1_exp.bln)
    
        * step_betw_station : is offset between two stations .
        * self.prof_loc_east : Easting value (longitude(1,))
        * self.prof_loc_north : Nothing value (latitude (1,))
        * self.elev_value     : Elevation  value (Elevation(1,))
    
    Methods : 
    =========
    add_topo : return num_elev_points,  np.ndarray(2, num_elev_points) 
                row1 : relative coordinates 
                row2 : elevation points 
    add_topo_model : return elevation_model , np.ndarray(3, elevation_model)
                row1 : easting coordinates 
                row2 : northing coordinates 
                row3 : elevation points 
    """
    encodage='utf-8'
    
    def __init__(self, path_elev="K1_exp.bln", step_betw_station=50, ): 
        self.path_elev=path_elev 
        self.step_sta=step_betw_station
        
        self.prof_loc_east=None
        self.prof_loc_north=None
        self.elev_value=None
        
    
    def add_topo(self) : 
        """
        fonction  to add elevation from topography 
        x is the realtive offset 0  to n*step_betw_station 
        """
    
        file_elev=os.path.basename(self.path_elev)
        
        try: 
            
            ofi=open(file_elev,'r', encoding=self.encodage)
            
        except: 
            print("*** No such file in the directory***")
            return 
        exten_file=file_elev[-3:].upper()
        if file_elev[-3:] != exten_file: 
            print("** {} is a wrong format , should be *.bln or .BLN **".format(file_elev))

        head=ofi.readline() # read of headfile , not necessary.
        repere="east_DMS north_DMS UTM_zone   east_utm"\
            "      north_utm   azim_m  elev_m "
        if repere not in head:
            print(" * Please select the appropriate *.bln-file *")
            return
        
        prof_loc_list_east, prof_loc_list_north, elev_value_list=[],[],[]
        
        while 1 : 
            ligne=ofi.readline()
            ligne=ligne[:-1] # or ligne = ligne.strip()
            if ligne=='' or ligne == '\n': 
                break
            ligne=ligne[:-1] 
            list_line=ligne.split()
            prof_loc_list_east.append(float(list_line[6]))
            prof_loc_list_north.append(float(list_line[7]))
            elev_value_list.append(float(list_line[9]))
        ofi.close()
        
        # building numpy array (2, num_elev_points)
        rela_coord_list=[ii for ii in range(0,self.step_sta*len(elev_value_list),self.step_sta)] 
        
        prof_loc_east=np.asarray(prof_loc_list_east)
        prof_loc_north=np.asarray(prof_loc_list_north)
        elev_value=np.asarray(elev_value_list)
        rela_coord=np.asarray(rela_coord_list)
        
        rela_coord=rela_coord.reshape((1,rela_coord.shape[0]))
        
        self.prof_loc_east= prof_loc_east.reshape((1, prof_loc_east.shape[0]))

        self.prof_loc_north= prof_loc_north.reshape((1, prof_loc_north.shape[0]))

        self.elev_value=elev_value.reshape((1,elev_value.shape[0]))
        
        num_elev_points=np.concatenate((rela_coord,self.elev_value), axis=0)
        
        min_,max_=min(elev_value_list), max(elev_value_list)
        indmin, indmax=elev_value_list.index(min_),elev_value_list.index(max_)
        
        #Affichage 
        print("="*55)
        print("Line{0:^55}".format(file_elev[:-8]).upper())
        print("="*55)
        print("Nombre de stations: {:>5}".format(len(prof_loc_list_east)))
        print("Elevation mininale: station {0:>4} ={1:>6} {2:>2}".format("S"+str(indmin),min_,"m." ))
        print("Elevation maximale: station {0:>4} ={1:>6} {2:>2}".format("S"+str(indmin) ,max_,'m.'))
        print("{:^55}".format('------'))
        return num_elev_points
        
    def add_topo_model(self): 
        """
        generate a elevation model ndarray(3, num_elev_points)
        easting , northing , and elevation """
        
        list_models=[self.prof_loc_east, self.prof_loc_north, self.elev_value]
        elevationModel=list_models[0]
        for array in list_models[1:]: 
            elevationModel=np.concatenate((elevationModel,array), axis=0)
            
        return elevationModel
    

# TEST : 

if __name__=='__main__': 
    # path= r"C:/Users/kouao/OneDrive/paper2/process_paper2/"+\
    #     "data_stn/stnToDat_exportfiles/K1_cor.DAT"
    # dir_= r"C:/Users/kouao/OneDrive/paper2/process_paper2/"+\
    #     "data_stn/stnToDat_exportfiles"
    dir_=r'E:/OneDrive/Python/CodesExercices/'+\
        'ex_avgfiles/data/process_paper2/data_edifiles/K1_edi/K1_exp.BLN' # change the edipath
    
    os.chdir(r'E:/OneDrive/Python/CodesExercices/'+\
    'ex_avgfiles/data/process_paper2/data_edifiles/K1_edi')    
    
    top=Topo(path_elev=dir_, 
             step_betw_station=50)
    num_=top.add_topo()
    mod_=top.add_topo_model()

