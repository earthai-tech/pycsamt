# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 15:49:10 2020

@author: Administrator

"""
import os 
import glob
import shutil
import csamtpy.pyCS.deprecatedmodules.pyconvert as mconv

#=============================


class Gui_Ediwrite:
    """Create Graphical User Interface from Tkinter """
    
    def __init__(self,edipath, 
                 bln_coordinates_filepath, coordinateFile,savepath):

        self.edipath=edipath
        self.blnpath=bln_coordinates_filepath
        self.coordinateFile=coordinateFile
        self.savepath=savepath
        
        if self.edipath is None :
            raise FileExistsError
        
        if os.path.exists(self.blnpath):
            self.path_to_blnfile=os.path.join(self.blnpath,self.coordinateFile)
            
        if self.savepath is None :
            self.savepath=os.getcwd()
            
        self.rewrite_edifiles()
            
    def rewrite_edifiles(self):
        "methode to rewrite edifiles"
        os.chdir(self.edipath)
        edir=mconv.EdiRewrite(Path_to_edifiles=self.edipath,
            coordinates_bln_file=self.path_to_blnfile)
        edir.check_edi()
        edir.ex_DMS()
        edir.ediRewrite()
        self.info=edir.print_stats()
        # movefile 
        try : 
            os.mkdir('new_EDI')    
        except OSError as e :
            print(os.strerror(e.errno))
        self.savepath=self.edipath+'/new_EDI'
        
        slst=[file for file in os.listdir(self.edipath) if (file[0]=='S'and file[-4:]==".edi")]
        for file in slst: 
            shutil.move(file,self.savepath)

                
class Gui_AvgToJformat:
    """create graphical user from Thinker """
    
    def __init__(self,avgpath,profileFormat, 
                 mode,rotate_angle=180,azimuth=True ):
        self.dir_=avgpath 
        self.type=profileFormat
        self.azimut=azimuth
        self.mode=mode 
        self.rotate_angle=rotate_angle
        
        if self.dir_==None :
            raise FileNotFoundError
        
        self._writeToJformat()
        
    def _writeToJformat(self):
        '''write to jformat file'''
        os.chdir(self.dir_)
        filenames=glob.glob("*.AVG")  # reprend tous les fichiers *avg 
        # profile type files :
        # type1,type2='*.csv','*.txt'
        #change the type of file : 
        ptypefile=glob.glob(self.type)
        comp=0
        #change azimuth to False if azimuth=0 else Default is True   
        # azimut=False 
        
        #mode of survey : Transverse =TM=RYX 
        # mode ='RXY'
        
        for avgfile in filenames : 
        
            if filenames !='':
                fold="{}".format(avgfile[:-4])
                
            os.mkdir(self.dir_+'/{}'.format(fold))# creer un nouveau dossier .
                 
            path_pfile=ptypefile[comp]
            comp+=1
        
            
            # mode ='RXY'
            
            
            # main executor program ---
            ex_f =mconv.AvgToStations(avg_path=avgfile,
                                     mode=self.mode,
                                     rotate_angle=self.rotate_angle)
            ex_f.rawAvgToArray()
            ex_f.arrayToDf()
            ex_f.shrunk_df()
            ex_f.jchainf()
            pf_= mconv.Profile(path_pfile,azimuth=self.azimut)
            pf_.ex_profile()
            pf_.azimuth()
            pf_.build_jfile()
            ex_f.print_()
            
            save_path=os.path.join(self.dir_,'{}/'.format(fold) )
            
            # save_path=r'F:/OneDrive/Python/CodesExercices/' +\
            #             'ex_avgfiles/data/data_avg_control2/{}/'.format(fold)
            savefiles=glob.glob('*.dat')
            for file in savefiles :
        
                shutil.move(file,save_path)
                    
            

class Gui_StnToDat():
    """ Generate graphical user Interface from Tkinter """

    def __init__(self,stnpath=None, normalisation=False ,
                 utm_X_cor=300238.702, utm_Y_cor=2369.252 ):
        
        self.dir_=stnpath
        self.normalised=normalisation
        self.utm_X_cor=utm_X_cor
        self.utm_Y_cor=utm_Y_cor
        
        if self.utm_X_cor==None :
            self.utm_X_cor=0.
        if self.utm_Y_cor is None :
            self.utm_Y_cor=0.
            
        if self.dir_ is None :
            self.dir_=os.getcwd()
        else :
            raise FileNotFoundError
            
        self._rewrite_coordinate_file()
        
    def _rewrite_coordinate_file(self):
  
        os.chdir(self.dir_)
        
        try : 
            os.mkdir('stnToDat_exportfiles')    
        except OSError as e :
            print(os.strerror(e.errno))
        savepath=self.dir_+'/stnToDat_exportfiles'
        filenames=glob.glob('*.stn')
        for stn_file in filenames :
            #Fonctions d'arborescence    
            inst=mconv.StnToDat(filename=stn_file,
                                normalized=self.normalised,
                                utmVar_X=self.utm_X_cor,
                                utmVar_Y=self.utm_Y_cor)
            inst.convertToDat()
        datfile=glob.glob('*.dat')
        for file in datfile: 
            shutil.move(file,savepath)
        

class Gui_ITER2OASIS:
    """Generate GUI from Tkinter """
    
    def __init__(self, path_to_stn,file_format,
                 number_of_stations=1,step_btw_station=50, 
                 file_type=False,choix=False,savefig=False):
        self.dir_=path_to_stn
        self.number_of_stations=number_of_stations
        self.step_btw_station=step_btw_station
        self.file_format=file_format
        self.choix=choix
        self.file_type=file_type
        self.savefig=savefig
        
        self._write_oasfile()
        
    def _write_oasfile(self):
    
        os.chdir(self.dir_)
        
        #-------Main Program---------
        
        iter_filenames=glob.glob('*iter.dat')
        topo_cor_filenames=glob.glob("*_cor.dat")
        bln_filenames=glob.glob("*.bln")
        
        for file in iter_filenames :
            name=file[:-9]
            topofile=name+"_cor.dat"
            if topofile in topo_cor_filenames:
                blnfile=name+".bln"
                if blnfile in bln_filenames: 
                    # Fonction d(arborescence :)
                    exp_=mconv.iter2dToOasis(iterdatfiles=file, 
                              stndatfiles=topofile, 
                              blnfiles=blnfile, 
                              File_format=self.fileformat, step= self.step_btw_station, 
                              NumberOfStation=self.number_of_stations, 
                              FileType=self.file_type,
                              savefig=self.savefig)
                    exp_.arrangeIter2d()
                    exp_.extra_offDepthRho()
                    exp_.reoder()
                    exp_.useChoice(choice=self.choix)
            # Move to savepath
            #choice==True : Dont need to run this code:
        if self.choix==False :
            try : 
                os.mkdir('oasis_exportfiles')    
            except OSError as e :
                print(os.strerror(e.errno))
            savepath=self.dir_+'/oasis_exportfiles'
            ch=['*oas.csv','*oas.xlsx']
            for extend in ch :  
                datfile=glob.glob('*oas.csv')
                for file in datfile: 
                    shutil.move(file,savepath)
        #----End program-----

    
# if __name__=="__main__":
    
#     from tkinter import *
#     fen=Tk()
#     fen.title('AVG_App')
#     fen.geometry("500x50")
#     frame=Frame(fen)
#     frame.grid(row =4, column =4)
#     Button(frame, text="EDIREWRITE",relief="groove", command=Gui_Ediwrite).grid(row=1, column=1)
#     Button(frame, text='AVG-->JFORMAT', relief="groove",command=Gui_AvgToJformat).grid(row=1,column=2)
#     Button(frame, text="EXIT", command=fen.quit).grid(row=2, column=1,sticky=W)
#     fen.mainloop()
#     fen.destroy()
#     # ff.mainloop()

    