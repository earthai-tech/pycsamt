# -*- coding: utf-8 -*-
"""
Created on Sat Sep  5 10:46:37 2020

@author: @Daniel03
"""

import os 
from tkinter import *
from csamtpy.pyCS.gui.gui_avgapp import Gui_Ediwrite, Gui_StnToDat,Gui_AvgToJformat,Gui_ITER2OASIS

#=============================


class Press_EDIR():
    
    def __init__(self):
        self.root=Tk()
        self.root.title('EDIREWRITE')
        self.root.geometry("1216x450-0+20")
        self.root.minsize("1216","450")
        self.root.iconbitmap("ag2dg-xripg-001.ico")
        
        Label(self.root,text="Path to edifiles:").grid(row=0, column=0, sticky=W)
        self.input_path=Entry(self.root, width=70)
        self.input_path.grid(row=0,column=1,padx=5, pady=5)
        
        Label(self.root,text="Path to Coordinates File:").grid(row=0, column=2,sticky=W)
        self.input_blnpath=Entry(self.root,width=70)
        self.input_blnpath.grid(row=0, column=3,padx=10, pady=5)
        
        
        Label(self.root,text="Coordinates Filename:").grid(row=1, column=0,sticky=W)
        self.input_blnFile=Entry(self.root,width=70)
        self.input_blnFile.grid(row=1, column=1,padx=10, pady=5)
        
        Label(self.root,text="savepath:").grid(row=1, column=2,sticky=W)
        self.input_savepath=Entry(self.root,width=70)
        self.input_savepath.grid(row=1, column=3,padx=10, pady=5)
        

        can = Canvas(self.root, bg ="white", width =1000, height =300)
        can.grid(row =2, column =1, columnspan =3)
        
        #create Frame 
        bou1=Button(self.root, text="REWRITE", command=self.rewrite_edi)
        bou1.grid(row=3, column=1, columnspan=2, sticky=W, padx=10, pady=5)
        bou2=Button(self.root,text="QUIT", command=self.root.quit)
        bou2.grid(row=3 ,column=3, sticky=E, padx=10, pady=5)
        
        self.root.mainloop()
        self.root.destroy()


    def rewrite_edi(self):
        Gui_Ediwrite(edipath=self.input_path.get(),
                          bln_coordinates_filepath=self.input_blnpath.get(),
                          coordinateFile=self.input_blnFile.get(),
                          savepath=self.input_savepath.get())


class Press_STND:
    
    def __init__(self):

        self.root=Tk()
        self.root.title('STN-->DATCOR')
        self.root.geometry("800x150-0+20")
        self.root.minsize("800","150")
        self.root.iconbitmap("ag2dg-xripg-001.ico")
        
        Label(self.root,text="Path to coordinates file:").grid(row=0, column=0,sticky=W)
        self.input_coordinatepath=Entry(self.root,width=95)
        self.input_coordinatepath.grid(row=0, column=1,columnspan=2,padx=10, pady=5)

        Label(self.root,text="diff_UTM_X:").grid(row=1, column=0,sticky=E)
        self.input_utmX=Entry(self.root,width=35)
        self.input_utmX.grid(row=1, column=1,padx=10, pady=5)
        
        Label(self.root,text="Diff_UTM_Y:").grid(row=2, column=0,sticky=E)
        self.input_utmY=Entry(self.root,width=35)
        self.input_utmY.grid(row=2, column=1,padx=10, pady=5)  
        
        self.v1chk=BooleanVar()
        Checkbutton(self.root,text="Normalised coordinates ", variable=self.v1chk,
                    command=self._set_variable).\
            grid(row=1,column=2, padx=10, pady=5,sticky=E)        
        
        bou1=Button(self.root, text="STN-->DATCOR", command=self.stn_to_datcor)
        bou1.grid(row=3, column=1,columnspan=2, sticky=W, padx=10, pady=5)
        
        bou2=Button(self.root,text="QUIT", command=self.root.quit)
        bou2.grid(row=3 ,column=2, sticky=E, padx=10, pady=5)
        
        self.root.mainloop()
        self.root.destroy()


    def stn_to_datcor(self):
        
        if self._set_variable() ==True :
            normalization=True 
        else :
            normalization=False
            
        Gui_StnToDat(stnpath=self.input_coordinatepath.get(),
                          normalisation=normalization,
                          utm_X_cor=float(self.input_utmX.get()),
                          utm_Y_cor=float(self.input_utmY.get()))
        
    
    def _set_variable(self):
        if self.v1chk is True :
            return True 
        else :
            return False 
        
        
class Press_AVG:
    
    def __init__(self):
        
        self.root=Tk()
        self.root.title('AVG-->JFORMAT')
        self.root.geometry("1000x220-0+20")
        self.root.minsize("1000","220")
        self.root.iconbitmap("ag2dg-xripg-001.ico") 
        
        Label(self.root,text="Path to avgfiles:").grid(row=0, column=0,sticky=W)
        self.input_avgpath=Entry(self.root,width=125)
        self.input_avgpath.grid(row=0, column=1,columnspan=4,padx=10, pady=5)
        
        Label(self.root,text="Include Azimuth:").grid(row=1, column=0, sticky=W)
        self.azchk=BooleanVar()
        Checkbutton(self.root,text="YES/NO", variable=self.azchk,
                     command=self._set_azim).\
             grid(row=1,column=1, padx=5 ,sticky=W)
             
             
        Label(self.root,text="Choose Profile Extension:").grid(row=2, column=0,sticky=E)
        self.csvchk=BooleanVar()
        Checkbutton(self.root,text="*.CSV", variable=self.csvchk,
                     command=self._set_type_format).\
             grid(row=2,column=1, padx=5,sticky=W)
             
        self.txtchk=BooleanVar()
        Checkbutton(self.root,text="*.TXT", variable=self.txtchk,
                     command=self._set_type_format).\
             grid(row=2,column=2, padx=5,sticky=W)
             
        self.datchk=BooleanVar()
        Checkbutton(self.root,text="*.DAT", variable=self.datchk,
                     command=self._set_type_format).\
             grid(row=2,column=3, padx=5,sticky=W)
             
        Label(self.root,text="others(*.extension):").grid(row=3, column=1,sticky=E)
        self.other=Entry(self.root,width=20)
        self.other.grid(row=3, column=2,padx=5,sticky=W)
             
        Label(self.root,text="Phase rotation angle:").grid(row=4, column=0,sticky=E)
        self.input_phase_rot=Entry(self.root,width=20)
        self.input_phase_rot.grid(row=4, column=1,padx=5, sticky=W)
        
           
        Label(self.root,text="Select Mode:").grid(row=5, column=0,sticky=E)
        self.modechk=BooleanVar()
        Checkbutton(self.root,text="TE MODE<-->RXY(Ex/Hy)", variable=self.modechk,
                     command=self._set_mode).\
             grid(row=5,column=1, padx=10, pady=5,sticky=E)
             
        self.mode2chk=BooleanVar()
        Checkbutton(self.root,text="TM MODE<-->RYX(Ey/Hx)", variable=self.mode2chk,
                     command=self._set_mode).\
             grid(row=5,column=2, padx=5,sticky=E)
             
         
         # set command Button   
        bou1=Button(self.root, text="AVG-->JFORMAT",fg="blue", command=self._write_jformatFile)
        bou1.grid(row=6, column=1,columnspan=2, sticky=W, padx=10, pady=5)
        bou2=Button(self.root,text="QUIT", fg="red", command=self.root.quit)
        bou2.grid(row=6 ,column=3, columnspan=2, sticky=E, padx=10, pady=5)    
          
        
        self.root.mainloop()
        self.root.destroy()
        
        
    def _set_type_format (self):
        if self.csvchk is True :
            return "*.csv"
        elif self.txtchk is True :
            return "*.txt"
        elif self.datchk is True :
            return "*.dat"
        elif self.other is not None :
            self.other.get()
            if "." in self.other.get():
                return "*{0}".format(self.other.get())
            elif "*." in self.other.get():
                return self.other.get()
            else :
                return "*.{0}".format(self.other.get())
        else :
            raise TypeError
            
    def _set_azim(self):
        if self.azchk is True:
            return True 
        else :
            return False 
    
    def _set_mode(self):
        if self.modechk is True :
            self.mode="RXY"
        elif self.mode2chk is True :
            self.mode="RYX"
        elif (self.mode2chk==True and self.modechk==True) or \
            (self.mode2chk is False and self.modechk is False) :
            self.mode="RXY"

        return self.mode

    
    def _write_jformatFile(self):
        
        Gui_AvgToJformat(avgpath=self.input_avgpath.get(),
                         profileFormat=self._set_type_format(), 
                         mode=self._set_mode,
                         rotate_angle=int(self.input_phase_rot.get()),
                         azimuth=self._set_azim)
        
        
class Press_OAS:
    
    def __init__(self):
        self.root=Tk()
        self.root.title("OASIS")
        self.root.geometry("980x220-0+20")
        self.root.minsize("980","250")
        self.root.iconbitmap("ag2dg-xripg-001.ico")
        
        for ii, val in enumerate(["Path to files:","Written Files","Choose the stations step:",
                "Plot Data:"]):
            if ii>=2:
                
                Label(self.root,text=val).grid(row=ii+1, column=0,sticky=W)
            else :
                Label(self.root,text=val).grid(row=ii, column=0,sticky=W)
            

        # Label(self.root,text="Path to files:").grid(row=0, column=0,sticky=W)
        self.input_pathFiles=Entry(self.root,width=125)
        self.input_pathFiles.grid(row=0, column=1,columnspan=4,padx=10, pady=5)
        
        # Label(self.root,text="Written Files:").grid(row=1, column=0,sticky=W)
        self.styleFormat=["*.csv","*.xlsx"]
        self.valueFormat=["*.csv", "*.xlsx"]
        self.rad1=StringVar()
        self.rad1.set(self.styleFormat[1])
        
        for ii in range(2):
            boutr=Radiobutton(self.root, text=self.styleFormat[ii],variable=self.rad1,
                    value=self.valueFormat[ii], 
                    command=self.writeFormat)

            boutr.grid(row=1, column=ii+1, sticky=W)         
        
        # Label(self.root,text="Choose the stations step:").grid(row=3, column=0,sticky=W)
        sca=Scale(self.root, length=700,orient=HORIZONTAL, sliderlength=40,
              label='Step in meter(m):',from_=0, to=300, tickinterval=50, 
              resolution=0.25, showvalue=50, command=self.stepCount)
        sca.grid(row=3, column=1, columnspan=2)
        
        
        # Label(self.root,text="Plot Data:").grid(row=4, column=0, sticky=W)
        self.choix_var=BooleanVar()
        Checkbutton(self.root,text="YES/NO", variable=self.choix_var,
                     command=self._set_plot_choice).\
             grid(row=4,column=1, padx=5,sticky=W )
             
        
        for idx, val in enumerate (["Number of the station:","save figure:"]) :
            Label(self.root, text=val).grid(row=idx+4, column=2,sticky=W)
            
        # Label(self.root,text="Number of the station:").grid(row=4, column=1, )
        
        self.input_numberOfStations=Entry(self.root,width=30)
        self.input_numberOfStations.grid(row=4, column=2,columnspan=2,padx=10, pady=5)
        
        # Label(self.root,text="save figure:").grid(row=5, column=1, sticky=W)
        self.plotFig=BooleanVar()
        Checkbutton(self.root,text="YES/NO", variable=self.plotFig,
                    command=self._keep_savefig).\
                        grid(row=5,column=2, padx=5 )    

        Label(self.root,text=" Choose Your Write Excel sheet Type :").\
            grid(row=2, column=1, sticky= E)
        self.boxchoice=['Individual worksheet', 'Many sheets']
        self.value_box=["worksheet", "workbook"]
        self.boxchoice_var=StringVar()
        self.boxchoice_var.set(self.boxchoice[1])
        
        for ii in range(2):
            bou=Radiobutton(self.root,text=self.boxchoice[ii],
                            variable=self.boxchoice_var, value=self.value_box[ii],
                            command=self._written_type_file)
            if ii==0:
                bou.grid(row=2, column=ii+2, padx=50, sticky=E)
            else :
                bou.grid(row=2, column=ii+2, sticky=W)
             
    
        # set command Button   
        bou1=Button(self.root, text="ITER-->OASIS",fg="dark green", command=self._write_oas_file)
        bou1.grid(row=6, column=1,columnspan=2, sticky=W, padx=10, pady=5)
        bou2=Button(self.root,text="QUIT", fg="red", command=self.root.quit)
        bou2.grid(row=6 ,column=3, columnspan=2, sticky=E, padx=10, pady=5) 
        
        self.root.mainloop()
        self.root.destroy()
             
    def stepCount(self, step):
        
        self._step=int(float(step))

    
    def _keep_savefig(self):
        if self.plotFig.get():
            return True 
    
    def _set_plot_choice(self):
        
        if self.choix_var.get()==True:
            return int(float(self.input_numberOfStations.get()))
        else :
            return False 
         
    def writeFormat(self):
        if self.rad1.get()==self.valueFormat[0]:
            return True

        elif self.rad1.get()==self.valueFormat[1]:
            return False
            
                
    def _written_type_file(self):
        
        if self.boxchoice_var.get()==self.value_box[0]:
            return False 
        else :
            return True
            
 
    def _write_oas_file(self):

        Gui_ITER2OASIS(path_to_stn=self.input_pathFiles.get(),
                       file_type=self._written_type_file(),
                       file_format=self.writeFormat(),
                       choix=self._set_plot_choice(),
                       number_of_stations=self._set_plot_choice(),
                       step_btw_station=self._step,
                       savefig=self._keep_savefig())
  

if __name__=="__main__":
    # ff=Press_EDIR()
    fn=Press_OAS()






    