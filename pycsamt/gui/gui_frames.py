# -*- coding: utf-8 -*-
"""
Created on Sat Sep  5 10:46:37 2020

@author: Daniel 
    Module is deprecated !
"""

import os 
from tkinter import *
from pycsamt.ff.gui.gui_avgapp import Gui_Ediwrite, Gui_StnToDat,Gui_AvgToJformat,Gui_ITER2OASIS

#=============================


class Press_EDIR(Frame):
    
    def __init__(self, boss=None ):
        Frame.__init__(self, borderwidth=.5, width=900,height=300,relief=GROOVE)
        
        # Fram1=Frame(self,relief=GROOVE, width =700, height =300, borderwidth=2).pack(side=LEFT)
        
        self.haut, self.larg=150,300
        
        Label(self,text="Path to edifiles:").grid(row=0, column=0,sticky=E)
        self.input_path=Entry(self, width=60)
        self.input_path.grid(row=0, column=1,padx=7,pady=5,sticky=W)
        
        Label(self,text="Path to Coordinates File:").grid(row=1, column=0,sticky=E)
        self.input_blnpath=Entry(self,width=60)
        self.input_blnpath.grid(row=1, column=1, padx=7,pady=5,sticky=W)
        
        
        Label(self,text="Coordinates Filename:").grid(row=2, column=0,sticky=E)
        self.input_blnFile=Entry(self,width=35)
        self.input_blnFile.grid(row=3, column=1, padx=7, pady=5,sticky=W)
        
        Label(self,text="savepath:").grid(row=3, column=0,sticky=E)
        self.input_savepath=Entry(self,width=60)
        self.input_savepath.grid(row=2, column=1, padx=7, pady=5, sticky=W)
        
        # create Canvas and Scroolbar 

        self.can = Canvas(self, bg ="black", width =self.larg, height=self.haut, borderwidth=2, 
                     scrollregion=(0,0,2*self.larg,2*self.haut),relief=SUNKEN)
        self.can.grid(row =0, column =2,  rowspan=4)
        
        # trace les axes dans le canvas 
        self.traceAxes()        
        #scroolbars 
        hScroll=Scrollbar(self, orient=HORIZONTAL, command=self.can.xview)
        hScroll.grid(row=4,column=2,columnspan=2, sticky="ew")
        
        vScroll=Scrollbar(self,orient=VERTICAL, command=self.can.yview)
        vScroll.grid(row=0, column=3, rowspan=4, sticky="ns")
        self.can.configure(xscrollcommand=hScroll.set,
                      yscrollcommand=vScroll.set)
        # self.can.configure(scrollregion=self.can.bbox(ALL))
        
        #create Frame 
        bou1=Button(self, text="REWRITE EDI", command=self.rewrite_edi, width=15,height=1, bg="orange")
        bou1.grid(row=4, column=0, columnspan=2, padx=10, pady=5,sticky=W)
        # bou2=Button(boss,text="QUIT", command=boss.quit)
        # bou2.grid(row=3 ,column=3, sticky=E, padx=10, pady=5)
        
        
    def traceAxes(self):
        "Méthode traçant les axes de référence (pourra être surchargée)."
        # axes horizontal (X) et vertical (Y) :
        self.can.create_line(0, 2*self.haut-25, 2*self.larg, 2*(self.haut)-25, 
                             width=.5,fill="white", arrow=LAST)
        self.can.create_line(10, 2*self.haut,10, 0, width=.5,
                             fill="white", arrow=LAST)
        # indication des grandeurs physiques aux extrémités des axes :
        self.can.create_text(30,10, anchor =CENTER,
                             fill="white", text = "Ω.m")
        self.can.create_text(2*(self.larg)-10, 2*(self.haut)-10, fill="white", anchor=CENTER, text="m")       


    def rewrite_edi(self):
        Gui_Ediwrite(edipath=self.input_path.get(),
                          bln_coordinates_filepath=self.input_blnpath.get(),
                          coordinateFile=self.input_blnFile.get(),
                          savepath=self.input_savepath.get())


class Press_STND(Frame):
    
    def __init__(self, boss=None ):
        Frame.__init__(self,borderwidth=.5, width=900,height=300)
        Label(self,text="Path to coordinates file:").grid(row=0, column=0,sticky=W)
        self.input_coordinatepath=Entry(self,width=95)
        self.input_coordinatepath.grid(row=0, column=1,columnspan=2,padx=10, pady=5)

        Label(self,text="diff_UTM_X:").grid(row=1, column=0,sticky=E)
        self.input_utmX=Entry(self,width=35)
        self.input_utmX.grid(row=1, column=1,padx=10, pady=5)
        
        Label(self,text="Diff_UTM_Y:").grid(row=2, column=0,sticky=E)
        self.input_utmY=Entry(self,width=35)
        self.input_utmY.grid(row=2, column=1,padx=10, pady=5)  
        
        self.v1chk=BooleanVar()
        Checkbutton(self,text="Normalised coordinates ", variable=self.v1chk,
                    command=self._set_variable).\
            grid(row=1,column=2, padx=10, pady=5,sticky=E)        
        
        bou1=Button(self, text="STN-->DATCOR", width=15,height=1, bg="orange",command=self.stn_to_datcor)
        bou1.grid(row=3, column=2, columnspan=2, padx=10, pady=5,sticky=E)
        
        # bou2=Button(boss,text="QUIT", command=boss.quit)
        # bou2.grid(row=3 ,column=2, sticky=E, padx=10, pady=5)
        
        # self.root.mainloop()
        # self.root.destroy()


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
        
        
class Press_AVG(Frame):
    
    def __init__(self, boss=None ):
        Frame.__init__(self,borderwidth=.5, width=900,height=300,relief=GROOVE)

        Label(self,text="Path to avgfiles:").grid(row=0, column=0,sticky=E)
        self.input_avgpath=Entry(self,width=100)
        self.input_avgpath.grid(row=0, column=1,columnspan=4,padx=10, pady=5)
        
        Label(self,text="Include Azimuth:").grid(row=1, column=0, sticky=E)
        self.azchk=BooleanVar()
        Checkbutton(self,text="YES/NO", variable=self.azchk,
                     command=self._set_azim).\
             grid(row=1,column=1, padx=5 ,sticky=W)
             
             
        Label(self,text="Choose Profile Extension:").grid(row=2, column=0,sticky=E)
        self.csvchk=BooleanVar()
        Checkbutton(self,text="*.CSV", variable=self.csvchk,
                     command=self._set_type_format).\
             grid(row=2,column=1, padx=5,sticky=W)
             
        self.txtchk=BooleanVar()
        Checkbutton(self,text="*.TXT", variable=self.txtchk,
                     command=self._set_type_format).\
             grid(row=2,column=2, padx=5,sticky=W)
             
        self.datchk=BooleanVar()
        Checkbutton(self,text="*.DAT", variable=self.datchk,
                     command=self._set_type_format).\
             grid(row=2,column=3, padx=5,sticky=W)
             
        Label(self,text="others(*.extension):").grid(row=3, column=1,sticky=E)
        self.other=Entry(self,width=20)
        self.other.grid(row=3, column=2,padx=5,sticky=W)
             
        Label(self,text="Phase rotation angle:").grid(row=4, column=0,sticky=E)
        self.input_phase_rot=Entry(self,width=20)
        self.input_phase_rot.grid(row=4, column=1,padx=5, sticky=W)
        
           
        Label(self,text="Select Mode:").grid(row=5, column=0,sticky=E)
        self.modechk=BooleanVar()
        Checkbutton(self,text="TE MODE<-->RXY(Ex/Hy)", variable=self.modechk,
                     command=self._set_mode).\
             grid(row=5,column=1, padx=10, pady=5,sticky=E)
             
        self.mode2chk=BooleanVar()
        Checkbutton(self,text="TM MODE<-->RYX(Ey/Hx)", variable=self.mode2chk,
                     command=self._set_mode).\
             grid(row=5,column=2, padx=5,sticky=E)
             
         
         # set command Button   
        bou1=Button(self, text="AVG-->JFORMAT",fg="blue",width=15,height=1, bg="orange", command=self._write_jformatFile)
        bou1.grid(row=6, column=3, columnspan=2, padx=10, pady=5,sticky=W)
        # bou2=Button(boss,text="QUIT", fg="red", command=boss.quit)
        # bou2.grid(row=6 ,column=3, columnspan=2, sticky=E, padx=10, pady=5)    
          
        
        # self.root.mainloop()
        # self.root.destroy()
        
        
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
        
        
class Press_OAS(Frame):
    
    def __init__(self, boss=None ):
        Frame.__init__(self,borderwidth=.5, width=900,height=300,relief=GROOVE)
        
        for ii, val in enumerate(["Path to files:","Written Files:","Choose the stations step:",
                "Plot Data:"]):
            if ii>=2:
                
                Label(self,text=val).grid(row=ii+1, column=0,padx=10,sticky=E)
            else :
                Label(self,text=val).grid(row=ii, column=0,padx=10, sticky=E)
            

        # Label(self.root,text="Path to files:").grid(row=0, column=0,sticky=W)
        self.input_pathFiles=Entry(self,width=100)
        self.input_pathFiles.grid(row=0, column=1,columnspan=4,padx=10, pady=5)
        
        # Label(self.root,text="Written Files:").grid(row=1, column=0,sticky=W)
        self.styleFormat=["*.csv","*.xlsx"]
        self.valueFormat=["*.csv", "*.xlsx"]
        self.rad1=StringVar()
        self.rad1.set(self.styleFormat[1])
        
        for ii in range(2):
            boutr=Radiobutton(self, text=self.styleFormat[ii],variable=self.rad1,
                    value=self.valueFormat[ii], 
                    command=self.writeFormat)

            boutr.grid(row=1, column=ii+1, sticky=W)         
        
        # Label(self.root,text="Choose the stations step:").grid(row=3, column=0,sticky=W)
        sca=Scale(self, length=400,orient=HORIZONTAL, sliderlength=40,
              label='Step in meter(m):',from_=0, to=250, tickinterval=50, 
              resolution=0.25, showvalue=50, command=self.stepCount)
        sca.grid(row=3, column=1, columnspan=2)
        
        
        # Label(self.root,text="Plot Data:").grid(row=4, column=0, sticky=W)
        self.choix_var=BooleanVar()
        Checkbutton(self,text="YES/NO", variable=self.choix_var,
                     command=self._set_plot_choice).\
             grid(row=4,column=1, padx=5,sticky=W )
             
        
        for idx, val in enumerate (["Number of the station:","save figure:"]) :
            Label(self, text=val).grid(row=idx+4, column=2,sticky=W)
            
        # Label(self.root,text="Number of the station:").grid(row=4, column=1, )
        
        self.input_numberOfStations=Entry(self,width=15)
        self.input_numberOfStations.grid(row=4, column=2,padx=5, pady=5,sticky=E)
        
        # Label(self.root,text="save figure:").grid(row=5, column=1, sticky=W)
        self.plotFig=BooleanVar()
        Checkbutton(self,text="YES/NO", variable=self.plotFig,
                    command=self._keep_savefig).\
                        grid(row=5,column=2, padx=5 )    

        Label(self,text=" Choose Your Write Excel sheet Type :").\
            grid(row=2, column=1, sticky= E)
        self.boxchoice=['Individual worksheet', 'Many sheets']
        self.value_box=["worksheet", "workbook"]
        self.boxchoice_var=StringVar()
        self.boxchoice_var.set(self.boxchoice[1])
        
        for ii in range(2):
            bou=Radiobutton(self,text=self.boxchoice[ii],
                            variable=self.boxchoice_var, value=self.value_box[ii],
                            command=self._written_type_file)
            if ii==0:
                bou.grid(row=2, column=ii+2, padx=50, sticky=E)
            else :
                bou.grid(row=2, column=ii+2, sticky=W)
             
    
        # set command Button   
        bou1=Button(self, text="ITER-->OASIS",fg="dark green",width=15,height=1, bg="orange", command=self._write_oas_file)
        bou1.grid(row=6, column=0, columnspan=2, padx=10, pady=5,sticky=W)
        # bou2=Button(self,text="QUIT", fg="red", command=self.root.quit)
        # bou2.grid(row=6 ,column=3, columnspan=2, sticky=E, padx=10, pady=5) 
        
        # self.root.mainloop()
        # self.root.destroy()
             
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
    root=Tk()
    root.title('AVGAPP')
    root.geometry("1000x1000+01+01")
    root.minsize("1000","1000")
    root.iconbitmap("ag2dg-xripg-001.ico")
    # include EDIRERWITE
    main_frame=Frame(root,height=100, width=800,relief=RAISED)
    main_frame.grid(row=0, column=0)
    
    frame_Edir=Frame(main_frame, relief=RAISED, borderwidth=3)
    frame_Edir.grid(row=0,column=0)
    Label(frame_Edir, text="EDIREWRITE").grid(row=0,column=0)
    # Label(root, text="-"*190, fg="grey", borderwidth=2 ).grid(row=1, column=0)
    ff=Press_EDIR(frame_Edir)
    ff.configure(borderwidth=2)
    ff.grid(row=1, column=0)
    
    # INCLUDE STNDATCOR 
            
    frame_stn=Frame(main_frame, relief=RAISED,borderwidth=4)
    frame_stn.grid(row=1, column=0,pady=5)
    Label(frame_stn, text="-"*190, fg="grey", borderwidth=2 ).grid(row=0, column=0)
    Label(frame_stn, text="STN2DATCOR").grid(row=1, column=0)
    fn=Press_STND(frame_stn)
    fn.configure(width=1250)
    fn.grid(row=2, column=0 )
    
    # INCLUDE AVG
    
       
    frame_avg=Frame(main_frame, relief=SUNKEN, borderwidth=3)
    frame_avg.grid(row=2, column=0, pady=5)
    Label(main_frame, text="-"*190, fg="grey", borderwidth=2 ).grid(row=0, column=0) 
    Label(frame_avg, text="AVG2JFORMAT").grid(row=1, column=0)
    fa=Press_AVG(frame_avg)
    fa.grid(row=2, column=0)
    
    #INCLUDE OAS 
         
    frame_oas=Frame(main_frame, relief=SUNKEN, borderwidth=3)
    frame_oas.grid(row=3, column=0, padx=5)
    Label(frame_oas, text="-"*190, fg="grey", borderwidth=2 ).grid(row=0, column=0) 
    Label(frame_oas, text="ITER2OASIS").grid(row=1, column=0)
    fo=Press_OAS(frame_oas)
    fo.grid(row=3, column=0)
    
    Button(root, text="EXIT",  bg="red",command=root.destroy).grid(row=4, column=0, sticky=W)
    # main_frame.config(borderwidth=3)
    
    #add scrolbarr 
    
    # scrollbar=Scrollbar(root,relief=RAISED).grid(column=1, rowspan=3)
    # mylist=Listbox(root,yscrollcommand=scrollbar.set)
    # for i in range(3):
    # # for i in [frame_Edir, frame_stn,frame_avg]:
    #     # listbox.insert(END,str(i))
    #     mylist.insert(END,str(i))
    # mylist.pack(side=LEFT, fill=BOTH)
    # scrollbar.configure(command=mylist.yview)
    
    root.mainloop()
    root.destroy()
    






    