# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 12:01:17 2020

@author: Administrator
"""
import os 
from tkinter import *
from gui_avgapp import Gui_Ediwrite, Gui_StnToDat,Gui_AvgToJformat
from gui_edirewrite import Press_EDIR,Press_STND, Press_AVG, Press_OAS

#=============================

class MenuBar_EDIR(Frame):
    
    def __init__(self, boss=None):
        Frame.__init__(self, borderwidth=2,relief=GROOVE )
        #scrollingMenu Bar 
        fileMenu=Menubutton(self,text="EDIREWRITE")
        # fileMenu.pack(side=LEFT)
        fileMenu.grid(row=0, column=0, padx=20, pady=7, sticky=W )
        ##### Menu EDIREWRIE######
        me1=Menu(fileMenu)
        me1.add_command(label="Open", underline=0,
                        command=Press_EDIR)
        # integration of file Menu 
        fileMenu.configure(menu=me1)
        
        
        #### add AVG-->JFORMAT to Menu Bar ######
        self.avg_jformat=Menubutton(self, text="AVG-->JFORMAT")
        # self.avg_jformat.pack(side=LEFT, padx=10)
        self.avg_jformat.grid(row=0, column=1, padx=20,pady=7, sticky=W )
        # scrolling of Menu<AVG-->JFORMAT> 
        me1=Menu(self.avg_jformat)
        me1.add_command(label="Open", underline=0,
                        command=Press_AVG)

        
        self.avg_jformat.configure(menu=me1)
        
        #### add STN-->DATCOR to Menu Bar ######
        self.stn=Menubutton(self, text="STN-->DATCOR")#,font=('Times',12,"bold"))
        # self.stn.pack(side=LEFT, padx=10)
        self.stn.grid(row=0, column=2, padx=20,pady=7, sticky=W )
        # scrolling of Menu<STN-->DATCOR> 
        me1=Menu(self.stn)
        me1.add_command(label="Open", underline=0,
                        command=Press_STND)

        
        self.stn.configure(menu=me1)
        
        #### add ITER-->OASISto Menu Bar ######
        self.iteroas=Menubutton(self, text="ITER-->OASIS")
        # self.iteroas.pack(side=LEFT, padx=10)
        self.iteroas.grid(row=0, column=3, padx=20,pady=7, sticky=W )
        # scrolling of Menu<ITER-->OASIS> 
        me1=Menu(self.iteroas)
        me1.add_command(label="Open", underline=0,
                        command=Press_OAS)

        self.iteroas.configure(menu=me1)
        
        # Frame.pack(self,expand=YES)
        
                
class Application(Frame):
    
    def __init__(self,boss=None):
        Frame.__init__(self)
        self.master.title("AVG APP")
        self.master.geometry("650x250-0+20")
        self.master.minsize("650","250")
        self.master.iconbitmap("ag2dg-xripg-001.ico")
        
        mBar=MenuBar_EDIR(self)
        mBar.pack()
        # mBar.grid(row=0, column=0, sticky=W)
        self.can = Canvas(self.master, bg ="light grey", height =180,width =520, borderwidth=2,
                          relief=RAISED)
        self.can.pack()
        # self.can.grid(row=1, column=0)
        self.can.create_text(20,70,anchor=NW, text="Welcome to AVG-APPLICATION", 
                             font=('Comic sans MS',20,"italic"),fill="red")
        self.can.create_text(430,165,anchor=NW, text="Version V.00.01", 
                             font=('Comic sans MS',10,"italic"),fill="red")
        
        self.bou=Button(self.master, text="QUIT", font=("bold",10), command=self.master.quit)
        self.bou.pack(side=RIGHT)#height=90, width=500)
        
        
    def clean_can(self):
        self.can.delete(ALL)
    
if __name__=="__main__":
    Application().mainloop()

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    