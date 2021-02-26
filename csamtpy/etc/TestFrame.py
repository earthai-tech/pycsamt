#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 19:15:17 2020

@author: Daniel03
"""

import os 
from tkinter import *
from gui_avgapp import Gui_Ediwrite

#=============================


class Press_EDIR(Frame):
    
    def __init__(self, master =None, coul='grey'):
        Frame.__init__(self)
        
        Label(self,text="Path to edifiles:").grid(row=0, column=0, sticky=W)
        self.input_path=Entry(self, width=70)
        self.input_path.grid(row=0,column=1,padx=5, pady=5)
        
        Label(self,text="Path to Coordinates File:").grid(row=0, column=2,sticky=W)
        self.input_blnpath=Entry(self,width=70)
        self.input_blnpath.grid(row=0, column=3,padx=10, pady=5)
        
        
        Label(self,text="Coordinates Filename:").grid(row=1, column=0,sticky=W)
        self.input_blnFile=Entry(self,width=70)
        self.input_blnFile.grid(row=1, column=1,padx=10, pady=5)
        
        Label(self,text="savepath:").grid(row=1, column=2,sticky=W)
        self.input_savepath=Entry(self,width=70)
        self.input_savepath.grid(row=1, column=3,padx=10, pady=5)
        
        
        can = Canvas(self, bg ="white", width =1000, height =300)
        can.grid(row =2, column =1, columnspan =3)
        
        #create Frame 
        bou1=Button(self, text="REWRITE", command=self.rewrite_edi)
        bou1.grid(row=3, column=1, columnspan=2, sticky=W, padx=10, pady=5)
        # bou2=Button(self,text="QUIT", command=self.root.quit)
        # bou2.grid(row=3 ,column=3, sticky=E, padx=10, pady=5)
        
        # self.root.mainloop()
        # self.root.destroy()
        self.master.bind("<control-Z>", self.rewrite_edi)
        self.master.title("EDIREWRITE")
        self.pack()


    def rewrite_edi(self):
        Gui_Ediwrite(edipath=self.input_path.get(),
                          bln_coordinates_filepath=self.input_blnpath.get(),
                          coordinateFile=self.input_blnFile.get(),
                          savepath=self.input_savepath.get())

if __name__=="__main__":
    fen=Tk()
    ff=Press_EDIR(master=fen)
    ff.mainloop()

    # fn=Press_STND()






    