# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 19:08:08 2020

@author: Daniel03

Main Avg Application Using Tkinter graphical Interface 

"""

#++++++++++++++++++++++++++
from tkinter import *
import gui_avgapp as gui
from gui_edirewrite import Press_EDIR, Press_STND, Press_AVG,Press_OAS

#++++++++++++++++++++++++++

fen=Tk()
fen.title('MyApplication')
fen.geometry("500x70")
fen.minsize("500","70")
fen.iconbitmap("ag2dg-xripg-001.ico")
# frame=LabelFrame(fen,text="Important controls")
frame=Frame(fen)
frame.grid(row =4, column =4)
Button(frame, text="EDIREWRITE",relief="groove", 
       command=Press_EDIR).grid(row=0, column=0,pady=5)
Button(frame, text='AVG-->JFORMAT', relief="groove",
       command= Press_AVG).grid(row=0,column=1,padx=10, pady=5)
Button(frame, text='STN-->DATCOR', relief="groove",
       command=Press_STND).grid(row=0,column=2,padx=10, pady=5)
Button(frame, text='ITER-->OASIS', relief="groove",
       command=Press_OAS).grid(row=0,column=3,padx=10, pady=5)
Button(frame, text="EXIT", command=fen.quit).\
    grid(row=1, column=0, sticky= W)

fen.mainloop()
fen.destroy()


