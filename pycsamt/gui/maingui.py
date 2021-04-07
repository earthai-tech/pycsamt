# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 19:08:08 2020

    !!!Deprecated module !!!

@author: Daniel03

Main Avg Application Using Tkinter graphical Interface 

"""

#++++++++++++++++++++++++++
import os 
from tkinter import *
import pycsamt.ff.gui.gui_avgapp as gui
from pycsamt.ff.gui.gui_frames import Press_EDIR, Press_STND, Press_AVG,Press_OAS

#++++++++++++++++++++++++++
def create_rect(cannev, x, y, color="white"):
    canev.create_rectangle(x,y,x+r,y+r, fill=color)


fen=Tk()
fen.title('pyCS')
fen.geometry("900x920")
fen.minsize("900","920")
imgpath=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                     "ag2dg-xripg-001.ico")

# fen.iconbitmap("ag2dg-xripg-001.ico")
# os.chdir(os.path.dirname(imgpath))
fen.iconbitmap(os.path.basename(imgpath))
# frame=LabelFrame(fen,text="Important controls")

can=Canvas(fen, relief=SUNKEN, borderwidth=2, scrollregion=(0,0,2000,1500),
           bg="grey",confine=False)
can.grid(row=0, column=0,padx=5, pady=5)

#scroolbars 
hScroll=Scrollbar(fen, orient=HORIZONTAL, command=can.xview)
hScroll.grid(row=1,column=0, sticky="ew")

vScroll=Scrollbar(fen,orient=VERTICAL, command=can.yview)
vScroll.grid(row=0, column=1, sticky="ns")
can.configure(xscrollcommand=hScroll.set,
              yscrollcommand=vScroll.set)
# can.configure(scrollregion=can.bbox(ALL))

main_frame=Frame(can,relief=RAISED, borderwidth=4)
main_frame.grid(in_=can, row=0,column=0, pady=10, padx=10)

#"STN2DATCOR":Press_STND
  

#make a loop"STN2DATCOR":Press_STND
dico_func={ "AVG2JFORMAT": Press_AVG,
            "EDIREWRITE":Press_EDIR ,
            "ITER2OASIS":Press_OAS}

for ii, val in enumerate(list(dico_func)):
    comp=0
    if ii%2 !=0 :
        comp=1
    ff=Frame(main_frame)
    Label(ff,text=val.upper(), relief=RAISED,fg="red").grid(in_=ff, row=comp, column=0)
    fs=dico_func[val](boss=ff)
    fs.grid(in_=ff, row=comp+1,column=0, pady=2)
    fs.configure(borderwidth=.5, width=850,height=250)
    fs.grid_propagate(0)
    ff.grid(row=ii, column=0, pady=5)
    
fen.mainloop()



