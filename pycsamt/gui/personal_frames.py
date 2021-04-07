# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 13:08:02 2020
Deprecated Module !

@author: Daniel03

"""
from tkinter import *


class ButtonPerso(Canvas):
    
    def __init__(self, boss=None, **kwargs):
        Canvas.__init__(self,**kwargs)
        self.boss=boss
        self.border=kwargs.pop("borderwidth",3)
        self.height=kwargs.pop("height",600)
        self.width=kwargs.pop("width",300)
        self.command=kwargs.pop("command",None)
        for key, value in kwargs.items():
            setattr(self, key, value)
        
class CanvasPerso(Canvas):
    
    def __init__(self, boss=None, **kwargs):
        Canvas.__init__(self,borderwidth=borderwith, 
                        height=height, width=width, **kwargs)
        self.boss=boss
        self.border=kwargs.pop("borderwidth",3)
        self.height=kwargs.pop("height",600)
        self.width=kwargs.pop("width",300)
        
        for key, value in kwargs.items():
            setattr(self, key, value)
        
class FramePerso(Frame):
    
    def __init__(self, boss=None, **kwargs):
        Frame.__init__(self, boss,borderwidth=borderwidth,
                       height=height, width=width, **kwargs)
        self.boss=boss
        self.border=kwargs.pop("borderwidth",2)
        self.height=kwargs.pop("height",300)
        self.width=kwargs.pop("width",150)
        
        for key, value in kwargs.items():
            setattr(self, key, value)
        
        
        
if __name__=="__main__":
    
    root=Tk()
    root.title('AVGAPP')
    root.geometry("550x350+01+01")
    root.minsize("550","350")
    # root.iconbitmap("ag2dg-xripg-001.ico")
    
    fs=CanvasPerso(root,borderwidth=3,height=200,
                   width=300,relief="sunken" )
    fs.pack(expand=YES)
        
    
        





