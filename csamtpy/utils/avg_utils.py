# -*- coding: utf-8 -*-
"""
===============================================================================
    Copyright © 2021  Kouadio K.Laurent
    
    This file is part of pyCSAMT.
    
    pyCSAMT is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    pyCSAMT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
    
    You should have received a copy of the GNU Lesser General Public License
    along with pyCSAMT.  If not, see <https://www.gnu.org/licenses/>.

===============================================================================  

Created on Tue Aug  4 16:03:30 2020
    module will be deprecated soon !

@author: @Daniel03

-------
Classes :
    
    ReadFile : Extract informations of files with index 
    ++++++++
    
        methods :
        ---------
            * read_file : Readile file for extractiong data 
                
            * assert_index_value : assert index for sclinf infrmations on file 

Function utils : 
    
    * ll2DMS: 
        fonction convert degree decimal to DD:MM:SS
    * roundVar : 
        round value with exactitude 
    * convert_phase_from_pdSeries(df, columns=None)
        return dataframe 
    * ll_to_utm2 : 
        return utm_zone, utm_easting, utm_northing
    *convert_rotate_phase : 
        return exactitude float (phase value ) according to geometry principal phase value
    *convert_phase_from_pdSeries : 
        return  a specified pd.Series :df[columns]
   
"""
#==============================================================================
import os 

import csamtpy.utils.gis_tools as gis 
import numpy as np 
import pandas as pd 

from scipy.interpolate import interp1d

#==============================================================================


def ll2DMS(position_val): 
    """Conversion function : Convert degree decimals (long and lat) 
    value on string DD:MM:SS.
    
     Argument : 
    ------------
    **longitude or latitude value in degree decimals  ***
    
    returns : 
    ---------
    ***Value converted in  string format : DD:MM:SS
    """
    # verifie la condition si la valeur est bien 
    mm,ss="00","00"
    if position_val in [None, 'None']:
        return None
    if type(position_val) is str:
        raise ValueError("{0} not correct format, should be {{Degree decimals}}".format(position_val))
        return 
    else : 
        deg=int(position_val)
        decim=position_val- deg
        if decim==0 :
            return "{0}:{1}:{1}:".format(deg,mm,ss)
        elif decim !=0 : 
            minute=decim*60 
            dec_1=minute-int(minute)
            if int(minute)<10: 
                mm="0{}".format(int(minute))
            else : 
                mm="{}".format(int(minute))
            if dec_1==0 :
                return"{0}:{1}:{2}".format(deg,mm,ss)
            else : 
                sec=dec_1*60
                seconde=roundVar(sec)
                if int(seconde)<10: 
                    ss="0{}".format(int(seconde))
                elif int(seconde)==60:
                    min_=int(minute)+1
                    if min_<10:
                        return "{0}:0{1}:{2}".format(deg,min_,ss)
                    elif min_==60: 
                        deg=deg+1
                        mm="00"
                        return "{0}:{1}:{2}".format(deg,mm,ss)
                    else:
                        return "{0}:{1}:{2}".format(deg,min_,ss)
                else : 
                    ss="{}".format(int(seconde))
                return "{0}:{1}:{2}".format(deg,mm,ss)
    
def roundVar(value): 
    """
    Fonction rounds value with eatitude 
    Parameters : Float number 
    ----------
    Returns
    -------
    return value rounded.
    """     
    decim=value-int(value)
    if decim<=.5:
        val=int(value)
    elif decim>.5 : 
        new_value=value+.5
        val=int(new_value)
    return val

def ll_to_utm2( lat, lon,reference_ellipsoid=23):
    """
    Function redirect convertion lat, long to Easting , northing, 
    
    Parameters
    ----------
    lat : float
        latitude in degree decimals.
    lon : TYPE
        longitude in degree decimals 
    reference_ellipsoid : integer, Check module gis_tools:
         The default is 23.matching the "WGS-84"

    Returns
    -------
    utm_zone : string
        Zone of survey .
    utm_easting : float
        cordinates of the area survey in meters(UTM) .
    utm_northing : float
        cordinates of the of the area survey in meter (UTM) .

    """
    utm_zone, utm_easting, utm_northing = gis.ll_to_utm(reference_ellipsoid,lat,lon)
    
    return utm_zone, utm_easting, utm_northing

def convert_rotate_phase(value,rot_phase=180):
    """
    Parameters
    ----------
    value : float or int 
        DESCRIPTION.
    Returns
    -------
    value : float value 
        Geometry excatitude phase value.
    
    Examples :
    ---------
        import numpy as np
        np.random.seed(0)
        ss=np.array([convert_rotate_phase(float(ii)) for\
                 ii in np.linspace(1,525,50) if float(ii)<=382] )
        print(ss,"\n")
        print(np.linspace(1,525,50))
    """
    rot_phase=int(rot_phase)

    assert type(value)==float or type(value)==int,\
        'Unacceptable value:"{}", must be a float or int.'.format(value)

    if value > 0:
        if value >=rot_phase:
          return float(value%rot_phase)  

    elif value < 0: 
        if abs(value) >=rot_phase :
          return float(value%(-rot_phase))  

    return float(value) 


def convert_phase_from_pdSeries(df, columns=None,rotate=180) :
    """
    Parameters
    ----------
    df: object 
        ***df ***           must be a pandas dataframe 
        ***columns ***.     A string : a specified columns of the dataframe 

    Raises
    ------
    ValueError
        Not a float number

    Returns
    -------
    df[columns]: PandaSeries
        returns a specified pd.Series.
    """

    if columns==None : 
        return None 
    
    try :
        df.reset_index(drop=True)
    except:
        raise TypeError('{} must be a pd.Series'.format(df[columns]))

    for index, val in enumerate (df[columns]):
        b=convert_rotate_phase(val,rot_phase=rotate)
        df[columns].iloc[index]=b    
    return df[columns]
        
class ReadFile:
    """
    Class to read simple file : 
        it uses file wich contains head 
        
    Attributes :
    ============
    
    *FilePath : str or bytes
                path of work directory 
        
    * file_to_read : str 
                    file object can be a *.csv, .txt , or .bln, .dat .
        
    * sep: NoneType
                    string for spliting data  : ex : "," or "." or ":" or ";"
                    defaut is None used conventional method object.split()
        
    * ex_index : int or tuple
                    extraction index for the data in the file (columns).
                    set tuple for many index extraction :
                        ex: ex_index=9 or ex_index = -1
                            ex_index=(9,2,0) or (-10,-2,3)
            
    * headline : boolean 
                    the head of the file . if the file contains no head , 
                    set to FALSE , if not default is TRUE.
                    
    * head_range : int
                    the range of the head in the file
        
     *slicing : tuple or tuple of list 
                    if headline =True , default slice is (1, None), it removes 
                    the first headline
                     sets option : 
                         slicing must be a list : 
                             slicing=([1,4,9],[2,7,None] )
                             first index range :[1,4,9] is the start 
                             second index range :[2,7,None] is the end 
                        ---> [1:2],[4:7],[9:]
        
    methods :
    =========
    
    * read_file : return a dictionary for the data extract : 
                    key is the ilist_{number of ex_index}:
                    value : list of ex_index  contains string elements. 
                    ex (ex_index=9): {ilist_09:['573.40', '566.80', '553.95',...,
                                    '424.35', '411.50', '423.75']}
                    ex(ex_index=(9,2,0)) : {index_09:[573.40', '566.80'...],
                                            index_02:['110.48647', '110.48690',...,
                                            index_00:['S00','S01','S03'',...]}
                                              
    * assert_index_value : return int value for indexing 
                    np.random.seed(0)
                    su=(None,None,None)
                    print(ff.assert_index_value(ta))
                    ch=[round(15*ii,3) for ii in np.random.rand(15)]
                    ch2=tuple(ch)
                    ss=ff.assert_index_value(ch2)
                    ts= [ff.assert_index_value(ii) for ii in su ]
                    print(ch,"\n",ss,'\n',ts)
            
    * search_file : return datafile 
                    if FilePath is None ,search the Path of matching file
                        on the current directory. 
                    if file_to_read is None , search file on current directory. 
            
         
    Test: 
        if __name__=="__main__":
        >>>os.chdir(os.getcwd())
        >>>ff=ReadFile(FilePath =None, 
                       file_to_read="K1_exp.bln",
                       head_range=None,
                       ex_index=(8,-1),
                       slicing=([1,4,9],[2,7,None]))
        >>>  dd=ff.read_file()
        >>>  d1=[float(ss) for ss in dd['index_08']]
        >>>  d2=[float(ss) for ss in dd['index_09']]
        >>>  d1,d2=np.array(d1),np.array(d2)
        >>>  ary=np.concatenate((d1.reshape((1,d1.shape[0])),\
                                d2.reshape((1,d2.shape[0]))),axis=0)
        >>>  print(ary)
        
    """
    
    
    def __init__(self, **kwargs):
        
        self.path=kwargs.pop("FilePath",None)
        self.data_index=kwargs.pop("ex_index",None)
        self.separate=kwargs.pop("sep",None)
        self.headline=kwargs.pop("headline",True)
        self.range_head=kwargs.pop("head_range",None)
        self.slicing=kwargs.pop("slicing",(None,None))
        
        self._flag,self._comp,self._emp=-1,-1,[]
        self._option, self._start, self._end=-1,None,None
                  
        if self.path is not None :
            if os.path.exists(self.path):
                self.file=kwargs.pop("file_to_read",None)
                if self.file is not None :
                    assert os.path.isfile(self.file),"{} don't exist"\
                        " in this directory.Please change your directory".\
                            format(self.file)
                    self.file=self.file
                elif self.file is None :
                    self.file=self.search_file_(self.path)
            else :
                raise FileExistsError('No such file in such directory')
                
        elif self.path is None :
            self.file=self.file=kwargs.pop("file_to_read",None)
            if self.file is None :
                self.file=self.search_file_(self.path)
                
        self.data_index=np.array(self.data_index)
        
        if np.iterable(self.data_index):
            self.indexl=self.data_index.tolist()
            for i,value in enumerate (self.indexl) :
                if type(value)==str :
                    self.indexl[i]=int(value)
            self.max_index=abs(np.asarray(self.indexl)).max()    
            self.indexl=[int(ii) for ii in self.indexl \
                    if (type(ii)==int) or (type(ii)==float)]
            self._flag=1

        if np.iterable(self.data_index)==False :
            self.indexl=self.data_index
            if type(self.indexl) == float:
                self.indexl=int(self._index) 
            self._flag=0
            
        self._start=self.assert_index_value(self.slicing[0])
        self._end=self.assert_index_value(self.slicing[1])
        
        if self.headline == True:
            if self.range_head is None :
                self.range_head=1
            elif self.range_head is not None :
                assert type(self.range_head)==int or type(self.range_head)==float,\
                    "Unacceptable value {}, must be integer".format(self.range_head)
                self.range_head=self.assert_index_value(self.range_head)
                
        if np.iterable(self._start) or np.iterable(self._end):
            
            assert len(self._start)==len(self._end),"slicing argument {0},"\
                " start= {1} and end= {2} must be the same length".\
                    format(self.slicing, self._start, self._end)  
                    
            if type(self._start)==tuple or type(self._end)==tuple :
                self._start, self._end=list(self._start),list(self._end)
            self._option=1

        elif np.iterable(self._start)==False or np.iterable(self._end)==False :
            self._option=0
            
    @property
    def separate(self):
        """check property attribute"""
        return self._separate
    
    @separate.setter
    def separate(self,value):
        if type(value)==float or type(value)==int :
            raise TypeError('Wrong attribute value: {0}.'\
                            ' Must be a string. Separate file Options are :'\
                                '";" or "." or ":" or "."'.format(value))
        self._separate= value
        

    def __str__(self):
        print("Welcome To Readline Class:")
        ReadFile.__doc__

    def read_file(self):
        """
        Methode to extract data in the file from index 
        
        parameters :
        -----------
            * headline : Boolean 
                if headline =True , default slice is (1, None), 
                it removes the first headline
                 sets option slicing : 
                     slicing must be a list : 
                         slicing=([1,4,9],[2,7,None] )
                         fisrt index :[1,4,9] is the start 
                         second index :[2,7,None] is the end 
                    ---> [1:2],[4:7],[9:]
        Returns
        -------
        dictionary of index data
        
        Examples:
        --------
           import os 
           import numpy as np
           path=os.getcwd()
           dir =os.path.basename(path)
           >>> /K1_exp.bln
           >>> dic_ex=read_file(dir,(0,-9,3,-2,6),
                            headline=True,sep=" ")
           
          >>> sta=dic_ex['index_00']
          >>> azim=np.array([float(ss) for ss in dic_ex['index_08']])
          >>> print(azim)
        
        """
        
        # if self.file ==None :
        #     raise FileExistsError("No file in this directory: {}".format(os.getcwd()))
        if self.indexl== None :
            raise IndexError("No index for extraction! Index must be int, float , or tuple")
        
        try:
            with open (self.file,'r') as f :
                listfiles=f.readlines()
        
        except :
            FileNotFoundError('No such file in this directory: {0}.'.format(os.getcwd()))
        
        if self.headline:
            listfiles=listfiles[self.range_head:]
        
        # Option of slicing
        if self._option==0 :

            if self._start==None and self._end !=None :
                listfiles=listfiles[:self._end]
            elif self._start !=None and self._end ==None :
                listfiles =listfiles[self._start:]
            elif self._start != None and self._end != None \
                and type(self.slicing) != list :
                    assert np.array(self._start).min() < \
                np.array(self._end).max()<= len(listfiles)-1," Index {0} is out of the range."\
                    " Maximum index is {1}".format(np.array(self._end).max(),len(listfiles)-1)

                    listfiles=listfiles[self._start:self._end]
                    
        if self._option ==1:

            #assert slicing  value
            max_temp_start,max_temp_end=[],[]
            for ii in list(self._start) :
                if ii==None :
                    continue 
                max_temp_start.append(ii)
            for ii in list(self._end):
                if ii==None :
                    continue 
                max_temp_end.append(ii)
            
            assert np.array(max_temp_start).min() < \
                np.array(max_temp_end).max()<= len(listfiles)-1," Index {0} is out of the range."\
                    " Maximum index is {1}".format(np.array(max_temp_end).max(),len(listfiles)-1)
               
            tempfiles=listfiles[self._start[0]:self._end[0]]  
            for ii,val in enumerate(self._end) :
                if ii==0:
                    pass
                if ii !=0 :
                    if val == None :
                        assert (ii==len(self._start)-1) or \
                            (ii==len(self._end)-1), "Wrong mode of slicing !'None'"\
                                " nust be the first index or the last index of slicing"
                        temp= listfiles[self._start[ii]:]

                        for ss in temp :
                            tempfiles.append(ss)
                        temp=[]
                    else :
                        temp=listfiles[self._start[ii]:self._end[ii]]
                        for ss in temp :
                            tempfiles.append(ss)
                        temp=[]

            listfiles=tempfiles
         
        for ii, val in enumerate(listfiles) :
            self._comp=self._comp+1         # Avoid dictionnary reinitialization 
            if self.separate==None: 
                line=val.split()
            else : 
                line=val.split(self.separate)
            head=["index_{:02}".format(ind) for ind in range( len(line))]
    
            # extract multiple value from index 
            if self._flag==1 :
                assert 0 <= self.max_index < len(head),\
                    "Argument[ex_index={0:}]Input index {1:^3} is out of the range."\
                        "|Maximum index :{2:^3}|".\
                        format(tuple(self.indexl), self.max_index,len(head)-1)
                
                if self._comp==0:
                    headtemp=[head[ss] for ss in self.indexl]
                    none=[[]for ii in range(len(headtemp))]
                    dico={headtemp:none for headtemp,none in zip (headtemp,none)}
                for idx in self.indexl:
                    for key in dico: 
                        if 'index_{:02}'.format(idx) in key:
                            dico[key].append(line[idx])
                    if idx < 0 and (abs(idx)<=len(head)):
                        reverse=head[idx]
                        for key in dico :
                            if reverse in key :
                                dico[key].append(line[idx])
                                
            # extract single value from  index               
            if self._flag==0:
                if self.indexl<0 :
                    assert 0> self.indexl >= -len(head),"Index {0:^4} out of the range."\
                        " Minimum index egal to {1:^3}.".format(self.data_index,-len(head))
                        
                    reverse=head[self.indexl]
                    if self._comp==0:
                        dico={head[self.indexl]:[]}
                    for key in dico:
                        if reverse in key :
                            dico[key].append(line[self.indexl])
                if self.indexl>=0:
                    assert 0<= self.indexl < len(head),"Index {0:^4} out of the range."\
                        " Maximum index egal to {1:^3}.".format(self.data_index,len(head)-1)
    
                    for idx, car in enumerate( line) :
                        if idx==self.indexl:
                            self._emp.append(car)
                            dico={head[idx]:self._emp}
        return dico
    
    def assert_index_value(self, value):
        """
        Method to assert index value
        
        parameter:
        ----------
            * value : int, list , tuple, float or string 

        Returns
        -------
        value of index 
        
        Example :
        --------
            ff=ReadFile(file_to_read="K1_exp.bln", 
            ex_index=(8,-1),
            slicing=((2,15),(10,None)))
            np.random.seed(0)
            su=[12,[2],(3)]
            ch=[round(15*ii,3) for ii in np.random.rand(15)]
            ch2=tuple(ch)
            ss=ff.assert_index_value(ch2)
            ts= [ff.assert_index_value(ii) for ii in su ]
            print(ch,"\n",ss,'\n',ts)

        """
        # f=-1
        temp=[]
        
        if np.iterable(value)==False:
            
            if value == None:
                value=None
            else :
                value=int(value)
                
        if np.iterable (value):
            typi=type(value)
    
            if None in value :
                f=0
            elif None not in value :
                f=1
                
            if f==0:
                value=list(value)
                for ii in value : 
                    if ii==None :
                        # pass 
                        temp.append(ii)
                    else :
                        temp.append(int(ii))
                    
                if typi==tuple:
                    
                    value=tuple(temp)
                else :
                    value=temp
                    
            elif f==1 :
                value=[int(ii) for ii in value]
                if typi==tuple :
                    value=tuple(value)
        return value
    
    def search_file_(self, Pathfile):
        """
        Method to seach the matching file. 
        
        Returns
        -------
        *select_file * : str or basestring
            Path is None:
                search the list of the file present in the current work directory.
            file is None :
                input the matching file :
                users must choose the data file he want to extract. 
        Ex :
        ---
            ff=ReadFine()
            Path=ff.path
            >>> None 
            ds=ff.search_file(Pathfile=Path)
            >>> ----------------------------------------------------------
                             list of current files                              
                Current directory : F:\OneDrive\Python\
                    CodesExercices\ex_avgfiles\modules
                Number of files found :   20
                filenames :(...)
                Enter the data filename:K1_exp.bln
                -----------------End of searching file.-----------------
        """
        temp=[]
        if Pathfile is None :
            Pathfile=os.getcwd()
            
        listfiles=os.listdir(Pathfile)
        for ss in listfiles :
            if os.path.isfile(ss):
                temp.append(ss)
            else  :
                continue
            
        filenames =temp
        print("-"*80)
        print('{:^80}\n'.format("list of current files"))
        print("Current directory : {:<}".format(Pathfile))
        print("Number of files found : {:>4}".format(len(filenames)))
        print("filenames :")
        
        for ii in filenames :
            print("{:<30}".format(ii),end='') 
        print("-"*80)
        select_file=input("Enter the data filename:")
        print('{:-^80}\n'.format("End of searching file."))
        
        if select_file in filenames :
            return select_file
        else :
            raise FileNotFoundError("No such file in this directory.")
            
            
def interpol_Scipy (x_value, y_value,x_new, kind="linear",
                    fill="extrapolate"):
    
    f=interp1d(x_value, y_value, kind=kind,fill_value=fill)
    y_new=f(x_new)
    
    return y_new
                
#Test : 
if __name__=='__main__':
    dir_=r"C:\Users\kouao\OneDrive\Python\CodesExercices\ex_avgfiles\modules"
    # chdir(dir_)
    # p=ll2DMS(26.04987)
    # print(p)
    # s,east,north=ll_to_utm2(26.051521962922,110.486260190638)
    # print(s,east,north)
    # np.random.seed(0)
    # ser=np.array([ii*152 for ii in np.random.randn(100)])
    # df= pd.DataFrame(data=ser, columns=["TEST"])
    # print(df.shape)
    # print(ss)
    # print(ss%-180)

 