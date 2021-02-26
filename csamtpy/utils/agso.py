# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 20:30:41 2020

@author: KouaoLaurent alias  @Daniel03

Class : 
    **Agso ** . 
            Data of Geological Welllogs 

"""
import os 
# import numpy as np
import warnings
from csamtpy.utils._csamtpylog import csamtpylog
import pandas as pd


_logger=csamtpylog.get_csamtpy_logger(__name__)
# _logger.basicConfig(filename='LogTest')


class Agso (object):
    """
    Read Agso as pandas Series
    Geological conventional rocks and structurals handling. 
    
    staticsmethods :
    ----------------
        * _agso_ * : func 
            return a dataframe from AGS0 file 
        
        * _agso_on_dict_: func
            return a special dictionnay of element of AGSO file 
        
        e.g. :: 
            >>>     import os 
                    path=r'F:\OneDrive\Python\CodesExercices\ex_avgfiles\modules\_geocodes'
                    os.chdir(path)
                    
                    #TEST with AGSO.csv file
                    geo=Agso._agso_on_dict_(set_agsoDataFrame=False, return_orientation="SERIES")
                    
                    # print(agf)
                    # agf_struc=Agso._agso_(configuration_agso_filename=os.path.basename(os.path.join(path,"AGSO_STCODES.csv")))
                    
                    #test with "AGSO_STCODES
                    # geo=Agso._agso_on_dict_(set_agsoDataFrame=True, return_orientation="series", 
                    #                         agso_codefile=os.path.basename(os.path.join(path,"AGSO_STCODES.csv")))
                    
                    #EXTRACT ELEMENT WITH THEIR CORRESPONDANCE 
                    dico_test={key:value for key , value in zip (geo["CODE"],geo['__DESCRIPTION'])}
                    dico_lst=sorted(dico_test.items())
                    # print(dico_lst)
        
    """
    @staticmethod
    def _agso_(configuration_agso_filename =None):
        """
        Returns
        -------
        df_ : object 
            Read pandas series of Agsofiles.

        """
        
        if configuration_agso_filename is None :

            # path =os.path.normpath(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
            #                                     'geodrill','_geocodes','AGSO.csv'))
            path = os.path.join(os.environ ['pyCSAMT'], 'geodrill', '_geocodes', 'AGSO.csv')
            # os.chdir(os.path.dirname(path))
            # fileagso =os.path.basename(path)
            # os.chdir(path)
            df_=pd.read_csv(path)
        elif configuration_agso_filename is not None :
            df_=pd.read_csv(configuration_agso_filename)
        
        return df_
    
    @staticmethod
    def _agso_on_dict_(**kwargs):#single_series_extraction=False,
        """

        Parameters
        ----------
        *set_agsoDataFrame* : Boolean  
                if False it will take the agso default filename 'AGSO.csv'
                .
            
        * single_series_extraction* : Boolean  (deprecatedWarning)
                way to creat a dictionnary 
                Each column match its values  no including index 
                keys are columns 
                if False it will return a dictionnary with index 
                
        * return_orientation * : str 
                specify the way to extract the dictionary 
                
        * agso_codefile : str , 
                agsoFilename input, it will call _agso_ staticmethod.
                and straighforwardly output the dictionnay . No need to call _agso_ methods
                default is AGSO.csv file.

        Returns
        -------
        dictionary of agso elements
        
        
        """
        
        set_agsoDataFrame=kwargs.pop("set_agsoDataFrame", False)
        agso_codefile=kwargs.pop("agso_codefile", None)
        return_orientation=kwargs.pop("return_orientation", None)
        
        orients=["series", "index",
                                 "records", "split"]
        
        # files_contents=[ii for ii in os.listdir(os.path.getcwd) if os.path.isfile(ii)]
        dico={}
        
        if set_agsoDataFrame is False:
            agsoDataFrame=Agso._agso_()
        elif set_agsoDataFrame is True:
            if agso_codefile is not None :
                assert type(agso_codefile)==str , "Wrong Type ! Must be a string format"
                if not agso_codefile.endswith(".csv"):
                    _logger.error("Wrong Format of agsoFile ! Check your  _geocodes Folder !")
                    warnings.warn("Wrong agso Format ! must  choose "\
                                  "file amongs files on _geocode Folder!",category="error")#, \
                                  #format(",".join([ii for ii in os.listdir(os.getcwd()) if os.path.isfile(ii)])))
    
                agsoDataFrame=Agso._agso_(configuration_agso_filename=agso_codefile)
            elif agso_codefile is None :
                raise Exception('Must input at least the name of geocodefile ! '\
                                'or set argument "set_agsoDataFrame" to False'\
                                    ' to select the default agsoFile [AGSO.csv]. ')
                
        df=agsoDataFrame.copy()
        
        if return_orientation is None : # return the defaut orientation 
            dico=df.to_dict()
        else :
            try :
                if return_orientation.lower() in orients :
                    dico=df.to_dict(return_orientation.lower())
                elif return_orientation not in orients :
                    _logger.error("Must input the pandas series return_orientation."\
                                  " Some orientation are {0}".\
                                      format(','.join(["\{0}".format(jj) for jj in orients])))
                    print("Wrong orientation ! List of return_orientations are :[{0}]".\
                          format(','.join(["'{0}'".format(ii) for ii in orients])) )    
            except TypeError :
                pass
        ###################
        # this particular loop below are not necessary :
        # warnings.info(text=, )
        # if single_series_extraction is True :
        #     columnsNames=list(df.columns)
        #     for ii , val in enumerate(columnsNames):
        #         values=df[val]
        #         dico[val]=values
        # else :
            
        return dico
        
     
if __name__=="__main__":

    path=r'F:\OneDrive\Python\CodesExercices\ex_avgfiles\modules\_geocodes'
    os.chdir(path)
    
    #TEST with AGSO.csv file
    geo=Agso._agso_on_dict_(set_agsoDataFrame=False, return_orientation="SERies")
    
    # print(agf)
    # agf_struc=Agso._agso_(configuration_agso_filename=os.path.basename(os.path.join(path,"AGSO_STCODES.csv")))
    
    #test with "AGSO_STCODES
    # geo=Agso._agso_on_dict_(set_agsoDataFrame=True, return_orientation="series", 
    #                         agso_codefile=os.path.basename(os.path.join(path,"AGSO_STCODES.csv")))
    
    #EXTRACT ELEMENT WITH THEIR CORRESPONDANCE 
    dico_test={key:value for key , value in zip (geo["CODE"],geo['__DESCRIPTION'])}
    dico_lst=sorted(dico_test.items())
    # print(dico_lst)
    print(geo
          )
 
        
        
        
        
        
        
        
        
        