# -*- coding: utf-8 -*-
#       Author: Kouadio K.Laurent<etanoyau@gmail.com>
#       Licence: LGPL

import os
import  numpy as np
import pandas as pd

from pycsamt.utils.ml_utils import neuron
from pycsamt.geodrill.requestmanager import ManageDB
from pycsamt.utils._csamtpylog import csamtpylog 
_logger =csamtpylog().get_csamtpy_logger(__name__)


class ML: 
    """
    Implement ANN to predict layer name from the SQL database. Each database 
    columns is a database properties 
    Each column of Γ  is a database property considered as a vector γ 
    with N_Γ- number of rocks. The column of electrical property of rocks 
    (E_props) is composed of the representative chart of  Palacky (1988)
    and  the rock and mineral property classification of Slichter and 
    Telkes (1942). The E_props column can be considered as a vector
    γ_(E_props ). Moreover, the powerful tool of SQLite  allows to build a 
    sub-table of E_props containing the minimal (v-) and the maximum (v+)
    resistivity values of each rock straightforwardly linked to the main column
    E_props of  Γ.  For instance, according to (Palacky, 1988), the rock of 
    massive sulfides ranges between  0.01 Ω.m to 1 Ω.m, the saprolite  
    between 2 and 200 Ω. m and  according to  (Slichter and Telkes, 1942), 
    the gneiss ranges between 10 000 to 〖10〗^7 Ω.m and so on. 
      https://pubs.usgs.gov/tm/2006/11A02/
      https://pubs.usgs.gov/of/2003/ofr-03-420/ofr-03-420.html

    
    reference: 
        * Slichter, L. B., and M. Telkes, 1942, Electrical Properties of Rocks 
            and Minerals, Geological. (F. B. C. Spicer, J. F. Schairer, 
             and S. H. Cecil, eds.,): © 1942 Geological Society of America.
            https://doi.org/10.1130/SPE36-p299
        * Palacky, G. J., 1988, Resistivity Characteristics of Geologic 
            Targets: Geophysics, 3, 
            https://doi.org/52–129.10.1190/1.9781560802631.ch3 
    .. note:: 
        The `ML`class is separated from the whole package to be an 
        online additional module. It development is still ongoing and 
        it does not affect the software to work properly. Indeed, 
        its development consists to fetching the better parameters from the 
        online database trained using the artificial neural network to 
        predict the layer name according to the given attributes (layer
        porosity, permeability, density,mineral composition, etc.).
        The better way is to keep the parameters updated rather than those 
        computed using only a case study previously done. Huger the data from 
        the different geological areas should be, better should be 
        the parameters for the layer name prediction.
    """
    database_ = os.path.join(os.path.realpath ('.'),
                             'pycsamt/geodrill/_geomemo')
    
    def __init__(self, **kwargs): 
        self._logging = csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        self.DB =None
        self._df = None
        
        for key in list(kwargs.keys()): 
            self.__setattr__(key, kwargs[key])
            
    @property 
    def df(self):
        """Fetch dataframe from memory."""
        self.DB =ManageDB(db_host= self.database_, 
                          db_name = '__m.sq3'
                          )
        self._df = pd.read_sql (
            'SELECT * from GEODB', con = self.DB.connexDB)
        self.DB.closeDB()
        return  self._df
    

    
            
# if __name__=='__main__': 
#     # db_dirsql ='pycsamt/geodrill/_geomemo/memory.sq3'
#     mlObj = ML()
#     print(mlObj.df)
    
