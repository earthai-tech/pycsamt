# -*- coding: utf-8 -*-
#       Copyright © 2021  Kouadio K.Laurent
#       Author:  ~Daniel03 <etanoyau@gmail.com>
#       Created:on Wed Oct 14 13:38:13 2020
#       Licence: LGPL
"""
.. _module-GeoDataBase::`pycsamt.geodrill.geodatabase`
 
    :synopsis: Special class to manage outputs-input requests from-into 
                SQL  database .
                ...
.. warnings:: Editing this module presume that you are aware of what you  
              are doing. The module is a core of Geodrill subpackages.
              However the the way the dataBase is arranged can be enhanced  
              and adapted for better convenient or other suitable purposes.
"""
import os
import sys
import numpy as np  
import pandas as pd 
import warnings
import datetime

import  pycsamt.utils.exceptions as CSex 
from pycsamt.geodrill.structural import Geo_formation 
from pycsamt.geodrill._dictapp import Glob
from pycsamt.geodrill.requestmanager import ManageDB
from  pycsamt.utils.decorator import redirect_cls_or_func
from pycsamt.utils._csamtpylog import csamtpylog 
#set LogInfos
try :
    _logger=csamtpylog.get_csamtpy_logger(__name__)
    # _logger.setLevel(logging.DEBUG)
    filename=os.path.basename(__file__)
    # lineno=__file__.__code__.co_firstlineno + 1
except Exception as error :
    _logger.warning("No basic configuration found. Try to set logging.basicConfig()"\
                    "file as .json file or .yml file",error )
    
    warnings.warn_explicit('you may configure logFile', category=DeprecationWarning,
                           filename=filename, lineno=1)
    pass 

# let set the systeme path find memory dataBase
 
sys.path.insert(0, os.path.abspath('.'))  
sys.path.insert(0, os.path.abspath('..')) 
sys.path.insert(0, os.path.abspath('../..'))  # for consistency 
sys.path.insert(0, os.path.abspath('pycsamt/geodrill/_geomemo'))   

# =============================================================================

class GeoDataBase (object): 
    """
    Geodatabase class . Currently we do not create the specific pattern 
    for each geostructures. DataBase is built is built following the   
       codef `code`,   `label`,`__description`,`pattern`, `pat_size`,`pat_density`,
       `pat_thickness`,`RGBA`, `electrical_props`, `hatch`, `colorMPL`, `FGDC` .

    Arguments
    -----------
        **geo_structure_name** : str 
                Name of geological rocks , strata or layer.
                
    .. seealso:: FGDC-Digital cartographic Standard for Geological 
                Map Symbolisation. 
    
    """
    #  FGDC is not set yet , we use the  matplotlib pattern symbol makers 
    make_pattern_symbol =["/", "\\", "|", '-', '+', 'x', 'o', 'O', '.', 
                          '*', '\-', '\+', '\o', '\O', '\.', '\*'] 
    #use '\\' rather than '\'.
    # latter , i will be deprecated to FGDC geological map symbolisation . 
    
    codef = ['code','label','__description','pattern', 'pat_size',	
             'pat_density','pat_thickness','rgb','electrical_props', 
                 'hatch', 'colorMPL', 'FGDC' ]
    # geoDataBase=os.path.join(os.environ['pyCSAMT'],'pycsamt',
    #                'geodrill', 'geoDB','sql_utils', 'sql_DB', 'memory.sq3') 
    # locate the geodataBase
    geoDataBase = os.path.join(
        os.path.abspath('pycsamt/geodrill/_geomemo'),'memory.sq3')
 
    # :memory: is faster we chose this options :geoDataBase.sq3 
    #in sql_DB contains drill holes and wells Tables 
    # will develop in the future extensions 

    def __init__(self, geo_structure_name=None) :
        self._logging =csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        self.geo_structure_name = geo_structure_name
        self.dateTime= datetime.datetime.now().utcnow()   # Get the date time now  
        self.comment =None                                # initialise comment  text
        
  
        self._mplcolor =None 
        self._rgb =None 
        self._electrical_props=None 
        self._pattern=None

        try : 
        # to connect geodataBse 
            self.manage_geoDataBase =ManageDB(
                db_host=os.path.dirname(self.geoDataBase), 
                                      db_name ='memory.sq3')
        except : 
            mess =''.join(['Connection to geoDataBase failed! Sorry we can ',
                           'not give a suitable reply for your request!', 
                    'It would  be process with "geological structural class." !'])
            
            warnings.warn(mess)
            self.success= 0
 
        else : self.success = 1 

    def _avoid_injection (self): 
        """
        For secure , we do not firstly introduce directly the request. We will 
        check whether the object `request`   exists effectively  in our
        GeoDatabase . if not , request will be redirect to structural and 
        strata class issue from  module `structural`
        
        """
        # self.manage_geoDataBase.executeReq(" select __description  from AGS0")
        
        self.geo_structure_name=self.geo_structure_name.lower() # for consistency 
        manage = self.manage_geoDataBase.curs.execute(
            " select __description  from AGS0")
        __description = [ geoform[0] for geoform in  list(manage)]


        if self.geo_structure_name in __description : 
            self.geo_structure_exists =True 
            
        else :
            mess ='Structure <%s> does not exist  in'\
                ' our GeoDataBase yet.' % self.geo_structure_name
            
            self.geo_structure_exists =False
            warnings.warn(mess)
            self._logging.debug ('Could not find a {0} in our geoDataBase.'\
                                 ' It would be redirect to _strata and '
                                 '_structural classes for suitable processing.'.
                                 format(self.geo_structure_name))
            
            # self.manage_geoDataBase.closeDB() # close the database
            
        # if self.geo_structure_exists : 
        #     self._get_geo_structure()
        
    def _retreive_databasecolumns (self, columns): 
        """ Retreive data from database columns
        
        :param columns: Columns name is `str`. To retreive data of many columns 
                      please put the columns name on list.
        :returns: 
            list of data of each columns.
        
        :Exemple:
            >>> from pycsamt.geodrill.geodatabase import GeoDataBase 
            >>> dbObj = GeoDataBase()
            >>>  values = dbObj._retreive_databasecolumns(
                    ['__description', 'electrical_props'])    
        
        """
        
        if isinstance(columns, str): 
            columns =[columns]
            
        new_columns = []
        for obj in columns : 
            if obj =='name'or obj.find('name')>=0: 
                obj = '__description'
            if obj not in self.codef :
                self._logging.debug(
                    f'Object `{obj}` not found in {self.codef}!'
                    'Please provide the right column name.')
                warnings.warn(f'Object `{obj}` not found in {self.codef}!'
                    'Please provide the right column name.')
            else:
                new_columns.append(obj)
            
        if len(new_columns) ==0 : # Object not found in database 
            self._logging.error('None object found in the database !')
            return 
        
        _l=[]  
        
        for obj in new_columns : 
            manage = self.manage_geoDataBase.curs.execute(
                " select %s  from AGS0"% obj)
            valuesdb = [ geoform[0] for geoform in  list(manage)]
            if len(new_columns)==1: 
                _l= valuesdb
            else: _l.append(valuesdb)
            
        self.manage_geoDataBase.closeDB() # close the database
        return _l
            
        
    def _get_geo_structure(self, structure_name=None):
        """
        After checking wether the name of structures exists , 
        let find the geoformation properties   from geodatabase  . 
      
        :param struture_name:   name of geological rock or layer 
        :type struture_name: str
        
        """
        if structure_name is not None :
            self.geo_structure_name = structure_name.lower() 
        
        if self.geo_structure_name is None : 
            warnings.warn('No name is inputted as geological formation. Sorry ,'
                          ' your request is aborted !')
            raise CSex.pyCSAMTError_SQL_manager(
                'No name is given as geological formation. Sorry ,'\
                 ' your request is aborted !')
                
        if self.geo_structure_name is not None :
            self._avoid_injection()
            
        if self.geo_structure_exists : 
            __geonames = list(self.manage_geoDataBase.curs.execute(
                "Select * from AGS0 where __description = '{0}'".\
                  format(self.geo_structure_name.lower())))[0]  
     
            # once value is get then set attribute 
            #for each column of geoDataBase
            for ii, codec in enumerate(self.codef) :
                if self.codef [10] == codec : 
                    self.colorMPL= __geonames[ii]
                    # if '#' in __geonames[ii] : # we assume that value
                    #colorMPL is hexadecimal value eg : #0000ff
                    self.__setattr__(codec, self.colorMPL) # color MPL format
                    # to  string tuple  like '(1.0, 0.5, 0.23)'
                    # else :# value in RGB color (0-1 --> 0 to 255 bits)

                elif self.codef [8]== codec : # set electrical properties
                # (min,max) line (1e-5,5.2e0 )
                    self.electrical_props = __geonames[ii] 
                    self.__setattr__(codec , self.electrical_props) # get 
                    #color matplotlib from property attributes

                else :  self.__setattr__(codec , __geonames[ii])
                
    def _reminder_geo_recorder(self, geo_structure_name ): 
        """
        To have reminder of geological formation into the geodatabase ,
        this method allow to output information 
        if the structure does not exist, An error will occurs. 
        
        :param geo_structure_name: name of geological formation 
        :type geo_structure_name: str 
        
        """
        mess ='--->"{}" Querry  successfully executed ! '  
        
        if geo_structure_name is not None :
            self.geo_structure_name = geo_structure_name.lower() 
        
        if self.geoDataBase is not None : 
            self._avoid_injection() 
        if self.geo_structure_exists is True : 
                # keep only the tuple values 
            __geonames = list(self.manage_geoDataBase.curs.execute(
                "Select * from AGS0 where __description = '{0}'".\
                          format(self.geo_structure_name)))[0] 
            

            print(mess.format(self.geo_structure_name))
            print(__geonames)
            
        
        self.manage_geoDataBase.closeDB() # close the database 
        
                
    def _update_geo_structure (self, geo_formation_name =None, **kws): 
        """
        Update _indormation into geoDataBase . 
        
        Remember that the geodatabase is build following this table
        codef 'code','label','__description','pattern', 'pat_size',
        'pat_density','pat_thickness','rgb','electrical_props',  'hatch',
        'colorMPL', 'FGDC'. 
        
        :param geo_formation_name:  name of formation be sure the formation
                             already exists in the geoDataBase 
                             if not an error occurs 
        :type geo_formation_name: str
        
        - Update the electrical property of basement rocks = [1e99, 1e6 ]
        
        :Example:
            
            >>> from pycsamt.geodrill.geodatabase import GeoDataBase 
            >>> GeoDataBase()._update_geo_structure(
                **{'__description':'basement rocks', 
                    'electrical_props':[1e99, 1e6 ]})
        """

        # find geological rocks name if is in keywords dict
        if geo_formation_name is None : 
            if '__description' in list(kws.keys()) : 
                geo_formation_name=str(kws['__description'])
            elif 'name' in list(kws.keys()) : 
                geo_formation_name=str(kws['name'])
            else : 
                raise CSex.pyCSAMTError_SQL_update_geoinformation(
                    ' ! Unable to find a new geological structure name.')
    
        if not isinstance(geo_formation_name, str) : 
            raise CSex.pyCSAMTError_SQL_update_geoinformation(
                'Unacceptable rock/layer name ={0}.'\
                    ' Please provide a right rock/layer name.')
        geo_formation_name=str(geo_formation_name) # for consistency 
        
        if geo_formation_name is not None :
            self.geo_structure_name = geo_formation_name.lower() 
        
        # build new_dictionnary without the keyname and key value of dictionnay 
        tem_geodict ={geokey:geovalue for geokey , geovalue in kws.items() 
                      if not (geokey =='__description' or geokey =='name') 
                      }
        
        if self.geo_structure_name is not None :
            self._avoid_injection()
        if self.geo_structure_exists is False : 
            mess ="".join([
                ' Actually geological formation name =  <{0}> can '.format(
                    self.geo_structure_name),
                 'not be updated because', 
                ' it doesnt  exist in our DataBase. To set new geological ',
                'formation with their corresponding values ', 
                ' see { _add_geo_structure } method .'])
            self._logging.warn(mess)
            raise CSex.pyCSAMTError_SQL_update_geoinformation(
                'Update name= {0}  failed ! it doesnt not exist in '
                'geoDataBase'.format(self.geo_structure_name))
                    
            
        elif self.geo_structure_exists : # even the geostructure exists ,
        #let check wether the key provided 
            for geo_key in list (tem_geodict.keys()) : # is among the geocodes keys
                if geo_key not in self.codef : #if key provided not in geodatable keys 
                    
                    mess =''.join([
                        "Sorry the key = {0} is wrong! key doesnt".format(geo_key),
                        " exist in geoDataBase.Please provide a right ", 
                        " keys among = {0}". format(tuple(self.codef[2:]))])
                    self._logging.error(mess)
                    raise CSex.pyCSAMTError_SQL_manager(mess)
        
                elif geo_key in self.codef :
                    if geo_key.find('pat') >= 0 : 
                        try : # keep value to real
                            update_geo_values = float(kws[geo_key]) 
                        except : 
                            msg =''.join([
                                'update failed ! Could not convert',
                                ' value = {0} to float.'.format(kws[geo_key]), 
                                'Please try again later.'])
                            self._logging.error(msg)
                            raise CSex.pyCSAMTError_SQL_update_geoinformation(msg)
                            
                    # let get the formatage of all values  (properties values )  
                    elif geo_key.find('colorMPL') >=0 : 
                        self.colorMPL = kws[geo_key] 
                        update_geo_values = self.colorMPL
                    # let fill automatically the "rgb" and the colorMPL                
                    elif geo_key.find('rgb')>=0 :  
                        # keep the rgb value g : R125G90B29 and compute the colorMPL
                        self.rgb = kws[geo_key]
                        update_geo_values = self.rgb
                    elif geo_key .find('hatch') >=0 : 
                        self.hatch = kws[geo_key]
                        update_geo_values = self.hatch
                    elif geo_key.find('electrical_props') >=0 : 
                        self.electrical_props =kws[geo_key]
                        update_geo_values = self.electrical_props

                    else : update_geo_values =str (kws[geo_key])
                    # now must be put on the data base
                    # fill automaticall colorMPL when rgb is provided
                    if geo_key.find('rgb') >= 0 :  
                    
                        self.manage_geoDataBase.curs.execute(
                            "update AGS0 set rgb = '{0}'  where __description ='{1}'".
                              format( update_geo_values[0] , 
                                     self.geo_structure_name))
                        self.manage_geoDataBase.curs.execute(
                            "update AGS0 set colorMPL= '{0}'  where __description ='{1}'".
                              format(update_geo_values[1] , 
                                     self.geo_structure_name))
                    
                    else :
                        __oldvalues = list(
                            self.manage_geoDataBase.curs.execute(
                                "Select * from AGS0 where __description ='{0}'".
                                   format(self.geo_structure_name)))[0]
                        self.manage_geoDataBase.curs.execute(
                            "update AGS0 set {0}= '{1}'  where __description ='{2}'".
                               format(geo_key, update_geo_values ,
                                      self.geo_structure_name))
                        __newvalues = list(
                            self.manage_geoDataBase.curs.execute(
                                "Select * from AGS0 where __description ='{0}'".\
                                   format(self.geo_structure_name)))[0]
                    # inpout new info to database from cursor         
                    self.manage_geoDataBase.commit()    
                    
                    if geo_key.find('rgb') >=0 : 
                        print('---> {0} colors was successfully set to '
                              'rgb = {1} & matplotlib rgba = {2}.'.
                              format(self.geo_structure_name,
                                     update_geo_values[0],
                                     update_geo_values[1]))
                    else : 
                        fmt_mess = '---> {0} was successfully set to'\
                            ' geoDataBase.\n ** Old value = {1} \n is '\
                                '**updated to \n New value = {2}'
                        print(fmt_mess.format(
                            self.geo_structure_name,__oldvalues, __newvalues ))
                        
            self.manage_geoDataBase.closeDB() # close the database
            
    @property 
    def hatch (self):
        return self._hatch 
    @hatch.setter 
    def hatch (self, mpl_hatch):
        mm=0  # counter  of hach symbol present on the value provided 
        mpl_hatch =str(mpl_hatch) # for consitency put on string value 
        # removed the tuple sign  "( and )" if provided 
        if '('  in mpl_hatch  : mpl_hatch =mpl_hatch.replace('(', '')
        if   ')' in mpl_hatch: mpl_hatch =mpl_hatch.replace(')', '')
        
        # chech whether the value provided is right 
        if mpl_hatch == 'none' : self._hatch =mpl_hatch 
        else :
            for mpstr  in mpl_hatch : 
                if mpstr in self.make_pattern_symbol  : 
                    mm +=1 
            if len(mpl_hatch) == mm :  #all value are symbols  and put the  
                self._hatch = '(' + mpl_hatch +')'
            else : self._hatch ='none' # abandon value and initialise to None 


    @staticmethod
    def _add_geo_structure( new_geological_rock_name=None , **kws) : 
        """
        Add new _geological information  into geodatabase .
        
        DataBase properties:
        --------------------
            - code
            - label
            - __description
            - pattern 
            - pat_size	
            - pat_density
            - pat_thickness
            - rgb'
            - electrical_props
            - hatch
            - colorMPL
            - FGDC 
            
        .. note:: `__description` could be replaced by `name`.
                    `code` , `label` and `FGDC` dont need to be fill. 
                    Values are rejected if given.
    
        :param new_geological_rock_name: new name of geological formation to add 
        :type new_geological_rock_name: str 
        
        :param informations: 
            dict , must be on keyward keys  when keywords keys 
                are provided , program will check whether all keys are
                effectively the right keys. if not will aborted the process.
        :type informations: dict
        
        :Example: 
            
            >>> from pycsamt.geodrill.geodatabase import GeoDataBase 
            >>> geodatabase_obj= GeoDataBase._add_geo_structure( **{
            ...                                     'name': 'massive sulfure', 
            ...                                     'pattern': 218., 
            ...                                     'pat_size': 250., 
            ...                                     'pat_density': 0.75, 
            ...                                     'pat_thickness': 2., 
            ...                                     'rgb': 'R128B28',
            ...                                     'hatch': '+.+.o.+', 
            ...                                     'electrical_props':[1e0 , 1e-2],
            ...                                     } )
        """
        
        def __generate_structure_code (__description , __geocodeList) : 
            """
            Will geological description will create a code and label 

            :param __description: name of geological formation 
            :type __description: str 
  
            :returns: geological formation code 
            :rtype: str 
            
            """
            def _rev_func_code (code, CODE): 
                """
                generate code and check thin new code  doesnt not exist in 
                amongs the code of geoDataBase .

                :param code:  new_generate code 
                :type code: str 
                :param CODE: codes already exists in dataBase 
                :type CODE: str 
                
                """
                # actually the lencode is > than 3 
                mm=0
                while code  in CODE :
                    if code not in CODE : 
                        break 
                    if mm > len(code): 
                        mm=0
                    code = code + code[mm]
                    mm=mm+1
            
                return code
                
            # fisrtly code is marked by three ,main letters 
            if len(__description) == 3 :code=  __description.upper() 
            elif len(__description) > 3 : code =__description[:4].upper() 
            if len(__description) < 3 :
                nadd = 0
                code =__description # loop thin you find a code =3 
                while nadd < 2 :
                    if len(code )== 3 : 
                        break 
                    code += code[nadd]
                    nadd +=1 
                code =code.upper()
            # then check whether the new code generate exist or not.
            
            for cof in __geocodeList  : 
                if code  not in __geocodeList : return code 
                if code in __geocodeList : 
                    code =_rev_func_code(code =code , CODE=__geocodeList)  

            return code             # return new code that not exist in geocodes 
            
        _logger.info (
            'Starting process  new geological information into GeoDatabase')
        
        # find geological rocks name if is in keywords dict
        if new_geological_rock_name is None : 
            if '__description' in list(kws.keys()) : 
                new_geological_rock_name=str(kws['__description'])
            elif 'name' in list(kws.keys()) : 
                new_geological_rock_name=str(kws['name'])
            else : 
                raise CSex.pyCSAMTError_SQL_update_geoinformation(
                    ' ! Unable to find a new geo_logical structure name.')
    
        if not isinstance(new_geological_rock_name, str) : 
            raise CSex.pyCSAMTError_SQL_update_geoinformation(
                'Unacceptable rocks names ={0}.'
                 ' Please provide a right rock name.')
        new_geological_rock_name=str(new_geological_rock_name) # 
        
        
        # ---------------------------call Geodatabse --------------------------
        
        geoDataBase_obj =  GeoDataBase(new_geological_rock_name)
        
        # initialise to 'none' value 
        mmgeo={geokey : 'none' for geokey in geoDataBase_obj.codef } 
        
        if geoDataBase_obj.success ==1 :  geoDataBase_obj._avoid_injection()
        elif geoDataBase_obj.success  ==0:
            mess = "Connection to SQL geoDataBase failed ! Try again later." 
            warnings.warn(mess)
            raise CSex.pyCSAMTError_SQL_geoDataBase(mess)
            
        #----------------------------------------------------------------------
       
        if geoDataBase_obj.geo_structure_exists : 
            mess ='! Name {0} already exists in our GeoDataBase. Could not add '\
                'again as new geostructure. Use "_update_geo_geostructure method"'\
                    ' to update infos if you need !'.format(
                        geoDataBase_obj.geo_structure_name)
            warnings.warn(mess)
            _logger.error(mess)
            raise CSex.pyCSAMTError_SQL_update_geoinformation(mess)
            
        if geoDataBase_obj.geo_structure_exists is False : 
  
            # make an copy of codef useful in the case where 
            # user provided "name" as key instead of "__description" 
            import copy 
            new_codef = copy.deepcopy(geoDataBase_obj.codef)
            #  set the first value of keys to fill 
            mmgeo ['__description'] = str(new_geological_rock_name) 
            if 'name' in list(kws.keys()) : 
                 new_codef[2]= 'name'
            
            # get the list of geo_code values in DataBase 
            geoDataBase_obj.manage_geoDataBase.curs.execute(
                'Select code  from AGS0 ')
            
            __geocode =[codegeo[0] for codegeo in  
                        list(geoDataBase_obj.manage_geoDataBase.curs)]
            
            
            for key in list(kws.keys()) : 

                if key not in new_codef  : 
                    raise CSex.pyCSAMTError_SQL_update_geoinformation(
                        'Process aborted ! wrong <{0}> key!'
                         ' could not add new informations. '
                         'Please check your key !'.format(key))
                #  set code and labels from geo_description name 
                elif key  in new_codef:  
                
                    try : # check if name is provided intead of  __description 
                    # name (generaly code and label are the same)
                        mmgeo['code'] = __generate_structure_code (
                            new_geological_rock_name, __geocode)
                        mmgeo['label']= __generate_structure_code (
                            new_geological_rock_name,  __geocode)
 
                    except : # user can provide name instead of __description 
                        mmgeo['code'] = __generate_structure_code (
                            new_geological_rock_name, __geocode)
                        mmgeo['label']= __generate_structure_code (
                            new_geological_rock_name, __geocode )
 
                    
                    if key.find('pat')>= 0 :
                        geoDataBase_obj.pattern = kws[key]
                
                        for kvalue  in ['pattern', 'pat_size',
                                        't_density','pat_thickness']: 
                            if key == kvalue : 
                                mmgeo[kvalue]= kws[key]
                                
                    # set RGB value and MPL colors eg : R128G128B --.(0.50, 0.5, 1.0)
                    if key =='rgb' :
                        geoDataBase_obj.rgb= kws[key]
                        # set at the same time rgb value and color MPL 
                        mmgeo['rgb']=  geoDataBase_obj.rgb[0] 
                        mmgeo['colorMPL']=  geoDataBase_obj.rgb[1]
                        
                    # set Matplotlib color whether the rgb is not provided .
                    # if provided will skip 
                    if key=='colorMPL': 
                        if mmgeo['colorMPL'] =='none' : 
                            geoDataBase_obj.colorMPL= kws[key]
                            mmgeo['colorMPL']= geoDataBase_obj.colorMPL
                    # optional keys 
                    
                    if key == 'electrical_props' : 
                        geoDataBase_obj.electrical_props= kws[key]
                        mmgeo['electrical_props']= geoDataBase_obj.electrical_props
                    if key == 'hatch': 
                        
                        mmgeo['hatch']= str(kws[key])
                    if key == 'FGDC': 
                        mmgeo['FGDC']= str(kws[key])
                        
            # print(geoDataBase_obj.success)
            # now build insert all info 
            mm_sql=[]
            for codk in new_codef:  # build info in order and input to GeodataBase
                if codk == 'name' : codk = '__description'
                mm_sql.append(mmgeo[codk])

        
            reqSQL = 'insert into AGS0 ({0}) values ({1})'.format(
                ','.join(['{0}'.format(key) for key in geoDataBase_obj.codef ]),
                ','.join(['?' for ii in range(len(geoDataBase_obj.codef))]))
 
            try : 
                # geoDataBase_obj.manage_geoDataBase.curs.execute(reqSQL , mm_sql )
                geoDataBase_obj.manage_geoDataBase.curs.execute(
                    reqSQL , mm_sql )
                
            except : 
                mess='Process to set {0} infos failed!  Try again '\
                    'later! '.format(new_geological_rock_name)
                warnings.warn (mess)
                _logger.error(mess)
                
                raise CSex.pyCSAMTError_SQL_update_geoinformation(mess)
                
            else :
   
                geoDataBase_obj.manage_geoDataBase.commit()     
                print('---> new data ={} was successfully set into'
                      ' GeoataBase ! '.format(new_geological_rock_name)) 
            
        geoDataBase_obj.manage_geoDataBase.closeDB() # close the database 
        
                  
    @property 
    def pattern (self):
        "return geopattern"
        return self._pattern 
        
    @pattern.setter 
    def pattern (self, pattern_value):
        "configure geopattern"
        try : 
             float(pattern_value)
        except : 
            mes ='Process aborted ! Could not convert'\
                f' {pattern_value} to float number.'
            self._logging.warning(mes)
            raise CSex.pyCSAMTError_SQL_update_geoinformation(mes)

        else : 
            self._pattern = float(pattern_value)
    
    @property 
    def colorMPL(self): 
        "return geocolorMPL"
        return self._mplcolor 
    
    @colorMPL.setter 
    def colorMPL (self, mpl_color): 
        """
        configure geocolorMPL
        to set matplotlib _color in rgb value 
        value is range (0 to 1) coding to 0 to 255 bits.
        """
 
        if isinstance(mpl_color, str): # get value from database 
   
            if mpl_color.find('(')>=0 and mpl_color.find(')')>=0 :
                # build the tuple of mpl colors 
                self._mplcolor = tuple([ float(ss) for ss in 
                                         mpl_color.replace(
                                             '(', '').replace(')',
                                                              '').split(',')])
            # we assume that value colorMPL is hexadecimal value eg : #0000ff
            elif '#' in mpl_color : 
                 # assume the colorMPL is in hexadecimal 
                self._mplcolor =  str(mpl_color).lower() 
            elif 'none' in mpl_color : # initilisation  value 
                self._mplcolor =  'none' # keep the value on the dataBase 
            else : 
                import matplotlib as mpl 
                try : # try to convert color to rgba
                
                    self._mplcolor = mpl.colors.to_rgb(str ( mpl_color))
                except : 
                    raise  CSex.pyCSAMTError_SQL_manager(
                        ' Unsupported {0} color!'.format(mpl_color))  
                else :  # keep only R, G, B and abandon alpha . 
                        #Matplotlib give tuple of 4 values 
                        # as (R, G, B, alpha)
                    self._mplcolor =self._mplcolor[:3] 
                    
             # set value to database way 
        elif  isinstance(mpl_color, (list, tuple, np.ndarray)): 
            if 3 <len(mpl_color) < 3 : 
                msg =''.join(['update failed ! value = {0} '.format(mpl_color),
                              'must be a tuple of 3 values= (Red, Green, Blue)',
                             'values. Please provided a right number', 
                             '  again later.'])
                
                self._logging.error(msg)
                raise CSex.pyCSAMTError_SQL_update_geoinformation(msg)
            # let check whether the value provided can be converted to float    
            if len(mpl_color)==3 :  
                try : 
                     self._mplcolor= tuple( [float(ss) for
                                             ss in list( mpl_color)])
                except : 
                     msg =''.join(['update failed ! Could not convert value ',
                                   '= {0} to float.'.format(mpl_color), 
                                   'Please try again later.'])
                     self._logging.error(msg)
                     raise CSex.pyCSAMTError_SQL_update_geoinformation(msg)
                else : 
                    # try to check if value is under 1.
                    # because color is encoding to 1 to 255 bits 
                    for ival in  self._mplcolor: 
                        if 1 < ival <0  : 
                            if ival > 1 : fmt ='greater than 1' 
                            elif  ival <0 :
                                fmt= 'less than 0'
                            msg = ''.join([
                                'update failed ! Value provided  =',
                                f' `{ival}` is UNacceptable value ! Input ',
                                f' value is {fmt}. It must be encoding from ',
                                '1 to 255 bits as MPL colors.'])
                            raise CSex.pyCSAMTError_SQL_update_geoinformation(msg)
                            
                self._mplcolor=str( self._mplcolor) # put on str for consistency 
                
    @property 
    def rgb(self):
        "return georgb"
        return self._rgb 
    @rgb.setter 
    def rgb(self, litteral_rgb): 
        """
        configure georgb
        
        .. note:: return the rgb value and the convert rgb palette value: 
                keep the rgb value eg `R125G90B29` and compute the colorMPL
                let fill automatically the "rgb" and the colorMPL 
        """
        from pycsamt.geodrill.geoCore.structural import get_color_palette

        self._rgb=(litteral_rgb, str(
            get_color_palette(RGB_color_palette=litteral_rgb)))
        
    @property 
    def electrical_props(self): 
        "return electrical property"
        return self._electrical_props
    
    @electrical_props.setter 
    def electrical_props(self,range_of_rocks_resvalues):
        """
        configure electrical property
        
        .. note:: Electrical_property of rocks must a tuple of resisvity ,
                 max and min bounds  eg : [2.36e-13, 2.36e-3]
        """
        # electrical props were initialised by float 0. 
        if isinstance(range_of_rocks_resvalues , str) : 
            if '(' in range_of_rocks_resvalues  : 
                self._electrical_props = tuple([ 
                    float(ss) for ss in # build the tuple of mpl colors 
                    range_of_rocks_resvalues .replace('(',
                                                '').replace(')',
                                                            '').split(',')]) 
            elif 'none' in range_of_rocks_resvalues : 
                self._electrical_props =.0
                
        elif isinstance(range_of_rocks_resvalues,(list,tuple, np.ndarray)): 
            if len(range_of_rocks_resvalues) ==2  : 
                try : 
                    self._electrical_props =[float(res) 
                                       for res in range_of_rocks_resvalues]
                except : 
                    raise CSex.pyCSAMTError_SQL_update_geoinformation(
                        ' !Could not convert input values to float.')
                else :# range the values to min to max 
                    self._electrical_props =sorted(self._electrical_props) 
                    self._electrical_props = str(tuple(self._electrical_props))
            else : 
                # force program to format value to 0. float 
                self._electrical_props = .0   
                    
        elif not isinstance(range_of_rocks_resvalues,
                            (list,tuple, np.ndarray)) or len(
                                range_of_rocks_resvalues)!=2:
            try : # 0 at initialization
                range_of_rocks_resvalues = float(range_of_rocks_resvalues)  
            except  : 
                if len(range_of_rocks_resvalues) > 1: fmt ='are'
                else :fmt ='is'
    
                mess = ''.join([
                    'Unable to set electrical property of rocks.', 
                    ' We need only minimum and maximum resistivities',
                    ' bounds. {0} {1} given'.format(
                        len(range_of_rocks_resvalues), fmt)])
                        
                self._logging.error(mess)
                warnings.warn(mess)
                raise CSex.pyCSAMTError_SQL_update_geoinformation(mess)
            else : 
                
                self._electrical_props = .0 #mean value  initialised 
                

    @property
    def _setGeoDatabase(self): 
        """
        .. note:: property of  GeoDataBase -create the GeoDataBase
                 Setting geoDataBase table
                 No Need to reset the DataBase at least you dropped the table,
                 avoid to do that if you are not sure of what you are doing. 
        """
        import  pycsamt.utils.func_utils as func 
        #set other attributes
        # create connection 
        # try : 
        #     # call manage DB  so to connect geodataBse 
        #     manage_geoDataBase =ManageDB(db_host=os.path.dirname(self.geoDataBase), 
        #                                  db_name ='memory.sq3')
        # except : pass 
        # if success is True  : 
        # createa tABLE IN THE DATABASE
        req  = ''.join(['create table AGS0', 
                      ' (', 
                      '{0} TEXT,', 
                       ' {1} TEXT,', 
                       ' {2} TEXT,', 
                       ' {3} REAL,', 
                       ' {4} REAL,', 
                       ' {5} REAL,', 
                       ' {6} REAL,', 
                       ' {7} TEXT,', 
                       ' {8} REAL,', 
                       ' {9} TEXT,', 
                       ' {10} TEXT,', 
                       ' {11} TEXT', 
                       ')'
                       ]) 
        
        # create Request to insert value into table 
        mes ='insert into AGS0 ('+ ','.join(['{}'.format(icode)
                                             for icode in self.codef])
        enter_req = mes+ ')'
        

        # create  Geoformation objets 
        geo_formation_obj =Geo_formation()      # set oBject of geostructures 
        
        new_codef ,  new_codef[-1], new_codef[0]= self.codef [:8], 'color', 'codes'
        # get attribute from geoformation and #build new geo_formations _codes 
        geo_form_codes =[ getattr(geo_formation_obj, codehead) 
                         for codehead in new_codef ]
        geo_form_codes =func.concat_array_from_list(geo_form_codes,
                                                    concat_axis=1)
        # generate other main columns of tables to fill laters 
        geo_add_form_codes = func.concat_array_from_list(
            list_of_array = [ np.zeros((len(geo_formation_obj.codes),)), 
                            np.full((len(geo_formation_obj.codes),), 'none'), 
                            np.array([str (clsp) 
                                      for clsp in geo_formation_obj.mpl_colorsp]),
                            np.full((len(geo_formation_obj.codes),),'none'),
                                                                      ],
                                                         concat_axis=1 )
        # create Table resquest 
        
        req = req.format(*self.codef)              # generate a request for table creation 
        enter_req = enter_req.format(*self.codef)  # generate interrequest 

        GDB_DATA = np.concatenate((geo_form_codes,geo_add_form_codes), axis =1)
        
        # generate values Host string so to avoid injection 
        # print(req)
        values_str = 'values (' + ','.join(['?' 
                                for itg in range( GDB_DATA.shape[1])]) +')'
        insert_request = ''.join([enter_req , values_str])

        # create Table
        
        try : 
            self.manage_geoDataBase.executeReq(query=req ) 
            
        except : 
            warnings.warn('Could not create {AGS0} Table !')
            raise CSex.pyCSAMTError_SQL_geoDataBase(
                'Table AGS0 already exists !')

        if self.manage_geoDataBase.success ==1: 
            # enter the record 
            for ii, row_geoDataBase in enumerate(GDB_DATA ): 
                row_geoDataBase =tuple(row_geoDataBase)
                self.manage_geoDataBase.executeReq(
                    query=insert_request , param =row_geoDataBase ) 
        
        self.manage_geoDataBase.commit()
        self.manage_geoDataBase.closeDB()


    
class Recorder_sql(object):
    """
    Class to record data from file or pd.core.DataFrame and to tranfer 
    into SQL database. 
    
    Arguments 
    -----------
        **database** : str , 
                name of sql database 
                
        **table** : str , 
                name of table in dict_app
    
        **Glob.dicoT** : dict,
                dicoT is from dict_app module , Global class for sql variable 
                encapsulated on particular dictionnary.
                
    =========================  ==============================================
    Methods                    Description 
    =========================  ==============================================
    transferdata_to_sqlDB      transfer Data to SQL Database 
    keepDataInfos              keep informations from datafile '*csv'
    arrangeData_for_dictapp    to arrange data, acording the dict_app 
                               arangement
    =========================  ==============================================

    .. deprecated::deprecated methods 
                - recordData --> to  keepDataInfos (staticmethod)
                - set_on_dict_app --> to tarrangeData_for_dictapp(staticmethod)
    
    :Example: 
        
        >>> from pycsamt.pycsamt.geodrill.geodatabase import Recorder_sql
        >>> realpath=os.path.dirname(os.path.realpath(__file__)) 
        >>> #where 'the file'sql_recorder is located'
        >>> print(realpath)
        >>> path_to_files =os.path.normpath(os.path.join(realpath,'sql_utils',
        ...                                                 'sql_DB'))
        >>> os.chdir(path_to_files)
        >>>  filename='nofacies_data_2.csv'
        >>> memory_DB='memory.sq3'
        >>> path_to_memory=os.path.dirname(os.path.realpath(memory_DB))
        >>> Rec=Recorder_sql(database=memory_DB,table=None)
        >>> recordList1=Recorder_sql.recordData(
            data=filename,new_tablename='TEST2',sep=',')
        >>> sdico_app=Recorder_sql.set_on_dict_app( datalist=recordList0)
        >>>  sdico_app=Recorder_sql.arrangeData_for_dictapp( 
            datalist=recordList0)
        >>> trand_sql= Rec.transferdata_to_sqlDB( 
            filename=filename,record_list= None,
        ...                 table_name='exam_zju_2',comments=None,
        ...                 visualize_table_creating_query=True,
        ...                 path_to_sqlDataBase=path_to_files,
        ...                 Drop_DB_Tables='none', fetchall=True, 
        ...                 ready_to_transfer='n',
        ...                 close_connexion =False)
    """
    
    def __init__(self, database , table=None,**kwargs):
        
        self.sqlmanager=database
        self.table=table
        self.GlobDicoT=Glob.dicoT
        
        # self.descriptif =Glob.dicoT[table]
        
        
    def transferdata_to_sqlDB(self, record_list=None,
                              filename=None, table_name=None, **kwargs):
        """
        Function to transfer Data from Dict_app to SQL DataBase. Users 
        can use this function byincluding several arguments.
        The function will build the data , arrange  it and put it in 
        the dataBase by commit the dataBase. Use only this func is benefit.
        It is better to revise arguments of that function.
        
        Parameters
        ----------
            * record_list : dict, optional
                    Dictionnay build according the dict_app model. 
                    The *default* is None.
                
            *  filename : str, optional
                    file must be on ".csv" format. The default is None.
            *  table_name : str, optional
                    Name of DataBase Table. The default is None.
            
            *  comments : str  
                    little comment to identify your database table. 
                    
            *  path_to_sqlDataBase : str  
                    path where the SQL dataBase is located . 
                
            *  visualize_table_creating_query : bool 
                   If the connexion to server is unlikable  
                   set to True to see whether query entered is right or wrong.
                   
            *  Drop_DB_Tables : str, 
                   way to drop table in SQL Database . set litteral arguments like 
                  the name of database user want to drop or 
                  [ no "*" or all to drop all tables. 
                                                             
            *  Ready_to_transfer : str 
                   process to commit Database , the curso tranfered
                   the dataBase to SQL connexion.
                   set litteral 'no' or 'yes' to do.
              
            *  close_connexion : bool , 
                  set True when transfer is done . it seems connexion.close()

        Raises
        ------
            Exception occurs when Table Name is not set on dict_app.

        .. note:: The process of organization is full request of  PostgreSQL 
        
        """
        
        comments=kwargs.pop('comments', None)
        sql_DB_loc =kwargs.pop('path_to_sqlDataBase',None)
        visualize_req=kwargs.pop('visualize_table_creating_query',False)
        drop_table=kwargs.pop('Drop_DB_Tables','none')
        fetch_data=kwargs.pop('fetchall', False)
        ready_to_transfer=kwargs.pop('Ready_to_transfer', 'n')
        close_connex=kwargs.pop('close_connexion',False)
        
        if record_list == None :
            record_list=Recorder_sql.keepDataInfos(
                data=filename,new_tablename= table_name)
            
        sdico_app=Recorder_sql.arrangeData_for_dictapp(
            datalist=record_list, comments=comments)
        
        for keys , values  in sdico_app.items():
            # or keys=[keys for keys in sdico_app]--> keys =keys[0]
            #or keys = [*sdico_app]---> keys=keys[0]
            valapp,insertvalues=values
            self.GlobDicoT.__setitem__(keys, valapp)
        # name_of_sql_DATABASE=self.sqlmanager
        manDB=ManageDB(db_name=self.sqlmanager, 
                               db_host=sql_DB_loc)
        table_req=manDB.dicT_sqlDB(dicTables=self.GlobDicoT,
                         visualize_request=visualize_req)
        
        #build request for values : 

        if self.table is None : 
            self.table =table_name
            
        if self.table not in list(self.GlobDicoT.keys()):
            raise  f'No such Table <{self.table}> '\
                'in the dataBase{self.sqlmanager} ! '\
                'You may enter the right Table Name.'
            
        self.descriptif =self.GlobDicoT[self.table]
        #"INSERT INTO oeuvres(comp,  titre,  duree, interpr) VALUES(?,?,?,?)"
        value_req=table_req.replace('CREATE TABLE', 'INSERT INTO')
        
        beacons="("         # balises= "beacons" in english
        # fieldValue=[]
        # for ii, rowline in enumerate(self.descriptif):
        #     fieldValue.append(rowline[5])   # rowline of  value in dicoT.
        
        len_values =len(insertvalues[0])
        beacons = beacons+'?,'*len_values 
        beacons=beacons [:-1] +")"   # sustract the ',' and close the parenthesis 
        value_req= value_req + ' VALUES %s' % beacons
        
        # for itemvalues in insertvalues: 
        manDB.executeReq( query =value_req, param=insertvalues)
        # manDB.executeReq( query =value_req, param=None)        
        if fetch_data is True :
            manDB.print_query()
        

        if drop_table ==False or drop_table.lower() in [
                'none',"no", 'not yet','wait','not ready']:
            pass
        elif drop_table ==True or drop_table.lower() in [
                'total', 'all', 'all tables','*']:
            
            manDB.drop_TableDB( dicTables=self.GlobDicoT,drop_all=True)
            
        elif drop_table ==True and drop_table.lower() in list[
                self.GlobDicoT.keys()]:
            
            manDB.drop_TableDB( dicTables=self.GlobDicoT,
                               drop_table_name=drop_table )
        else :
            pass 
        
        if ready_to_transfer ==True or ready_to_transfer.lower() in [
                'y','alright,''yes','ok','fine',"right"] :
            manDB.commit()
                     
        if close_connex is True :
            manDB.closeDB()
            
    @staticmethod
    def keepDataInfos(data, new_tablename=None , **kwargs ):
        """
        Function to KeepData from file infos . the function 
        is otherwritten fo RecordData. The difference between two function 
        is that function organise data from each row of columns 
        
        Parameters
        ----------
            * data : str, np.array, list, or pd.core.DataFrame object
                    Data ca, be on the format above or filename of data 
                    if the argument "data" is a filename, we must be convert 
                    on ".csv" format.
                
            * new_tablename : str, optional
                    Name of database. if name is not given , the 
                    function return only list . The default is None.
 
        Raises
        ------
            IndexError.
                 if lengh of number of columns like heads of data 
                does not match the data.shape[0], then errors will occurs.  

        Returns
        -------
            list  
                list of value in the case of no name is providen for tablename.
                else  return dict if name of datatable is providen.
            
        """
        
        nump_columns=kwargs.pop('columns_of array', None)
        separate_colist=kwargs.pop('seperate_columns_list', ',')
        seperate_csvfile=kwargs.pop('sep',',')
        
        tempkeylenvalue,temp,index_item =[],[],[]
        if nump_columns is not None :
            if type (nump_columns) ==np.ndarray :
                nump_columns=nump_columns.tolist()
            elif type (nump_columns) == str :
                nump_columns=nump_columns.split(separate_colist)
            else :
                if separate_colist not in nump_columns:
                    nump_columns=[nump_columns]  # put on a list 
    
        if os.path.isfile(data) and data.endswith('.csv'): # read csv file .
            data=pd.read_csv(data, sep=seperate_csvfile,
                             delimiter=None,header='infer',
                             index_col=False)
        
        if type (data)==pd.core.frame.DataFrame:
            columns=data.columns        # get the dataframe columns
            columns=[*columns]          # unpack pd.Index on list 
            datavalues=data.to_numpy()  # recover numpy array
            
            for ii, rowlines in enumerate( datavalues) : 
                temp.append(tuple(rowlines.tolist()))
                # if ii ==0 : 
                #     for jj, item in rowlines.tolist(): 
                #         type_of_value.append(type(item))
            type_of_value=[type(item) for  item in datavalues[0,:].tolist()] 
            index_item=[columns.index(item)+1 for item in columns]
            num_of_val=[datavalues[:,ii].shape[0]
                        for ii, item in enumerate(columns)]
            tempkeylenvalue=(columns, type_of_value,
                             index_item,num_of_val,temp)
            
            if new_tablename is not None : 
                return {new_tablename:tempkeylenvalue}

            return tempkeylenvalue
            
        elif type(data) ==np.ndarray:
            
            if nump_columns is None :
                raise Exception (
                    'The primary keys is needed ! You may provide at'
                    ' least {0} fields to fill the database '
                    'like  dictionnary values first key.'.format(
                        data.shape[1]))
            if np.dim(data)==1 :
                data=data.reshape((1,data.shape[0]))
            asize=data.shape[1]         # check the number of columns
            len_nump=len(nump_columns)
            assert asize == len_nump ,'The shape at index 1 must be the same'\
                ' like the lenght << {0}>> of data you provide'.format(
                    len(nump_columns))
            # if as_key_for_DB ==None : 
            #     as_key_for_DB=['_'.join([val[:-2]for val in nump_columns])] 
                
            if asize >1 :
                for ii, rowlines in enumerate( data) : 
                    temp.append(tuple(rowlines.tolist()))
                    
                type_of_value=[type(item) for  item in data[0,:].tolist()] 
                index_item=[nump_columns.index(item)+1 
                            for item in nump_columns]
                num_of_val=[datavalues[:,ii].shape[0] 
                            for ii, item in enumerate(nump_columns)]
                tempkeylenvalue=(nump_columns,
                                 type_of_value,index_item,num_of_val,temp)
                if new_tablename is not None : 
                    return {new_tablename:tempkeylenvalue}
    
                return tempkeylenvalue
                
            elif asize ==1 :
                # type_of_value=type(data.tolist()[0])
                type_of_value=type(data.tolist())       
                if new_tablename is not None : 
                    return {new_tablename:(nump_columns[0],
                                           type_of_value,asize, data.tolist())}
                return (nump_columns[0],type_of_value,asize,
                        data.shape[1], data.shape[0],data.tolist())
            
            
    @staticmethod    
    def arrangeData_for_dictapp( datalist, **kwargs):
        """
        Function overwritten from "set_on_dictapp func". Reorganise 
        data to dict_app model.
        
        Parameters
        ----------
            * datalist : list, dict
                list of value providen for fill the dict_app.

        Raises
        ------
             pyCSAMTError_SQL_manager
                None dataname detected

        Returns
        -------
            dict
                datalist, Data arranged according to  dict_app arrangement.
        """
        
        comments=kwargs.pop('comments', None)
        dataname=kwargs.pop('name_of_table', None )
        
        build_tup=[]
        dicosqlDB ={}
        # we can populate later 
        field_dico_inv ={int:'i',bool:'b',
                         np.ndarray:'l',list:'t',
                         float:'f',str:'t',None:'n',
                         dict:'h'} # date:'date',buffer:'by'
        try : 
            from datetime import datetime , timezone
        except ImportError:
            raise 'Importation datetime module failed ! please try again'
        
        if type (datalist) is list :
            if dataname is None :
                raise 'Please set a table name !'
            elif dataname is not None :
                datalist={dataname:datalist}

        if type(datalist) is dict :

            for key, values in datalist.items():
                main_f , type_f, posi , num_of_val, value= values
                for ii, item  in enumerate(main_f) :
                    if type_f[ii] in field_dico_inv.keys():
                        type_f[ii]='{0}'.format(field_dico_inv[type_f[ii]])
                    else :
                        type_f[ii]='{0}'.format('k')
                    posi[ii]="{0:<3}".format(posi[ii])
                    num_of_val[ii]='{0:<7}'.format(num_of_val[ii])

                    # find other param 
                    id_='{0:<10}'.format(id(main_f[ii]))
                    dtime='{0:} - {1:}'.format(datetime.now(),timezone.utc)
                    
                    idkey='{0:<}'.format(key[:-2]+main_f[ii][:-2]+id_[0]+\
                                         num_of_val[ii][0]+main_f[ii][-1])
                    if comments is None :
                        comments ='{0:}{1:<7}{2}'.format(key[:int(len(key)/2)],
                                               str(id(datalist))[-7:],main_f[ii])
                    
                    build_tup.append((item, type_f[ii], id_, posi[ii],
                                     num_of_val[ii],dtime, idkey, comments))

                # build_tup.append(value) # better to use append instead of extend 
 
                    
                dicosqlDB ={key:(build_tup,value)} 
                
        return dicosqlDB 
    
    @redirect_cls_or_func( keepDataInfos,'func"recordData" was redirected '
                          'to func "keepDataInfos" on called it.')
    @staticmethod    
    def recordData (data, new_tablename=None , **kwargs ):
        """
        Func created to record data for easily set on dictionnary for sql arrangement. 
        Indeed the func displays values along each colummns and determine the type of 
        each column. This sort of arrangement may not be use for PostgreSQL.Reason why 
        the func is deprecated . Once call , it redirected to onother function above. 
        
        Parameters
        ----------
            * data : str, np.array, list,  pd.core.DataFrame object
                    Data ca, be on the format above or filename of data 
                    if the argument "data" is a filename, we must be 
                    convert on ".csv" format.
                
            * new_tablename : str, optional
                    Name of database. if name is not given , the 
                function return only list . The default is None.

        Raises
        ------
            IndexError
                 if lengh of number of columns like heads of data 
                 does not match the data.shape[0], then errors will occurs.  

        Returns
        -------
            list , dict
                list of value in the case of no name is providen for tablename.
                else return dict if name of datatable is providen.

        :Example:
            
            >>> from pycsamt.geodrill.geodatabase import recordData
            >>> filename='nofacies_data.csv'
            >>> recordList0=recordData(data=filename,
                                    new_tablename='example')
        """
        nump_columns=kwargs.pop('columns_of array', None)
        separate_colist=kwargs.pop('seperate_columns_list', ',')
        seperate_csvfile=kwargs.pop('sep',',')
        # table_id =kwargs.pop('id_',None)
        # comment=''
        

        tempkeylenvalue =[]
        if nump_columns is not None :
            if type (nump_columns) ==np.ndarray :
                nump_columns=nump_columns.tolist()
            elif type (nump_columns) == str :
                nump_columns=nump_columns.split(separate_colist)
            else :
                if separate_colist not in nump_columns:
                    nump_columns=[nump_columns]  # put on a list 
    
        if os.path.isfile(data) and data.endswith('.csv'): # read csv file .
            data=pd.read_csv(data, sep=seperate_csvfile,
                             delimiter=None,header='infer',
                             index_col=False)
        
        if type (data)==pd.core.frame.DataFrame:
            columns=data.columns        # get the dataframe columns 
            datavalues=data.to_numpy()  # recover numpy array
            
            for ii, kkey in enumerate( columns) : 
                # type_of_value=type(datavalues[:,ii].tolist()[0])
                type_of_value=type(datavalues[:,ii].tolist())                
                tempkeylenvalue.append((kkey,type_of_value,ii+1, 
                                        datavalues[:,ii].shape[0],
                                        datavalues[:,ii].tolist()))
            if new_tablename is not None : 
                return {new_tablename:tempkeylenvalue}

            return tempkeylenvalue
            
        elif type(data) ==np.ndarray:
            
            if nump_columns is None :
                raise Exception (
                    'The primary keys is needed ! You may provide at least'
                    ' {0} fields to fill the database '\
                                 'like  dictionnary values first key.'.format(
                                     data.shape[1]))
            if np.dim(data)==1 :
                data=data.reshape((1,data.shape[0]))
            asize=data.shape[1]         # check the number of columns
            len_nump=len(nump_columns)
            assert asize == len_nump ,'The shape at index 1 must be the same like the '\
                'lenght << {0}>> of data you provide'.format(
                    len(nump_columns))
            # if as_key_for_DB ==None : 
            #     as_key_for_DB=['_'.join([val[:-2]for val in nump_columns])] 
                
            if asize >1 :
                for ii, kkeys in enumerate(nump_columns):
                    datatolist=datavalues[:,ii].tolist()
                    # type_of_value=type(datatolist[0])
                    type_of_value=type(datatolist)                    
                    tempkeylenvalue.append((kkeys,type_of_value,
                                            ii+1,datavalues[:,ii].shape[0],
                                            datatolist))
                
                if new_tablename is not None : 
                    return {new_tablename:tempkeylenvalue}
                
                return tempkeylenvalue
                
            elif asize ==1 :
                # type_of_value=type(data.tolist()[0])
                type_of_value=type(data.tolist())       
                if new_tablename is not None : 
                    return {new_tablename:(nump_columns[0],
                                           type_of_value,asize, data.shape[1],
                                            data.tolist())}
                return (nump_columns[0],type_of_value,asize,
                        data.shape[1], data.tolist())
    
    @redirect_cls_or_func(arrangeData_for_dictapp,
                          'set_on_dict will redirect to '
                          '"arrangeData_for_dictapp once called." ')
    @staticmethod    
    def set_on_dict_app( datalist, **kwargs):
        """
        function to arrange value according to dict_app  arangement 
        The func can call python datatype to SQLdatatype . function include 
        the datetime , and the id of each data to fullfill dictionnary.
        NB : func is redirected because it's linked to RecordData whom it also
        redirected . if you use this function , it will redirect you to 
        "arrangeData_for_dict_app".

        Parameters
        ----------
            * datalist : list, dict
                list of value providen for fill the dict_app.

        Raises
        ------
            pyCSAMTError_SQL_manager 
                DataBase no found

        Returns
        -------
            dict
                datalist , Data arranged according to  dict_app arrangement.
        """
        
        comments=kwargs.pop('comments', None)
        dataname=kwargs.pop('name_of_table', None )

        # we can populate later 
        field_dico_inv ={int:'i',bool:'b',
                         np.ndarray:'l',list:'t',
                         float:'f',str:'t',None:'n',
                         dict:'h'} # date:'date',buffer:'by'
        try : 
            from datetime import datetime , timezone
        except ImportError:
            raise 'Importation datetime module failed ! please try again'
        
        if type (datalist) is list :
            if dataname is None :
                raise 'Please set a table name !'
            elif dataname is not None :
                datalist={dataname:datalist}

        if type(datalist) is dict :

            for key, values in datalist.items():
                for ii, rowval in enumerate(values) :
                    main_f, type_f, posi, num_of_val,value = rowval
                    main_f='{0}'.format(main_f)
                    if type_f in field_dico_inv.keys():
                        type_f='{0}'.format(field_dico_inv[type_f])
                    else :
                        type_f='{0}'.format('k')
                    posi="{0:<3}".format(posi)
                    num_of_val='{0:<7}'.format(num_of_val)
                    if type(value[0])==str :
                        value=['{0:<17}'.format(elm) for elm in value]
                    # find other param 
                    id_='{0:<10}'.format(id(main_f))
                    dtime='{0:} - {1:}'.format(datetime.now(),timezone.utc)
                    
                    idkey='{0:<}'.format(key[:-2]+main_f[:-2]+\
                                         id_[0]+num_of_val[0]+main_f[-1])
                    if comments is None :
                        comments ='{0:}{1:<7}{2}'.format(key[:int(len(key)/2)],
                                               str(id(datalist))[-7:],main_f)
                    valueT=(main_f, type_f, id_, posi,
                               num_of_val,value,dtime,idkey,comments)
                    values[ii]=valueT
                    
                datalist[key]=values
                
        return datalist     

if __name__=='__main__'   : 

    dbObj = GeoDataBase()    
    values = dbObj._retreive_databasecolumns (['__description', 'electrical_props'])
    print(values)
    


            
                

                
            
            
            
