# -*- coding: utf-8 -*-
#       Created on Tue Oct 13 15:28:57 2020
#       Author: Kouadio K.Laurent<etanoyau@gmail.com>
#       Licence: LGPL
"""
.. _module-requestmanager::`pycsamt.geodrill.requestmanager`
 
    :synopsis: Specially designed to Manage SQL 
"""
import os
import shutil 
import warnings
import pandas as pd 
# from pg8000 import DBAPI
import sqlite3 as sq3 

try : 
    from pycsamt.utils._csamtpylog import csamtpylog
    _logger=csamtpylog.get_csamtpy_logger(__name__)
except :
    pass


class ManageDB(object) : 
    """
    build a datable postgre Sql  from dict_app.py  simple way to make a transit
    between two objects One object dict_app to populate DataBase
        
    Arguments
    ------------
           **db_name** : str  
               name of dataBase 
           **db_host** : st 
                path to database 

    ====================  ==============  ==================================== 
    Attributes              Type            Explanation 
    ====================  ==============  ==================================== 
    connex                  object            DataBase connection 
    curs                    object            Database cursor
    ====================  ==============  ==================================== 

    ==========================  ===============================================
     Methods                     Explanation 
    ==========================  ===============================================
    dicT_sqlDB                  send infos as  dictionnary to dataBase 
    execute req                 execute a sql_request
    drop_TableDB                drop all Tables in sql memory DB or single Table 
    closeDB                     close after requests the connection and the cursor
    commit                      transfer the data to DataBase. if not the data 
                                will still in the cursor and  not in the dataBase 
    print_last_Query            print the last operating system 
    export_req                  export the request on datasheet like excelsheet . 
    ==========================  ===============================================
        
        
    :Example: 
        
        >>> from pycsamt.geodrill.requestmanager import ManageDB
        >>> path= os.getcwd()
        >>> nameofDB='memory.sq3'
        >>> manDB=ManageDB(db_name=nameofDB, 
        ...                   db_host=path)
        ... print(SqlQ.sql_req[-1])
        ... manDB.executeReq(SqlQ.sql_req[2])
        ... ss=manDB.print_last_Query()
        ... print(ss)
        ... manDB.export_req(SqlQ.sql_req[-1],
                             export_type='.csv')
        ... manDB.dicT_sqlDB(dictTables=Glob.dicoT, 
                           visualize_request=False)
        """
    def __init__(self, db_name =None, db_host=None): 
        self._logging=csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        
        self.db_host=db_host 
        self.db_name=db_name
        
        if self.db_name is not None  : 
            self.connect_DB()
            
            
    def connect_DB(self, db_host=None , db_name=None): 
        """
        Create sqqlite Database 
        
        :param db_host:   DataBase location path
        :type db_host: str
        
        :param db_name: str , DataBase name 
        :type db_name: str 
        """
        if db_host is not None :
            self.db_host = db_host 
        if db_name is not None : 
            self.db_name = db_name 
        
        mess= ''
        if self.db_name is None  :
            mess ='Could not create a DataBase ! Need to input the DataBase name.'
            
        if self.db_host is None : 
            mess ='Could not create a DataBase : No "{0}" Database path detected.'\
                ' Need to input  the path for Database location.'.format(self.db_name)
        
        if mess !='': 
            warnings.warn(mess)
            self._logging.error(mess)
            # raise CSex.pyCSAMTError_SQL(mess )
         
        # try to connect to de dataBase 
        if self.db_host is not None :
            try : 
                self.connexDB=sq3.connect(os.path.join(self.db_host, self.db_name))
            except :
                warnings.warn("Connection to SQL %s failed !." %self.db_name)
                # raise CSex.pyCSAMTError_SQL_manager('Connection to SQL dataBase failed ,Please try again.') 
                
                self.success=0
            else : 
                self.curs =self.connexDB.cursor()
                self.success=1
        
    def dicT_sqlDB(self, dicTables, **kwargs): 
        """
        Method to create Table for sqlDataBase . 
        Enter Data to DataBase from dictionnary. 
        Interface objet : Database _Dictionnary  
        to see how dicTable is arranged , may consult dict_app module 
                
        Parameters
        ----------
        * dictTables : dict
               Rely on dict_app.py module. it populates the  datababse 
               from dictionnay app 
               
        Returns
        ---------
          str 
            execute queries from dict_app 
            
        :Example :
            
            >>>  mDB=GestionDB (dbname='memory.sq3, 
            ...                   db_host =os.getcwd()')
            >>> mDB.dicT_sqlDB(dicTables=Glob.dicoT,
            ...                  visualize_request=False)
            >>> ss=mB.print_last_query()
            >>> print(ss)
        """
        visual_req=kwargs.pop('visualize_request', False)
        
        field_dico ={'i':'INTEGER',"t":"TEXT",'b':'BOOL',
                     'd': 'HSTORE',"k": "SERIAL", 'n':'NULL',
                     "f": 'REAL','s':'VARCHAR','date':'date',
                     'by':'BYTHEA','l':'ARRAY',
                     }

        for table in dicTables: 
            req="CREATE TABLE %s (" % table
            pk=""
            for ii, descrip in enumerate(dicTables[table]): 
                field=descrip[0]
                tfield=descrip[1] # Type of field 
                # for keys in  field_dico.keys():
                if tfield in field_dico.keys():
                    # if tfield == keys : 
                    typefield=field_dico[tfield]
                else :
                    # sql vriable nom :'s':'VARCHAR'
                    typefield='VARCHAR(%s)'%tfield
                        
                req= req+'%s %s, ' %(field, typefield)
                

            if pk=='': 
                req=req [:-2] + ")" # delete the last ',' on req.
            else : 
                req =req +"CONSTRAINT %s_pk PRIMARY KEYS(%s))" %(pk,pk)
                
            if visual_req is True : # print the request built .
                print(req)                
            try : 
                self.executeReq(req)
            except : 
                pass # the case where the table already exists. 
            
        return req

    
    def executeReq(self, query, param=None):
        """
        Execute request  of dataBase  with detection of error. 
        
        Parameters  
        -----------
        * query : str  
                sql_query 
        * param : str 
                Default is None . 
        
        raise : 
        -------
            layout the wrong sql queries . 
            
        return : 
        -------
         int  
             1, the request has been successuful run .
              
        :Example: 
            
                >>>  for keys in Glob.dicoT.keys(): 
                ...        reqst='select * from %s'%keys
                >>>  ss=manageDB.executeReq(query=reqst)
                >>>  print(ss)
        """
        
        try :
            if param is None :
                self.curs.execute(query)# param)
            else :
                self.curs.execute(query, param)
                
        except: 
            warnings.warn(f'Request SQL {query} failed. May trouble of SQL  '
                          'server connexion. Please try again later ')
            # raise (f'Request SQL {query}executed failed',err)
            return 0
        else : 
            return 1
    
    def drop_TableDB(self, dicTables, drop_table_name=None ,
                     drop_all=False): 
        """
        Drop the name of table on dataBase or all databases.

        Parameters
        ----------
        * dicTables : dict
                application dictionnary. Normally provide from 
                dict_app.py module 
            
        * drop_table_name : str, optional
                field name of dictionnay (Table Name). 
                The default is None.
            
        * drop_all : Bool, optional
                Must select if you need to drop all table. 
                The default is False.

        Raises
        ------
            Exception : Errors occurs ! .

        """
        
        if drop_all is False and drop_table_name is None : 
            raise 'Must be input at least one name contained of keys in the dicT_app'
        elif drop_all is True : 
            for keys in dicTables.keys():
                req="DROP TABLE %s" %keys
                self.executeReq(req)
        elif drop_table_name is not None : 
            if drop_table_name in dicTables.keys(): 
                req="DROP TABLE %s" % drop_table_name
            else : 
                raise'No such name in the dictionnary application Table!'\
                    'Dict_app keys Tables Names are : ¨{0}'.format(dicTables.keys())
                    
            self.executeReq(req)
            
        self.connexDB.commit()
    
    def closeDB(self):
        """
        simple method to close Database. 
        """
        
        if self.connexDB : 
            # self.curs.close()
            self.connexDB.close()
            
    def commit(self) :
        """
        special commit method for the database when cursor  and connexion 
        are still open.
        """
        if self.connexDB : 
            self.connexDB.commit()
    
    def print_query(self, column_name=None ) : 
        """
        return the result of the previous query.
        
        Parameters :
        ------------
            * query_table_nam : str
                name of table to fetch colounm data .
        """
        if  column_name is not None :
            return self.curs.fetchone()
        else :
            return self.curs.fetchall()
    
    def export_req(self, query =None ,
                   export_type='.csv', **kwargs):
        """
        method to export data  from DataBase 

        Parameters
        ----------
        * query : str, optional
            Sql requests. You may consult sql_request files. 
            The default is None.
        * export_type : Str, optional
            file extension. if None , it will export on simple file. 
            The default is '.csv'.
        * kwargs : str
            Others parameters.

        Raises
        ------
        Exception
            Print wrong sqlrequests.

        :Example:
            
            >>> from sqlrequests import SqlQ
            >>> manageDB.executeReq(SqlQ.sql_req[2])
            >>> ss=manageDB.print_last_Query()
            >>> print(ss)
            >>> manageDB.export_req(SqlQ.sql_req[-1],
                                export_type='.csv',
                                )
        """
        
        exportfolder=kwargs.pop('savefolder','savefiles')
        filename =kwargs.pop('filename','req_file')
        indexfile=kwargs.pop('index',False)
        headerfile=kwargs.pop("header",True)
        
        if query is None : 
            raise Exception ("SQL requests (%s) no found ! Please Try to put "\
                             " your sql_requests"% query)
        elif query is not None : 
            
            df_sql=pd.read_sql_query(query,self.connexDB)
            
        if filename.endswith('.csv'):
            export_type ='.csv'
            filename=filename[:-4]
        elif filename.endswith(('.xlxm', '.xlsx', '.xlm')): 
            export_type='.xlsx'
            filename=filename[:-5]
        else : 
              assert export_type  is not None , 'Must input the type to export file.'\
                  ' it maybe ".csv" or ".xlsx"'

        #-----export to excel sheet  
        if export_type in ['csv','.csv', 'comma delimited',
                             'comma-separated-value','comma sperated value',
                                       'comsepval']:
            # export to excelsheet:  
            df_sql.to_csv(filename+'.csv', header=headerfile,
                  index =indexfile,sep=',', encoding='utf8')
            
            sql_write =1 
        elif export_type in ['xlsx','.xlsx', 'excell',
                             'Excell','excel','Excel','*.xlsx']: 
            df_sql.to_excel(filename+'.xlsx',sheet_name=filename[:-3],
                index =indexfile)
            sql_write =0

        #wite a new folder
        if exportfolder is not None : 
            try : 
                os.mkdir('{0}'.format(exportfolder))    
            except OSError as e :
                print(os.strerror(e.errno))
        sql_path=os.getcwd()
        savepath=sql_path+'/{0}'.format(exportfolder)
        
        if sql_write ==1 : 
            shutil.move(filename +'.csv', savepath)
            print('---> outputDB_file <{0}.csv> has been written.'.format(filename))
        elif sql_write ==0 :  
            shutil.move(filename +'.xlsx',savepath)
            print('---> outputDB_file <{0}.xlsx> has been written.'.format(filename))
        
        


    
    
    
                  
                        
                    
            
        



 
