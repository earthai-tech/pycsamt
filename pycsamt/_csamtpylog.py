# -*- coding: utf-8 -*-
#       Author: Kouadio K.Laurent<etanoyau@gmail.com>
#       Licence: LGPL
"""
    .looging file was tst only on .yaml file.  
    . see load_configure_set_logfile (only for .yaml file )
    
Created on Sun Sep 27 19:29:59 2020

A more Pythonic way of logging:
Define a class MtPyLog to wrap the python logging module;
Use a (optional) configuration file (yaml, ini, json) to configure the logging,
It will return a logger object with the user-provided config setting.
see also: http://www.cdotson.com/2015/11/python-logging-best-practices/


"""
import os 
import yaml
import logging 
import logging.config
import inspect 


class csamtpylog:
    """
    Field to logs csamtpy module Files  in order to tracks all 
    exceptions
    
    """
    
    @staticmethod
    def load_configure (path2configure =None, OwnloggerBaseConf=False) :
        """
        configure/setup the logging according to the input configfile

        :param configfile: .yml, .ini, .conf, .json, .yaml.
        Its default is the logging.yml located in the same dir as this module.
        It can be modified to use env variables to search for a log config file.
        """
        
        configfile=path2configure
        
        if configfile is None or configfile == "":
            if OwnloggerBaseConf ==False :
                logging.basicConfig()
            else :
                csamtpylog().set_logger_output()
            
        elif configfile.endswith(".yaml") or configfile.endswith(".yml") :
            this_module_path=os.path.abspath(__file__)
            # print(this_module_path)
            
            logging.info ("this module is : %s", this_module_path)
            yaml_path=os.path.join(os.path.dirname(this_module_path),
                                   configfile)
            
            logging.info('Effective yaml configuration file %s', yaml_path)

            if os.path.exists(yaml_path) :
                with open (yaml_path,"rt") as f :
                    config=yaml.safe_load(f.read())
                logging.config.dictConfig(config)
            else :
                logging.exception(
                    "the config yaml file %s does not exist?", yaml_path)
                
        elif configfile.endswith(".conf") or configfile.endswith(".ini") :
            logging.config.fileConfig(configfile,
                                     disable_existing_loggers=False)
            
        elif configfile.endswith(".json") :
            pass 
        else :
            logging.exception("logging configuration file %s is not supported" %
                configfile)
            
    
    @staticmethod        
    def get_csamtpy_logger(loggername=''):
        """
        create a named logger (try different)
        :param loggername: the name (key) of the logger object in this Python interpreter.
        :return:
        """
        logger =logging.getLogger(loggername)
        csamtpylog.load_configure() #set configuration 

        return logger
    
    
    @staticmethod
    def set_logger_output(logfilename="myapplog.log", 
                      date_format='%m %d %Y %I:%M%S %p',
                      filemode="w",format_="%(asctime)s - %(name)s - %(levelname)s - %(message)s",#'%s(asctime)s - %(message)s',
                       level=logging.DEBUG ):
        """
        
        Parameters
        ----------
        * logfilename : str, optional
                        logfilename output. The default is "LogFileTest.log".
        * date_format :str, 
                        format of date . consult module time. The default is '%m %d %Y %I:%M%S %p'.
        * filemode : str, 
                        mode to write the out[ut file. Must be ["a","w", "wt',"wt"].
                                                    The default is "a".
        * format_ : str,
                        format of date and time. Must consult the module of datetime.
                        The default is '%s(asctime)s %(message)s'.
        * level : object, optional
                    logging object.  Must be C.E.W.I.D :
                        >>> logging.CRITICAL
                        >>> Logging.ERROR
                        >>> logging.WARNING
                        >>> logging.INFO
                        >>> logging.DEBUG
                    The default is logging.DEBUG.

        Returns : NoneType object
        -------
        Procedure staticFunctions

        """
        # #create logger 
        logger=csamtpylog.get_csamtpy_logger(__name__) # logger=logging.getLogger(__name__)
        logger.setLevel(level)
        
        # #create console handler  and set the level to DEBUG 
        ch=logging.StreamHandler()
        ch.setLevel(level) # level=logging.DEBUG
        
        #create formatter :
        formatter =logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        
        #add Formatter to ch 
        ch.setFormatter(formatter)
        
        #add ch to logger 
        logger.addHandler(ch)
        
        
        # logging.basicConfig(filename="test.log", level=logging.DEBUG)
 
        logging.basicConfig(filename=logfilename, format=format_,
                            datefmt=date_format,
                            filemode=filemode, level=level)
        logging.debug("debug message")
        logging.info("info message ")
        logging.warning("warn message")
        logging.error("error message")
        logging.critical("critical message")
        
    @staticmethod
    def load_configure_set_logfile (path2configfile=None): # loggername =None, setLevel=Name
        """
        configure/setup the logging according to the input configure .yaml file.

        :param configfile: .yml, or add ownUsersConfigYamlfile (*.yml) 
        Its default is the logging.yml located in logfiles folder 
        It can be modified to use env variables to search for a log config file.
        
        """
        
        ownUserLogger="main_logging_configfile.yml"
        if path2configfile is None :
            env_var=os.environ['pyCSAMT']
            path2configfile =os.path.join( env_var, 'csamtpy','geoDrill','_logfiles',
                ownUserLogger)
            

        elif path2configfile is not None :
            if os.path.isdir(os.path.dirname(path2configfile)):
                if path2configfile.endswith('.yml') or path2configfile.endswith('.yaml'):
                    
                    logging.info('Effective yaml configuration file :%s', path2configfile)
                else :
                    logging.exception('File provided {%s}, is not a .yaml config file !'%os.path.basename(path2configfile))
            else :
                
                logging.exception ('Wrong path to .yaml config file.')
        
        yaml_path=path2configfile
        os.chdir(os.path.dirname(yaml_path))
        if os.path.exists(yaml_path) :
            with open (yaml_path,"rt") as f :
                config=yaml.safe_load(f.read())
            logging.config.dictConfig(config)
        else :
            logging.exception(
                "the config yaml file %s does not exist?", yaml_path) 
                
                    

def test_yaml_configfile(yamlfile="logging.yml"):
    
    this_func_name = inspect.getframeinfo(inspect.currentframe())[2]
    # print("this_func_name:\n",this_func_name)
    # print(("starting", this_func_name))
    
    UsersOwnConfigFile = yamlfile
    csamtpylog.load_configure(UsersOwnConfigFile)
    # print(myobj)
    logger = csamtpylog.get_csamtpy_logger(__name__)
    
    print((logger, id(logger), logger.name, logger.level, logger.handlers))
    
    # 4 use the named-logger
    logger.debug(this_func_name + ' __debug message')
    logger.info(this_func_name + ' __info message')
    logger.warn(this_func_name + ' __warn message')
    logger.error(this_func_name + ' __error message')
    logger.critical(this_func_name + ' __critical message')

    print(("End of: ", this_func_name))
    
def test_config_configfile(configfile="logging.conf"): 
    
    this_func_name=inspect.getframeinfo(inspect.currentframe())[1]
    
    print(('starting', this_func_name))
    UsersOwnConfigFile =configfile
    
    csamtpylog.load_configure(UsersOwnConfigFile)
    
    logger=csamtpylog.get_csamtpy_logger(__name__)
    
    print((logger, id(logger), logger.name, logger.level, logger.handlers))
    
    logger.debug(this_func_name + ' __debug message')
    logger.info(this_func_name + ' __info message')
    logger.warn(this_func_name + ' __warn message')
    logger.error(this_func_name + ' __error message')
    logger.critical(this_func_name + ' __critical message')
    
    print(("End of :", this_func_name))
    
    
def test_none_configfile():
    
    this_func_name =inspect.getframeinfo(inspect.currentframe())[2]
    
    print(('Starting of :', this_func_name))
    
    UsersOwnConfigFile =''
    
    csamtpylog.load_configure(UsersOwnConfigFile)
    
    logger=csamtpylog.get_csamtpy_logger(__name__)
    
    print((logger, id(logger), logger.name, logger.level, logger.handlers))
    
    logger.debug(this_func_name + ' __debug message')
    logger.info(this_func_name + ' __info message')
    logger.warning(this_func_name + ' __warn message')
    logger.error(this_func_name + ' __error message')
    logger.critical(this_func_name + ' __critical message')
    
    print(("End of :", this_func_name))
    
    
        
def test_jsaon_configfile():
    pass 



if __name__=="__main__":

    csamtpylog.load_configure_set_logfile()
    logger=csamtpylog.get_csamtpy_logger(__name__)
    # this_func_name=inspect.getframeinfo(inspect.currentframe())[1]
    # test_yaml_configfile(ownUserLogger)
    logger.debug(' __debug message')
    logger.info(' __info message')
    logger.warning(' __warn message')
    logger.error(' __error message')
    logger.critical(' __critical message')
    
   
    
    
    
    
    
    
    
    
            
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
