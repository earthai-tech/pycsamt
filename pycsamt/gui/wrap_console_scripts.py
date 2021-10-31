# -*- coding: utf-8 -*-
#       Copyright © 2021  Kouadio K.Laurent, Licence: LGPL
#       @author: KouaoLaurent alias @Daniel03 <etanoyau@gmail.con>
#       Created on Tue Apr  6 12:58:45 2021
"""
    .Console scripts. 
    write pseudo-wrapper  to interact with python shell environnement. 
    .. note:: Try to use dict_kwargs to collect input arguments values, 
            as well as defaults values .
"""

from pycsamt.utils._csamtpylog import csamtpylog 
_logger =csamtpylog.get_csamtpy_logger(__name__)


class wrap_cscripts(object): 
    """
    get input arguments from user , sanitize it and convoy towards required 
    modules 
    :param default_kwargs: dictionnary of kwargs arguments with the defaults 
        values 
    :type default_kwargs: dict 
    
    .. note:: to sanitize  input arguments , provided `*` for arguments which
        must be converte into `integer` and double `**` for those  who should be 
        convert into `float`. If conversion is not possible, will
        take the default arguments
        
    :Example :
        >>> from pycsamt.gui import wrap_console_scripts as wrs
        >>>  indico ={'Path to avg file':os.path.join(os.path.abspath('.') ,
        ...                                              'data', 'avg') 
        ...     'number of skin depth': 9,'Reference frequency':1024. }
        >>> indico = None 
        >>> cs_obj = wrs.wrap_cscripts(default_kwargs = default_kwargs, 
        ...                   input_kwargs =indico)
        >>> print(cs_obj.sanitize_kwargs)
    """
    
    def __init__(self, default_kwargs, **kws) :
        
        self._logging = csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        
        self.default_kwargs =default_kwargs
  
        self.input_kwargs =kws.pop('input_kwargs', {})
        if self.input_kwargs is None : self.input_kwargs={}
        
        self.sanitize_kwargs ={}

        if self.default_kwargs is not None : 
            self._console_scripts()
        
        
    def _console_scripts(self, default_kwargs=None, **kws):
        """
        Console scripts methods used to collect input values for sanitizeed to 
        the corresponding functions before calling main function. 

        :param default_kwargs: default values acceptables for the func
        :type default_kwargs: dict
        :param input_kwargs: dictionnary of input values . If values of dict is 
            given,values should be considered as consoles input  values and will
            not request again. 
        :type input_kwargs: dict 
        
        :return: sanitized input request and values 
        :rtype: dict

        """
        self._logging.info('wraps input arguments value from user ! ')

        if default_kwargs is not None : 
            self.default_kwargs =default_kwargs 
            
        input_kwargs =kws.pop('input_kwargs', None)
        
        if input_kwargs is not None :
            if not isinstance(input_kwargs , dict):
                self.input_kwargs ={}
            else : 
                self.input_kwargs = input_kwargs 
      
        for kk, (knames, values) in enumerate(list(self.default_kwargs.items())):

            if len(list(self.input_kwargs.keys())) !=0 : 
                tem_= list(self.input_kwargs.keys())
                for jj, iname in enumerate(tem_): 
                    if knames.lower().find(iname.lower()) >=0 : #value found 
                        user_response = self.input_kwargs[iname]
                        break 
                    else : user_response = None  # value no found then ask the user 
                if user_response is None : #  request the user   
                    user_response = input('> {0} [{1}]: '.format(knames, values))
                    if user_response in ['','\n' , ' ']:
                        user_response = values
                        print('--> set to default = {}'.format(values))

            elif len(list(self.input_kwargs.keys())) ==0 :

                user_response = input('> {0} [{1}]: '.format(knames, values))
                if user_response in ['','\n' , ' ']:
                    user_response = values
                    print('--> set to default = {}'.format(values))
                
            if user_response =='None': user_response =None 
            litteral_txt = knames + ' [{0}]'.format(values)
            
            if '*' or '**' in litteral_txt : 
                litteral_text , user_response = self._sanitize_user_response(
                    litteral_text =litteral_txt, user_response=user_response,
                    default_value = values)
                
            self.sanitize_kwargs[litteral_text]= user_response 
            
    @staticmethod   
    def _sanitize_user_response(litteral_text, user_response,  default_value):
        """
        Sanitize user response with default types. In fact recheck string 
        input values from user . Then figure out  value with its 
        corresponding type using  jokers `**` and `*`. 
            * `**` means ---> convert into float values 
            * `*` means convert into interger values .
            
        """
        _logger.info(' pipeline to sanitize input argument !')
    
        dc =None 
        if  litteral_text.find('*')>=0 : dc = '*'
        elif litteral_text.find('**')>=0 : dc ='**'
        
        if dc is not None : 
            try : 
                if dc =='*' : user_response = int(user_response) 
                elif dc=='**' :user_response = float(user_response)
            
            except : 
                print(' Wrong given value ! '
                          ' should set to default = {0}'.format(default_value))
                user_response =default_value 
                
            litteral_text= litteral_text.replace(dc, '')     

        return litteral_text , user_response 
