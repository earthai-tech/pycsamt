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
.. _module-Structural::`geodrill.geoCore.structural` 
    
    :synopsis:class for  geological structural analysis 
               contains some conventional structure can populate the data
    
Created on Sat Nov 28 21:19:13 2020

@author: @Daniel03

"""


import os, warnings
import numpy as np 
from csamtpy.ff.processing.callffunc import set_stratum_on_dict as strato
from csamtpy.utils  import exceptions as CSex
from csamtpy.utils._csamtpylog import csamtpylog
_logger=csamtpylog.get_csamtpy_logger(__name__)

#---- End import modules --------

class geo_pattern: 
    """
    Singleton class to deal with geopattern  with other modules. It is and exhaustive pattern
    dict, can be add and change. This pattern  will be depreacted later , to create for pyCSAMT,
    its owwn geological pattern in coformity with the conventional geological swatches .
     deal with USGS(US Geological Survey ) swatches- references and FGDC (Digital cartographic 
    Standard for Geological  Map Symbolisation):
         
    .. FGDCgeostdTM11A2_A-37-01cs2.eps :: 
        make _pattern:{'/', '\', '|', '-', '+', 'x', 'o', 'O', '.', '*'}
                /   - diagonal hatching
                \   - back diagonal
                |   - vertical
                -   - horizontal
                +   - crossed
                x   - crossed diagonal
                o   - small circle
                O   - large circle
                .   - dots
                *   - stars
    """
    pattern={
                        "basement rocks" :      ['.+++++.', (.25, .5, .5)],
                        "igneous rocks":        ['.o.o.', (1., 1., 1.)], 
                        "duricrust"   :         ['+.+',(1., .2, .36)],
                        "gravel" :              ['oO',(.75,.86,.12)],
                        "sand":                 ['....',(.23, .36, .45)],
                        "conglomerate"    :     ['.O.', (.55, 0., .36)],
                        "dolomite" :            ['.-.', (0., .75, .23)],
                        "limestone" :           ['//.',(.52, .23, .125)],
                       "permafrost"  :          ['o.', (.2, .26, .75)],
                        "metamorphic rocks" :   ['*o.', (.2, .2, .3)],
                        "tills"  :              ['-.', (.7, .6, .9)],
                        "standstone ":          ['..', (.5, .6, .9)],
                        "lignite coal":         ['+/.',(.5, .5, .4)],
                        "coal":                 ['*.', (.8, .9, 0.)],
                        "shale"   :             ['=', (0., 0., 0.7)],
                        "clay"   :              ['=.',(.9, .8, 0.8)],
                        "saprolite" :           ['*/',(.3, 1.2, .4)],
                        "sedimentary rocks":    ['...',(.25, 0., .25)],
                        "fresh water"  :        ['.-.',(0., 1.,.2)],
                        "salt water"   :        ['o.-',(.2, 1., .2)],
                        "massive sulphide" :     ['.+O',(1.,.5, .5 )],
                        "sea water"     :       ['.--',(.0, 1., 0.)],
                        "ore minerals"  :       ['--|',(.8, .2, .2)],
                        "graphite"    :         ['.++.',(.2, .7, .7)],
                        }
            
class Structure :
    """
    Class for typical geological strutural conventions 
    for AGSO_STCODES .  All geological structural informations are
    geostructral object.
                  
    Holds the following information:
        
    ==========================  ===============  ===============================
    Attributes                  Type             Explanation
    ==========================  ===============  ===============================
    boudin_axis                 geos_obj         boudin    
    fold_axial_plane            geos_obj         axial plam of structural fold.
    banding_gneissosity         geos_obj         gneissossity of boudin plan  
    s_fabric                    geos_obj         fabric plan
    fault_plane                 geos_obj         fault plan 
    fracture_joint_set          geos_obj         fracture joint 
    undifferentiated_plane      geos_obj         unnamed geological structure 
    sharp_contact               geos_obj         sharp contact `various discrepancy` 
                                                 contact `stratigraphy discrepancy`  
                                                 fracture or fault discrepancies
    ==========================  ===============  ================================

    More attributes can be added by inputing a key word dictionary

    :Example: 
        
        >>> from geodrill.geoCore.structural import Structure
        >>> structural=Structure()
        >>> boudin=boudin_axis()
        >>> print(boudin.code) 
        >>> print(structural.boudin_axis.name)
        >>> print(structural.boudin_axis.color)
    """  
    _logger.info('Set Structural main geological informations. ')
    
    def __init__(self, **kwargs):
        
        
        self.boudin_axis =boudin_axis()
        self.fold_axial_plane=fold_axial_plane()
        self.banding_gneissosity=banding_gneissosity()
        self.s_fabric=s_fabric()
        self.fault_plane=fault_plane()
        self.fracture_joint_set=fracture_joint_set()
        self.undifferentiated_plane=undifferentiated_plane()
        self.sharp_contact=sharp_contact()
        
        
        for keys in list(kwargs.keys()):
            setattr(self, keys, kwargs[keys])

#====================================================
#   geological structural objects 
#=====================================================
class boudin_axis(object):
    """
    Special class for boudins_axis
    
    Holds the following information:
    
    =================  =============  ========================================
    Attributes         Type             Explanation
    =================  =============  ========================================
    code                str             conventional code    
    label               str             named label
    size                str             drawing size  
    pattern             str             drawing pattern 
    density             str             elmts density 
    thickness           str             drawing thickness 
    color               str             color set
    =================  =============  ========================================

    More attributes can be added by inputing a key word dictionary
    """
    def __init__(self, **kwargs):
        
        self.code =None 
        self.label=None 
        self.name=None 
        self.pattern=None 
        self.size =None 
        self.density=None 
        self.thickness=None
        self.color =None
        
        for keys in list(kwargs.keys()):
            self.__setattr__(keys, kwargs[keys])
            
        self._set_boudin_axis()   
        
    def _set_boudin_axis(self):
        """
        methode to populates attributes

        """
        baxis=strato()[1]
        for keys, values in baxis.items():
            if keys == self.__class__.__name__:
                self.code =values.__getitem__('code')
                self.label=values.__getitem__('label')
                self.name=values.__getitem__('name')
                self.pattern=values.__getitem__('pattern')
                self.size=values.__getitem__('size')
                self.color=values.__getitem__('color')
                self.density=values.__getitem__('density')
                self.thickess=values.__getitem__('thickness')

class fold_axial_plane(object):
    """
    Special class for fold_axial_plane
    
    Holds the following information:

    =================  =============  ========================================
    Attributes              Type          Explanation
    =================  =============  ========================================
    code                str             conventional code    
    label               str             named label
    size                str             drawing size  
    pattern             str             drawing pattern 
    density             str             elmts density 
    thickness           str             drawing thickness 
    color               str             color set
    =================  =============  ========================================

    More attributes can be added by inputing a key word dictionary
    """
    def __init__(self, **kwargs):
        
        self.code =None 
        self.label=None 
        self.name=None 
        self.pattern=None 
        self.size =None 
        self.density=None 
        self.thickness=None
        self.color =None
        
        self._set_fold_axial_plane()
        
        for keys in list(kwargs.keys()):
            self.__setattr__(keys, kwargs[keys])
            
    def _set_fold_axial_plane(self):
        """
        methode to populates attributes

        """
        baxis=strato()[1]
        for keys, values in baxis.items():
            if keys == self.__class__.__name__:
                self.code =values.__getitem__('code')
                self.label=values.__getitem__('label')
                self.name=values.__getitem__('name')
                self.pattern=values.__getitem__('pattern')
                self.size=values.__getitem__('size')
                self.color=values.__getitem__('color')
                self.density=values.__getitem__('density')
                self.thickess=values.__getitem__('thickness')

class banding_gneissosity(object):
    """
    Special class for banding_gneissosity
    
    Holds the following information:

    =================  =============  ========================================
    Attributes              Type          Explanation
    =================  =============  ========================================
    code                str             conventional code    
    label               str             named label
    size                str             drawing size  
    pattern             str             drawing pattern 
    density             str             elmts density 
    thickness           str             drawing thickness 
    color               str             color set
    =================  =============  ========================================

    More attributes can be added by inputing a key word dictionary
    """
    def __init__(self, **kwargs):
        
        self.code =None 
        self.label=None 
        self.name=None 
        self.pattern=None 
        self.size =None 
        self.density=None 
        self.thickness=None
        self.color =None
        
        self._set_banding_gneissosity()
        
        for keys in list(kwargs.keys()):
            self.__setattr__(keys, kwargs[keys])
            
    def _set_banding_gneissosity(self):
        """
        methode to populates attributes

        """
        baxis=strato()[1]
        # print(baxis)
        for keys, values in baxis.items():
            if keys == self.__class__.__name__:
                self.code =values.__getitem__('code')
                self.label=values.__getitem__('label')
                self.name=values.__getitem__('name')
                self.pattern=values.__getitem__('pattern')
                self.size=values.__getitem__('size')
                self.color=values.__getitem__('color')
                self.density=values.__getitem__('density')
                self.thickess=values.__getitem__('thickness')             

class s_fabric(object):
    """
    Special class for s_fabric
    
    Holds the following information:

    =================  =============  ========================================
    Attributes              Type          Explanation
    =================  =============  ========================================
    code                str             conventional code    
    label               str             named label
    size                str             drawing size  
    pattern             str             drawing pattern 
    density             str             elmts density 
    thickness           str             drawing thickness 
    color               str             color set
    =================  =============  ========================================

    More attributes can be added by inputing a key word dictionary
    """
    def __init__(self, **kwargs):
        
        self.code =None 
        self.label=None 
        self.name=None 
        self.pattern=None 
        self.size =None 
        self.density=None 
        self.thickness=None
        self.color =None
        
        self._set_s_fabric()
        
        for keys in list(kwargs.keys()):
            self.__setattr__(keys, kwargs[keys])
            
    def _set_s_fabric(self):
        """methode to populates attributes"""
        baxis=strato()[1]
        # print(baxis)
        for keys, values in baxis.items():
            if keys == self.__class__.__name__:
                self.code =values.__getitem__('code')
                self.label=values.__getitem__('label')
                self.name=values.__getitem__('name')
                self.pattern=values.__getitem__('pattern')
                self.size=values.__getitem__('size')
                self.color=values.__getitem__('color')
                self.density=values.__getitem__('density')
                self.thickess=values.__getitem__('thickness')  

class fault_plane(object):
    """
    Special class for fault_plane
    
    Holds the following information:

    =================  =============  ========================================
    Attributes          Type            Explanation
    =================  =============  ========================================
    code                str             conventional code    
    label               str             named label
    size                str             drawing size  
    pattern             str             drawing pattern 
    density             str             elmts density 
    thickness           str             drawing thickness 
    color               str             color set
    =================  =============  ========================================

    More attributes can be added by inputing a key word dictionary
    """
    def __init__(self, **kwargs):
        
        self.code =None 
        self.label=None 
        self.name=None 
        self.pattern=None 
        self.size =None 
        self.density=None 
        self.thickness=None
        self.color =None
        
        self._set_fault_plane()
        
        for keys in list(kwargs.keys()):
            self.__setattr__(keys, kwargs[keys])
            
    def _set_fault_plane(self):
        """
        methode to populates attributes

        """
        baxis=strato()[1]
        # print(baxis)
        for keys, values in baxis.items():
            if keys == self.__class__.__name__:
                self.code =values.__getitem__('code')
                self.label=values.__getitem__('label')
                self.name=values.__getitem__('name')
                self.pattern=values.__getitem__('pattern')
                self.size=values.__getitem__('size')
                self.color=values.__getitem__('color')
                self.density=values.__getitem__('density')
                self.thickess=values.__getitem__('thickness')
                
class fracture_joint_set(object):
    """
    Special class for fracture_joint_set
    
    Holds the following information:

    =================  =============  ========================================
    Attributes         Type           Explanation
    =================  =============  ========================================
    code                str             conventional code    
    label               str             named label
    size                str             drawing size  
    pattern             str             drawing pattern 
    density             str             elmts density 
    thickness           str             drawing thickness 
    color               str             color set
    =================  =============  ========================================

    More attributes can be added by inputing a key word dictionary
    """
    def __init__(self, **kwargs):
        
        self.code =None 
        self.label=None 
        self.name=None 
        self.pattern=None 
        self.size =None 
        self.density=None 
        self.thickness=None
        self.color =None
        
        self._set_fracture_joint_set()
        
        for keys in list(kwargs.keys()):
            self.__setattr__(keys, kwargs[keys])
            
    def _set_fracture_joint_set(self):
        """methode to populates attributes
        """
        baxis=strato()[1]
        # print(baxis)
        for keys, values in baxis.items():
            if keys == self.__class__.__name__:
                self.code =values.__getitem__('code')
                self.label=values.__getitem__('label')
                self.name=values.__getitem__('name')
                self.pattern=values.__getitem__('pattern')
                self.size=values.__getitem__('size')
                self.color=values.__getitem__('color')
                self.density=values.__getitem__('density')
                self.thickess=values.__getitem__('thickness')
                
class undifferentiated_plane(object):
    """
    Special class for undifferentiated_plane
    
    Holds the following information:

    =================  =============  ========================================
    Attributes         Type             Explanation
    =================  =============  ========================================
    code                str             conventional code    
    label               str             named label
    size                str             drawing size  
    pattern             str             drawing pattern 
    density             str             elmts density 
    thickness           str             drawing thickness 
    color               str             color set
    =================  =============  ========================================

    More attributes can be added by inputing a key word dictionary
    """
    def __init__(self, **kwargs):
        
        self.code =None 
        self.label=None 
        self.name=None 
        self.pattern=None 
        self.size =None 
        self.density=None 
        self.thickness=None
        self.color =None
        
        self._set_undifferentiated_plane()
        
        for keys in list(kwargs.keys()):
            self.__setattr__(keys, kwargs[keys])
            
    def _set_undifferentiated_plane(self):
        """
        methode to populates attributes

        """
        baxis=strato()[1]
        # print(baxis)
        for keys, values in baxis.items():
            if keys == self.__class__.__name__:
                self.code =values.__getitem__('code')
                self.label=values.__getitem__('label')
                self.name=values.__getitem__('name')
                self.pattern=values.__getitem__('pattern')
                self.size=values.__getitem__('size')
                self.color=values.__getitem__('color')
                self.density=values.__getitem__('density')
                self.thickess=values.__getitem__('thickness')             
                
class sharp_contact(object):
    """
    Special class for sharp_contact
    
    Holds the following information:

    =================  =============  ========================================
    Attributes          Type           Explanation
    =================  =============  ========================================
    code                str             conventional code    
    label               str             named label
    size                str             drawing size  
    pattern             str             drawing pattern 
    density             str             elmts density 
    thickness           str             drawing thickness 
    color               str             color set
    ==================  =============  ========================================
    
    More attributes can be added by inputing a key word dictionary
    """
    def __init__(self, **kwargs):
        
        self.code =None 
        self.label=None 
        self.name=None 
        self.pattern=None 
        self.size =None 
        self.density=None 
        self.thickness=None
        self.color =None
        
        self._set_sharp_contact()
        
        for keys in list(kwargs.keys()):
            self.__setattr__(keys, kwargs[keys])
            
    def _set_sharp_contact(self):
        """
        methode to populates attributes
        """
        baxis=strato()[1]
        # print(baxis)
        for keys, values in baxis.items():
            if keys == self.__class__.__name__:
                self.code =values.__getitem__('code')
                self.label=values.__getitem__('label')
                self.name=values.__getitem__('name')
                self.pattern=values.__getitem__('pattern')
                self.size=values.__getitem__('size')
                self.color=values.__getitem__('color')
                self.density=values.__getitem__('density')
                self.thickess=values.__getitem__('thickness')             
                
class Geo_formation (object): 
    """
    This class is an axilliary class to supplement geodatabase , 
    if the GeodataBase doesnt reply  to SQL request  , then use this class
    to secah information about structures .  If SQL is done as well ,
    program won't call this class as rescure . 
    Containers of more than  150 geological strutures.
        
   
    
    .. note :: replace in attributes param "**" by  the *name of struture*
    
    ==================  ============  =========================================
    Attributes          Type           Explanation
    ==================  ============  =========================================
    names               array_like      names of all geological strutures 
    codes               array_like      names of all geological codes 
    **code              str             code of specific geological structure 
    **label             str             label of specific structure
    **name              str             label of specific structure
    **pattern           str             pattern of specific structure  
    **pat_size          str             pattern size  of specific structure
    **pat_density       str             pattern density l of specific structure
    **pat_thickness     str             pttern thickess of specific structure
    **color             str             color of specific structure
    ==================  ============  =========================================

    1.  To see the names of strutures , write the script below 
    
    :Example:
        
        >>> from geodrill.geoCore.strutural import Geo_formation as gf 
        >>> geo_structure = gf()
        >>> geo_structure.names
        
    2.  To extract color  and to get the code of structure  like amphibolite 
    
    :Example:
        
        >>> from geodrill.geoCore.strutural import Geo_formation as gf 
        >>> geo_structure = gf() 
        >>> geo_structure.amphibolite['color'] 
        >>> geo_structure.amphibolite['code'] 
        >>> geoformation_obj.AMP['color']
        ... 'R128G128'
    """ 
    codef =['code','label','__description','pattern', 'pat_size',	'pat_density',
            'pat_thickness','color']

    def __init__(self, agso_file =None , **kwargs) :
        self._logging = csamtpylog.get_csamtpy_logger(self.__class__.__name__)
        
        self._agso_fn =agso_file
        
        for key in list(kwargs.keys()): 
            self.__setattr__(key, kwargs[key])
            
        if self.agso_fn is not None : 
            self._set_geo_formation()

    @property 
    def agso_fn (self): 
        if self._agso_fn is None : 
             self._agso_fn = os.path.join(os.environ['pyCSAMT'], 'geodrill', '_geocodes',
                                         'agso'.upper() + '.csv' )
        if self._agso_fn is not None : 
            if os.path.isfile(self._agso_fn) : 
                with open(self._agso_fn, 'r', encoding ='utf8') as _fags0: 
                    agso_lines = _fags0.readlines()
                    _secags0 = agso_lines[0].strip().split(',')
                    _secags0 =[_code.lower() for _code in _secags0]
                    self._AGS0_DATA =agso_lines[1:]
                    
                    if _secags0 != self.codef :
                        warnings.warn(' !Problem to decode geostructures structures!')
                        raise CSex.pyCSAMTError_structural('Geostructures files provided is wrong !')
            else : 
                
                warnings.warn('No Geostructure file detected .{_structures unfound}. It seems the file moved from its host folder.')
                raise CSex.pyCSAMTError_structural('pyCSAMT inner geocodes property is not in host folder {_geocodes_}.'\
                                                   ' Please provide a right structures codes.')
            return self._AGS0_DATA 
    
    def _set_geo_formation(self, agso_fn =None): 
        """
        Read and set  attributes and decode agso geostructures .
        
        :param agso_fn:  full path geological structural  file.
        :type agso_fn: str 
        
        :Example: 
            
            >>> geoformation_obj =Geo_formation()
            >>> DATA = geoformation_obj._AGS0_DATA
            >>> print( geoformation_obj.argillite['name'])
            >>> print( geoformation_obj.argillite['color'])
            >>> print( geoformation_obj.wood['color'])
            >>> print( geoformation_obj.amphibolite['code'])
            >>> print( geoformation_obj.names)
        """
        self._logging.info('Read &  and decodes geostructures files . ')
        if agso_fn is not None : 
            self.agso_fn= agso_fn 
            
        #split the code files 
        
        self._AGS0_DATA = [ geos.strip().split(',') for ii, geos in enumerate(self._AGS0_DATA)]
        self._AGS0_DATA= np.array (self._AGS0_DATA)

        # sanitize the geonames . 
        self.__setattr__('names',
                         np.array( [name.replace('"', '') for name in self._AGS0_DATA[:, 2]] ))
        
        self.__setattr__('codes', self._AGS0_DATA[:, 0])
        
        # set for all values in geofomations codes 
        for ii, codeff in enumerate(self.codef [1:],1) :  # count starting from one .
            if codeff.find('pat') >=0 : 
                for jj, pp in enumerate(self._AGS0_DATA[:, ii]): 
                    if pp == '' or pp ==' ' :       # repalec all None value by 0. later will filled it 
                        self._AGS0_DATA[:, ii][jj]=0.
                        
                self.__setattr__(codeff, 
                                 np.array([float(pp) for pp in self._AGS0_DATA[:, ii]]))
            else : 
               self.__setattr__(codeff,
                         np.array( [name.replace('"', '') for name in self._AGS0_DATA[:, ii]] )) 
        
        DECODE ={}
        # change RGBA color palette into Matplotlib color 
        self.mpl_colorsp= [ get_color_palette(rgb) for rgb in self._AGS0_DATA[:,-1]]

        for jj,  decode  in enumerate(self._AGS0_DATA):
           for ii, codec  in enumerate( self.codef) : 
               if codec =='__description' : codec = 'name'
               if codec =='color': 
                   DECODE[codec]= self.mpl_colorsp[jj]
               else : DECODE [codec]= decode[ii]
               self.__setattr__(decode[2].strip().replace('"', '').replace(' ', '_'),
                                DECODE)
               self.__setattr__(decode[0].strip().replace('"', '').replace(' ', '_'),
                                DECODE)
               
           DECODE ={}  
        
        
def get_color_palette (RGB_color_palette): 
    """
    Convert RGB color into matplotlib color palette. In the RGB color system two bits
    of data are used for each color, red, green, and blue. That means that each color runs
    on a scale from 0 to 255. Black would be 00,00,00, while white would be 255,255,255.
    Matplotlib has lots of pre-defined colormaps for us . They are all normalized to 255,
    so they run from 0 to 1. So you need only normalize data, then we can manually 
    select colors from a color map  

    :param RGB_color_palette: str value of RGB value 
    :type RGB_color_palette: str 
        
    :returns: rgba, tuple of (R, G, B)
    :rtype: tuple
     
    :Example: 
        
        >>> from geodrill.geoCore.structural import get_color_palette 
        >>> get_color_palette (RGB_color_palette ='R128B128')
    """   
    def ascertain_cp (cp): 
        if cp >255. : 
            warnings.warn(' !RGB value is range 0 to 255 pixels , '
                          'not beyond !. Your input values is = {0}.'.format(cp))
            raise CSex.pyCSAMTError_parameter_number('Error color RGBA value ! '
                                                     'RGB value  provided is = {0}.'
                                                     ' It is larger than 255 pixels.'.format(cp))
        return cp
    if isinstance(RGB_color_palette,(float, int, str)): 
        try : 
            float(RGB_color_palette)
        except : 
              RGB_color_palette= RGB_color_palette.lower()
             
        else : return ascertain_cp(float(RGB_color_palette))/255.
    
    rgba = np.zeros((3,))
    
    if 'r' in RGB_color_palette : 
        knae = RGB_color_palette .replace('r', '').replace('g', '/').replace('b', '/').split('/')
        try :
            _knae = ascertain_cp(float(knae[0]))
        except : 
            rgba[0]=1.
        else : rgba [0] = _knae /255.
        
    if 'g' in RGB_color_palette : 
        knae = RGB_color_palette .replace('g', '/').replace('b', '/').split('/')
        try : 
            _knae =ascertain_cp(float(knae[1]))
        except : 
            rgba [1]=1.
            
        else :rgba[1]= _knae /255.
    if 'b' in RGB_color_palette : 
        knae = knae = RGB_color_palette .replace('g', '/').split('/')
        try : 
            _knae =ascertain_cp(float(knae[1]))
        except :
            rgba[2]=1.
        else :rgba[2]= _knae /255.
        

    return tuple(rgba)
    

    
                
if __name__=="__main__":
    
    # structural=Structure()
    # boudin=boudin_axis()
    # # print(boudin.code) 
    # print(structural.boudin_axis.name)
    # print(structural.boudin_axis.color)
    
    geoformation_obj =Geo_formation()
    DATA = geoformation_obj._AGS0_DATA
    rgb=  geoformation_obj.amphibolite['color']
    print(rgb)


    
                
    
    