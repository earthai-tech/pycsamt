# -*- coding: utf-8 -*-
#   Licence:BSD 3-Clause
#   Author: LKouadio <etanoyau@gmail.com>
#   Created date: Fri Apr 10 08:46:56 2022 

""" 
`WATex`_ Type variables
======================== 

.. |ERP| replace:: Electrical resistivity profiling 

.. _WATex: https://github.com/WEgeophysics/watex/
.. _pandas DataFrame: https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html
.. _Series: https://pandas.pydata.org/docs/reference/api/pandas.Series.html

Some customized type variables  need to be explained for easy understanding 
in the whole package. Indeed, customized type hints is used to define the 
type of arguments. 

**M**: Suppose to be the interger variable `IntVar` to denote the number of 
    rows in the ``Array``. 
    
**N**: Like the ``M``, *N* means the number of column in the ``Array``. It 
    is bound with  integer variable. 
    
**T**: Is known as generic type standing for `Any` type of variable. We keep 
    it unchanged. 

**U**: Unlike `T`, `U` stands for nothing. Use to sepcify the one dimentional 
    array. For instance:: 
        
        >>> import numpy as np 
        >>> array = np.arange(4).shape 
        ... (4, )
        
**S**: Indicates the `Shape` status. It is bound by `M`, `U`, `N`. 'U' stands
    for nothing for one dimensional array. While, the common shape expects 
    for one of two dimensional arrays, it is possible to extend array for 
    more than one dimensional. The class object :class:`AddShape` is 
    created to grand all the remaining value of integers shape. 
    
**D**: Stands for  dtype object. It is bound with  :class:`DType`.

**Array**: Defined for  one dimensional array and `DType` can be specify. For 
    instance, we generated two arrays (`arr1`and `arr2`) for different types:: 
        
        >>> import numpy as np
        >>> from watex.typing import TypeVar, Array, DType
        >>> T = TypeVar ('T', float) 
        >>> A = TypeVar ('A', str, bytes )
        >>> arr1:Array[T, DType[T]] = np.arange(21) # dtype ='float'
        >>> arr2: Array[A, DType[A]] = arr1.astype ('str') # dtype ='str'
        
**NDArray**: Stands for multi-dimensional arrays i.e more than two. Here, the 
    difference between the one dimensional type variable ``Array`` is that 
    while the latter accepts the ``DType`` argument  as the second parameter. 
    It could be turn to the number of multidimentional rows including the 
    `Array as first argument and specify the DType as the second argument 
    like this:: 
        
        >>> import numpy as np 
        >>> from watex.typing import TypeVar, Array, NDarray, DType 
        >>> T =TypeVar ('T', int)
        >>> U = TypeVar ('U')
        >>> multidarray = np.arange(7, 7).astype (np.int32)
        >>> def accept_multid(
                arrays: NDArray[Array[T, U], DType [T]]= multidarray
                ):
            ''' asserted with MyPy and work-fine.'''
                ...
                
**Sub**: Stands for subset. Indeed, the class is created to define the 
    conductive zone. It is a subset ``Sub`` of ``Array``. For example, we first 
    build an array secondly extract the conductive zone from |ERP| line.
    Finally, we checked the type hint to assert whether the extracted zone 
    is a subset of the whole |ERP| line. The demo is given below:: 
        
        >>> import numpy as np 
        >>> from watex.typing import TypeVar, DType, Array , Sub
        >>> from watex.utils.exmath import _define_conductive_zone
        >>> T= TypeVar ('T', float)
        >>> erp_array: Array[T, DType[T]] = np.random.randn (21) # whole line 
        >>> select_zone, _ = _define_conductive_zone (erp = erp_array , auto =True)
        >>> select_zone: Array[T, DType[T]]
        >>> def check_cz (select_zone: Sub[Array]): 
                ''' assert with MyPy and return ``True`` as it works fine. '''
                ... 
                
**SP**: Stands for Station positions. The unit of position may vary, however, 
    we keep for :mod:`watex.method.electrical.ElectricalResistivityProfiling`
    the default unit in ``meters`` by starting at position 0. Typically,
    positions are recording according to the dipole length. For the example, 
    we can generated a position values for ``121 stations`` with dipole 
    length equals to ``50m`` i.e the length of the survey line is ``6 km``. 
    Here we go: 
        
        * Import required modules and generate the whole survey line::
            
            >>> import numpy as np 
            >>> from watex.typing import TypeVar, DType, SP, Sub 
            >>> T =TypeVar ('T', bound =int)
            >>> surveyL:SP = np.arange(0, 50 *121 , 50.).astype (np.int32)
            ... (work fine with MyPy )
            
        * Let's verify whether the extract data from surveyL is also a subset 
            of station positions:
                
            -  We use the following fonction to to extract the specific
                part of whole survey line `surveyL`:: 
                    
                    >>> from watex.utils.exmath import define_conductive_zone
                    >>> subpos,_ = define_conductive_zone (surveyL, s='S10') 
                    
            -  Now, we check the instance value `subpos` as subset array of 
                of `SP`. Note that the station 'S10' is included in the 
                extracted locations and is extented for seven points. For 
                further details, refer to `define_conductive_zone.__doc__`:: 
                
                    >>> def checksup_type (sp: Sub[SP[T, DType[T]]] = subpos ): 
                            ''' SP is an array of positions argument `sp`  
                            shoud be asserted as a subestof the whole line.'''
                            ... 
                    ... (test verified. subpos is a subset of `SP`) 
                    
**Series**: Stands for `pandas Series`_ object rather than using the specific 
    ``pandas.Series`` everywhere in the package. 
    
**DataFrame**: Likewise the ``Series`` generic type hint, it stands for 
    ``pandas DataFrame`_ object. It used to replace ``pandas.DataFrame`` object
    to identify the callable arguments in the whole packages. 
    Both can be instanciated as below:: 
        
        >>> import numpy as np 
        >>> import pandas pd 
        >>> from watex.typing import TypeVar , Any, DType , Series, DataFrame
        >>> T  =TypeVar('T')
        >>> seriesStr = pd.Series ([f'obx{s}' for s in range(21)],
                                 name ='stringobj')
        >>> seriesFloat = pd.Series (np.arange(7).astype(np.float32),
                                 name =floatobj)
        >>> SERs = Series [DType[str]] # pass 
        >>> SERf =Series [DType [float]] # pass 
    
        ..
    
        >>> dfStr= pd.DataFrame {'ser1':seriesStr , 
                            'obj2': [f'none' for i in range (21)]}
        >>> dfFloat= pd.DataFrame {'ser1':seriesFloat , 
                            'obj2': np.linspace (3, 28 , 7)}
        >>> dfAny= pd.DataFrame {'ser1':seriesStr, 
                            'ser2':seriesFloat}
        >>> DFs  = DataFrame [SERs] | DataFrame [DType[str]]
        >>> DFf  = DataFrame [SERf] | DataFrame [DType[float]]
        >>> DFa =  DataFrame [Series[Any]] | DataFrame [DType[T]]
 
**EDIO**: Stands for Electrical Data Interchange (EDI) Object. It is an object 
    built from `pycsamt`_ or `MTpy`_ packages. It holds 'T' or str as type 
    variables supposed to be an object created from the aforementioned packages.
    
---

Additional definition for common arguments 
=========================================== 

To better construct a good API, an explanation of some arguments is useful 
to let the user aware when meeting such argument as a type variable in fromt 
a callable function sperated by a colon. 

**erp** : Stand for Electrical Resistivity Profiling. Typically, the type hint 
    for |ERP| is ``Array[float, DType [float]]`` or ``List[float]``. Its
    array is supposed to hold the apparent resistivy values  collected 
    during the survey. 
    
**p**: Typically mean position but by preference means station location
    positions. The type hint used to defined the `p` is ``
    ``Array[int, DType [int]]`` or ``List[int]``. Indeed, the position 
    supposed to be on integer array and the given values enven in float 
    should be casted to integers. 
     
**cz**: Stands for Conductive Zone. It is a subset of |ERP| so they share the 
    same type hint. However, for better demarcation, ``Sub`` is convenient to 
    use to avoid any confusion about the full |ERP| and the extracted  
    conductive as demontrated in the example above in ``Sub`` type hint
    definition.
        
"""
from __future__ import annotations 
# from typing import TYPE_CHECKING

# if TYPE_CHECKING:
from typing import (
    List,
    Tuple,
    Sequence, 
    Dict, 
    Iterable, 
    Callable, 
    Union, 
    Any , 
    Generic,
    Optional,
    Type , 
    Mapping,
    Text,
    TypeVar, 
    Iterator,
    SupportsInt,

)

__all__=[ 
    "List",
    "Tuple",
    "Sequence", 
    "Dict", 
    "Iterable", 
    "Callable", 
    "Any" , 
    "Generic",
    "Optional",
    "Union",
    "Type" , 
    "Mapping",
    "Text",
    "Shape", 
    "DType", 
    "NDArray", 
    "ArrayLike", 
    "EDIO", 
    "Sub", 
    "SP", 
    "F",
    "T", 
    "V", 
    "Series", 
    "Iterator",
    "SupportsInt",
    ]

T = TypeVar('T')
V = TypeVar('V')
K = TypeVar('K')
M =TypeVar ('M', bound= int ) 
N= TypeVar('N',  bound =int )
U= TypeVar('U')
D =TypeVar ('D', bound ='DType')
S = TypeVar('S', bound='Shape')

class AddShape (Generic [S]): 
    """ Suppose to be an extra bound to top the `Shape` for dimensional 
    more than two. 
    
    Example 
    ------- 
    >>> import numpy as np 
    >>> np.random.randn(7, 3, 3) 
    >>> def check_valid_type (
        array: NDArray [Array[float], Shape[M, AddShape[N]]]): 
        ... 
    
    """
class Shape (Generic[M, S], AddShape[S]): 
    """ Generic to construct a tuple shape for NDarray. `Shape` has is 
    written wait for two dimensional arrays with M-row and N-columns. However 
    for three dimensional,`Optional` Type could be: 
        
    :Example: 
        >>> import numpy as np 
        >>> # For 1D array 
        >>> np
        >>> np.random.rand(7)
        >>> def check_array1d( 
            array: Array[float, Shape[M, None]])
        >>> np.random.rand (7, 7).astype('>U12'):
        >>> def check_array2d_type (
            array: NDArray[Array[str], Shape [M, N], DType ['>U12']])
        
    """
    def __getitem__ (self, M, N) -> S: 
        """ Get the type of rown and type of columns 
        and return Tuple of ``M`` and ``N``. """
        ... 
    
class DType (Generic [T]): 
    """ DType can be Any Type so it holds 'T' type variable. """
    def __getitem__  (self, T) -> T: 
        """ Get Generic Type object and return Type Variable"""
        ...  
       
class ArrayLike(Generic[T, D]): 
    """ Array Type here means the 1D array i.e singular column. 
    For multi-dimensional array we used NDArray instead. 
    """
    
    def __getitem__ (self, T) -> Union ['ArrayLike', T]: 
        """ Return Type of the given Type variable. """ 
        ... 
    
    
class NDArray(ArrayLike[T, DType [T]], Generic [T, D ]) :
    """NDarray has ``M``rows, ``N`` -columns, `Shape` and `DType` object. 
    and Dtype. `Shape` is unbound for this class since it does not make sense
    to specify more integers. However, `DType` seems useful to provide. 
    
    :Example: 
        >>> import numpy as np 
        >>> T= TypeVar (T, str , float) # Dtype here is gone to be "str" 
        >>> array = np.c_[np.arange(7), np.arange(7).astype ('str')]
        >>> def test_array (array: NDArray[T, DType [T]]):...
    """
    def __getitem__ (self,T ) -> T: 
        """ Return type variable. Truly the ``NDArray``"""
        ... 
    
class F (Generic [T]): 
    """ Generic class dedicated for functions, methods and class and 
    return the given types i.e callable object with arguments or `Any`. 
    
    :Example: 
        >>> import functools 
        >>> def decorator (appender ='get only the documention and pass.'):
                @functools.wraps(func):
                def wrapper(*args, **kwds)
                    func.__doc__ = appender + func.__doc__
                    return func (*args, **kwds) 
                return wrapper 
        >>> @decorator  # do_nothing = decorator (anyway)
            def anyway(*args, **kwds):
                ''' Im here to '''
                ...
        >>> def check_F(anyway:F): 
                pass 
    """
    def __getitem__ (self, item: Callable [...,T]
                     ) -> Union ['F', Callable[..., T], T, Any]:
        """ Accept any type of variable supposing to be a callable object 
        functions, methods or even classes and return the given type 
        object or another callable object  with its own or different specific 
        parameters or itself or Any."""
        return self 
    
class Sub (Generic [T]): 
    """ Return subset of whatever Array"""
    ... 
     
class SP(Generic [T, D]): 
    """ Station position arrays hold integer values of the survey location.
    Most likely, the station position is given according to the dipole length.
    Assume the dipole length is ``10 meters`` and survey is carried out on 
    21 stations. The station position array  should be an array of interger 
    values from 0. to 200 meters. as like:: 
        
        >>> import numpy as np 
        >>> positions: SP = np.arange(0, 21 * 10, 10.
                                     ).astype (np.int32) # integer values 
    """
    ... 
    
class Series (DType[T], Generic [T]): 
    """ To reference the pandas `Series`_ object. 
    
    .. _Series: https://pandas.pydata.org/docs/reference/api/pandas.Series.html
    
    :Example:
        >>> import numpy as np
        >>> import pandas as pd 
        >>> from watex.typing import DType, Series  
        >>> ser = pd.Series (np.arange (21), name ='nothing')
        
    .. code: Python 
        
        def check_type (serObj:Series): 
            ''' pass anyway'''
            ... 
        check_type (seObj: Series[DType[str]]=ser ) 
    
    """
    def __getitem__ (self, item: T) -> 'Series': 
        """ Get the type variable of item T and return `Series`_ object."""
        return self 
    
class EDIO(Generic [T]): 
    """ EDIO stand for Electrical Data Interchange (EDI) Object. It is an
    EDI object  built from  :class:`watex.edi.Edi` or  from `pycsamt`_ 
    or `MTpy`_. 
    
    It holds 'T' or str as type variables supposed to be an object created from 
    the aforementioned packages. Indeed the `str` is related to the EDI-file or
    the path-like object to where the EDI file is located. If given,  it assumed
    to be read under the hood and return an EDI object. 
    
    """
    def __getitem__  (self, T: str | T ) -> object: 
        """ Get Generic Type object and return an object presumed to be an
        EDI Object `EDIO`."""
        ...      
        
class ZO(Generic [T]): 
    """ ZO stand for Impedance tensor Object. It is an Impendance tensor object 
    built from :class:`watex.extenals.z.Z` or  :class:`pycsamt.core.z.Z` or 
    :class:`mtpy.core.z.Z`. It is a tridimensional data with 
    dimension equals to (n_freq, 2, 2) where `n_freq` equals to the number of 
    collected frequency and 2x2 matrix referred to components
    xx (0, 0), xy (0, 1), yx (1, 0) and yy (1, 1). 

    """
    def __getitem__  (self, T: str | T ) -> object: 
        """ Get Generic Type object and return an object presumed to be an
        Z Object `ZO`."""
        ...      
class DataFrame (Series[T], Generic[T]): 
    """ Type hint variable to illutsrate the `pandas DataFrame`_ object. 
    
    .. _pandas DataFrame: https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html
    .. _Series: https://pandas.pydata.org/docs/reference/api/pandas.Series.html
    
    Indeed, `pandas DataFrame`_ can be considered as an aggregation of `Series`_, 
    thus, the generic type hint variable is supposed to hold a `Series`_
    object. 
    
    :Example:
        
        >>> import numpy as np
        >>> import pandas as pd 
        >>> from watex.typing import DType, DataFrame 
        
    .. code: Python 
         
        df =pd.DataFrame ({serie1: np.arange(7), 
                           serie2: np.linspace (0, 1000, 7), 
                           serie3: [f'0b{i} for i in range(7)]
                                    })
        def check_type (dfObj:DataFrame): 
            ... 
        ckeck_type (dfObj: DataFrame [DType [object]] =df)
    
    """
    
    def __getitem__(self, item: T)->'DataFrame':
        """ Get the type hint variable of `pandas DataFrame`_ and return the 
        object type variable."""
        
        return self     
    
if __name__=='__main__': 
    def test (array:Sub[SP[ArrayLike[int, DType[int]], DType [int]]]):... 
    def test2 (array:Sub[SP[ArrayLike, DType [int]]]):... 
    
    DFSTR  = DataFrame [Series[DType[str]]]
    DF = DataFrame [DType [object]]
    



























