
from pycsamt.utils._set_gdal import set_GDAL 

HAS_GDAL , EPSG_DICT  = set_GDAL()

__all__=['HAS_GDAL', 'EPSG_DICT']