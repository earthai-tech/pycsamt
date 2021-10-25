# -*- coding: utf-8 -*-
import os 
from pycsamt.modeling.occam2d import getMisfit 


resPath= r'data/inversionFiles'
pathresp =os.path.join(resPath,'K1.resp' )
path_data =os.path.join(resPath,'K1.dat' )
getMisfit(data_fn = path_data, resp_fn = pathresp, 
                                    kind='phase')

# display figure in non-interactive mode
# see  http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.show
import matplotlib as mpl
mpl.pyplot.show() 
