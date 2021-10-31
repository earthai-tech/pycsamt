# -*- coding: utf-8 -*-
import os 
from pycsamt.modeling.occam2d import getMisfit 


resPath= r'data/occam2D'
pathresp =os.path.join(resPath,'RESP17.resp' )
path_data =os.path.join(resPath,'OccamDataFile.dat' )
getMisfit(data_fn = path_data, resp_fn = pathresp, 
                                    kind='phase')

# display figure in non-interactive mode
# see  http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.show
import matplotlib as mpl
mpl.pyplot.show() 
