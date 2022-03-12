# -*- coding: utf-8 -*-
#       author:K.L ~ @Daniel03


import os 
from pycsamt.viewer.plot import plot_dataAndFits

#path to avg files 
path ='data/avg'

pathData = [os.path.join(path, file) 
            for file in ['K1.AVG','K2.AVG',
                          ]]

plot_dataAndFits(data_fn = pathData, stations=['S00', 'S04', 's08', 'S12']
                                      )

# display figure in non-interactive mode
# see  http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.show
import matplotlib as mpl
mpl.pyplot.show() 