# -*- coding: utf-8 -*-


from pycsamt.modeling.occam2d import plotResponse 

# path to your Occam inversion inversion files  composed of 
resPath='data/inversionFiles'
plotResponse(data_fn =resPath,
                stations = ['S00', 'S04'],# 's08', 'S12'],  # sites to visualize 
                 rms =['1.013', '1.451'],# '1.00', '1.069'], # rms of each line
                  error_type ='resi' )

# display figure in non-interactive mode
# see  http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.show
import matplotlib as mpl
mpl.pyplot.show() 