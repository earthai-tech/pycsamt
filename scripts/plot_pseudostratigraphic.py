# -*- coding: utf-8 -*-
"""
Created on Sat Oct 30 21:38:24 2021
    Script to plot the pseudostratigraphic  log
    Use the script `plot_geostratigraphy_model.py` to run your model first 
    then come to visualize the pseudostratigraphic log at each station.
    
    Run the script `plot_geostratigraphy_model.py` only one time 
     is sufficient to visualize the station of all block model.
"""

from pycsamt.geodrill.geocore import GeoStratigraphy

zoom =[60, 120]#'25' #  '25%' or [60, 700] --> [top bottom] list value
station ='S00'
annotate_kws = {'fontsize':12}
GeoStratigraphy.plotPseudostratigraphic(station =station, zoom =zoom, 
                                        annotate_kws =annotate_kws )


