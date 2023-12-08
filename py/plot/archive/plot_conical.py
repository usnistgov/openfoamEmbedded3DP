#!/usr/bin/env python
'''Functions for plotting standard layouts for conical nozzle simulations'''

# external packages
import sys
import os
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging
import traceback

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from plot.var_plots import varPlots
from plot.folder_plots import folderPlots
from plot.xs_plot import XSPlot
from plot.txt_plot import txtPlot
from plot.time_plot import timePlot
from plot.rate_plot import ratePlot
from plot.measurement_plot import measurementPlot
from tools.config import cfg

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib', 'imageio', 'IPython', 'PIL']:
    logging.getLogger(s).setLevel(logging.WARNING)


#-------------------------------------------

def qualityPlots0(topFolder:str, exportFolder) -> varPlots:
    return superSummaryPlot(os.path.join(cfg.path.server, 'conical')
                , os.path.join(cfg.path.fig, 'conical')
                 , 'nozzle_angle', ['arean', 'vertdispn', 'aspectratio', 'speeddecay']
                            , 2.5, 8, xunits='niw'
                 , splitxvar = 'yvar'
               , **kwargs)
