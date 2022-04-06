#!/usr/bin/env python
'''Functions for putting figure labels on figures'''

# external packages
import matplotlib.pyplot as plt
from typing import List, Dict, Tuple, Union, Any
import logging
import string

# local packages


# plot settings
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams.update({'font.size': 10})
plt.rcParams['svg.fonttype'] = 'none'

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# info
__author__ = "Leanne Friedrich"
__copyright__ = "This data is publicly available according to the NIST statements of copyright, fair use and licensing; see https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software"
__credits__ = ["Leanne Friedrich"]
__license__ = "NIST"
__version__ = "1.0.0"
__maintainer__ = "Leanne Friedrich"
__email__ = "Leanne.Friedrich@nist.gov"
__status__ = "Development"

#----------------------------------------------

    
def subFigureLabel(ax, label:str, inside:bool=True) -> None:
    '''add a subfigure label to the top left corner'''
    if inside:
        x=0.05
        y = 0.95
        ha = 'left'
        va = 'top'
    else:
        x = -0.05
        y = 1.05
        ha = 'right'
        va = 'bottom'
    ax.text(x, y, label, fontsize=12, transform=ax.transAxes, horizontalalignment=ha, verticalalignment=va)
    
def subFigureLabels(axs, horiz:bool=True, inside:bool=True) -> None:
    '''add subfigure labels to all axes'''
    alphabet_string = string.ascii_uppercase
    alphabet_list = list(alphabet_string)
    if len(axs.shape)==1:
        # single row
        for ax in axs:
            subFigureLabel(ax, alphabet_list.pop(0), inside=inside)
    else:
        if horiz:
            # 2d array
            for axrow in axs:
                for ax in axrow:
                    subFigureLabel(ax, alphabet_list.pop(0), inside=inside)
        else:
            w = len(axs[0])
            h = len(axs)
            for i in range(w):
                for j in range(h):
                    subFigureLabel(axs[j][i], alphabet_list.pop(0), inside=inside)