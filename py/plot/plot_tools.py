#!/usr/bin/env python
'''General plotting tools'''

# external packages
import sys
import os
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
from typing import List, Dict, Tuple, Union, Any
import logging
import traceback

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)



#################################################################

def plotCircle(ax:plt.Axes, x0:float, y0:float, radius:float, caption:str, color, sigmapos:int=0) -> None:
    '''plotCircle plots a circle, with optional caption.
    ax is axis to plot on
    x0 is the x position
    y0 is the y position
    radius is the radius of the circle, in plot coords
    caption is the label. use '' to have no caption
    color is the color of the circle
    sigmapos is the position in the sigma list for this circle. Useful for timeplots, so labels don't stack on top of each other. If we're using this to plot an ideal cross-section or some other sigma-unaffiliated value, sigmapos=0 will put the label inside the circle or just above it.'''
    
    circle = plt.Circle([x0,y0], radius, color=color, fill=False) # create the circle
    ax.add_artist(circle)                                         # put the circle on the plot

    if len(caption)>0:
        if radius>0.3:
            # if the circle is big, put the label inside
            txtx = x0                         
            txty = y0+0.2*sigmapos
            ax.text(txtx, txty, caption, horizontalalignment='center', verticalalignment='center', color=color)
        else:
            # if the circle is small, put the label outside
            angle=(90-30*sigmapos)*(2*np.pi/360)  # angle to put the label at, in rad                            
            dar = 0.2                             # length of arrow
            arrowx = x0+radius*np.cos(angle)      # where the arrow points   
            arrowy = y0+radius*np.sin(angle)
            txtx = arrowx+dar*np.cos(angle)       # where the label is
            txty = arrowy+dar*np.sin(angle)
            ax.annotate(caption, (arrowx, arrowy), color=color, xytext=(txtx, txty), ha='center', arrowprops={'arrowstyle':'->', 'color':color})

            
            
def setSquare(ax):      
    '''make the axis square'''
    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')