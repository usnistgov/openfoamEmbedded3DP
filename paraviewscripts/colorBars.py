#!/usr/bin/env pvpython
'''Functions for generating images of filaments from vtk files.'''

# external packages
import os
import numpy as np
import csv
import re
import time
from datetime import datetime
from typing import List, Dict, Tuple, Union, Any, TextIO
from vtkmodules.vtkCommonCore import vtkLogger
vtkLogger.SetStderrVerbosity(vtkLogger.VERBOSITY_OFF)
import logging
import traceback
from paraview.simple import * # import the simple module from the paraview

# local packages
from paraview_general import *
from folder_stats import folderStats

# logging
logger = logging.getLogger(__name__)


#################################################################

FONT = 24

##############################################
######## SCREENSHOTS #########################

def positionCB(ColorBar) -> None:
    '''put the color bar in the bottom'''
    ColorBar.AutoOrient = 0
    ColorBar.Orientation = 'Horizontal'
    ColorBar.WindowLocation = 'LowerCenter'
    ColorBar.ScalarBarLength = 0.7
    ColorBar.UseCustomLabels = 1
    ColorBar.TitleFontSize = FONT
    ColorBar.LabelFontSize = FONT
    ColorBar.ComponentTitle = ''


def uColorBar(renderView1, display, umax:float=0.01, umin:float=0, name:str='|Velocity|') -> None:
    '''velocity magnitude color bar'''
    # set scalar coloring
    #ColorBy(display, ('POINTS', 'U', 'Magnitude'))

    # show color bar/color legend
    display.SetScalarBarVisibility(renderView1, True)
    uLUT = GetColorTransferFunction('U')
    uPWF = GetOpacityTransferFunction('U')

    uLUT.RescaleTransferFunction(umin, umax)
    uPWF.RescaleTransferFunction(umin, umax)
    if (umax-umin)>0.015:
        du = 0.005
    elif (umax-umin)>0.006:
        du = 0.002
    elif (umax-umin)>0.003:
        du = 0.001
    elif (umax-umin)>=0.001:
        du = 0.0003
    else:
        du = 0.0001

    uLUTColorBar = GetScalarBar(uLUT, renderView1)
    positionCB(uLUTColorBar)
    uLUTColorBar.CustomLabels = np.arange(umin, umax+du, du)
    uLUTColorBar.Title = name+' (m/s)'
    uLUTColorBar.RangeLabelFormat = '%-#0.4f'
    return uLUTColorBar

def pColorBar(renderView1, display, pmin:float=-50000, pmax:float=50000, name:str='p_rgh') -> None:
    '''velocity magnitude color bar'''
    # set scalar coloring
    #ColorBy(display, ('POINTS', 'U', 'Magnitude'))

    # show color bar/color legend
    display.SetScalarBarVisibility(renderView1, True)
    pLUT = GetColorTransferFunction(name)
    pPWF = GetOpacityTransferFunction(name)

    pLUT.RescaleTransferFunction(pmin, pmax)
    pPWF.RescaleTransferFunction(pmin, pmax)
    if pmax==50000 and pmin==-50000:
        dp = 25000
    if pmax==100 and pmin==-100:
        dp = 50
    else:
        dp = round((pmax-pmin)/4)

    pLUTColorBar = GetScalarBar(pLUT, renderView1)
    positionCB(pLUTColorBar)
    if name=='p_rgh':
        pLUTColorBar.Title = 'Adjusted pressure (Pa)'
    else:
        pLUTColorBar.Title = 'Pressure (Pa)'
    pLUTColorBar.CustomLabels = np.arange(pmin, pmax+dp, dp)
    return pLUTColorBar

def nuColorBar(renderView1, display, color) -> None:
    '''viscosity color bar'''
    # show color bar/color legend
    display.SetScalarBarVisibility(renderView1, True)
    uLUT = GetColorTransferFunction(color)
    uPWF = GetOpacityTransferFunction(color)

    # rescale to full range of newt, HB experiments
#     uLUT.RescaleTransferFunction(10**-6, 10**2)
#     uPWF.RescaleTransferFunction(10**-6, 10**2)
    emin = -3
    emax = 1
    uLUT.RescaleTransferFunction(10**emin, 10**emax)
    uPWF.RescaleTransferFunction(10**emin, 10**emax)
 
    # convert to log space. If we've already converted to log space during creation of another image, skip this
    if not uLUT.UseLogScale==1:
        uLUT.MapControlPointsToLogSpace()
        uLUT.UseLogScale = 1

    uLUTColorBar = GetScalarBar(uLUT, renderView1)
    positionCB(uLUTColorBar)
    # uLUTColorBar.CustomLabels = [10**i for i in range(-6, 3)]
    uLUTColorBar.CustomLabels = [10**i for i in range(emin, emax+2)]
    uLUTColorBar.Title = 'Viscosity (m^2/s)' # if the density is 1000 kg/m^3, the viscosity is in  MPa.s
    uLUTColorBar.RangeLabelFormat = '%-#2.0e'
    return uLUTColorBar

def viscColor(display, visc:str) -> None:
    '''This is a preset list of colors for viscosity, where the kinematic viscosity range is [10**i for i in range(-6, 3)] m^2/s'''
    # colordict = {0.000001:[59, 76, 192],\
    #              0.00001:[97, 128, 233],\
    #              0.0001:[140, 175, 254],\
    #              0.001:[183, 207, 249], \
    #              0.01:[220, 221, 221], \
    #              0.1:[244, 197, 173], \
    #              1:[244, 153, 122],\
    #              10:[222, 97, 77],\
    #              100:[180, 4, 38]}
    colordict = {10**-3:[59, 76, 192],\
                 10**-2:[140, 175, 254],\
                 10**-1:[220, 221, 221], \
                 10**0:[244, 153, 122],\
                 10:[180, 4, 38]}
    try:
        c = colordict[float(visc)]
    except:
        raise Exception('Viscosity not in list')
    c = [i/255 for i in c]
    display.AmbientColor = c
    display.DiffuseColor = c

def shearRateColorBar(renderView1, display, rmin:float=0.1, rmax:float=1000, name:str='ScalarGradient') -> None:
    '''shear rate color bar'''
    # show color bar/color legend
    display.SetScalarBarVisibility(renderView1, True)
    gLUT = GetColorTransferFunction(name)
    gPWF = GetOpacityTransferFunction(name)

    # rescale to full range of newt, HB experiments
    gLUT.RescaleTransferFunction(rmin, rmax)
    gPWF.RescaleTransferFunction(rmin, rmax)
    
    # convert to log space. If we've already converted to log space during creation of another image, skip this
    if not gLUT.UseLogScale==1:
        gLUT.MapControlPointsToLogSpace()
        gLUT.UseLogScale = 1

    gLUTColorBar = GetScalarBar(gLUT, renderView1)
    positionCB(gLUTColorBar)
    gLUTColorBar.CustomLabels = [10**i for i in range(int(np.ceil(np.log10(rmin))), int(np.floor(np.log10(rmax))))]
    gLUTColorBar.Title = 'Shear rate (1/s)' # if the density is 1000 kg/m^3, the viscosity is in  MPa.s
    gLUTColorBar.RangeLabelFormat = '%-#2.0e'
    return gLUTColorBar

def shearStressColorBar(renderView1, display, rmin:float=1, rmax:float=50, name:str='shearStress') -> None: # RG
    '''shear stress color bar'''
    # show color bar/color legend
    display.SetScalarBarVisibility(renderView1, True)
    gLUT = GetColorTransferFunction(name)
    gPWF = GetOpacityTransferFunction(name)

    # rescale to full range of newt, HB experiments
    gLUT.RescaleTransferFunction(rmin, rmax)
    gPWF.RescaleTransferFunction(rmin, rmax)
    
    # convert to log space. If we've already converted to log space during creation of another image, skip this
    if not gLUT.UseLogScale==1:
        gLUT.MapControlPointsToLogSpace()
        gLUT.UseLogScale = 1

    gLUTColorBar = GetScalarBar(gLUT, renderView1)
    positionCB(gLUTColorBar)
    gLUTColorBar.CustomLabels = [1, 5.0, 10.0, 50.0, 100]
    gLUTColorBar.Title = 'Shear stress (Pa)'
    gLUTColorBar.RangeLabelFormat = '%-#2.0f'
    return gLUTColorBar
    

def alphaColorBar(renderView1, display) -> None:
    '''volume fraction color bar'''
    # set scalar coloring
    ColorBy(display, ('CELLS', 'alpha.ink'))
    ColorBy(display, ('POINTS', 'alpha.ink'))

    # rescale color and/or opacity maps used to include current data range
    display.RescaleTransferFunctionToDataRange(True, False)
    display.SetScalarBarVisibility(renderView1, True)
    alphainkLUT = GetColorTransferFunction('alphaink')
    alphainkLUT.RescaleTransferFunction(0.0, 1)
    alphainkPWF = GetOpacityTransferFunction('alphaink')

    aLUTColorBar = GetScalarBar(alphainkLUT, renderView1)
    positionCB(aLUTColorBar)
    aLUTColorBar.CustomLabels = np.arange(0, 1.25, 0.25)
    aLUTColorBar.Title = 'Volume fraction ink'
    aLUTColorBar.RangeLabelFormat = '%-#0.2f'
    return aLUTColorBar