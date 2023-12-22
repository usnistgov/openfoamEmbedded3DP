#!/usr/bin/env python
'''Functions for plotting standard layouts for adjacent filament simulations'''

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
from plot.super_summary_plot import superSummaryPlot
from plot.trace_plots import tracePlot
from plot.pic_plot import picPlot
from plot.convergence_plot import convergencePlot
from plot.cells_plot import cellsPlot
from plot.colors import sigmaVelocityFunc
from tools.config import cfg

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib', 'imageio', 'IPython', 'PIL']:
    logging.getLogger(s).setLevel(logging.WARNING)


#-------------------------------------------



def sigmaVelocityDict(x2:Any='') -> dict:
    '''get a color dictionary that represents sigma,ink_velocity combinations'''
    if x2=='':
        return dict([[f'{sigma}, {velocity}', sigmaVelocityFunc(sigma,velocity)] for sigma in [0,40] for velocity in [0,10]])
    else:
        return dict([[f'{sigma}, {velocity}, {x2}', sigmaVelocityFunc(sigma,velocity)] for sigma in [0,40] for velocity in [0,10]])

def topFolder(tf:str='adjacent', **kwargs):
    if tf=='adjacent':
        return os.path.join(cfg.path.server, 'adjacent')
    elif tf=='tight spacing':
        return os.path.join(cfg.path.server, 'adjacent', 'tight spacing')
    elif tf=='combined':
        return [os.path.join(cfg.path.server, 'adjacent')
                  , os.path.join(cfg.path.server, 'adjacent', 'tight spacing')]
    else:
        raise ValueError('tf must be adjacent or tight spacing')
        
def exportFolder(ef:str='images', **kwargs):
    if ef=='images':
        return cfg.path.fig
    elif ef=='vidframes':
        return os.path.join(cfg.path.fig, 'vidframes')


def adjXSPlot(*args, x2:float=8, **kwargs) -> folderPlots:
    '''plot cross-sections at t=2.5 and positions ahead and behind the nozzle'''
    return XSPlot(topFolder(**kwargs), exportFolder(**kwargs)
                 , 2.5, [-3, x2]
                 , xunits='niw', plotType='paper'
                  , cvar='sigma,ink_velocity,xbehind'
                , colorDict=sigmaVelocityDict(x2)
               , **kwargs)

def adjTxtPlot(**kwargs) -> folderPlots:
    '''plot folder names'''
    return txtPlot(topFolder(**kwargs), exportFolder(**kwargs)
                , plotType='paper'
               , **kwargs)

def adjPicPlot(tag:str, **kwargs) -> folderPlots:
    '''plot paraview screenshots'''
    return picPlot(topFolder(**kwargs), exportFolder(**kwargs)
                   , tag
                , plotType='paper'
               , **kwargs)

def adjTimePlot(**kwargs) -> folderPlots:
    '''plot simulation progress'''
    return timePlot(topFolder(**kwargs), exportFolder(**kwargs)
                , plotType='paper'
               , **kwargs)

def adjRatePlot(**kwargs) -> folderPlots:
    '''plot simulation rates'''
    return ratePlot(topFolder(**kwargs), exportFolder(**kwargs)
                , plotType='paper'
               , **kwargs)

def adjMetricPlot(var:str, **kwargs) -> folderPlots:
    '''plot a specific metric as a gradient plot'''
    return measurementPlot(topFolder(**kwargs), exportFolder(**kwargs)
                    , var, 2.5, 8, xunits='niw'
                , plotType='paper'
               , **kwargs)

def adjCellsPlot(**kwargs) -> varPlots:
    '''plot the number of cells at the end'''
    return cellsPlot(topFolder(**kwargs), exportFolder(**kwargs), plotType='paper', **kwargs)


def adjSSPlot(vvar, **kwargs) -> varPlots:
    '''plot a specific metric against spacing as a scatter plot'''
    return superSummaryPlot(topFolder(**kwargs)
                , exportFolder(**kwargs)
                , os.path.join(cfg.path.fig, 'adjacent')
                , 'spacing', vvar, 2.5, 8, xunits='niw', cvar='sigma,ink_velocity', mvar='ink_velocity'
                , xticks=[0.5, 0.75, 1, 1.25]
                , splitxvar = 'ink_transportModel,sup_transportModel,adjacent_filament_orientation'
                , splityvar='yvar'
                , restrictions={'ink_transportModel,sup_transportModel':['HerschelBulkley, HerschelBulkley', 'Newtonian, Newtonian']}
                , plotType='paper'
                , markerDict={0:'o', 10:'$8$'}
                , lineDict = {0:'dotted', 10:'solid'}
                , colorDict={'0, 0':'#ad5555', '0, 10':'#852113', '40, 0':'#52abcc', '40, 10':'#1a5f78'}
               , **kwargs)

def adjTracePlot(vvar, **kwargs) -> varPlots:
    '''plot a specific metric against xbehind and spacing'''
    return tracePlot(topFolder(**kwargs), exportFolder(**kwargs), vvar, **kwargs)

def adjConvergencePlot(**kwargs) -> varPlots:
    '''plot the residuals over time'''
    convergencePlot(topFolder(**kwargs), exportFolder(**kwargs), 'ralpha', scaling=10**-8, **kwargs)
    
                                       

#----------------------------------------------------------   

def spacingColorPlot(func, *args, **kwargs) -> varPlots:
    '''plot a specific metric colored by spacing, for matched transportmodels'''
    return func(*args
                       , cvar='spacing'
                       , splitxvar='ink_transportModel,sup_transportModel,adjacent_filament_orientation'
                       , splityvar='sigma,ink_velocity'
                       , restrictions={'ink_transportModel,sup_transportModel':['HerschelBulkley, HerschelBulkley', 'Newtonian, Newtonian']
                                      }
                       , plotType='paper'
                       , **kwargs)

def timeColorPlot(func, sigma, *args, **kwargs) -> varPlots:
    '''plot a specific metric colored by time'''
    return func(*args, cvar='time'
                       , splitxvar='ink_transportModel,sup_transportModel'
                       , splityvar='adjacent_filament_orientation,ink_velocity'
                       , sigma_list=[sigma]
                       , cname='coolwarm'
                       , spacing_list=[0.875]
                       , plotType='paper'
                       , **kwargs)

def spacingSplitPlot(func, model:str, direction:str, *args, **kwargs) -> varPlots:
    '''plot each folder on a different axis, and split the axes by spacing and sigma,ink_velocity'''
    return func(*args
                , cvar='sigma,ink_velocity'
                , colorDict=sigmaVelocityDict()
                , makeLegend=False
                   , splitxvar='spacing'
                   , splityvar='sigma,ink_velocity'
                   , restrictions={'ink_transportModel,sup_transportModel':[f'{model}, {model}']
                                   ,'adjacent_filament_orientation':[direction]
                                  }
                   , plotType='paper'
                   , **kwargs)

def spacingSplitPlots(func, *args, **kwargs) -> List[varPlots]:
    return [spacingSplitPlot(func, model, direction, *args, **kwargs) 
            for model in ['Newtonian', 'HerschelBulkley'] 
            for direction in ['y', 'z']]

def mismatchModelPlot(func, *args, **kwargs) -> varPlots:
    '''a plot where the ink and support model are not the same'''
    return func(*args
                , cvar='sigma,ink_velocity'
                , colorDict=sigmaVelocityDict()
                , makeLegend=False
                   , splitxvar='ink_transportModel,sup_transportModel,adjacent_filament_orientation'
                   , splityvar='sigma,ink_velocity'
                   , restrictions={'ink_transportModel,sup_transportModel':['Newtonian, HerschelBulkley', 'HerschelBulkley, Newtonian']
                                  }
                   , plotType='paper'
                   , **kwargs)
            


def NNPlot(func, *args, **kwargs) -> folderPlots:
    '''Newtonian-Newtonian simulations'''
    return func(*args
                 , xvar='spacing', yvar='sigma,ink_velocity'
                 , ink_transportModel_list=['Newtonian']
                 , sup_transportModel_list=['Newtonian']
                 , splitxvar='adjacent_filament_orientation'
                , makeLegend=False
               , **kwargs)
    
def HBHBPlot(func, *args, **kwargs) -> folderPlots:
    '''HerschelBulkley-HerschelBulkley simulations'''
    return func(*args
                 , xvar='spacing', yvar='sigma,ink_velocity'
                 , ink_transportModel_list=['HerschelBulkley']
                 , sup_transportModel_list=['HerschelBulkley']
                 , splitxvar='adjacent_filament_orientation'
                , makeLegend=False
               , **kwargs)

def NNPlotz(func, *args, **kwargs) -> folderPlots:
    '''Newtonian-Newtonian simulations'''
    return func(*args
                 , xvar='spacing', yvar='sigma,ink_velocity'
                 , ink_transportModel_list=['Newtonian']
                 , sup_transportModel_list=['Newtonian']
                 , adjacent_filament_orientation_list=['z']
                , makeLegend=False
               , **kwargs)

def NNPloty(func, *args, **kwargs) -> folderPlots:
    '''Newtonian-Newtonian simulations'''
    return func(*args
                 , xvar='spacing', yvar='sigma,ink_velocity'
                 , ink_transportModel_list=['Newtonian']
                 , sup_transportModel_list=['Newtonian']
                 , adjacent_filament_orientation_list=['y']
                , makeLegend=False
               , **kwargs)
    
def HBHBPloty(func, *args, **kwargs) -> folderPlots:
    '''HerschelBulkley-HerschelBulkley simulations'''
    return func(*args
                 , xvar='spacing', yvar='sigma,ink_velocity'
                 , ink_transportModel_list=['HerschelBulkley']
                 , sup_transportModel_list=['HerschelBulkley']
                 , adjacent_filament_orientation_list=['y']
                , makeLegend=False
               , **kwargs)

def HBHBPlotz(func, *args, **kwargs) -> folderPlots:
    '''HerschelBulkley-HerschelBulkley simulations'''
    return func(*args
                 , xvar='spacing', yvar='sigma,ink_velocity'
                 , ink_transportModel_list=['HerschelBulkley']
                 , sup_transportModel_list=['HerschelBulkley']
                 , adjacent_filament_orientation_list=['z']
                , makeLegend=False
               , **kwargs)
    
def rhePlot(func, *args, **kwargs) -> folderPlots:
    '''varied rheology simulations'''
    return func(*args
               , xvar='ink_transportModel,sup_transportModel'
               , yvar='sigma,ink_velocity'
                 , spacing_list=[0.875]
                , makeLegend=False
               , **kwargs)

def orientationPlot(func, *args, **kwargs) -> folderPlots:
    '''just one image for disturb/print and for in-plane-out-of-plane'''
    return func(*args
                , xvar='ink_velocity'
                , yvar='adjacent_filament_orientation'
                , spacing_list=[0.875]
                , ink_transportModel_list=['Newtonian']
                , sup_transportModel_list=['Newtonian']
                , **kwargs
               )

def rhePlotWide(func, *args, **kwargs) -> folderPlots:
    '''varied rheology simulations, where orientations are split into different plots horizontally'''
    return rhePlot(func, *args, splitxvar='adjacent_filament_orientation', **kwargs)

def rhePlotTall(func, *args, **kwargs) -> folderPlots:
    '''varied rheology simulations, where orientations are stacked vertically'''
    return rhePlot(func, *args, splityvar='adjacent_filament_orientation', **kwargs)

def doublePlot(func, *args, **kwargs) -> folderPlots:
    return func(*args
                 , xvar='spacing', yvar='sigma,ink_velocity'
                 , restrictions={'ink_transportModel,sup_transportModel':['HerschelBulkley, HerschelBulkley', 'Newtonian, Newtonian']}
                 , splitxvar='adjacent_filament_orientation'
                , splityvar='ink_transportModel,sup_transportModel'
                , makeLegend=False
               , **kwargs)

def triplePlot(func, *args, **kwargs) -> folderPlots:
    '''plot all of the simulations on one plot'''
    return func(*args
                 , xvar='spacing'
                 , yvar='sigma,ink_velocity'
                 , cvar='ink_transportModel,sup_transportModel'
                 , splitxvar='adjacent_filament_orientation'
               , **kwargs)

def tripleD(func, d, *args, **kwargs) -> folderPlots:
    '''plot all of the simulations on one plot for a single displacement direction'''
    return func(*args
                 , xvar='spacing'
                 , yvar='sigma,ink_velocity'
                 , cvar='ink_transportModel,sup_transportModel'
                 , adjacent_filament_orientation_list=[d]
               , **kwargs)

def vidSeries(outerFunc, innerFunc, *args, **kwargs) -> None:
    '''iterate through times and export frames for a video series'''
    for t in [round(0.1*i,1) for i in range(26)]:
        outerFunc(innerFunc, *args, time=t, ef='vidframes', png=True, svg=False, export=True, display=False, overwrite=True, **kwargs)

def tripleDList(func, *args, **kwargs) -> List[folderPlots]:
    '''plot all rheologies on one plot, one plot each for y and z displacement'''
    return [tripleD(func, d, *args, **kwargs) for d in ['y', 'z']]

def threePlots(func, *args, **kwargs) -> List[folderPlots]:
    '''plot the newt-newt, HBHB and rhePlot using the same function and inputs'''
    return [f(func, *args, **kwargs) for f in [NNPlot, HBHBPlot, rhePlotWide]]
    
def fivePlots(func, *args, **kwargs) -> List[folderPlots]:
    '''plot the newt-newt, HBHB and rhePlot using the same funcion and inputs'''
    return [f(func, *args, **kwargs) for f in [NNPlotz, NNPloty, HBHBPlotz, HBHBPloty, rhePlotWide]]
