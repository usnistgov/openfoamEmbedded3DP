#!/usr/bin/env python
'''Functions for plotting images of filaments and baths'''

import sys
import os
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging

import interfacemetrics as intm
from plot_general import *

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.size'] = 10

__author__ = "Leanne Friedrich"
__copyright__ = "This data is publicly available according to the NIST statements of copyright, fair use and licensing; see https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software"
__credits__ = ["Leanne Friedrich"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Leanne Friedrich"
__email__ = "Leanne.Friedrich@nist.gov"
__status__ = "Production"

#-------------------------------------------

    
def timestr(t:float) -> str:
    '''convert the time to a 3 digit string'''
    tstr = str(t)
    if t<10:
        tstr = "00" + tstr
    elif t<100:
        tstr = "0" + tstr
    return tstr
    

def picFromFolder(folder:str, t:float, tag:str='y_umag'):
    '''gets one picture from a folder
    returns the picture
    t is the time since extrusion started
    tag is the name of the image type, e.g. 'y_umag'. Used to find images. '''
    imfile = os.path.join(folder, 'images', 't'+timestr(int(round(10*t)))+'_'+tag+'.png')
    if os.path.exists(imfile):
        im = Image.open(imfile)
        return im
    else:
        raise ValueError


def picPlot(folder:str, cp:comboPlot, t:float, dx:float, cropx:int, cropy:int, tag:str='y_umag') -> None:
    '''plots picture from just one folder. 
    folder is the full path name
    cp is the comboPlot object that stores the plot
    t is the time since extrusion started
    dx is the spacing between images in plot space, e.g. 0.7
    cropx is the size to crop from the left and right edges of the piture in px
    cropy is the size to crop from top and bottom in px
    tag is the name of the image type, e.g. 'y_umag'. Used to find images. '''
    try:
        im = picFromFolder(folder, t, tag=tag)
        width, height = im.size 
        im = im.crop((cropx, cropy, width-cropx, height-cropy))
    except:
        return
    width, height = im.size 
    dy = dx*(height/width)
    try:
        color, x0, y0, sigmapos = vvplot(folder, cp)
    except:
        return
    cp.axs[0].imshow(im, extent=[x0-dx/2, x0+dx/2, y0-dy/2, y0+dy/2])
        

def picPlotLegend(folder:str, cp:comboPlot, t:float, tag):
    '''Crops the legend from the image and puts it on the plot.
    folder is a full path name
    cp holds the plot
    t is the time for the picture to use'''
    try:
        im = picFromFolder(folder, t, tag=tag)
        im = im.crop((100, 1100, 1100, 1200))
    except:
        return
    x0 = np.mean(cp.xlistreal)
    y0 = max(cp.ylistreal)+cp.dy
    dx = cp.dx
    width, height = im.size 
    dy = dx*(height/width)
    cp.axs[0].imshow(im, extent=[x0-dx/2, x0+dx/2, y0-dy/2, y0+dy/2])


def picPlots(folderList:List[str], cp:comboPlot, t:float, dx:float, cropx:int, cropy:int, tag:str='y_umag') -> None:
    '''plot all pictures for simulations in a folder
    folderList is a list of paths
    cp holds the plot
    t is the time in s
    dx is the spacing between images in plot space, e.g. 0.7
    cropx is the size to crop from the left and right edges of the piture in px
    cropy is the size to crop from top and bottom in px
    tag is the name of the image type, e.g. 'y_umag'. Used to find images. '''
    for folder in folderList:
        picPlot(folder, cp, t, dx, cropx, cropy, tag=tag)
    cp.figtitle = 't = '+str(t)+' s'
    cp.clean()
    picPlotLegend(folderList[0], cp, t, tag=tag)


def picPlots0(topFolder:str, exportFolder:str, time:float, sigma:float, tag:str='y_umag', overwrite:bool=False, **kwargs) -> None:
    '''plot all pictures for simulations in a folder, but use automatic settings for cropping and spacing and export the result
    topFolder is the folder that holds the simulations
    exportFolder is the folder to export the images to
    time is the time in s since flow started
    sigma is the surface tension
    tag is the name of the image type, e.g. 'y_umag'. Used to find images.
    other kwargs can be used to style the plot
    '''
    
    label = 'pics_'+tag+'_t'+str(time)+'_sigma_'+str(sigma)
    fn = intm.imFn(exportFolder, label, topFolder, **kwargs)
    if not overwrite and os.path.exists(fn+'.png'):
        return
    
    flist = intm.fp.caseFolders(topFolder)
    if len(flist)==0:
        return
    
    dx = 0.7
    cp = comboPlot(topFolder, [-dx, dx], [-dx, dx], 6.5, gridlines=False, sigmalist=[sigma], **kwargs)
    cp.legendList()
    cp.addLegend()
    if tag.startswith('y'):
        cropx = 100
        cropy = 120
        dx = dx*2
    elif tag.startswith('x'):
        cropx = 375
        cropy = 375
        dx = dx*2
    else:
        cropx = 100
        cropy = 120
    picPlots(flist, cp, time, dx, cropx, cropy, tag=tag)
    intm.exportIm(fn, cp.fig)
    