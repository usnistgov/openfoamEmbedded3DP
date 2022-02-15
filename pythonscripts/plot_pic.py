#!/usr/bin/env python
'''Functions for plotting images of filaments and baths'''

# external packages
import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from PIL import Image
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging
import traceback

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import interfacemetrics as intm
from plot_general import *

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib', 'imageio', 'IPython', 'PIL']:
    logging.getLogger(s).setLevel(logging.WARNING)

# plotting
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rc('font', family='Arial')
matplotlib.rc('font', size='10.0')

# info
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
        tstr = f"00{tstr}"
    elif t<100:
        tstr = f"0{tstr}"
    return tstr
    

def picFromFolder(folder:str, t:float, tag:str='y_umag'):
    '''gets one picture from a folder
    returns the picture
    t is the time since extrusion started
    tag is the name of the image type, e.g. 'y_umag'. Used to find images. '''
    ts = timestr(int(10*t))
    imfile = os.path.join(folder, 'images', f't{ts}_{tag}.png')
    if os.path.exists(imfile):
        im = Image.open(imfile)
        return im
    else:
        raise ValueError


def picPlot(folder:str, cp:comboPlot, t:float, dx:float, crops, tag:str='y_umag', **kwargs) -> None:
    '''plots picture from just one folder. 
    folder is the full path name
    cp is the comboPlot object that stores the plot
    t is the time since extrusion started
    dx is the spacing between images in plot space, e.g. 0.7
    crops is a dictionary holding the % of the image to crop to
    tag is the name of the image type, e.g. 'y_umag'. Used to find images. '''
    try:
        im = picFromFolder(folder, t, tag=tag)
        width, height = im.size 
        im = im.crop((int(crops['cropxl']*width), int(height*crops['cropyb']), int(width*crops['cropxr']), int(height*crops['cropyt'])))
        width, height = im.size 
    except:
        return
    ar = height/width
    dy = dx*ar
    try:
        color, x0, y0, sigmapos = vvplot(folder, cp)
    except Exception as e:
        return
    
    cp.axs[0].imshow(im, extent=[x0-dx/2, x0+dx/2, y0-dy/2, y0+dy/2])
        

def picPlotLegend(folder:str, cp:comboPlot, t:float, tag:str):
    '''Crops the legend from the image and puts it on the plot.
    folder is a full path name
    cp holds the plot
    t is the time for the picture to use'''
    try:
        im = picFromFolder(folder, t, tag=tag)
        im = im.crop((im.size[0]*0.1, im.size[0]*0.9, im.size[1]*0.9, im.size[1])) # RG
    except:
        return

    dx = cp.dx*2 # image size RG
    x0 = (max(cp.xmlist)+min(cp.xmlist))/2
    y0 = min(cp.ymlist)-cp.dy/2
    # dx = cp.dx
    width, height = im.size 

    dy = dx*(height/width)
#     x0 = cp.xmlist[-1]/2 # center the legend RG
#     y0 = -0.55 # put the legend at the bottom of the plot RG
    cp.axs[0].imshow(im, extent=[x0-dx/2, x0+dx/2, y0-dy/2, y0+dy/2])
    cp.indicesreal = cp.indicesreal.append({'x':0, 'y':-0.25}, ignore_index=True)


def picPlots(folderList:List[str], cp:comboPlot, t:float, dx:float, crops:dict, tag:str='y_umag',**kwargs) -> None:
    '''plot all pictures for simulations in a folder
    folderList is a list of paths
    cp holds the plot
    t is the time in s
    dx is the spacing between images in plot space, e.g. 0.7
    cropx is the size to crop from the left and right edges of the piture in px
    cropy is the size to crop from top and bottom in px
    tag is the name of the image type, e.g. 'y_umag'. Used to find images. '''
    for folder in folderList:
        picPlot(folder, cp, t, dx, crops, tag=tag, **kwargs)
    cp.figtitle = f't = {t} s'
    picPlotLegend(folderList[0], cp, t, tag=tag)
    cp.clean()


def picPlots0(topFolder:str, exportFolder:str, time:float, sigma:float, tag:str='y_umag', overwrite:bool=False, **kwargs) -> None:
    '''plot all pictures for simulations in a folder, but use automatic settings for cropping and spacing and export the result
    topFolder is the folder that holds the simulations
    exportFolder is the folder to export the images to
    time is the time in s since flow started
    sigma is the surface tension
    tag is the name of the image type, e.g. 'y_umag'. Used to find images.
    other kwargs can be used to style the plot
    '''
    
    label = f'pics_{tag}_t{time}_sigma_{sigma}'
    fn = intm.imFn(exportFolder, label, topFolder, **kwargs)
    if not overwrite and os.path.exists(f'{fn}.png'):
        return
    
    flist = intm.fp.caseFolders(topFolder)
    if len(flist)==0:
        return
    
    
    if tag=='y_shearStressy' or tag=='y_viscy':
        crops = {'cropxl':60/1216, 'cropxr':(1216-600)/1216, 'cropyt':(1216-120)/1216, 'cropyb':240/1216}
    else:
        if tag.startswith('y'):
            cropx = 60/1216 # RG
            cropy = 240/1216
        elif tag.startswith('x'):
            cropx = 350/1216 # RG
            cropy = 375/1216
        else:
            cropx = 100/1216
            cropy = 120/1216
        crops = {'cropxl':cropx, 'cropxr':1-cropx, 'cropyt':1-cropy, 'cropyb':cropy}
        
    dx = 0.5
    ar = (crops['cropyt']-crops['cropyb'])/(crops['cropxr']-crops['cropxl'])
    cp = comboPlot(topFolder, [-dx, dx], [-dx*ar, dx*ar], 6.5, gridlines=False, sigmalist=[sigma], **kwargs)
    dx = dx*2
    if len(cp.flist)==0:
        return
    cp.legendList()
    
    picPlots(flist, cp, time, dx*0.9, crops, tag=tag, **kwargs)
    
#     display(cp.fig)
    cp.addLegend()
    intm.exportIm(fn, cp.fig, **kwargs)
    return cp.fig
    