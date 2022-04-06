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
__credits__ = ["Leanne Friedrich", "Ross Gunther"]
__license__ = "NIST"
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
        im = im.crop((int(crops['cropxl']*width), int(height*crops['cropyt']), int(width*crops['cropxr']), int(height*crops['cropyb']))) # left, top, right, bottom
        width, height = im.size 
    except Exception as e:
        return 0,0
    ar = height/width
    dy = dx*ar
    try:
        color, x0, y0, sigmapos = vvplot(folder, cp)
    except Exception as e:
        return dx,dy

    cp.axs[0].imshow(im, extent=[x0-dx/2, x0+dx/2, y0-dy/2, y0+dy/2], aspect='auto')
    return dx,dy
        

def picPlotLegend(folder:str, cp:comboPlot, t:float, imdy:float, crops:dict, tag:str):
    '''Crops the legend from the image and puts it on the plot.
    folder is a full path name
    cp holds the plot
    t is the time for the picture to use'''
    try:
        im = picFromFolder(folder, t, tag=tag)
        im = im.crop((im.size[0]*0.1, im.size[0]*0.9, im.size[1]*0.9, im.size[1])) # RG
    except Exception as e:
        return

    dx = (max(cp.xmlist)-min(cp.xmlist))+cp.dx*0.9 # width in plot space
    width, height = im.size              # height in px
    dy = dx/(width/height)               # height in plot space
    
    x0 = (max(cp.xmlist)+min(cp.xmlist))/2 # x midpoint
    y0 = min(cp.ymlist)-imdy/2-dy/2-0.1    # y midpoint
    y0ind = y0/cp.dy       # index of y midpoint
    
    cp.legdy = dy

    cp.axs[0].imshow(im, extent=[x0-dx/2, x0+dx/2, y0-dy/2, y0+dy/2], aspect='auto')
    cp.indicesreal = cp.indicesreal.append({'x':0, 'y':y0ind}, ignore_index=True)


def picPlots(folderList:List[str], cp:comboPlot, t:float, dx:float, crops:dict, tag:str='y_umag', legend:bool=True, **kwargs) -> None:
    '''plot all pictures for simulations in a folder
    folderList is a list of paths
    cp holds the plot
    t is the time in s
    dx is the spacing between images in plot space, e.g. 0.7
    cropx is the size to crop from the left and right edges of the piture in px
    cropy is the size to crop from top and bottom in px
    tag is the name of the image type, e.g. 'y_umag'. Used to find images. '''
    imdy0 = 0
    for folder in folderList:
        imdx, imdy = picPlot(folder, cp, t, dx, crops, tag=tag, **kwargs)
        if not imdy==0:
            imdy0 = imdy
    cp.figtitle = f't = {t} s'
    if legend:
        picPlotLegend(folderList[0], cp, t, imdy0, crops, tag=tag)
    cp.clean()

def picCrops(tag:str, **kwargs) -> dict:
    '''get the crop range as a percentage of the image size. cropxl is crop size from left side, cropxr is on right side measured from the left, cropyt is on top measured from top, cropyb is from bottom
    xl = left
    yb = bottom
    xr = right
    yt = top
    '''
    if 'crops' in kwargs:
        crops = kwargs['crops']
    else:
        if tag=='y_shearStressy' or tag=='y_viscy':
            crops = {'cropxl':60/1216, 'cropxr':(1216-600)/1216, 'cropyb':(1216-120)/1216, 'cropyt':240/1216}
        elif tag=='a_stre':
            crops = {'cropxl':600/1216, 'cropxr':(1216-300)/1216, 'cropyb':(550)/1216, 'cropyt':240/1216}
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
            crops = {'cropxl':cropx, 'cropxr':1-cropx, 'cropyb':1-cropy, 'cropyt':cropy}
    return crops
    

def picPlots0(topFolder:str, exportFolder:str, time:float, sigma:float, tag:str='y_umag', overwrite:bool=False, imsize:float=6.5, **kwargs) -> None:
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
    
    
    crops = picCrops(tag, **kwargs)
    if 'crops' in kwargs:
        kwargs.pop('crops')

    ar = (crops['cropyb']-crops['cropyt'])/(crops['cropxr']-crops['cropxl'])
    if ar>1:
        dy = 0.5
        dx = dy/ar
    else:
        dx = 0.5
        dy = dx*ar
    cp = comboPlot(topFolder, [-dx, dx], [-dy, dy], imsize, gridlines=False, sigmalist=[sigma], **kwargs)
    dx = dx*2

    if len(cp.flist)==0:
        return
    cp.legendList()
    
    picPlots(flist, cp, time, dx*0.9, crops, tag=tag, **kwargs)
    
#     display(cp.fig)
    cp.addLegend()
    intm.exportIm(fn, cp.fig, **kwargs)
    return cp.fig
    