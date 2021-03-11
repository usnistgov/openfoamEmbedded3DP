import os
import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import seaborn as sns
import statistics as st
from PIL import Image
import interfacemetrics as intm
import scipy as sp
from shapely.geometry import Polygon
import re
import folderparser as fp
import random
import math
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.size'] = 10
from typing import List, Dict, Tuple, Union, Any, TextIO
from interfacemetricsplots import *

###### picture plots

    # picture plot, plots frames at time t, with spacing dx, dy
    
def timestr(t:float) -> str:
    tstr = str(t)
    if t<10:
        tstr = "00" + tstr
    elif t<100:
        tstr = "0" + tstr
    return tstr
    
# gets one picture from a folder
# returns the picture
def picFromFolder(folder:str, t:float, tag:str='y_umag'):
    imfile = os.path.join(folder, 'images', 't'+timestr(int(round(10*t)))+'_'+tag+'.png')
    if os.path.exists(imfile):
        im = Image.open(imfile)
        return im
    else:
        raise ValueError

# plots picture from just one folder
def picPlot(folder:str, cp:comboPlot, t:float, dx:float, cropx:int, cropy:int, tag:str='y_umag') -> None:
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
        
# plots surface velocity legend 
def picPlotLegend(folder:str, cp:comboPlot, t:float):
    try:
        im = picFromFolder(folder, t)
        im = im.crop((100, 1100, 1100, 1200))
    except:
        return
    x0 = np.mean(cp.xlistreal)
    y0 = max(cp.ylistreal)+cp.dy
    dx = cp.dx
    width, height = im.size 
    dy = dx*(height/width)
    cp.axs[0].imshow(im, extent=[x0-dx/2, x0+dx/2, y0-dy/2, y0+dy/2])

# plot all pictures for simulations in a folder
def picPlots(folderlist:List[str], cp:comboPlot, t:float, dx:float, cropx:int, cropy:int, tag:str='y_umag') -> None:
    for folder in folderlist:
        picPlot(folder, cp, t, dx, cropx, cropy, tag=tag)
    cp.figtitle = 't = '+str(t)+' s'
    cp.clean()
    picPlotLegend(folderlist[0], cp, t)

# plot all pictures for simulations in a folder, but use automatic settings for cropping and spacing and export the result
def picPlots0(topFolder:str, exportFolder:str, time:float, sigma:float, tag:str='y_umag', overwrite:bool=False, **kwargs) -> None:
    
    label = 'pics_'+tag+'_t'+str(time)+'_sigma_'+str(sigma)
    fn = intm.imFn(exportFolder, label, topFolder, **kwargs)
    if not overwrite and os.path.exists(fn+'.png'):
        return
    
    flist = intm.caseFolders(topFolder)
    if len(flist)==0:
        return
    
    dx = 0.7
    cp = comboPlot(topFolder, [-dx, dx], [-dx, dx], 6.5, gridlines=False, sigmalist=[sigma], **kwargs)
#     cp.sigmalist = [sigma]
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
    