#!/usr/bin/env python
'''Functions for plotting paraview screenshots of simulations'''

# external packages
import sys
import os
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging
import traceback
from PIL import Image

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from plot.combo_plot import comboPlot
from folder_stats import folderStats


# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib', 'imageio', 'IPython', 'PIL']:
    logging.getLogger(s).setLevel(logging.WARNING)


#-------------------------------------------

class picPlot(comboPlot):
    '''plot screenshots
    topFolder is the folder that holds all the files
    overwrite True to overwrite plots
    '''
    
    def __init__(self, topFolder:str
                 , exportFolder:str
                 , tag:str
                 , time:float = 2.5
                 , **kwargs):
        self.tag = tag
        self.time = time
        self.topFolder = topFolder
        self.getCrops(tag, **kwargs)   # get the area to crop
        self.getxryr()   # determine shape of plot
        super().__init__(topFolder, exportFolder=exportFolder, xr=self.xr, yr=self.yr, gridlines=False, insideLabels=False, makeFrame=False, **kwargs)
        self.getFN(titleVars={'time_list':[f'{time:0.1f} s']}, **kwargs)
        if not self.checkOverwrite(export=self.export, overwrite=self.overwrite):
            return
        
        self.legendPlotted = False
        self.getTimestr()   # get the time to search for
        # iterate through folders and plot data
        for i,row in self.filedf.iterrows():
            self.plotFolder(row)
            
        self.clean()
            
        if self.export:
            self.exportIm(**kwargs) 
            
        #-------------------
        
    def getFileLabel(self) -> str:
        '''get a label for exporting file names'''
        self.label = f'pic_{self.tag}_{self.time}'
        
    def getCrops(self, tag:str, **kwargs) -> dict:
        '''get the crop range as a fraction of the image size. cropxl is crop size from left side, cropxr is on right side measured from the left, cropyt is on top measured from top, cropyb is from bottom
        xl = left
        yb = bottom
        xr = right
        yt = top
        '''
        if 'crops' in kwargs:
            crops = kwargs.pop('crops')
        else:
            if (tag=='y_shearStressy' or tag=='y_viscy') and 'conical' in self.topFolder:
    #             crops = {'cropxl':60/1216, 'cropxr':(1216-600)/1216, 'cropyb':(1216-120)/1216, 'cropyt':240/1216}
                crops = {'cropxl':240/1216, 'cropxr':(1216-740)/1216, 'cropyb':(1216-600)/1216, 'cropyt':380/1216}
            elif tag=='a_stre':
                crops = {'cropxl':600/1216, 'cropxr':(1216-300)/1216, 'cropyb':(550)/1216, 'cropyt':240/1216}
            else:
                if tag.startswith('y'):
                    cropx = 60/1216 # RG
                    cropy = 350/1216
                elif tag.startswith('x'):
                    cropx = 350/1216 # RG
                    cropy = 375/1216
                elif tag.startswith('a'):
                    cropx = 300/1216
                    cropy = 220/1216
                else:
                    cropx = 100/1216
                    cropy = 120/1216
                crops = {'cropxl':cropx, 'cropxr':1-cropx, 'cropyb':1-cropy, 'cropyt':cropy}
        self.crops = crops
        
    def getxryr(self):
        '''get the xrange and yrange for initializing the plot'''
        ar = (self.crops['cropyb']-self.crops['cropyt'])/(self.crops['cropxr']-self.crops['cropxl'])
        if ar>1:
            dy = 0.5
            dx = dy/ar
        else:
            dx = 0.5
            dy = dx*ar
        self.dx = dx*2
        self.dy = dy*2
        self.xr = [-dx, dx]
        self.yr = [-dy, dy]
        
    #-----------------------------------------------------
        
    def getTimestr(self) -> str:
        '''convert the time to a 3 digit string'''
        tstr = str(int(self.time*10))
        tstr = tstr.zfill(3)
        self.tstr = tstr


    def picFromFolder(self, folder:str):
        '''gets one picture from a folder
        returns the picture
        t is the time since extrusion started
        tag is the name of the image type, e.g. 'y_umag'. Used to find images. '''
        imfile = os.path.join(folder, 'images', f't{self.tstr}_{self.tag}.png')
        if os.path.exists(imfile):
            im = Image.open(imfile)
            return im
        else:
            print(imfile)
            return []
            
    def cropIm(self, im):
        '''crop the image'''
        width, height = im.size 
        xl = self.crops['cropxl']
        yt = self.crops['cropyt']
        xr = self.crops['cropxr']
        yb = self.crops['cropyb']
        im = im.crop((int(xl*width), int(height*yt), int(width*xr), int(height*yb))) # left, top, right, bottom
        return im
    
    def plotLegend(self, im, ax):
        '''crop the legend from the image and put it at the bottom'''
        im = im.crop((im.size[0]*0.1, im.size[0]*0.9, im.size[1]*0.9, im.size[1])) # RG
        dx = (max(self.xmlist)-min(self.xmlist))+self.dx*0.9 # width in plot space
        width, height = im.size              # height in px
        dy = dx/(width/height)               # height in plot space

        x0 = (max(self.xmlist)+min(self.xmlist))/2 # x midpoint
        y0 = min(self.ymlist)-self.imdy/2-dy/2-0.1    # y midpoint
        y0ind = y0/self.dy       # index of y midpoint. can be fractional

        self.legdy = dy

        ax.imshow(im, extent=[x0-dx/2, x0+dx/2, y0-dy/2, y0+dy/2], aspect='auto')
        self.indicesreal = self.indicesreal.append({'x':0, 'y':y0ind}, ignore_index=True)

    def plotFolder(self, row) -> None:
        '''given a row in the pandas dataframe, plot the slices'''
        # identify if we need to add a new ideal plot
        folder = row['folder']
        im0 = self.picFromFolder(folder)
        if im0==[]:
            return
        
        im = self.cropIm(im0)
        width, height = im.size 
        ar = height/width
        dx = self.dx*0.95
        dy = dx*ar
        pos = self.getXYRow(row)
        xm = pos['x0']
        ym = pos['y0']
        pos['ax'].imshow(im, extent=[xm-dx/2, xm+dx/2, ym-dy/2, ym+dy/2], aspect='auto')
        self.imdy = dy
        
        if (not self.legendPlotted) and self.makeLegend:
            self.plotLegend(im0, pos['ax'])
        