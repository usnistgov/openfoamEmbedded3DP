#!/usr/bin/env python
'''Functions for combining images into videos'''

# external packages
import os, sys
import glob
from IPython.display import Image
import imageio
import subprocess, json
import re
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging
import cv2 as cv

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(currentdir)
sys.path.append(parentdir)
import file.file_handling as fh
import folder_scraper as fs

# plot settings
plt.rcParams.update({'font.size': 12})

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib', 'imageio', 'IPython', 'PIL']:
    logging.getLogger(s).setLevel(logging.WARNING)


#----------------------------------------------


def get_length(filename:str) -> float:
    '''Get the length of the video'''
    rl = ["ffprobe", "-v", "error", "-show_entries", "format=duration", "-of", "default=noprint_wrappers=1:nokey=1", filename]
    result = subprocess.run(rl, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    duration = float(result.stdout)
    return duration



def saveVid(folder:str, s:str, p:str, diag:bool=True) -> None:
    '''compile images into a video, where images contain the string {s}_{p}.png'''
    zfiles = glob.glob(os.path.join(folder,'images', f'*{s}_{p}.png'))
    if len(zfiles)>0:
        fn = os.path.join(folder,'images',  f'{s}_{p}.mp4')
        if os.path.exists(fn):
            fnlength = get_length(fn)
            times = [get_time(f) for f in zfiles]
            if fnlength < max(times)-0.1:
                if diag:
                    logging.info(f'remove {fn}, vid length = {fnlength}, sim time = {max(times)-0.1}')
                os.remove(fn)
        if not os.path.exists(fn):
            if diag:
                logging.info(f'write {fn}')
            mp4writer = imageio.get_writer(fn, fps=10)
            for filename in zfiles:
                image = imageio.imread(filename)[:,:,:3] # ignore alpha channel
                mp4writer.append_data(image)
            mp4writer.close()
            
            
class bigVideo:
    '''Compile all of the time series for all of the simulations into one big video. 
    folderList is a list of the folders to include. 
    filename is the name of the video to save. 
    titleLength is the time that the title cards are up, in s'''
    
    def __init__(self, folderList:str, filename:str, tags:List[str]=['y_umag','a_umag'], titleLength:float=1, fps:int=10, imShape:tuple=(1216, 1216, 3), diag:bool=True):
        self.folderList = folderList
        self.filename = filename
        self.tags=tags
        self.titleLength = titleLength
        self.diag = diag
        self.fps = fps
        self.imShape = imShape
        self.initializeWriter()
        try:
            self.writeFolders()
        except Exception as e:
            self.closeWriter()
            raise e
        if self.diag:
            logging.info(f'Done creating {self.filename}')
        self.closeWriter()
        
        
    def initializeWriter(self):
        self.mp4writer = imageio.get_writer(self.filename, fps=self.fps, quality=5, pixelformat='yuvj444p')
        
    def closeWriter(self):
        self.mp4writer.close()
        
    def addFrame(self, image):
        '''add a single frame to the video'''
        imshape = image.shape
        if not imshape==self.imShape:
            oldimshape = image.shape
            xfrac = self.imShape[0]/imshape[0]
            yfrac = self.imShape[1]/imshape[1]
            if xfrac>yfrac:
                image = cv.resize(image, (int(imshape[1]*yfrac), int(imshape[0]*yfrac)))
            else:
                image = cv.resize(image, (int(imshape[1]*xfrac), int(imshape[0]*xfrac)))
            imshape = image.shape
            dx = self.imShape[0]-imshape[0]
            dy = self.imShape[1]-imshape[1]
            color = [254, 254, 254]
            image = cv.copyMakeBorder(image, int(np.floor(dx/2)), int(np.ceil(dx/2)), int(np.floor(dy/2)), int(np.ceil(dy/2)), cv.BORDER_CONSTANT, value=color)
        self.mp4writer.append_data(image)
        
    def saveTitleCardToWriter(self, folder:str, subfolder:str='images') -> None:
        '''For a given folder, save the title card to the videowriter. n is number of frames of the title card to add'''
        n = int(round(self.titleLength*self.fps))
        if len(subfolder)>0:
            file = os.path.join(folder, 'images', 'titleCard.png')
        else:
            file = os.path.join(folder,  'titleCard.png')
        if os.path.exists(file):
            image = imageio.imread(file)[:,:,:3] # ignore alpha channel
            for i in range(n):
                self.addFrame(image)
        else:
            raise NameError(f'No title card in {folder}') 
            
    def get_time(self, filename:str) -> float:
        '''Get the time of the simulation timepoint. if there is a decimal, this time is already in seconds. Otherwise, it is in deciseconds.'''
    #     f = re.split('t|_', os.path.basename(filename))[1]
        spl = re.split('_', os.path.basename(filename))
        f = spl.pop(0)
        while not f[0]=='t' and len(spl)>0:
            f = spl.pop(0)
        if len(spl)==0:
            raise NameError(f'No time in file {filename}')
        f = f[1:]
        if '.' in f:
            return float(f)
        else:
            return float(f)/10
            
    def saveFramesToWriter(self, folder:str, tag:str, subfolder:str='images', debug:bool=False) -> None:
        '''For a given folder, save frames to an imageio writer
        tags should be a list for each folder, and within each list, a list of tags that should be in every file. For example, to get all sigma=0 and then all sigma=40, tags could be [['sigma_0', 'y_umag'],['sigma_40', 'y_umag']]. 
        subfolder is the folder inside of the simulation folder where the images are stored'''
        zfiles = []
        if len(subfolder)>0:
            searchFolder = os.path.join(folder, subfolder)
        else:
            searchFolder = folder
        for file in os.listdir(searchFolder):
            if 'png' in file and tag in file:
                f = os.path.join(searchFolder, file)
                zfiles.append({'file':f, 'time':self.get_time(f)})
        df = pd.DataFrame(zfiles)
        df.sort_values(by='time', inplace=True)
        df.reset_index(drop=True, inplace=True)
        
        for i, filename in enumerate(df.file):
            image = imageio.imread(filename)[:,:,:3] # ignore alpha channel
            self.addFrame(image)
    
    def writeFolders(self):
        '''write all the folders'''
        for folder in self.folderList:
            fhi = fh.folderHandler(folder)
            if self.diag:
                logging.info(f'Writing {fhi.shortName()}')
            self.saveTitleCardToWriter(folder)
            for tag in self.tags:
                self.saveFramesToWriter(folder, tag)

#------------------------

class bigVideoFiles(bigVideo):
    '''Compile all of the time series for all of the simulations into one big video. 
    folderList is a list of the folders to include. 
    filename is the name of the video to save. 
    titleLength is the time that the title cards are up, in s'''
    
    def __init__(self, folderList:str, filename:str, tags:List[str]=['y_umag','a_umag'], imShape:tuple=(1216, 1216, 3), **kwargs):
        super().__init__(folderList, filename, tags=tags, imShape=imShape, **kwargs)
    
    def get_time(self, filename:str) -> float:
        '''Get the time of the simulation timepoint. if there is a decimal, this time is already in seconds. Otherwise, it is in deciseconds.'''
    #     f = re.split('t|_', os.path.basename(filename))[1]
        spl = re.split('_', os.path.basename(filename))
        f = spl.pop(0)
        while not f[0]=='t' and len(spl)>0:
            f = spl.pop(0)
        if len(spl)==0:
            raise NameError(f'No time in file {filename}')
        f = f[1:]
        if '.' in f:
            return float(f)
        else:
            return float(f)/10
            
    def saveFramesToWriter(self, folder:str, tag:str, subfolder:str='images', debug:bool=False) -> None:
        '''For a given folder, save frames to an imageio writer
        tags should be a list for each folder, and within each list, a list of tags that should be in every file. For example, to get all sigma=0 and then all sigma=40, tags could be [['sigma_0', 'y_umag'],['sigma_40', 'y_umag']]. 
        subfolder is the folder inside of the simulation folder where the images are stored'''
        zfiles = []
        if len(subfolder)>0:
            searchFolder = os.path.join(folder, subfolder)
        else:
            searchFolder = folder
        for file in os.listdir(searchFolder):
            if 'png' in file and tag in file:
                f = os.path.join(searchFolder, file)
                zfiles.append({'file':f, 'time':self.get_time(f)})
        df = pd.DataFrame(zfiles)
        df.sort_values(by='time', inplace=True)
        df.reset_index(drop=True, inplace=True)
        
        for i, filename in enumerate(df.file):
            image = imageio.imread(filename)[:,:,:3] # ignore alpha channel
            self.addFrame(image)
    
    def writeFolders(self):
        '''write all the folders'''
        for folder in self.folderList:
            fhi = fh.folderHandler(folder)
            if self.diag:
                logging.info(f'Writing {fhi.shortName()}')
            self.saveTitleCardToWriter(folder)
            for tag in self.tags:
                self.saveFramesToWriter(folder, tag)             

#------------------------                
                
class bigVideoFigure(bigVideo):
    '''Compile time series figures into one video 
    folderList is a list of the folders to include. 
    filename is the name of the video to save. 
    titleLength is the time that the title cards are up, in s'''
    
    def __init__(self, folder:str, filename:str
                 , ftags:List[str]=['x_ink_tm,sup_tm_y_sigma,iv_spacing_0.875_splitx_Ddir'
                                   ,'x_spacing_y_sigma,iv_ink_tm_HB_sup_tm_HB_Ddir_y'
                                   ,'x_spacing_y_sigma,iv_ink_tm_HB_sup_tm_HB_Ddir_z'
                                   ,'x_spacing_y_sigma,iv_ink_tm_N_sup_tm_N_Ddir_y'
                                   ,'x_spacing_y_sigma,iv_ink_tm_N_sup_tm_N_Ddir_z']
                 , tags:List[str]=['y_umag','a_umag'], imShape:tuple=(1792,2000, 3), **kwargs):
        self.folder = folder
        self.ftags = ftags
        super().__init__([folder], filename, tags=tags, imShape=imShape, **kwargs)
            
    def get_time(self, filename:str) -> float:
        '''Get the time of the simulation timepoint. if there is a decimal, this time is already in seconds. Otherwise, it is in deciseconds.'''
    #     f = re.split('t|_', os.path.basename(filename))[1]
        spl = re.split('_', os.path.basename(filename))
        f = spl.pop(0)
        while not '.' in f and len(spl)>0:
            f = spl.pop(0)
        if len(spl)==0:
            raise NameError(f'No time in file {filename}')
        return float(f)
    
    def addFrameToList(self, file:str, zfiles:list) -> list:
        '''add the frame to the list if it belongs in our filters'''
        if not '.png' in file:
            return zfiles
        d = {'file':file, 'time':self.get_time(file)}
        for tag in self.tags:
            if tag in file:
                d['tag'] = tag
        if not 'tag' in d:
            return zfiles
        for ftag in self.ftags:
            if ftag in file:
                d['ftag'] = ftag
        if 'ftag' in d:
            zfiles.append(d)
        return zfiles
                
            
    def saveFramesToWriter(self) -> None:
        '''For a given folder, save frames to an imageio writer
        tags should be a list for each folder, and within each list, a list of tags that should be in every file. For example, to get all sigma=0 and then all sigma=40, tags could be [['sigma_0', 'y_umag'],['sigma_40', 'y_umag']]. 
        subfolder is the folder inside of the simulation folder where the images are stored'''
        zfiles = []
        for file in os.listdir(self.folder):
            zfiles = self.addFrameToList(os.path.join(self.folder, file), zfiles)
        df = pd.DataFrame(zfiles)
        df.sort_values(by=['tag', 'ftag', 'time'], inplace=True)
        df.reset_index(drop=True, inplace=True)
        for i, filename in enumerate(df.file):
            image = imageio.imread(filename)[:,:,:3] # ignore alpha channel
            self.addFrame(image)
    
    def writeFolders(self):
        '''write all the folders'''
        self.saveFramesToWriter()
    

    
            
#--------------------------------
            
def exp(n:float) -> str:
    '''Format a power of 10 in exponential notation'''
    return r'$10^{{{}}}$'.format(int(round(np.log10(n))))


#--------------------------------

def checkSimRate(folder:str, fix:bool=True) -> str:
    '''Check if the simulation rate makes sense. fix=True to rescrape the legend if the simulation rate doesn't make sense'''
    leg = fs.legendUnique(folder) # legend
    shortname = os.path.join(os.path.basename(os.path.dirname(folder)), os.path.basename(folder))
    # simrate = leg['simulation_rate_(hr/s)']
    simrate = leg['simulation_rate'] # RG
    warn = False
    try:
        simrate = '{:.1f}'.format(float(simrate))
    except:
        warn = True
        pass
    if simrate==0:
        warn = True
    if warn:
        logging.warning(f'{shortname}: Simulation rate is {simrate}')
        if fix:
            logging.info(f'{shortname}: Repopulating legend')
            fs.populate(folder)
            return checkSimRate(folder, fix=False)
    return simrate

def titleCard(folder:str, overwrite:bool=False, diag:bool=True) -> None: 
    '''Create and export a title card using the legend for the simulation'''
    fname = os.path.join(folder, 'images', 'titleCard.png')
    shortname = fh.folderHandler(folder).shortName()
    if os.path.exists(fname) and not overwrite:
        return
    
    # mydpi = 144
    mydpi = 96 # RG
    gammadot = [10**n for n in np.arange(-6, 6, 0.1)]
    fig,ax = plt.subplots(1,1, figsize=(3,3))
    ax.set_ylabel(r'Viscosity $\eta$ (Pa.s)', fontname="Arial")
    ax.set_xlabel(r'Strain rate $\dot{\gamma}$ (s$^{-1}$)', fontname="Arial")
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylim([10**-3, 10**6])
    ax.set_xlim([10**-6, 10**6])
    
    colors = ['#263e59', '#cc3d3d']
    leg = fs.legendUnique(folder) # legend
    simrate = checkSimRate(folder, True)
    rows = [[],[]]
    
    if os.path.basename(folder)[:2]=='aj': # RG
        ormap = {'y': 'Horizontal', 'z': 'Vertical'}
        plotvars = [['Simulation number', os.path.basename(shortname)],\
                    ['Corresponding simulation', leg['compare_to']],\
                    ['Nozzle offset', float(leg['adjacent_filament_offset'])/0.603],\
                    ['Line placement', ormap[leg['adjacent_filament_orientation'].strip()]],\
                    ['Surface tension (mJ/m^2)',round(1000*(float(leg['sigma'])))]]
    else:
        plotvars = [['Simulation number', shortname],\
                ['Simulation rate (hr/s)', simrate],\
                ['Interfacial tension (mJ/m^2)', round(1000*(float(leg['sigma'])))]]
                    
    
    for j,s in enumerate(['ink', 'sup']):
        model = leg[s+'_transportModel']
        density = float(leg[s+'_rho'])
        col = colors[j]
        plotvars.append([s+' density (g/mL)', density/1000])
        rows[j].append(len(plotvars)-1)
        plotvars.append([s+' transport model', model])
        rows[j].append(len(plotvars)-1)
        if model=='Newtonian':
            eta0 = float(leg[s+'_nu'])*density
            eta = [eta0 for i in gammadot]
            plotvars.append([s+' viscosity (Pa*s)', exp(eta0)])
            rows[j].append(len(plotvars)-1)
        else:
            eta = [0 for i in gammadot]
            eta0 = float(leg[s+'_nu0'])*density
            tau0 = float(leg[s+'_tau0'])*density
            k = float(leg[s+'_k'])*density
            n = float(leg[s+'_n'])
            plotvars.append([s+' plateau viscosity (Pa*s)', exp(eta0)])
            rows[j].append(len(plotvars)-1)
            plotvars.append([s+' yield stress (Pa)', tau0])
            rows[j].append(len(plotvars)-1)
            plotvars.append([s+' k (Pa*s^n)', k])
            rows[j].append(len(plotvars)-1)
            plotvars.append([s+' n', n])
            rows[j].append(len(plotvars)-1)
            for i, gammad in enumerate(gammadot):
                eta[i] = min(eta0, tau0/gammad + k*gammad**(n-1))
    
        ax.plot(gammadot, eta, label=s, color=col)
        
    ax.legend()
    t1 = ax.table(cellText=[[i[0]+'\n   '+str(i[1])] for i in plotvars], loc='right', edges='', cellLoc='left')
    for k, rows in enumerate(rows):
        for i in rows:
            for j in [0]:
                t1[(i, j)].get_text().set_color(colors[k])
    t1.auto_set_font_size(False)
    t1.set_fontsize(12)
    t1.auto_set_column_width(col=[0])
    t1.scale(1, 1)
    fig.set_size_inches(1216/mydpi, 1216/mydpi) # to account for a smaller monitor, was 1216 1216 RG
    fig.tight_layout()
    fig.savefig(fname, dpi=mydpi)
    plt.close()
    
    if diag:
        logging.info(f'Exported {shortname}\\titleCard.png')
        
        
        
def convertToFastGif(vid:str, factor:int=1, crop:dict={}):
    '''convert the mp4 to a fast gif, reducing frames by a factor of int'''
    video = imageio.get_reader(vid,  'ffmpeg')
    newName = vid.replace('.mp4', '_fast.gif')
    if os.path.exists(newName):
        os.remove(newName)
    dat = video.get_meta_data()
    fps = int(dat['fps'])
    result = imageio.get_writer(newName, fps=fps)
    skip = factor
    i = 0
    for num , im in enumerate(video):
        frame = video.get_data(num)
        skip = skip-1
        if skip==0:
            # Write the frame into the
            # file 'filename.avi'
            if 'x0' in crop:
                frame = frame[crop['y0']:crop['yf'], crop['x0']:crop['xf']]
            result.append_data(frame)
            skip = factor

      # When everything done, release 
    # the video capture and video 
    # write objects
    video.close()
    result.close()

    logging.info(f'Exported {newName}')
    return
