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

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(currentdir)
sys.path.append(parentdir)
import folderparser as fp
import folderscraper as fs

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

def get_time(filename:str) -> float:
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
            
def saveFramesToWriter(mp4writer, folder:str, tags:List[str], subfolder:str='images', debug:bool=False) -> None:
    '''For a given folder, save frames to an imageio writer
    tags should be a list for each folder, and within each list, a list of tags that should be in every file. For example, to get all sigma=0 and then all sigma=40, tags could be [['sigma_0', 'y_umag'],['sigma_40', 'y_umag']]. 
    subfolder is the folder inside of the simulation folder where the images are stored'''
    zfiles = []
    if len(subfolder)>0:
        searchFolder = os.path.join(folder, subfolder)
    else:
        searchFolder = folder
    for file in os.listdir(searchFolder):
        if 'png' in file:
            append = True
            for s in tags:
                if not s in file:
                    append = False
                    break
            if append:
                zfiles.append(os.path.join(searchFolder, file))
    times = [get_time(f) for f in zfiles]
    for i, filename in enumerate(zfiles):
        if i==0 or round(times[i]-times[i-1],3)==0.1:
            image = imageio.imread(filename)[:,:,:3] # ignore alpha channel
            mp4writer.append_data(image)
            
def saveTitleCardToWriter(mp4writer, folder:str, n:int, subfolder:str='images') -> None:
    '''For a given folder, save the title card to the videowriter. n is number of frames of the title card to add'''
    if len(subfolder)>0:
        file = os.path.join(folder, 'images', 'titleCard.png')
    else:
        file = os.path.join(folder,  'titleCard.png')
    if os.path.exists(file):
        image = imageio.imread(file)[:,:,:3] # ignore alpha channel
        for i in range(n):
            mp4writer.append_data(image)
    else:
        raise NameError(f'No title card in {folder}')
            
def saveBigVideo(folderList:str, filename:str, titleLength:float=1, diag:bool=True) -> None:
    '''Compile all of the time series for all of the simulations into one big video. 
    folderList is a list of the folders to include. 
    filename is the name of the video to save. 
    titleLength is the time that the title cards are up, in s'''
    fps = 10
    mp4writer = imageio.get_writer(filename, fps=fps)
    try:
        for folder in folderList:
            if diag:
                logging.info(f'Writing {fp.shortName(folder)}')
            saveTitleCardToWriter(mp4writer, folder, int(round(titleLength*fps)))
            for s in ['y', 'a']:
                saveFramesToWriter(mp4writer, folder, [s+'_umag'])
    except Exception as e:
        mp4writer.close()
        raise e
        
    if diag:
        logging.info('Done creating big video')
    
    mp4writer.close()
    
    
def saveFigureVideo(topfolderList:str, filename:str, tags:List[List[str]], titleLength:float=1, diag:bool=True) -> None:
    '''Compile a time series of the combined picture plots for all of the simulations into one big video. 
    folderList is a list of the topfolders to include.
    filename is the name of the video to save. 
    tags should be a list for each folder, and within each list, a list of tags that should be in every file. For example, to get all sigma=0 and then all sigma=40, tags could be [['sigma_0', 'y_umag'],['sigma_40', 'y_umag']]. 
    titleLength is the time that the title cards are up, in s'''
    fps = 10
    mp4writer = imageio.get_writer(filename, fps=fps)
    try:
        for folder in topfolderList:
            if diag:
                logging.info(f'Writing {fp.shortName(folder)}')
            saveTitleCardToWriter(mp4writer, folder, int(round(titleLength*fps)), subfolder='')
            for tag in tags:
                saveFramesToWriter(mp4writer, folder, tag, subfolder='', debug=diag)
    except Exception as e:
        mp4writer.close()
        raise e
        
    if diag:
        logging.info(f'Done creating video {filename}')
    
    mp4writer.close()
    
def tryFloat(s:str) -> Union[str, float]:
    '''try to convert the string to a float'''
    try:
        sout = float(s)
    except:
        return s
    else:
        return sout
    
def findFiles(topfolders:str, keyList:List[str]) -> List[str]:
    '''Find all of the simulation folders in topfolders, and put them in order based on their fp.legendUnique keys, given in priority order in keyList'''
    folders = []
    for topfolder in topfolders:
        for cf in fp.caseFolders(topfolder):
            leg = fp.legendUnique(cf)
            keylocal = [value for value in keyList if value in leg]
            folders.append(dict([['file', cf]] + [[s, tryFloat(leg[s])] for s in keylocal]))
    folders = pd.DataFrame(folders)
    keyList = [value for value in keyList if value in folders]
    folders = folders.sort_values(by=keyList)
    return folders
#     return list(folders['file'])
    
            
#--------------------------------
            
def exp(n:float) -> str:
    '''Format a power of 10 in exponential notation'''
    return r'$10^{{{}}}$'.format(int(round(np.log10(n))))


#--------------------------------

def checkSimRate(folder:str, fix:bool=True) -> str:
    '''Check if the simulation rate makes sense. fix=True to rescrape the legend if the simulation rate doesn't make sense'''
    leg = fp.legendUnique(folder) # legend
    shortname = os.path.join(os.path.basename(os.path.dirname(folder)), os.path.basename(folder))
    simrate = leg['simulation_rate_(hr/s)']
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
    shortname = fp.shortName(folder)
    if os.path.exists(fname) and not overwrite:
        return
    
    mydpi = 144
    gammadot = [10**n for n in np.arange(-6, 6, 0.1)]
    fig,ax = plt.subplots(1,1, figsize=(3,3))
    ax.set_ylabel(r'Viscosity $\eta$ (Pa.s)', fontname="Arial")
    ax.set_xlabel(r'Strain rate $\dot{\gamma}$ (s$^{-1}$)', fontname="Arial")
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylim([10**-3, 10**6])
    ax.set_xlim([10**-6, 10**6])
    
    colors = ['#263e59', '#cc3d3d']
    leg = fp.legendUnique(folder) # legend
    simrate = checkSimRate(folder, True)
    plotvars = [['Simulation number', shortname],\
                ['Simulation rate (hr/s)', simrate],\
                ['Surface tension (mJ/m^2)', round(1000*(float(leg['sigma'])))]]
    rows = [[],[]]
    
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
    fig.set_size_inches(1216/mydpi, 1181/mydpi) # to account for a smaller monitor, was 1216 1216 RG
    fig.tight_layout()
    fig.savefig(fname, dpi=mydpi)
    plt.close()
    
    if diag:
        logging.info(f'Exported {shortname}\\titleCard.png')
