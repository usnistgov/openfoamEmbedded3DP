#!/usr/bin/env python
'''Functions for combining images into videos'''

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

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import folderparser as fp

plt.rcParams.update({'font.size': 12})
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

__author__ = "Leanne Friedrich"
__copyright__ = "This data is publicly available according to the NIST statements of copyright, fair use and licensing; see https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software"
__credits__ = ["Leanne Friedrich"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Leanne Friedrich"
__email__ = "Leanne.Friedrich@nist.gov"
__status__ = "Production"

#----------------------------------------------


def get_length(filename:str) -> float:
    '''Get the length of the video'''
    rl = ["ffprobe", "-v", "error", "-show_entries", "format=duration", "-of", "default=noprint_wrappers=1:nokey=1", filename]
    result = subprocess.run(rl, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    duration =float(result.stdout)
    return duration

def get_time(filename:str) -> float:
    '''Get the time of the simulation timepoint'''
    return float(re.split('t|_', os.path.basename(filename))[1])/10

def saveVid(folder:str, s:str, p:str, diag:bool=True) -> None:
    '''compile images into a video, where images contain the string {s}_{p}.png'''
    zfiles = glob.glob(os.path.join(folder,'images', '*'+s+'_'+p+'.png'))
    if len(zfiles)>0:
        fn = os.path.join(folder,'images',  s+'_'+p+'.mp4')
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
            
def saveFramesToWriter(mp4writer, folder:str, s:str, p:str) -> None:
    '''For a given folder, given an angle s (e.g. 'y', 'a') and a coloring (e.g. 'umag') to an imageio writer'''
    zfiles = glob.glob(os.path.join(folder,'images', '*'+s+'_'+p+'.png'))
    times = [get_time(f) for f in zfiles]
    for i, filename in enumerate(zfiles):
        if i==0 or round(times[i]-times[i-1],3)==0.1:
            image = imageio.imread(filename)[:,:,:3] # ignore alpha channel
            mp4writer.append_data(image)
            
def saveTitleCardToWriter(mp4writer, folder:str, n:int) -> None:
    '''For a given folder, save the title card to the videowriter. n is number of frames of the title card to add'''
    file = os.path.join(folder, 'images', 'titleCard.png')
    if os.path.exists(file):
        image = imageio.imread(file)[:,:,:3] # ignore alpha channel
        for i in range(n):
            mp4writer.append_data(image)
            
def saveBigVideo(folderList:str, filename:str, titleLength:float=1, diag:bool=True) -> None:
    '''Compile all of the time series for all of the simulations into one big video. folderList is a list of the folders to include. filename is the name of the video to save. titleLength is the time that the title cards are up, in s'''
    fps = 10
    mp4writer = imageio.get_writer(filename, fps=fps)
    try:
        for folder in folderList:
            if diag:
                logging.info(folder)
            saveTitleCardToWriter(mp4writer, folder, int(round(titleLength*fps)))
            for s in ['y', 'a']:
                saveFramesToWriter(mp4writer, folder, s, 'umag')
    except Exception as e:
        mp4writer.close()
        raise e
        
    if diag:
        logging.info('Done creating big video')
    
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

def titleCard(folder:str, overwrite:bool=False, diag:bool=True) -> None: 
    '''Create and export a title card using the legend for the simulation'''
    
    fname = os.path.join(folder, 'images', 'titleCard.png')
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
    
    colors = ['#263e59', '#b03810']
    leg = fp.legendUnique(folder) # legend
    
    simrate = leg['simulation_rate_(hr/s)']
    try:
        simrate = '{:.1f}'.format(float(simrate))
    except:
        pass
    
    plotvars = [['Simulation number', leg['folder']],\
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
    t1 = ax.table(cellText=plotvars, loc='right', edges='', cellLoc='left')
    for k, rows in enumerate(rows):
        for i in rows:
            for j in [0,1]:
                t1[(i, j)].get_text().set_color(colors[k])
    t1.auto_set_font_size(False)
    t1.set_fontsize(12)
    t1.auto_set_column_width(col=[0,1])
    t1.scale(1, 1)
    fig.set_size_inches(1216/mydpi, 1216/mydpi)
    fig.tight_layout()
    fig.savefig(fname, dpi=mydpi)
    plt.close()
    
    if diag:
        print(os.path.join(os.path.basename(os.path.dirname(folder)), os.path.basename(folder)))
