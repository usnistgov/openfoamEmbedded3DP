#!/usr/bin/env pvpython
'''Functions for importing vtk files of simulated filaments.'''

# external packages
import os
import sys
import csv
import re
import numpy as np # RG
from paraview.simple import * # import the simple module from the paraview
import time
from datetime import datetime
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(currentdir)
sys.path.append(parentdir)
sys.path.append(os.path.join(parentdir, 'py'))  # add python folder
#import folderparser as fp
import py.file.file_handling as fh
from py.folder_stats import folderStats

# logging
logger = logging.getLogger(__name__)




#################################################################

def isNum(s:str) -> bool:
    '''check if the character is a number'''
    if s in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
        return True
    else:
        return False

def simNum(f:str) -> int:
    '''extract the simulation number from the folder'''
    i = 0
    while i<len(f) and not isNum(f[i]):
        i+=1
    j = i
    while j<len(f) and isNum(f[j]):
        j+=1
    if i==j: # no number RG
        return -100
    return int(f[i:j])


def filterSimNums(topfolders:List[str], nlist:List[int]) -> List[str]:
    '''get a list of folders that have sim numbers in the list'''
    folders = []
    for topfolder in topfolders:
        for f in fh.simFolders(topfolder):
            if simNum(os.path.basename(f)) in nlist:
                folders.append(f)
    return folders

def extractCorNums(topfolders:List[str], dirs:List[str], nlist:List[str]) -> Union[List[str],List[str]]: # RG
    '''get a list of folders corresponding to sims in the list'''
    nb = []
    nlist = [simNum(os.path.basename(n)) for n in nlist]
    for topfolder in topfolders:
        for f in fh.simFolders(topfolder):
            if simNum(os.path.basename(f))!=-100 and simNum(os.path.basename(f)) in nlist:
                geo = os.path.join(f,'geometry.csv')
                with open(geo, "r") as g:
                    data = list(csv.reader(g))
                    for row in data:
                        if row[0]=='corresponding simulation':
                            nb.append(row[1][3:])
    nb = [int(i) for i in nb]
    folders = filterSimNums(dirs, nb)
    unabridged = []
    unabridged = [filterSimNums(dirs,[i])[0] for i in nb]
    return folders, unabridged


class stateVars():
    '''This object holds important information about the folder and paraview objects that are used across functions'''
    
    def __init__(self, folder):
        mksubdirs(folder) # make images and interfacePoints folders
        self.folder = folder        
        self.caseVTMSeries = ""
        self.renderView1 = ""
        self.animationScene1 = ""
        self.timeKeeper1 = ""
        self.timeStamp = ""
        self.times = ""
        self.slice = ""
        self.initialized = False
        self.colorBars = []
        self.fs = folderStats(self.folder)
        
   
    def readVisc(self):
        '''Get viscosity settings from legend. This is necessary for making viscosity maps involving Newtonian fluids.'''

        self.inkmodel = self.fs.ink.transportModel
        if self.inkmodel=='Newtonian':
            self.inkfunc = 'inknu'
            self.inknu = self.fs.ink.kinematic['nu']
        else:
            self.inkfunc = 'nu1'
            self.inknu = self.fs.ink.kinematic['nu0']
        self.supmodel = self.fs.sup.transportModel
        if self.supmodel=='Newtonian':
            self.supfunc = 'supnu'
            self.supnu = self.fs.sup.kinematic['nu']
        else:
            self.supfunc = 'nu2'
            self.supnu = self.fs.sup.kinematic['nu0']
        
        return 0   

    def hideAll(self) -> None:
        for cb in self.colorBars:
            try:
                HideScalarBarIfNotNeeded(cb, self.renderView1)   
            except Exception as e:
                logging.error(f'Error hiding color bar: {e}')
        hideAll()
        
        
    def stressFunc(self) -> str:
        '''get the stress function from the legend dictionary'''
        ink_rho = self.fs.ink.kinematic['rho']
        sup_rho = self.fs.sup.kinematic['rho']
        if self.fs.ink.transportModel=='Newtonian':
            nuink = self.fs.ink.kinematic['nu']
            if self.fs.sup.transportModel=='Newtonian':
                nusup = self.fs.sup.kinematic['nu']
                return f'({ink_rho}*{nuink}*"alpha.ink"+{sup_rho}*{nusup}*(1-"alpha.ink"))*ScalarGradient'
            else:
                return f'({ink_rho}*{nuink}*"alpha.ink"+{sup_rho}*nu2*(1-"alpha.ink"))*ScalarGradient'
        else:
            if self.fs.sup.transportModel=='Newtonian':
                nusup = self.fs.sup.kinematic['nu']
                return f'({ink_rho}*nu1*"alpha.ink"+{sup_rho}*{nusup}*(1-"alpha.ink"))*ScalarGradient'
            else:
                return f'({ink_rho}*nu1*"alpha.ink"+{sup_rho}*nu2*(1-"alpha.ink"))*ScalarGradient'
        
def mksubdirs(folder:str) -> None:
    '''make the images and interfacePoints folders  '''
    os.makedirs(os.path.join(folder, 'images'), exist_ok = True)
    os.makedirs(os.path.join(folder, 'interfacePoints'), exist_ok = True)
    os.makedirs(os.path.join(folder, 'nozzleSlicePoints'), exist_ok = True) # RG




##-----------------------
# paraview


def cleanSession() -> None:
    '''reset the paraview session'''
    Disconnect()
    Connect()

    

def hideAll() -> None:
    '''hide all existing displays'''
    ss = GetSources()
    for s in ss:
        Hide(ss[s])


def initSeries0(sv:stateVars) -> Any:
    '''import the vtk series file
    this returns the datareader object that refers to the series'''
    sfile = sv.fs.fh1.series()
    if not os.path.exists(sfile):
        raise NameError('vtm.series file does not exist')
    if 'vtm' in sfile:
        caseVTMSeries = XMLMultiBlockDataReader(FileName=sfile)
        caseVTMSeries.CellArrayStatus = ['U', 'alpha.ink', 'p_rgh']
        caseVTMSeries.PointArrayStatus = ['U', 'alpha.ink', 'p_rgh']
    else:
        bn = os.path.basename(sfile)
        caseVTMSeries = LegacyVTKReader(registrationName=bn, FileNames=[sfile])
    sv.times = sv.fs.fh1.parseVTKSeries() # read times from the VTK.series file
    return caseVTMSeries
    


def setTime(time:float, sv:stateVars) -> None:
    '''Go to a time. If the requested time is not in the list, but you expect it to be there, try remaking the VTK.series file and look again.'''
    if time in sv.times:
        sv.animationScene1.AnimationTime = time
        sv.timeKeeper1.Time = time
    else:
        if max(sv.times)>0.1*len(sv.times): # this time should be in the list, but it is not
            sv.fs.fh1.redoVTKSeriesNoLog() # remake the VTK.series file
            sv.times = sv.fs.fh1.parseVTKSeries() # read times from the new VTK.series file
            if time in sv.times:
                sv.animationScene1.AnimationTime = time
                sv.timeKeeper1.Time = time
            else:
                logging.warning(f'time {time} is not in list {sv.times}')
                sv.animationScene1.AnimationTime = sv.times[-1]
                sv.timeKeeper1.Time = sv.times[-1]
           
        

def initializeP(sv:stateVars) -> stateVars:  
    '''initialize paraview'''
    ResetSession()     
    LoadPalette(paletteName='WhiteBackground')
        # make background white
    paraview.simple._DisableFirstRenderCameraReset()
        # disable automatic camera reset on 'Show'  
    hideAll()
        # hide all existing sources
    sv.animationScene1 = GetAnimationScene()
        # get animation scene   
    sv.timeKeeper1 = GetTimeKeeper()
        # get the time-keeper  
    sv.animationScene1.UpdateAnimationUsingDataTimeSteps()
        # update animation scene based on data timesteps
    sv.renderView1 = GetActiveViewOrCreate('RenderView')
        # get active view
    sv.renderView1.ViewSize = [1216, 1216]
        # set view size
    sv.renderView1.OrientationAxesVisibility = 0
        # hide orientation axes
    sv.renderView1.CameraParallelProjection = 1
        # turn off perspective
    return sv


def computeShearRate(sv:stateVars):
    '''get an object that represents the shear rate'''
    calculator = Calculator(Input=sv.caseVTMSeries) # RG
    calculator.ResultArrayName = 'magU'
    calculator.Function = 'mag(U)'
    
    computeDerivatives = ComputeDerivatives(Input=calculator)
    computeDerivatives.Scalars = ['POINTS', 'magU'] # RG
    computeDerivatives.Vectors = ['POINTS', 'U']
    computeDerivatives.OutputVectorType = 'Scalar Gradient'
    computeDerivatives.OutputTensorType = 'Vector Gradient'
    # show data in view
    # computeDerivativesDisplay = Show(computeDerivatives, renderView1, 'GeometryRepresentation')
    return computeDerivatives