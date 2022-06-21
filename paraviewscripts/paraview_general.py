#!/usr/bin/env pvpython
'''Functions for importing vtk files of simulated filaments.'''

# external packages
import os
import sys
import csv
import re
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
import folderparser as fp

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
    return int(f[i:j])


def filterSimNums(topfolders:List[str], nlist:List[int]) -> List[str]:
    '''get a list of folders that have sim numbers in the list'''
    folders = []
    for topfolder in topfolders:
        for f in fp.caseFolders(topfolder):
            if simNum(os.path.basename(f)) in nlist:
                folders.append(f)
    return folders


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
        
   
    def readVisc(self):
        '''Get viscosity settings from legend. This is necessary for making viscosity maps involving Newtonian fluids.'''
        # le = fp.importIf(self.folder)
        le = fp.legendUnique(self.folder)
        if len(le)==0:
            return 1

        self.inkmodel = le['ink_transportModel']
        if self.inkmodel=='Newtonian':
            self.inkfunc = 'inknu'
            self.inknu = le['ink_nu']
        else:
            self.inkfunc = 'nu1'
            self.inknu = le['ink_nu0']
        self.supmodel = le['sup_transportModel']
        if self.supmodel=='Newtonian':
            self.supfunc = 'supnu'
            self.supnu = le['sup_nu']
        else:
            self.supfunc = 'nu2'
            self.supnu = le['sup_nu0']
        
        return 0   

    def hideAll(self) -> None:
        for cb in self.colorBars:
            try:
                HideScalarBarIfNotNeeded(cb, self.renderView1)   
            except Exception as e:
                logging.error('Error hiding color bar: '+str(e))
        hideAll()
        
def mksubdirs(folder:str) -> None:
    '''make the images and interfacePoints folders  '''
    # fp.mkdirif(os.path.join(folder, 'images'))   
    # fp.mkdirif(os.path.join(folder, 'interfacePoints'))
    os.makedirs(os.path.join(folder, 'images'), exist_ok = True)
    os.makedirs(os.path.join(folder, 'interfacePoints'), exist_ok = True)
    os.makedirs(os.path.join(folder, 'nozzleSlicePoints'), exist_ok = True) # RG

def stressFunc(le:dict) -> str:
    '''get the stress function from the legend dictionary'''
    ink_rho = le['ink_rho']
    sup_rho = le['sup_rho']
    if le['ink_transportModel']=='Newtonian':
        nuink = float(le['ink_nu'])
        if le['sup_transportModel']=='Newtonian':
            nusup = float(le['sup_nu'])
            return f'({ink_rho}*{nuink}*"alpha.ink"+{sup_rho}*{nusup}*(1-"alpha.ink"))*ScalarGradient'
        else:
            return f'({ink_rho}*{nuink}*"alpha.ink"+{sup_rho}*nu2*(1-"alpha.ink"))*ScalarGradient'
    else:
        if le['sup_transportModel']=='Newtonian':
            nusup = float(le['sup_nu'])
            return f'({ink_rho}*nu1*"alpha.ink"+{sup_rho}*{nusup}*(1-"alpha.ink"))*ScalarGradient'
        else:
            return f'({ink_rho}*nu1*"alpha.ink"+{sup_rho}*nu2*(1-"alpha.ink"))*ScalarGradient'


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
    sfile = fp.series(sv.folder)
    if not os.path.exists(sfile):
        raise NameError('vtm.series file does not exist')
    if 'vtm' in sfile:
        caseVTMSeries = XMLMultiBlockDataReader(FileName=sfile)
        caseVTMSeries.CellArrayStatus = ['U', 'alpha.ink', 'p_rgh']
        caseVTMSeries.PointArrayStatus = ['U', 'alpha.ink', 'p_rgh']
    else:
        bn = os.path.basename(sfile)
        caseVTMSeries = LegacyVTKReader(registrationName=bn, FileNames=[sfile])
    sv.times = fp.parseVTKSeries(sv.folder) # read times from the VTK.series file
    return caseVTMSeries
    


def setTime(time:float, sv:stateVars) -> None:
    '''Go to a time. If the requested time is not in the list, but you expect it to be there, try remaking the VTK.series file and look again.'''
    if time in sv.times:
        sv.animationScene1.AnimationTime = time
        sv.timeKeeper1.Time = time
    else:
        if max(sv.times)>0.1*len(sv.times): # this time should be in the list, but it is not
            fp.redoVTKSeriesNoLog(sv.folder) # remake the VTK.series file
            sv.times = fp.parseVTKSeries(sv.folder) # read times from the new VTK.series file
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