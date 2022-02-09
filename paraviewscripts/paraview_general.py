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
import folderparser as fp

# logging
logger = logging.getLogger(__name__)

# info
__author__ = "Leanne Friedrich"
__copyright__ = "This data is publicly available according to the NIST statements of copyright, fair use and licensing; see https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software"
__credits__ = ["Leanne Friedrich"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Leanne Friedrich"
__email__ = "Leanne.Friedrich@nist.gov"
__status__ = "Production"


#################################################################

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
            self.supfunc = 'nu1'
            self.supnu = le['sup_nu0']


        # while not le[i0][0]=='transportProperties':
        #     i0+=1 # find the index where the transport properties start
        # self.inkmodel = le[i0+2][1] # find the row that contains the ink model.
        # if self.inkmodel=='Newtonian':
        #     # if it's newtonian, the viscosity is stored 3 rows below 'transportProperties'
        #     nuinki = i0+3
        #     self.inkfunc = 'inknu'
        # else:
        #     # otherwise the zero shear viscosity is 4 rows below
        #     nuinki = i0+4
        #     self.inkfunc = 'nu1'
        # nuink = le[nuinki][1] # get the string viscosity in mPa s
        # self.inknu = nuink

        # # find the support transport properties
        # # in some legends, there are no Herschel Bulkley boxes, and in some there are
        # # we need to find where the support section starts
        # nusupi = nuinki 

        # # iterate through rows until we find the support transportmodel
        # while not le[nusupi][0]=='transportModel':
        #     nusupi+=1

        # self.supmodel = le[nusupi][1]
        # if self.supmodel=='Newtonian':
        #     nusupi = nusupi+1
        #     self.supfunc = 'supnu'
        # else:
        #     nusupi = nusupi+2
        #     self.supfunc = 'nu2'
        # nusup = le[nusupi][1]
        # self.supnu = nusup
        
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
    os.makedirs(os.path.join(folder, 'nozzlePoints'), exist_ok = True) # RG



        
# # find the path of the vtk series file
# def series(folder:str) -> str:
#     cf = fp.caseFolder(folder)
#     vtkfolder = os.path.join(cf, 'VTK')
#     if os.path.exists(vtkfolder):
#         for file in os.listdir(vtkfolder):
#             if '.series' in file:
#                 return os.path.join(vtkfolder, file)
#         # if there is a vtk folder but no series file, generate one
#         fp.redoVTKSeriesNoLog(folder)
#         return fp.parseVTKSeries(folder)
#     return ""


# def readTimes(folder) -> List[float]:
#     times = fp.parseVTKSeries(folder)
#     if len(times)==0 or len(times)<fp.vtkFiles(folder):
#         fp.redoVTKSeriesNoLog(folder)
#         times = fp.parseVTKSeries(folder)
#     return times



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