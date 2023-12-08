#!/usr/bin/env pvpython
'''Functions for generating images of filaments from vtk files.'''

# external packages
import os
import numpy as np
import csv
import re
import time
from datetime import datetime
from typing import List, Dict, Tuple, Union, Any, TextIO
from vtkmodules.vtkCommonCore import vtkLogger
vtkLogger.SetStderrVerbosity(vtkLogger.VERBOSITY_OFF)
import logging
import traceback
from paraview.simple import * # import the simple module from the paraview

# local packages
from paraview_general import *
from folder_stats import folderStats
from colorBars import *

# logging
logger = logging.getLogger(__name__)


#################################################################

FONT = 24

##############################################
######## SCREENSHOTS #########################



def getpfunc(f, kwargs:Dict):
    '''Get the pressure function, given a function f and dict kwargs'''
    kw = []
    for s in ['pmin', 'pmax', 'name', 'x']:
        if s in kwargs:
            kw = kw+[s, kwargs[s]]
    return lambda sv: f(sv, **dict(kw))

def getShearRatefunc(f, kwargs:Dict):
    '''Get the shear rate function, given a function f and dict kwargs'''
    kw = []
    for s in ['rmin', 'rmax', 'x']:
        if s in kwargs:
            kw = kw+[s, kwargs[s]]
    return lambda sv: f(sv, **dict(kw))

def getShearStressfunc(f, kwargs:dict): # RG
    '''Get the shear stress function, given a function f and dict kwargs'''
    kw = {}
    for s in ['rmin', 'rmax', 'x', 'yieldSurface']:
        if s in kwargs:
            kw[s] = kwargs[s]
    return lambda sv: f(sv, **kw)
        
def getxfunc(f, kwargs:Dict):
    '''Get the function for an x slice, given function f and dict kwargs'''
    if 'x' in kwargs:
        return lambda sv: f(sv, x=kwargs['x'])
    else:
        return lambda sv: f(sv, **kwargs) # RG

class ssVars:
    '''Holds input variables that tell us what kind of images to generate.
    tag is the type of image we want to collect. 
        'volumes': in 3D, look at the interface between the ink and support
            optional keyword variables:
            'coloring': the variable on which to color the surface. default is 'umag', which uses 'U' coloring
            'volViewList':List[str]: list of view angles. Allowed angles = 'x' (view from x direction), 'y' (view from y direction), 'z' (view from z direction), 'a' (view at angle from above/downstream/side), 'b' (view at angle from above/upstream/side)
        'meshes': take a slice out of the middle of the volume, at y=0, and show the mesh
        'vectors': take a slice out of the middle of the volume, at y=0, and show the flow field with arrows
        'tubes': in 3D, look at the interface between the ink and support, and add streamlines
            optional keyword variables:
            'tubeh': the absolute z position (in m) at which to collect streamlines. default is z=0.001 m
            'volViewList':List[str]: list of view angles. Allowed angles = 'x' (slice on x plane), 'y' (slice on y plane), 'a' (view at angle from above/downstream/side), 'b' (view at angle from above/upstream/side)
        'viscy', 'viscx': take a slice and show the viscosity.
            optional keyword variables:
            'x' (only for viscx): the absolute x position (in m) at which to take the slice. Default is x=-0.001 m for a 0.603 mm ID nozzle.
        'alphaSlicex': take a slice and color the ink white and support black
        'uslicey', 'uslicex': take a slice and show the velocity magnitude
            optional keyword variables:
            'x' (only for uslicex): the absolute x position (in m) at which to take the slice. Default is x=-0.001 m for a 0.603 mm ID nozzle.
        'uzslicey', 'uzslicex': take a slice and show the z velocity
            optional keyword variables:
            'x' (only for uzslicex): the absolute x position (in m) at which to take the slice. Default is x=-0.001 m for a 0.603 mm ID nozzle.
        'py', 'px': take a slice and show the adjusted pressure (pressure - hydrostatic pressure)
            optional keyword variables:
            'x' (only for px): the absolute x position (in m) at which to take the slice. Default is x=-0.001 m for a 0.603 mm ID nozzle.
            'pmin', 'pmax': min and max pressure on the legend
            'name': 'p' if you want pressure or 'p_rgh' if you want adjusted pressur
    tlist is the list of times. leave empty to get pictures for every time       
    '''
    
    
    def __init__(self, tag:str, tlist:List[float], **kwargs):
        if tag=='volumes':
            if 'volViewList' not in kwargs:
                self.volList = ['y']
            else:
                self.volList = kwargs['volViewList']
            if not 'coloring' in kwargs:
                self.coloring = 'umag'
                self.function = vol
            else:
                self.coloring = kwargs['coloring']
                self.function = lambda sv: vol(sv, coloring=self.coloring)
        elif tag=='meshes':
            self.coloring = 'mesh'
            self.function = mesh
            self.volList = ['y']
        elif tag=='vectors':
            self.coloring = 'vecs'
            self.function = sliceandyarrows
            self.volList = ['y']
        elif tag=='tubes':
            if 'volViewList' not in kwargs:
                self.volList = ['a']
            else:
                self.volList = kwargs['volViewList']
            if 'tubeh' not in kwargs:
                self.tubeh = 0.001
            else:
                self.tubeh = kwargs['tubeh']
            self.coloring = 'stre'
            self.function = surfandtube
        elif tag=='viscy':
            self.coloring='viscy'
            self.function=viscy
            self.volList = ['y']
        elif tag=='viscx':
            self.coloring = 'viscx'
            self.function = getxfunc(viscx, kwargs)
            self.volList = ['x']
        elif tag=='viscNozx':
            self.coloring = 'viscNozx'
            self.function = getxfunc(viscNozx, kwargs)
            self.volList = ['x']
        elif tag=='uslicey':
            self.coloring = 'uslicey'
            self.function = uslicey
            self.volList = ['y']
        elif tag=='uslicex':
            self.coloring = 'uslicex'
            self.function = getxfunc(uslicex, kwargs)
            self.volList = ['x']
        elif tag=='uzslicey':
            self.coloring = 'uzslicey'
            self.function = uzslicey
            self.volList = ['y']
        elif tag=='uzslicex':
            self.coloring = 'uzslicex'
            self.function = getxfunc(uzslicex, kwargs)
            self.volList = ['x']
        elif tag=='uzsliceNozx':
            self.coloring = 'uzsliceNozx'
            self.function = getxfunc(uzsliceNozx, kwargs)
            self.volList = ['x']
        elif tag=='uyslicey':
            self.coloring = 'uyslicey'
            self.function = uyslicey
            self.volList = ['y']
        elif tag=='uyslicex':
            self.coloring = 'uyslicex'
            self.function = getxfunc(uyslicex, kwargs)
            self.volList = ['x']
        elif tag=='uysliceNozx':
            self.coloring = 'uysliceNozx'
            self.function = getxfunc(uysliceNozx, kwargs)
            self.volList = ['x']
        elif tag=='alphaSlicex': # RG
            if not 'coloring' in kwargs:
                self.coloring = 'alphaSlice_-1'
            else:
                self.coloring = kwargs['coloring']
            self.function = getxfunc(alphaSlicex, kwargs)
            self.volList = ['x']
        elif tag=='py':
            self.coloring='py'
            self.function = getpfunc(pslicey, kwargs)
            self.volList = ['y']
        elif tag=='px':
            self.coloring = 'px'
            self.function = getpfunc(pslicex, kwargs)
            self.volList = ['x']
        elif tag=='pNozx':
            self.coloring = 'pNozx'
            self.function = getpfunc(psliceNozx, kwargs)
            self.volList = ['x']
        elif tag=='shearRatex':
            self.coloring='shearRatex'
            self.function = getShearRatefunc(shearRateSlicex, kwargs)
            self.volList = ['x']
        elif tag=='shearRateNozx':
            self.coloring='shearRateNozx'
            self.function = getShearRatefunc(shearRateSliceNozx, kwargs)
            self.volList = ['x']
        elif tag=='shearRatey':
            self.coloring='shearRatey'
            self.function = getShearRatefunc(shearRateSlicey, kwargs)
            self.volList = ['y']
        elif tag=='shearStressx': # RG
            self.coloring='shearStressx'
            self.function = getShearStressfunc(shearStressSlicex, kwargs)
            self.volList = ['x']
        elif tag=='shearStressNozx': 
            self.coloring='shearStressNozx'
            self.function = getShearStressfunc(shearStressSliceNozx, kwargs)
            self.volList = ['x']
        elif tag=='shearStressy': # RG
            self.coloring='shearStressy'
            self.function = getShearStressfunc(shearStressSlicey, kwargs)
            self.volList = ['y']
        self.tlist = tlist
        self.tag = tag
        
    def prnt(self) -> str:
        '''convert self to a string'''
        if len(self.tlist)>0:
            tprint = self.tlist
        else:
            tprint = 'all'
        s = f'\t{self.tag}: times={tprint}'
        if self.tag=='volumes':
            s = s+f', views:{self.volList}, coloring:{self.coloring}'
        elif self.tag=='tubes':
            s = s+f', views:{self.volList}, tube position:{self.tubeh}'
        return s

####### file operations ######
    

def timeStr(t:float) -> str:
    '''convert a time to a 3 digit string'''
    tstr = '{0:0=3d}'.format(t)
    return tstr


def fullFN(t:float, view:str, coloring:str, folder:str):
    '''get a full file name for an image'''
    fn = f't{timeStr(t)}_{view}_{coloring}.png'
    filename = os.path.join(folder, 'images', fn)
    return filename


######### show/hide/initializes ############
#### import the vtk series file


def initializeSV(sv:stateVars) -> None:
    '''Initialize paraview, import series file, and create a time stamp. Use this if you already have a stateVars object'''
    logging.info(f'Screenshots: Initializing paraview for {sv.folder}')
    ResetSession()
    sv =  initializeP(sv)
        # initialize Paraview
    sv =  initSeries(sv)
        # import the vtk series files
    sv = timeStamp(sv) 
        # add time stamp
    sv.initialized = True
    return 

def initSeries(sv:stateVars) -> stateVars:
    '''Import the series file and set up the display'''
    caseVTMSeries = initSeries0(sv)
    caseVTMSeriesDisplay = Show(caseVTMSeries, sv.renderView1, 'GeometryRepresentation')
    sv.renderView1.ResetCamera()
    sv.renderView1.Update()
    sv.caseVTMSeries = caseVTMSeries
    return sv    
 

def timeStamp(sv:stateVars) -> stateVars:
    '''create a timestamp to add to images later'''
    annotateTimeFilter1 = AnnotateTimeFilter(Input=sv.caseVTMSeries)
    annotateTimeFilter1.Format = '%2.1f s' # RG
    annotateTimeFilter1Display = Show(annotateTimeFilter1, sv.renderView1, 'TextSourceRepresentation')
    annotateTimeFilter1Display.FontSize = FONT
    sv.timeStamp = annotateTimeFilter1
    return sv


###############
# camera operations

def resetCamera(sv:stateVars) -> None:
    '''reset the view area'''
    s = sv.fs.geo.bath_width/9.668
    sv.renderView1.ResetCamera(-0.005*s, 0.005*s, -0.002*s, 0.002*s, -0.002*s, 0.002*s)


def resetCam(sv:stateVars, time:float) -> None:
    '''go to a specific view area, at a given time'''
    setTime(time, sv)
    resetCamera(sv)

def setView(st:str, sv:stateVars) -> None:
    '''go to a specific view point, out of [x,y,z,a,b]. a and b are views from an angle, where a is downstream of the nozzle, and b is upstream of the nozzle.'''
    if st=="z":
        sv.renderView1.CameraFocalPoint = [0, 0,0]
        sv.renderView1.CameraPosition = [0, 0, 10]
        sv.renderView1.CameraViewUp = [0,1,0]
    elif st=="y":
        sv.renderView1.CameraFocalPoint = [0, 0,0]
        sv.renderView1.CameraPosition = [0, -10, 0]
        sv.renderView1.CameraViewUp = [0,0,1]
    elif st=="x":
        sv.renderView1.CameraFocalPoint = [0, 0,0]
        sv.renderView1.CameraPosition = [10, 0, 0]
        sv.renderView1.CameraViewUp = [0,0,1]
    elif st=="a":
        di = sv.fs.geo.niw # change back RG
        sv.renderView1.CameraFocalPoint = [0.0006*di/0.6, -0.0002*di/0.6, 0.0005*di/0.6]
        sv.renderView1.CameraPosition = [2*di/0.6, 1*di/0.6, 2*di/0.6]
        sv.renderView1.CameraViewUp = [0,0,1]
    elif st=="b":
        di = sv.fs.geo.niw 
        sv.renderView1.CameraFocalPoint = [0.001*di/0.6, 0, 0.001*di/0.6]
        sv.renderView1.CameraPosition = [-1*di/0.6,-1*di/0.6,1*di/0.6]
        sv.renderView1.CameraViewUp = [0,0,1]


def setAndUpdate(st:str, sv:stateVars) -> None:
    '''go to a specific viewpoint and view area'''
    setView(st, sv)
    sv.renderView1.Update()
    resetCamera(sv)


#### annotations





################ types of surfaces #########################


def vol(sv:stateVars, coloring:str='U') -> None:
    '''filament surface with velocity coloring'''
    sv.hideAll()
    alphasurface(sv, coloring)


def plainalphasurface(sv:stateVars, **kwargs) -> None:
    '''filament surface with no coloring'''
    sv.hideAll()
    alphasurface(sv, 'None', **kwargs)
    
def setDisplayColor(display, color:str, *args):
    '''Set the whole display to a single color'''
    try:
        x = ['POINTS', color]
        display.SetScaleArray = x
        display.OpacityArray = x
        if len(args)==0:
            ColorBy(display, ('POINTS', color))
        else:
            ColorBy(display, ('POINTS', color, args[0]))
    except Exception as e:
        logging.error('setDisplayColor:'+str(e))
        pass

    

    
    
def inkClip(sv:stateVars, clipinput, colVar:str, invert:int, **kwargs):
    '''Just clip out the ink portion and color it by viscosity or velocity. 
    clipinput is the Paraview object that you're trying to clip, e.g. sv.caseVTMSeries. 
    colVar is the variable being colored ('U', 'nu1', 'nu2', 'inknu', or 'supnu'). 'nu1' and 'nu2' use local viscosities measured for non-Newtonian fluids. 'inknu' and 'supnu' are constant viscosities collected from the legend, for Newtonian fluids. 
    invert 0 to get ink, 1 to get support. 
    A clipVal in kwargs allows you to custom set the alpha value where you want to clip.'''
    
    clip = Clip(Input=clipinput)

    # clip out just the fluid we're looking at
    clip.Scalars = ['POINTS', 'alpha.ink']

    if colVar=='nu1' or colVar=='inknu':
        clip.Value = 0.9
    elif colVar=='nu2' or colVar=='supnu':
        clip.Value = 0.1
    elif 'clipVal' in kwargs:
        clip.Value = kwargs['clipVal']
    else:
        clip.Value = 0.5
    clip.ClipType = 'Scalar'
    clip.Invert = invert

    # show data in view
    clipDisplay = Show(clip, sv.renderView1, 'UnstructuredGridRepresentation')
    clipDisplay.Representation = 'Surface'
    
    if colVar=="U":
        setDisplayColor(clipDisplay, 'U')
    elif colVar[0]=='U':
        setDisplayColor(clipDisplay, 'U', colVar[1])
    elif colVar=='p_rgh':
        setDisplayColor(clipDisplay, 'p_rgh')
    elif colVar=='p':
        setDisplayColor(clipDisplay, 'p')
    elif colVar=='nu1' or colVar=='nu2':
        setDisplayColor(clipDisplay, colVar) 
        sv.colorBars.append(nuColorBar(sv.renderView1, clipDisplay, colVar))
    elif colVar=='inknu' or colVar=='supnu':
        if colVar=='inknu':
            viscColor(clipDisplay, sv.inknu)
        else:
            viscColor(clipDisplay, sv.supnu)
    elif colVar=='shearRate':
        ColorBy(clipDisplay, ('CELLS', 'ScalarGradient', 'Magnitude'))
    elif colVar=='shearStress': # RG
        ColorBy(clipDisplay, ('POINTS', 'shearStress', 'Magnitude'))
    else:
#         clipDisplay.SetScaleArray = 'None'
#         ColorBy(clipDisplay, None)
        c = [0.9, 0.9, 0.9]
        clipDisplay.AmbientColor = c
        clipDisplay.DiffuseColor = c
        
    Hide(clipinput, sv.renderView1)
    
    return clipDisplay, clip


def stressClip(sv, display, invert, clipVal):
    '''clip by stress'''
    try:
        clip = Clip(Input=display)
    except:
        traceback.print_exc()

    # clip out just the fluid we're looking at
    clip.Scalars = ['POINTS', 'shearStressMag']
    clip.Value = clipVal
    clip.ClipType = 'Scalar'
    clip.Invert = invert
    
    return clip
    

def yieldClip(sv, display, tau0:float):
    '''clip just the yield stress out and color it black'''
    tauclip0 = stressClip(sv, display, 0, clipVal = (10**-0.01)*tau0)
    tauclip = stressClip(sv, tauclip0, 1, clipVal = (10**0.01)*tau0)
        # turn off scalar coloring
    clipDisplay = Show(tauclip, sv.renderView1, 'UnstructuredGridRepresentation')
    clipDisplay.Representation = 'Surface'
    c = [0.0,0.0,0.0]
    clipDisplay.AmbientColor = c
    clipDisplay.DiffuseColor = c
    return tauclip



def alphasurface(sv:stateVars, colVar:str, opacity:float=1, **kwargs) -> None:
    '''filament surface
    colVar is the variable being colored ('U' or 'None')'''
    clipDisplay,_ = inkClip(sv, sv.caseVTMSeries, colVar, 0)
    clipDisplay.Opacity = opacity
    sv.colorBars.append(uColorBar(sv.renderView1, clipDisplay, umax=0.01))
    resetCam(sv, (sv.times[-1])) 
    setAndUpdate('y', sv)
    sv.renderView1.Update()
    return 

def sliceMake(sv:stateVars, origin:List[float], normal:List[float], **kwargs):
    '''Get the slice object, to be further manipulated'''
    if 'sinput' in kwargs:
        slice1 = Slice(Input=kwargs['sinput'])
    else:
        slice1 = Slice(Input=sv.caseVTMSeries)
    slice1.SliceType = 'Plane'
    slice1.SliceOffsetValues = [0.0]
    Hide3DWidgets(proxy=slice1.SliceType)
    slice1.SliceType.Origin = origin
    slice1.SliceType.Normal = normal
    slice1Display = Show(slice1, sv.renderView1, 'GeometryRepresentation')
    slice1Display.SetScaleArray = ['POINTS', 'U'] # this doesn't matter because we're going to recolor later
    slice1Display.OpacityArray = ['POINTS', 'U']
    Hide(sv.caseVTMSeries, sv.renderView1)
    return slice1
    
#--------------------------

def viscSlice(sv:stateVars, origin:List[float], normal:List[float], view:str, yieldSurface:bool=True) -> None:
    '''Plot the viscosity map. Segment out the ink and support separately and color separately, leaving some white space at the interface so you can see where the interface is.'''
    out = sv.readVisc()
    if out>0:
        return
    slice1 = sliceMake(sv, origin, normal)
    inkClip(sv, slice1, sv.inkfunc, 0, clipVal = 0.9)
    inkClip(sv, slice1, sv.supfunc, 1, clipVal = 0.1)
    # color bar added in inkClip
    if yieldSurface:
        addYieldSurface(sv, origin, normal)
    resetCam(sv, (sv.times[-1])) 
    setAndUpdate(view, sv)
    sv.renderView1.Update()
    
def viscy(sv:stateVars) -> None:
    '''Viscosity map, looking down the y axis, at (0,0,0)'''
    sv.hideAll()
    viscSlice(sv, [0,0,0], [0, -1, 0], 'y')

def scaleDim(sv:stateVars, x:float) -> float:
    '''scale the dimension by the nozzle inner diameter'''
    return x * sv.fs.geo.niw/0.6
    
def viscx(sv:stateVars, x:float=-0.001) -> None:
    '''Viscosity map, looking down the x axis, at (-1,0,0) mm'''
    sv.hideAll()
    
    viscSlice(sv, [scaleDim(sv, x), 0, 0], [1,0,0], 'x')

def viscNozx(sv:stateVars, xs:float=-0.002412, **kwargs) -> None: 
    '''Shear stress map, looking down the x axis, at (-1,0,0) mm, scaled by nozzle_inner_width'''
    sv.hideAll()
    viscSlice(sv, [scaleDim(sv, xs), 0, 0], [-1,0,0], 'x', **kwargs)  
    
#--------------------------

def uSlice(sv:stateVars, origin:List[float], normal:List[float], view:str, name:str='U', umin:float=0, umax:float=0.02) -> None:
    '''Plot the velocity map. Segment out the ink and support separately and color separately, leaving some white space at the interface so you can see where the interface is. 
    origin is an [x,y,z] point, where the slice should be taken
    normal is an [x,y,z] vector, indicating the normal to the plane of the slice
    view is the name of the viewpoint, e.g. y,z,a,b
    name is the name of the variable, e.g. U, UX, UY, UZ
    umin, umax are the legend range'''
    slice1 = sliceMake(sv, origin, normal)
    d1, _ = inkClip(sv, slice1, name, 0, clipVal = 0.9)
    d2, _ = inkClip(sv, slice1, name, 1, clipVal = 0.1)
    for d in [d1, d2]:
        if name=='U':
            dispname = '|Velocity|'
        else:
            dispname = name[1]+' velocity'
        sv.colorBars.append(uColorBar(sv.renderView1, d, umax=umax, umin=umin, name=dispname))
    resetCam(sv, (sv.times[-1])) 
    setAndUpdate(view, sv)
    sv.renderView1.Update()
    
def uslicey(sv:stateVars) -> None:
    '''Velocity map, looking down the y axis, at (0,0,0)'''
    sv.hideAll()
    uSlice(sv, [0,0,0], [0, -1, 0], 'y')
    
def uslicex(sv:stateVars, x:float=-0.001) -> None:
    '''Velocity map, looking down the x axis, at (-1,0,0) mm'''
    sv.hideAll()
    uSlice(sv, [scaleDim(sv, x), 0, 0], [-1,0,0], 'x')  

def uzslicey(sv:stateVars) -> None:
    '''Velocity map, looking down the y axis, at (0,0,0)'''
    sv.hideAll()
    uSlice(sv, [0,0,0], [0, -1, 0], 'y', name='UZ', umin=-0.002, umax=0.002)
    
def uzslicex(sv:stateVars, x:float=-0.001) -> None:
    '''Velocity map, looking down the x axis, at (-1,0,0) mm'''
    sv.hideAll()
    uSlice(sv, [scaleDim(sv, x), 0, 0], [-1,0,0], 'x', name='UZ', umin=-0.002, umax=0.002)  
    
def uzsliceNozx(sv:stateVars, x:float=-0.002412) -> None:
    '''Velocity map, looking down the x axis, at (-1,0,0) mm'''
    sv.hideAll()
    uSlice(sv, [scaleDim(sv, x), 0, 0], [-1,0,0], 'x', name='UZ', umin=-0.002, umax=0.002)   
    
def uyslicey(sv:stateVars) -> None:
    '''Velocity map, looking down the y axis, at (0,0,0)'''
    sv.hideAll()
    uSlice(sv, [0,0,0], [0, -1, 0], 'y', name='UY', umin=-0.0006, umax=0.0006)
    
def uyslicex(sv:stateVars, x:float=-0.001) -> None:
    '''Velocity map, looking down the x axis, at (-1,0,0) mm'''
    sv.hideAll()
    uSlice(sv, [scaleDim(sv, x), 0, 0], [-1,0,0], 'x', name='UY', umin=-0.0006, umax=0.0006)  
    
def uysliceNozx(sv:stateVars, x:float=-0.002412) -> None:
    '''Velocity map, looking down the x axis, at (-1,0,0) mm'''
    sv.hideAll()
    uSlice(sv, [scaleDim(sv, x), 0, 0], [-1,0,0], 'x', name='UY', umin=-0.0006, umax=0.0006)   
    
#--------------------------

def alphaSlice(sv:stateVars, origin:List[float], normal:List[float], view:str, name:str='a') -> None: # RG
    '''Plot the ink and support, ink in white and support in black.
    name is the name of the variable, e.g. U, UX, UY, UZ'''
    slice1 = sliceMake(sv, origin, normal)
    di = inkClip(sv, slice1, name, 0, clipVal = 0.5)
    ds = inkClip(sv, slice1, name, 1, clipVal = 0.5)
    di.AmbientColor = [1.0, 1.0, 1.0]
    di.DiffuseColor = [1.0, 1.0, 1.0]
    ds.AmbientColor = [0.0, 0.0, 0.0]
    ds.DiffuseColor = [0.0, 0.0, 0.0]
    resetCam(sv, (sv.times[-1])) 
    setAndUpdate(view, sv)
    sv.renderView1.Update()
    
def alphaSlicex(sv:stateVars, **kwargs) -> None: # RG
    '''Alpha, looking down the x axis, at distance speficied in coloring'''
    sv.hideAll()
    c = kwargs['coloring'] # backpacking onto the color kwarg to get different x distances
    d = float(c.split('_')[1]) # get the x distance
    alphaSlice(sv, [d/1000, 0, 0], [-1,0,0], 'x')

#--------------------------

def pSlice(sv:stateVars, origin:List[float], normal:List[float], view:str, name:str='p_rgh', pmin:float=-100, pmax:float=100) -> None:
    '''Plot the pressure map. Segment out the ink and support separately and color separately, leaving some white space at the interface so you can see where the interface is. 
    origin is an [x,y,z] point, where the slice should be taken
    normal is an [x,y,z] vector, indicating the normal to the plane of the slice
    view is the name of the viewpoint, e.g. y,z,a,b
    name is the name of the variable, e.g. p, p_rgh, where p_rgh has the hydrostatic pressure removed
    pmin, pmax are the legend range
    '''
    slice1 = sliceMake(sv, origin, normal)
    d1, _ = inkClip(sv, slice1, name, 0, clipVal = 0.9)
    d2, _ = inkClip(sv, slice1, name, 1, clipVal = 0.1)
    for d in [d1, d2]:
        sv.colorBars.append(pColorBar(sv.renderView1, d, pmax=pmax, pmin=pmin))
    resetCam(sv, (sv.times[-1])) 
    setAndUpdate(view, sv)
    sv.renderView1.Update()
    
def pslicey(sv:stateVars, **kwargs) -> None:
    '''Velocity map, looking down the y axis, at (0,0,0)'''
    sv.hideAll()
    pSlice(sv, [0,0,0], [0, -1, 0], 'y', **kwargs)
    
def pslicex(sv:stateVars, xp:float=-0.001, **kwargs) -> None:
    '''Velocity map, looking down the x axis, at (-1,0,0) mm'''
    sv.hideAll()
    pSlice(sv, [scaleDim(sv, xp), 0, 0], [-1,0,0], 'x', **kwargs)
    
def psliceNozx(sv:stateVars, xp:float=-0.002412, **kwargs) -> None:
    '''Velocity map, looking down the x axis, at (-1,0,0) mm'''
    sv.hideAll()
    pSlice(sv, [scaleDim(sv, xp), 0, 0], [-1,0,0], 'x', **kwargs)
    
#--------------------------
    
def shearRateSlice(sv:stateVars, origin:List[float], normal:List[float], view:str, name:str='shearRate', rmin:float=0.1, rmax:float=1000, yieldSurface:bool=True) -> None:
    '''Plot the shear rate map. Segment out the ink and support separately and color separately, leaving some white space at the interface so you can see where the interface is.
    origin is an [x,y,z] point, where the slice should be taken
    normal is an [x,y,z] vector, indicating the normal to the plane of the slice
    view is the name of the viewpoint, e.g. y,z,a,b
    name is the name of the variable, e.g. p, p_rgh, where p_rgh has the hydrostatic pressure removed
    rmin, rmax are the legend range'''
    shearRate = computeShearRate(sv) # in paraview_general
    slice1 = sliceMake(sv, origin, normal, sinput=shearRate) # take slice from shear rate map
    d1, _ = inkClip(sv, slice1, name, 0, clipVal = 0.9)
    d2, _ = inkClip(sv, slice1, name, 1, clipVal = 0.1)
    for d in [d1, d2]:
        sv.colorBars.append(shearRateColorBar(sv.renderView1, d, rmax=rmax, rmin=rmin))
    if yieldSurface:
        addYieldSurface(sv, origin, normal)
    resetCam(sv, (sv.times[-1]))
    setAndUpdate(view, sv)
    sv.renderView1.Update()

def shearRateSlicey(sv:stateVars, **kwargs) -> None:
    '''Shear rate map, looking down the y axis, at (0,0,0)'''
    sv.hideAll()
    shearRateSlice(sv, [0,0,0], [0, -1, 0], 'y', **kwargs)
    
def shearRateSlicex(sv:stateVars, xs:float=-0.001, **kwargs) -> None:
    '''Shear rate map, looking down the x axis, at (-1,0,0) mm'''
    sv.hideAll()
    shearRateSlice(sv, [scaleDim(sv, xs), 0, 0], [-1,0,0], 'x', **kwargs)   
    
def shearRateSliceNozx(sv:stateVars, xs:float=-0.002412, **kwargs) -> None:
    '''Shear rate map, looking down the x axis, at (-1,0,0) mm'''
    sv.hideAll()
    shearRateSlice(sv, [scaleDim(sv, xs), 0, 0], [-1,0,0], 'x', **kwargs)   

#--------------------------

def addYieldSurface(sv:stateVars, origin:List[float], normal:List[float], name:str='shearStress') -> None:
    '''add a yield surface to the slice'''
    computeDerivatives = computeShearRate(sv) # calculate shear rate
    cellDatatoPointData = CellDatatoPointData(Input=computeDerivatives)
    cellDatatoPointData.CellDataArraytoprocess = ['ScalarGradient', 'U', 'VectorGradient', 'alpha.ink', 'cellID', 'nu1', 'nu2', 'p', 'p_rgh', 'rAU']
    
    slice1 = sliceMake(sv, origin, normal, sinput=cellDatatoPointData)
    
    calculator = Calculator(Input=slice1)
    calculator.ResultArrayName = 'shearStress'

    calculator.Function = sv.stressFunc()
 #   calculator.Function = '(1000*nu1*alpha.ink+1000*nu2*(1-alpha.ink))*ScalarGradient' # multiply nu by shear rate to get shear stress
    Hide(slice1, sv.renderView1)
    
    calculator3 = Calculator(Input=calculator)
    calculator3.ResultArrayName = 'shearStressMag'
    calculator3.Function = 'mag(shearStress)'

    dink, clipink = inkClip(sv, calculator3, name, 0, clipVal = 0.9)   # ink
    dsup, clipsup = inkClip(sv, calculator3, name, 1, clipVal = 0.1)   # support

    if sv.fs.ink.transportModel=='HerschelBulkley':
        tauink = sv.fs.ink.dynamic['tau0']
        inkyield = yieldClip(sv, clipink, tauink)
    if sv.fs.sup.transportModel=='HerschelBulkley':
        tausup = sv.fs.sup.dynamic['tau0']
        supyield = yieldClip(sv, clipsup, tausup)
    Hide(clipink, sv.renderView1)
    Hide(clipsup, sv.renderView1)
    

def shearStressSlice(sv:stateVars, origin:List[float], normal:List[float], view:str, name:str='shearStress', yieldSurface:bool=True, **kwargs) -> None: # RG
    '''Plot the shear stress map. Segment out the ink and support separately and color separately, leaving some white space at the interface so you can see where the interface is.
    origin is an [x,y,z] point, where the slice should be taken
    normal is an [x,y,z] vector, indicating the normal to the plane of the slice
    view is the name of the viewpoint, e.g. y,z,a,b
    name is the name of the value to clip
    rmin and rmax are the range of values in the legend
    '''

    computeDerivatives = computeShearRate(sv) # calculate shear rate
    cellDatatoPointData = CellDatatoPointData(Input=computeDerivatives)
    cellDatatoPointData.CellDataArraytoprocess = ['ScalarGradient', 'U', 'VectorGradient', 'alpha.ink', 'cellID', 'nu1', 'nu2', 'p', 'p_rgh', 'rAU']
    
    slice1 = sliceMake(sv, origin, normal, sinput=cellDatatoPointData)
    
    calculator = Calculator(Input=slice1)
    calculator.ResultArrayName = 'shearStress'

    calculator.Function = sv.stressFunc()   # this gives stress in Pa
 #   calculator.Function = '(1000*nu1*alpha.ink+1000*nu2*(1-alpha.ink))*ScalarGradient' # multiply nu by shear rate to get shear stress
    sv.hideAll()
    
    calculator3 = Calculator(Input=calculator)
    calculator3.ResultArrayName = 'shearStressMag'
    calculator3.Function = 'mag(shearStress)'

    dink, clipink = inkClip(sv, calculator3, name, 0, clipVal = 0.9)   # ink
    dsup, clipsup = inkClip(sv, calculator3, name, 1, clipVal = 0.1)   # support

    if yieldSurface:
        if sv.fs.ink.transportModel=='HerschelBulkley':
            tauink = sv.fs.ink.dynamic['tau0']
            inkyield = yieldClip(sv, clipink, tauink)
        if sv.fs.sup.transportModel=='HerschelBulkley':
            tausup = sv.fs.sup.dynamic['tau0']
            supyield = yieldClip(sv, clipsup, tausup)
    
    for d in [dink, dsup]:
        sv.colorBars.append(shearStressColorBar(sv.renderView1, d))
    resetCam(sv, (sv.times[-1]))
    setAndUpdate(view, sv)
    sv.renderView1.Update()
    return calculator


def shearStressSlicey(sv:stateVars, **kwargs) -> None: # RG
    '''Shear stress map, looking down the y axis, at (0,0,0)'''
    sv.hideAll()
    shearStressSlice(sv, [0,0,0], [0, -1, 0], 'y', **kwargs)
    
def shearStressSlicex(sv:stateVars, xs:float=-0.001, **kwargs) -> None: # RG
    '''Shear stress map, looking down the x axis, at (-1,0,0) mm, scaled by nozzle_inner_width'''
    sv.hideAll()
    shearStressSlice(sv, [scaleDim(sv, xs), 0, 0], [1,0,0], 'x', **kwargs)  
    
def shearStressSliceNozx(sv:stateVars, xs:float=-0.002412, **kwargs) -> None: 
    '''Shear stress map, looking down the x axis, at (-1,0,0) mm, scaled by nozzle_inner_width'''
    sv.hideAll()
    shearStressSlice(sv, [scaleDim(sv, xs), 0, 0], [-1,0,0], 'x', **kwargs)  

#--------------------------

def mesh(sv:stateVars) -> None:
    '''cross-section on the y=0 plane with mesh overlay and volume fraction coloring'''
    sv.hideAll()
    slice1 = Slice(Input=sv.caseVTMSeries)
    slice1.SliceType = 'Plane'
    slice1.SliceOffsetValues = [0.0]
    Hide3DWidgets(proxy=slice1.SliceType)
    slice1.SliceType.Normal = [0, -1, 0] # look down the y axis

    slice1Display = Show(slice1, sv.renderView1, 'GeometryRepresentation')
    slice1Display.SetScaleArray = ['POINTS', 'U']
    slice1Display.OpacityArray = ['POINTS', 'U']

    # show mesh edges
    slice1Display.SetRepresentationType('Surface With Edges')
    setAndUpdate('y', sv)
    sv.colorBars.append(alphaColorBar(sv.renderView1, slice1Display))   
    return 
 

def plainslice(sv:stateVars):
    '''colorless slice at the y=0 plane'''
    # create a new 'Slice'
    slice2 = Slice(Input=sv.caseVTMSeries)
    slice2.SliceType = 'Plane'
    slice2.HyperTreeGridSlicer = 'Plane'
    slice2.SliceOffsetValues = [0.0]
    Hide3DWidgets(proxy=slice2.SliceType)

    # Properties modified on slice1.SliceType
    slice2.SliceType.Normal = [0.0, -1.0, 0.0]

    # show data in view
    slice2Display = Show(slice2, sv.renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    slice2Display.ColorArrayName = [None, '']
    # turn off scalar coloring
    ColorBy(slice2Display, None)

    # change solid color
    slice2Display.AmbientColor = [1.0, 1.0, 1.0]
    slice2Display.DiffuseColor = [1.0, 1.0, 1.0]
    
    setAndUpdate('y', sv)
    
    return slice2
    
def yarrows(slice2, sv):
    '''Velocity arrows along the y=0 plane'''
    
    # create a new 'Glyph'
    glyph1 = Glyph(Input=slice2, GlyphType='Arrow')
    glyph1.OrientationArray = ['POINTS', 'U']
    glyph1.ScaleArray = ['POINTS', 'U']
    glyph1.ScaleFactor = 0.0006029999814927578
    glyph1.GlyphTransform = 'Transform2'
    glyph1.MaximumNumberOfSamplePoints = 2000
    glyph1.ScaleFactor = 0.03

    # show data in view
    glyph1Display = Show(glyph1, sv.renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    sv.colorBars.append(uColorBar(sv.renderView1, glyph1Display, umax=0.02))
    setAndUpdate('y', sv)
    
def sliceandyarrows(sv):
    '''Make a plain slice at the y=0 plane so you can see where the nozzle is, and add arrows for the velocity field'''
    sv.hideAll()
    slice2 =  plainslice(sv)
    yarrows(slice2, sv)
    
def tube(sv):
    '''Streamlines at the given z position tubeh'''
    try:
        streamTracer1 = StreamTracer(registrationName='StreamTracer1', Input=sv.caseVTMSeries, SeedType='Line') # RG
    except:
        streamTracer1 = StreamTracer(registrationName='StreamTracer1', Input=sv.caseVTMSeries)
        # unclear why SeedType Line throws an exception in PV 5.10.0, but it may have to do with version control.
    streamTracer1.Vectors = ['POINTS', 'U']
    streamTracer1.MaximumStreamlineLength = 0.01 # was 0.006 RG
    blc = sv.fs.geo.bath_left_coord/1000
    brc = sv.fs.geo.bath_right_coord/1000
    bfc = sv.fs.geo.bath_front_coord/1000
    bbc = sv.fs.geo.bath_back_coord/1000
    di = sv.fs.geo.nozzle_inner_width
    streamTracer1.SeedType.Point1 = [blc, bfc, sv.tubeh*di/0.6]
    streamTracer1.SeedType.Point2 = [brc, bbc, sv.tubeh*di/0.6]
    streamTracer1.SeedType.Resolution = 100
    # toggle 3D widget visibility (only when running from the GUI)
    Hide3DWidgets(proxy=streamTracer1.SeedType)
    # show data in view
    streamTracer1Display = Show(streamTracer1, sv.renderView1, 'GeometryRepresentation')
    # trace defaults for the display properties.
    streamTracer1Display.ColorArrayName = [None, '']
    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    streamTracer1Display.ScaleTransferFunction.Points = [-3.2023908359124786, 0.0, 0.5, 0.0, 3.5858197468415374, 1.0, 0.5, 0.0]
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    streamTracer1Display.OpacityTransferFunction.Points = [-3.2023908359124786, 0.0, 0.5, 0.0, 3.5858197468415374, 1.0, 0.5, 0.0]
    # update the view to ensure updated data information
    sv.renderView1.Update()
    # create a new 'Tube'
    tube1 = Tube(Input=streamTracer1)
    tube1.Scalars = ['POINTS', 'AngularVelocity']
    tube1.Vectors = ['POINTS', 'Normals']
    tube1.Radius = 1e-05*di/0.6
    # show data in view
    tube1Display = Show(tube1, sv.renderView1, 'GeometryRepresentation')
    ColorBy(tube1Display, ('POINTS', 'U', 'Z'))

    Hide(streamTracer1, sv.renderView1)
    sv.renderView1.Update()

    sv.colorBars.append(uColorBar(sv.renderView1, tube1Display, umax=0.002, umin=-0.002, name='z velocity'))

    # update the view to ensure updated data information
    sv.renderView1.Update()
    
def surfandtube(sv):
    '''Filament surface, with streamlines'''
    sv.hideAll()
    plainalphasurface(sv, opacity=0.5)
    tube(sv)
    
    
################ scripting #########################


def runThrough(v:ssVars, sv:stateVars, overwrite:bool=False) -> None:
    '''For a given folder, generate all the images'''
    if v.tag=='tubes':
        sv.tubeh = v.tubeh
 
    if len(v.tlist)>0:
        tlist = [int(round(10*i)) for i in v.tlist if i in sv.times]
    else:
        tlist = [int(round(10*i)) for i in sv.times]
    
    tlist2 = []
    fnlist = []
    viewlist = []
    # get a list of files to create

    
    for view in v.volList:
        for t in tlist:
            fn = fullFN(t, view, v.coloring, sv.folder)
                # find the full file name for this set of circumstances
            if (not os.path.exists(fn) or overwrite) and t/10 in sv.times:
                # if the file does not exist, add this set to the list
                tlist2.append(t)
                fnlist.append(fn)
                viewlist.append(view)
                
    # if we already have all the files we need, don't run this
    if len(tlist2)==0:
        return sv

    if not sv.initialized:
        try:
            initializeSV(sv)
        except Exception as e:
            logging.error(f'InitializeSV exception in {sv.folder}: {e}')
            raise

    try:
        f = v.function
        f(sv)                                  # create the image
        Show(sv.timeStamp, sv.renderView1)     # add time stamp
    except Exception as e:
        logging.error(f'Create image exception in {sv.folder}: {e}')
        return sv

    setTime(sv.times[-1]/10, sv) # for some reason we need to set this twice, or it will position the image wrong
    
    # iterate through times and take snapshots of the surfaces

    for i,t in enumerate(tlist2):
        logging.info(f'--t {t/10}, view: {viewlist[i]}, file: {fnlist[i]}' )
        setTime(t/10, sv)
        setView(viewlist[i], sv)
        SaveScreenshot(fnlist[i], sv.renderView1)
 
    return sv
         
def folderScript(folder:str, ssvList:List[ssVars], overwrite:bool=False):
    '''Initialize the folder and generate all images.'''
    try:
        if not os.path.exists(folder):
            return
        sv = stateVars(folder)
        sv.times = sv.fs.fh1.times()     
        for ssv in ssvList:
            sv = runThrough(ssv, sv, overwrite=overwrite)
    except Exception as e:
        logging.error(f'{folder}: {e}')
        traceback.print_exc()
        return
    cleanSession()
    return    
