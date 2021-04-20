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
from paraview.simple import * # import the simple module from the paraview

# local packages
from paraview_general import *

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

FONT = 24

##############################################
######## SCREENSHOTS

####### file operations ######
    

def timeStr(t:float) -> str:
    '''convert a time to a 3 digit string'''
    tstr = str(t)
    if t<10:
        tstr = "00" + tstr
    elif t<100:
        tstr = "0" + tstr
    return tstr


def fullFN(t:float, view:str, coloring:str, folder:str):
    '''get a full file name for an image'''
    fn = 't'+timeStr(t)+'_'+view+'_'+coloring+'.png'
    filename = os.path.join(folder, 'images', fn)
    return filename


# def vidFN(view:str, coloring:str, folder:str):
#     '''get the full file name for the video'''
#     fn = view+'_'+coloring+'.mp4'
#     filename = os.path.join(folder, 'images', fn)
#     return filename



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
    annotateTimeFilter1.Format = '%2.1f s'
    annotateTimeFilter1Display = Show(annotateTimeFilter1, sv.renderView1, 'TextSourceRepresentation')
    annotateTimeFilter1Display.FontSize = FONT
    sv.timeStamp = annotateTimeFilter1
    return sv


###############
# camera operations


def resetCam(sv:stateVars, time:float) -> None:
    '''go to a specific view area, at a given time'''
    setTime(time, sv)
    sv.renderView1.ResetCamera(-0.005, 0.005, -0.002, 0.002, -0.002, 0.002)

def setView(st:str, sv:stateVars) -> None:
    '''go to a specific view point, out of [x,y,z,a,b]. a and b are views from an angle, where a is downstream of the nozzle, and b is upstream of the nozzle.'''
    if st=="z":
        sv.renderView1.CameraFocalPoint = [0.001, 0,0]
        sv.renderView1.CameraPosition = [0.001, 0, 10]
        sv.renderView1.CameraViewUp = [0,1,0]
    elif st=="y":
        sv.renderView1.CameraFocalPoint = [0.001, 0,0]
        sv.renderView1.CameraPosition = [0.001, -10, 0]
        sv.renderView1.CameraViewUp = [0,0,1]
    elif st=="x":
        sv.renderView1.CameraFocalPoint = [0.001, 0,0]
        sv.renderView1.CameraPosition = [-10, 0.001, 0]
        sv.renderView1.CameraViewUp = [0,0,1]
    elif st=="a":
        sv.renderView1.CameraFocalPoint = [0.0006, -0.0002, 0.0005]
        sv.renderView1.CameraPosition = [2, 1, 2]
        sv.renderView1.CameraViewUp = [0,0,1]
    elif st=="b":
        sv.renderView1.CameraFocalPoint = [0.001, 0, 0.001]
        sv.renderView1.CameraPosition = [-1,-1,1]
        sv.renderView1.CameraViewUp = [0,0,1]


def setAndUpdate(st:str, sv:stateVars) -> None:
    '''go to a specific viewpoint and view area'''
    setView(st, sv)
    sv.renderView1.Update()
    sv.renderView1.ResetCamera(-0.005, 0.005, -0.002, 0.002, -0.002, 0.002)


#### annotations


def positionCB(ColorBar) -> None:
    '''put the color bar in the bottom'''
    ColorBar.AutoOrient = 0
    ColorBar.Orientation = 'Horizontal'
    ColorBar.WindowLocation = 'LowerCenter'
    ColorBar.ScalarBarLength = 0.7
    ColorBar.UseCustomLabels = 1
    ColorBar.TitleFontSize = FONT
    ColorBar.LabelFontSize = FONT
    ColorBar.ComponentTitle = ''


def uColorBar(renderView1, display, umax:float=0.01, umin:float=0, name:str='|Velocity|') -> None:
    '''velocity magnitude color bar'''
    # set scalar coloring
    #ColorBy(display, ('POINTS', 'U', 'Magnitude'))

    # show color bar/color legend
    display.SetScalarBarVisibility(renderView1, True)
    uLUT = GetColorTransferFunction('U')
    uPWF = GetOpacityTransferFunction('U')

    uLUT.RescaleTransferFunction(umin, umax)
    uPWF.RescaleTransferFunction(umin, umax)
    if (umax-umin)>0.015:
        du = 0.005
    elif (umax-umin)>0.006:
        du = 0.002
    elif (umax-umin)>0.003:
        du = 0.001
    else:
        du = 0.0001

    uLUTColorBar = GetScalarBar(uLUT, renderView1)
    positionCB(uLUTColorBar)
    uLUTColorBar.CustomLabels = np.arange(umin, umax+du, du)
    uLUTColorBar.Title = name+' (m/s)'
    uLUTColorBar.RangeLabelFormat = '%-#0.3f'
    return uLUTColorBar

def pColorBar(renderView1, display, pmax:float=100, pmin:float=-100, name:str='p_rgh') -> None:
    '''velocity magnitude color bar'''
    # set scalar coloring
    #ColorBy(display, ('POINTS', 'U', 'Magnitude'))

    # show color bar/color legend
    display.SetScalarBarVisibility(renderView1, True)
    pLUT = GetColorTransferFunction(name)
    pPWF = GetOpacityTransferFunction(name)

    pLUT.RescaleTransferFunction(pmin, pmax)
    pPWF.RescaleTransferFunction(pmin, pmax)
    if pmax==100 and pmin==-100:
        dp = 50
    else:
        dp = round((pmax-pmin)/4)

    pLUTColorBar = GetScalarBar(pLUT, renderView1)
    positionCB(pLUTColorBar)
    if name=='p_rgh':
        pLUTColorBar.Title = 'Adjusted pressure (Pa)'
    else:
        pLUTColorBar.Title = 'Pressure (Pa)'
    pLUTColorBar.CustomLabels = np.arange(pmin, pmax+dp, dp)
    return pLUTColorBar

def nuColorBar(renderView1, display, color) -> None:
    '''viscosity color bar'''
    # show color bar/color legend
    display.SetScalarBarVisibility(renderView1, True)
    uLUT = GetColorTransferFunction(color)
    uPWF = GetOpacityTransferFunction(color)

    # rescale to full range of newt, HB experiments
    uLUT.RescaleTransferFunction(10**-6, 10**2)
    uPWF.RescaleTransferFunction(10**-6, 10**2)
    
    # convert to log space. If we've already converted to log space during creation of another image, skip this
    if not uLUT.UseLogScale==1:
        uLUT.MapControlPointsToLogSpace()
        uLUT.UseLogScale = 1

    uLUTColorBar = GetScalarBar(uLUT, renderView1)
    positionCB(uLUTColorBar)
    uLUTColorBar.CustomLabels = [10**i for i in range(-6, 3)]
    uLUTColorBar.Title = 'Viscosity (m^2/s)' # if the density is 1000 kg/m^3, the viscosity is in  MPa.s
 #     uLUTColorBar.RangeLabelFormat = '%-#2.0e'
    return uLUTColorBar
    

def alphaColorBar(renderView1, display) -> None:
    '''volume fraction color bar'''
    # set scalar coloring
    ColorBy(display, ('CELLS', 'alpha.ink'))
    ColorBy(display, ('POINTS', 'alpha.ink'))

    # rescale color and/or opacity maps used to include current data range
    display.RescaleTransferFunctionToDataRange(True, False)
    display.SetScalarBarVisibility(renderView1, True)
    alphainkLUT = GetColorTransferFunction('alphaink')
    alphainkLUT.RescaleTransferFunction(0.0, 1)
    alphainkPWF = GetOpacityTransferFunction('alphaink')

    aLUTColorBar = GetScalarBar(alphainkLUT, renderView1)
    positionCB(aLUTColorBar)
    aLUTColorBar.CustomLabels = np.arange(0, 1.25, 0.25)
    aLUTColorBar.Title = 'Volume fraction ink'
    aLUTColorBar.RangeLabelFormat = '%-#0.2f'
    return aLUTColorBar


################ types of surfaces #########################


def vol(sv:stateVars) -> None:
    '''filament surface with velocity coloring'''
    sv.hideAll()
    alphasurface(sv, 'U')


def plainalphasurface(sv:stateVars) -> None:
    '''filament surface with no coloring'''
    sv.hideAll()
    alphasurface(sv, 'None')
    
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
    
def viscColor(display, visc:str) -> None:
    '''This is a preset list of colors for viscosity, where the kinematic viscosity range is [10**i for i in range(-6, 3)] m^2/s'''
    colordict = {0.000001:[59, 76, 192],\
                 0.00001:[97, 128, 233],\
                 0.0001:[140, 175, 254],\
                 0.001:[183, 207, 249], \
                 0.01:[220, 221, 221], \
                 0.1:[244, 197, 173], \
                 1:[244, 153, 122],\
                 10:[222, 97, 77],\
                 100:[180, 4, 38]}
    try:
        c = colordict[float(visc)]
    except:
        raise Exception('Viscosity not in list')
    c = [i/255 for i in c]
    display.AmbientColor = c
    display.DiffuseColor = c
    
    
def inkClip(sv:stateVars, clipinput, colVar:str, invert:int, **kwargs) -> None:
    '''Just clip out the ink portion and color it by viscosity or velocity. clipinput is the Paraview object that you're trying to clip, e.g. sv.caseVTMSeries. colVar is the variable being colored ('U', 'nu1', 'nu2', 'inknu', or 'supnu'). 'nu1' and 'nu2' use local viscosities measured for non-Newtonian fluids. 'inknu' and 'supnu' are constant viscosities collected from the legend, for Newtonian fluids. invert 0 to get ink, 1 to get support. A clipVal in kwargs allows you to custom set the alpha value where you want to clip.'''
    
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
    else:
#         clipDisplay.SetScaleArray = 'None'
#         ColorBy(clipDisplay, None)
        c = [0.9, 0.9, 0.9]
        clipDisplay.AmbientColor = c
        clipDisplay.DiffuseColor = c
        
    Hide(clipinput, sv.renderView1)
    
    return clipDisplay


def alphasurface(sv:stateVars, colVar:str) -> None:
    '''filament surface
    colVar is the variable being colored ('U' or 'None')'''
    clipDisplay = inkClip(sv, sv.caseVTMSeries, colVar, 0)
    sv.colorBars.append(uColorBar(sv.renderView1, clipDisplay, umax=0.01))
    resetCam(sv, (sv.times[-1])) 
    setAndUpdate('y', sv)
    sv.renderView1.Update()
    return 

def sliceMake(sv:stateVars, origin:List[float], normal:List[float]):
    '''Get the slice object, to be further manipulated'''
    slice1 = Slice(Input=sv.caseVTMSeries)
    slice1.SliceType = 'Plane'
    slice1.SliceOffsetValues = [0.0]
    Hide3DWidgets(proxy=slice1.SliceType)
    slice1.SliceType.Origin = origin
    slice1.SliceType.Normal = normal # look down the y axis
    slice1Display = Show(slice1, sv.renderView1, 'GeometryRepresentation')
    slice1Display.SetScaleArray = ['POINTS', 'U'] # this doesn't matter because we're going to recolor later
    slice1Display.OpacityArray = ['POINTS', 'U']
    Hide(sv.caseVTMSeries, sv.renderView1)
    return slice1

def viscSlice(sv:stateVars, origin:List[float], normal:List[float], view:str) -> None:
    '''Plot the viscosity map. Segment out the ink and support separately and color separately, leaving some white space at the interface so you can see where the interface is.'''
    out = sv.readVisc()
    if out>0:
        return
    slice1 = sliceMake(sv, origin, normal)
    inkClip(sv, slice1, sv.inkfunc, 0)
    inkClip(sv, slice1, sv.supfunc, 1)

    resetCam(sv, (sv.times[-1])) 
    setAndUpdate(view, sv)
    sv.renderView1.Update()
    
def viscy(sv:stateVars) -> None:
    '''Viscosity map, looking down the y axis, at (0,0,0)'''
    sv.hideAll()
    viscSlice(sv, [0,0,0], [0, -1, 0], 'y')
    
def viscx(sv:stateVars) -> None:
    '''Viscosity map, looking down the x axis, at (-1,0,0) mm'''
    sv.hideAll()
    viscSlice(sv, [-0.001, 0, 0], [1,0,0], 'x')
    
def uSlice(sv:stateVars, origin:List[float], normal:List[float], view:str, name:str='U', umin:float=0, umax:float=0.02) -> None:
    '''Plot the viscosity map. Segment out the ink and support separately and color separately, leaving some white space at the interface so you can see where the interface is. name is the name of the variable, e.g. U, UX, UY, UZ'''
    slice1 = sliceMake(sv, origin, normal)
    d1 = inkClip(sv, slice1, name, 0, clipVal = 0.9)
    d2 = inkClip(sv, slice1, name, 1, clipVal = 0.1)
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
    
def uslicex(sv:stateVars) -> None:
    '''Velocity map, looking down the x axis, at (-1,0,0) mm'''
    sv.hideAll()
    uSlice(sv, [-0.001, 0, 0], [1,0,0], 'x')  

def uzslicey(sv:stateVars) -> None:
    '''Velocity map, looking down the y axis, at (0,0,0)'''
    sv.hideAll()
    uSlice(sv, [0,0,0], [0, -1, 0], 'y', name='UZ', umin=-0.002, umax=0.002)
    
def uzslicex(sv:stateVars) -> None:
    '''Velocity map, looking down the x axis, at (-1,0,0) mm'''
    sv.hideAll()
    uSlice(sv, [-0.001, 0, 0], [1,0,0], 'x', name='UZ', umin=-0.002, umax=0.002)  


def pSlice(sv:stateVars, origin:List[float], normal:List[float], view:str, name:str='p_rgh', pmin:float=-100, pmax:float=100) -> None:
    '''Plot the viscosity map. Segment out the ink and support separately and color separately, leaving some white space at the interface so you can see where the interface is. name is the name of the variable, e.g. p, p_rgh, where p_rgh has the hydrostatic pressure removed'''
    slice1 = sliceMake(sv, origin, normal)
    d1 = inkClip(sv, slice1, name, 0, clipVal = 0.9)
    d2 = inkClip(sv, slice1, name, 1, clipVal = 0.1)
    for d in [d1, d2]:
        sv.colorBars.append(pColorBar(sv.renderView1, d, pmax=pmax, pmin=pmin))
    resetCam(sv, (sv.times[-1])) 
    setAndUpdate(view, sv)
    sv.renderView1.Update()
    
def pslicey(sv:stateVars) -> None:
    '''Velocity map, looking down the y axis, at (0,0,0)'''
    sv.hideAll()
    pSlice(sv, [0,0,0], [0, -1, 0], 'y')
    
def pslicex(sv:stateVars) -> None:
    '''Velocity map, looking down the x axis, at (-1,0,0) mm'''
    sv.hideAll()
    pSlice(sv, [-0.001, 0, 0], [1,0,0], 'x')
 
    


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
    
    streamTracer1 = StreamTracer(Input=sv.caseVTMSeries, SeedType='High Resolution Line Source')
    streamTracer1.Vectors = ['POINTS', 'U']
    streamTracer1.MaximumStreamlineLength = 0.006
    streamTracer1.SeedType.Point1 = [-0.003014999907463789, -0.0021104998886585236, sv.tubeh]
    streamTracer1.SeedType.Point2 = [0.003014999907463789, 0.0021104998886585236, sv.tubeh]
    streamTracer1.SeedType.Resolution = 50
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
    tube1.Radius = 1e-05
    # show data in view
    tube1Display = Show(tube1, sv.renderView1, 'GeometryRepresentation')
    ColorBy(tube1Display, ('POINTS', 'U', 'Magnitude'))

    Hide(streamTracer1, sv.renderView1)
    sv.renderView1.Update()

    sv.colorBars.append(uColorBar(sv.renderView1, tube1Display, umax=0.01))

    # update the view to ensure updated data information
    sv.renderView1.Update()
    
def surfandtube(sv):
    '''Filament surface, with streamlines'''
    sv.hideAll()
    plainalphasurface(sv)
    tube(sv)
    
    
################ scripting #########################



class ssVars:
    '''Holds input variables that tell us what kind of images to generate'''
    
    
    def __init__(self, tag, tlist, **kwargs):
        if tag=='volumes':
            if 'volViewList' not in kwargs:
                self.volList = ['y']
            else:
                self.volList = kwargs['volViewList']
            self.coloring = 'umag'
            self.function = vol
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
            self.function = viscx
            self.volList = ['x']
        elif tag=='uslicey':
            self.coloring = 'uslicey'
            self.function = uslicey
            self.volList = ['y']
        elif tag=='uslicex':
            self.coloring = 'uslicex'
            self.function = uslicex
            self.volList = ['x']
        elif tag=='uzslicey':
            self.coloring = 'uzslicey'
            self.function = uzslicey
            self.volList = ['y']
        elif tag=='uzslicex':
            self.coloring = 'uzslicex'
            self.function = uzslicex
            self.volList = ['x']
        elif tag=='py':
            self.coloring='py'
            self.function=pslicey
            self.volList = ['y']
        elif tag=='px':
            self.coloring = 'px'
            self.function = pslicex
            self.volList = ['x']
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

def runThrough(v:ssVars, sv:stateVars) -> None:
    '''For a given folder, generate all the images'''
    if v.tag=='tubes':
        sv.tubeh = v.tubeh
    
    if len(v.tlist)>0:
        tlist = [int(round(10*i)) for i in v.tlist]
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
            if (not os.path.exists(fn)) and t/10 in sv.times:
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

    if tlist2[0]==sv.times[-1]:
        setTime(sv.times[0]/10, sv) # for some reason we need to set this twice, or it will position the imagewrong
    else:
        setTime(sv.times[-1]/10, sv) # for some reason we need to set this twice, or it will position the imagewrong
    
    # iterate through times and take snapshots of the surfaces

    for i,t in enumerate(tlist2):
        logging.info(f'--t {t/10}, view: {viewlist[i]}, file: {fnlist[i]}' )
        setTime(t/10, sv)
        setView(viewlist[i], sv)
        SaveScreenshot(fnlist[i], sv.renderView1)
 
    return sv
         
def folderScript(folder:str, ssvList:List[ssVars]):
    '''Initialize the folder and generate all images.'''
    try:
        if not os.path.exists(folder):
            return
        sv = stateVars(folder)
        sv.times = fp.times(folder)        
        for ssv in ssvList:
            sv = runThrough(ssv, sv)
    except Exception as e:
        logging.error(folder+': '+str(e))
        return
    cleanSession()
    return    
