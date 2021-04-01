# trace generated using paraview version 5.8.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

import os
import numpy as np
from paraview.simple import * # import the simple module from the paraview
import re

# t0 = 5
# tf = 55
# dt = 5
# t0 = 5
# tf = 250
# dt = 50
# tlist = range(t0, tf, dt)
tlist = [10, 25, 50] # seconds * 10
modes = [0,1] # 0 for frames, 1 for gif frames

volumes = True
volviewlist = ['y']
#volviewlist = ['a']

meshes = False

vectors = False

tubes = False
tubeh = 0.001

folders = []
# topfolder = 'E:\\Leanne\\OpenFOAM\\newtHBsweep'
# flist = range(344, 350)
# for f in flist:
#     folders.append(os.path.join(topfolder, 'nb'+str(f)))
topfolder = r'C:\Users\lmf1\Documents\OpenFOAM\HBHBsweep'
flist = range(347, 474)
for f in flist:
    folders.append(os.path.join(topfolder, 'nb'+str(f)))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                             FUNCTIONS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
####### file operations ######

def ss_mkdirif(path):
    try:
        os.mkdir(path)
    except OSError:
        return 0
    else:
        print ("Created directory %s" % path)
        
def ss_mksubdirs(folder):
    mkdirif(os.path.join(folder, 'images'))   
    mkdirif(os.path.join(folder, 'images', 'frames'))  
    mkdirif(os.path.join(folder, 'images', 'videos')) 
        
def ss_series(folder):
    cf = casefolder(folder)
    vtkfolder = os.path.join(cf, 'VTK')
    if os.path.exists(vtkfolder):
        for file in os.listdir(vtkfolder):
            if '.vtm.series' in file:
                return os.path.join(vtkfolder, file)
    return ""

def ss_casefolder(folder):
    # if there is a folder in this folder named 'case', return that folder
    casefold = os.path.join(folder, 'case')
    if os.path.exists(casefold):
        return casefold
    # if this folder contains a folder called 'constant', this is the case folder, so return this folder
    constfold = os.path.join(folder, 'constant')
    if os.path.exists(constfold):
        return folder
    vtkfold = os.path.join(folder, 'VTK')
    if os.path.exists(vtkfold):
        return folder
    intfold = os.path.join(folder, 'interfacePoints')
    if os.path.exists(intfold):
        return folder
    else:
        #print('no case folder in ', folder)
        return ''
    
def ss_timestr(t):
    tstr = str(t)
    if t<10:
        tstr = "00" + tstr
    elif t<100:
        tstr = "0" + tstr
    return tstr

def ss_fullFN(t, view, coloring, folder, subfolder):
    fn = 't'+timestr(t)+'_'+view+'_'+coloring+'.png'
    filename = os.path.join(folder, 'images', subfolder, fn)
    return filename

def ss_vidFN(view, coloring, folder):
    fn = view+'_'+coloring+'.mp4'
    filename = os.path.join(folder, 'images', 'videos', fn)
    return filename

######### show/hide/initializes ############

def ss_initializeAll(folder):
    print('Initializing paraview for ', folder)
    ResetSession()
    sv = ss_stateVars(folder)
        # make subdirectories
    sv =  ss_initializeP(sv) 
        # initialize Paraview
    try:
        sv =  ss_initSeries(sv)
        # import the vtk series files
    except NameError as err:
        raise err
    sv = ss_timeStamp(sv)
        # add time stamp
    sv.initialized = True
    return sv

def ss_initializeSV(sv):
    # use this if you already have a stateVars object
    print('Initializing paraview for ', sv.folder)
    ResetSession()
    sv =  ss_initializeP(sv) 
        # initialize Paraview
    sv =  ss_initSeries(sv)
        # import the vtk series files
    sv = ss_timeStamp(sv)
        # add time stamp
    sv.initialized = True
    return 


#### hide everything
def ss_hideAll():
    ss = GetSources()
    for s in ss:
        Hide(ss[s])

        
#### initialize all of paraview
def ss_initializeP(sv):  
    ResetSession()     
    LoadPalette(paletteName='WhiteBackground')
        # make background white
    paraview.simple._DisableFirstRenderCameraReset()
        # disable automatic camera reset on 'Show'  
    ss_hideAll()
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
 
    
#### import the vtk series file
def ss_initSeries(sv):
    sfile = ss_series(sv.folder)
    if not os.path.exists:
        raise NameError('vtm.series file does not exist')
    casevtmseries = XMLMultiBlockDataReader(FileName=series(sv.folder))
    casevtmseries.CellArrayStatus = ['U', 'alpha.ink', 'p_rgh']
    casevtmseries.PointArrayStatus = ['U', 'alpha.ink', 'p_rgh']
    
    casevtmseriesDisplay = Show(casevtmseries, sv.renderView1, 'GeometryRepresentation')
    casevtmseriesDisplay.Representation = 'Surface'
    casevtmseriesDisplay.ColorArrayName = [None, '']
    casevtmseriesDisplay.OSPRayScaleArray = 'U'
    casevtmseriesDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    casevtmseriesDisplay.SelectOrientationVectors = 'None'
    casevtmseriesDisplay.ScaleFactor = 0.0006029999814927578
    casevtmseriesDisplay.SelectScaleArray = 'None'
    casevtmseriesDisplay.GlyphType = 'Arrow'
    casevtmseriesDisplay.GlyphTableIndexArray = 'None'
    casevtmseriesDisplay.GaussianRadius = 3.0149999074637892e-05
    casevtmseriesDisplay.SetScaleArray = ['POINTS', 'U']
    casevtmseriesDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    casevtmseriesDisplay.OpacityArray = ['POINTS', 'U']
    casevtmseriesDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    casevtmseriesDisplay.DataAxesGrid = 'GridAxesRepresentation'
    casevtmseriesDisplay.PolarAxes = 'PolarAxesRepresentation'
    casevtmseriesDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.009999999776482582, 1.0, 0.5, 0.0]
    casevtmseriesDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.009999999776482582, 1.0, 0.5, 0.0]
    sv.renderView1.ResetCamera()
    sv.renderView1.Update()
    ColorBy(casevtmseriesDisplay, None)
    times = casevtmseries.TimestepValues
    times2 = []
    for t in times:
        i = int(10*t)
        if i not in times2:
            times2.append(i)
    print(times2)
    sv.casevtmseries = casevtmseries
    sv.times = times2
    return sv

def ss_resetCam(sv, time):
    setTime(time, sv)
    sv.renderView1.ResetCamera(-0.005, 0.005, -0.002, 0.002, -0.002, 0.002)

    
def ss_setTime(time, sv):
    if time*10 in sv.times:
        sv.animationScene1.AnimationTime = time
        sv.timeKeeper1.Time = time
    else:
        print(f'time {time} is not in list {sv.times}')
    
def ss_setView(st, sv):
    if st=="z":
        sv.renderView1.CameraFocalPoint = [0.001, 0,0]
        sv.renderView1.CameraPosition = [0.001, 0, 10]
        sv.renderView1.CameraViewUp = [0,1,0]
    elif st=="y":
        sv.renderView1.CameraFocalPoint = [0.001, 0,0]
        sv.renderView1.CameraPosition = [0.001, -10, 0]
        sv.renderView1.CameraViewUp = [0,0,1]
    elif st=="a":
        sv.renderView1.CameraFocalPoint = [0.0006, -0.0002, 0.0005]
        sv.renderView1.CameraPosition = [2, 1, 2]
        sv.renderView1.CameraViewUp = [0,0,1]
    elif st=="b":
        sv.renderView1.CameraFocalPoint = [0.001, 0, 0.001]
        sv.renderView1.CameraPosition = [-1,-1,1]
        sv.renderView1.CameraViewUp = [0,0,1]
        
def ss_setAndUpdate(st, sv):
    ss_setView(st, sv)
    sv.renderView1.Update()
    sv.renderView1.ResetCamera(-0.005, 0.005, -0.002, 0.002, -0.002, 0.002)


#### annotations
FONT = 12

def ss_positionCB(ColorBar):
    ColorBar.AutoOrient = 0
    ColorBar.Orientation = 'Horizontal'
    ColorBar.WindowLocation = 'LowerCenter'
    ColorBar.ScalarBarLength = 0.7
    ColorBar.UseCustomLabels = 1
    ColorBar.TitleFontSize = FONT
    ColorBar.LabelFontSize = FONT
    ColorBar.ComponentTitle = ''

# velocity magnitude color bar
def ss_uColorBar(renderView1, display, umax):
    # set scalar coloring
    ColorBy(display, ('POINTS', 'U', 'Magnitude'))

    # show color bar/color legend
    display.SetScalarBarVisibility(renderView1, True)
    uLUT = GetColorTransferFunction('U')
    uPWF = GetOpacityTransferFunction('U')

    uLUT.RescaleTransferFunction(0.0, umax)
    uPWF.RescaleTransferFunction(0.0, umax)
    if umax>0.015:
        du = 0.005
    else:
        du = 0.002

    uLUTColorBar = GetScalarBar(uLUT, renderView1)
    ss_positionCB(uLUTColorBar)
    uLUTColorBar.CustomLabels = np.arange(0, umax+du, du)
    uLUTColorBar.Title = '|Velocity| (m/s)'
    uLUTColorBar.RangeLabelFormat = '%-#0.3f'
    
def ss_alphaColorBar(renderView1, display):
    # set scalar coloring
    ColorBy(display, ('CELLS', 'alpha.ink'))

    # rescale color and/or opacity maps used to include current data range
    display.RescaleTransferFunctionToDataRange(True, False)
    display.SetScalarBarVisibility(renderView1, True)
    alphainkLUT = GetColorTransferFunction('alphaink')
    alphainkLUT.RescaleTransferFunction(0.0, 1)
    alphainkPWF = GetOpacityTransferFunction('alphaink')

    aLUTColorBar = GetScalarBar(alphainkLUT, renderView1)
    ss_positionCB(aLUTColorBar)
    aLUTColorBar.CustomLabels = np.arange(0, 1.25, 0.25)
    aLUTColorBar.Title = 'Volume fraction ink'
    aLUTColorBar.RangeLabelFormat = '%-#0.2f'
    
def ss_timeStamp(sv):
    annotateTimeFilter1 = AnnotateTimeFilter(Input=sv.casevtmseries)
    annotateTimeFilter1.Format = '%2.1f s'
    annotateTimeFilter1Display = Show(annotateTimeFilter1, sv.renderView1, 'TextSourceRepresentation')
    annotateTimeFilter1Display.FontSize = FONT
    sv.timeStamp = annotateTimeFilter1
    return sv

################ types of surfaces #########################
       
def ss_vol(sv):
    ss_alphasurface(sv, 'U')

def ss_plainalphasurface(sv):
    ss_alphasurface(sv, 'None')

def ss_alphasurface(sv, color):
    # sv = stateVars object
    # color = 'U' or 'None'
    clip1 = Clip(Input=sv.casevtmseries)
    clip1.ClipType = 'Plane'
    clip1.HyperTreeGridClipper = 'Plane'
    clip1.Scalars = ['POINTS', 'alpha.ink']
    clip1.Value = 0.5
    clip1.ClipType = 'Scalar'
    clip1.Invert = 0

    # show data in view
    clip1Display = Show(clip1, sv.renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    clip1Display.Representation = 'Surface'
    clip1Display.ColorArrayName = [None, '']
    clip1Display.OSPRayScaleArray = 'U'
    clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    clip1Display.SelectOrientationVectors = 'None'
    clip1Display.ScaleFactor = 0.00018132181139662862
    clip1Display.GlyphType = 'Arrow'
    clip1Display.GlyphTableIndexArray = 'None'
    clip1Display.GaussianRadius = 9.066090569831431e-06
    
    if color=="U":
        clip1Display.SetScaleArray = ['POINTS', 'U']
        clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
        clip1Display.OpacityArray = ['POINTS', 'U']
        clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
        clip1Display.DataAxesGrid = 'GridAxesRepresentation'
        clip1Display.PolarAxes = 'PolarAxesRepresentation'
        clip1Display.ScalarOpacityUnitDistance = 5.84412949766512e-05
        clip1Display.ExtractedBlockIndex = 1
        clip1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]
        clip1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]
        ss_uColorBar(sv.renderView1, clip1Display, 0.01)
    else:
        clip1Display.SetScaleArray = 'None'
        ColorBy(clip1Display, None)
        
    Hide(sv.casevtmseries, sv.renderView1)

    ss_resetCam(sv, (sv.times[-1])/10) 
    ss_setAndUpdate('y', sv)
    sv.renderView1.Update()
    return [clip1, clip1Display]
    
    
def ss_mesh(sv):
    ss_hideAll()
    slice1 = Slice(Input=sv.casevtmseries)
    slice1.SliceType = 'Plane'
    slice1.SliceOffsetValues = [0.0]
    Hide3DWidgets(proxy=slice1.SliceType)
    slice1.SliceType.Normal = [0, -1, 0] # look down the y axis

    slice1Display = Show(slice1, sv.renderView1, 'GeometryRepresentation')
    slice1Display.Representation = 'Surface'
    slice1Display.ColorArrayName = [None, '']
    slice1Display.OSPRayScaleArray = 'U'
    slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    slice1Display.SelectOrientationVectors = 'None'
    slice1Display.ScaleFactor = 0.0006029999814927578
    slice1Display.SelectScaleArray = 'None'
    slice1Display.GlyphType = 'Arrow'
    slice1Display.GlyphTableIndexArray = 'None'
    slice1Display.GaussianRadius = 3.0149999074637892e-05
    slice1Display.SetScaleArray = ['POINTS', 'U']
    slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
    slice1Display.OpacityArray = ['POINTS', 'U']
    slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
    slice1Display.DataAxesGrid = 'GridAxesRepresentation'
    slice1Display.PolarAxes = 'PolarAxesRepresentation'
    slice1Display.ScaleTransferFunction.Points = [-0.003019144991412759, 0.0, 0.5, 0.0, 0.010859496891498566, 1.0, 0.5, 0.0]
    slice1Display.OpacityTransferFunction.Points = [-0.003019144991412759, 0.0, 0.5, 0.0, 0.010859496891498566, 1.0, 0.5, 0.0]

    # show mesh edges
    slice1Display.SetRepresentationType('Surface With Edges')
    ss_setAndUpdate('y', sv)
    ss_alphaColorBar(sv.renderView1, slice1Display)    
    return [slice1, slice1Display]
    
def ss_plainslice(sv):
    # create a new 'Slice'
    slice2 = Slice(Input=sv.casevtmseries)
    slice2.SliceType = 'Plane'
    slice2.HyperTreeGridSlicer = 'Plane'
    slice2.SliceOffsetValues = [0.0]
    Hide3DWidgets(proxy=slice2.SliceType)

    # Properties modified on slice1.SliceType
    slice2.SliceType.Normal = [0.0, -1.0, 0.0]

    # show data in view
    slice2Display = Show(slice2, sv.renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    slice2Display.Representation = 'Surface'
    slice2Display.ColorArrayName = [None, '']
    slice2Display.OSPRayScaleArray = 'U'
    slice2Display.OSPRayScaleFunction = 'PiecewiseFunction'
    slice2Display.SelectOrientationVectors = 'None'
    slice2Display.ScaleFactor = 0.0006029999814927578
    slice2Display.SelectScaleArray = 'None'
    slice2Display.GlyphType = 'Arrow'
    slice2Display.GlyphTableIndexArray = 'None'
    slice2Display.GaussianRadius = 3.0149999074637892e-05
    slice2Display.SetScaleArray = ['POINTS', 'U']
    slice2Display.ScaleTransferFunction = 'PiecewiseFunction'
    slice2Display.OpacityArray = ['POINTS', 'U']
    slice2Display.OpacityTransferFunction = 'PiecewiseFunction'
    slice2Display.DataAxesGrid = 'GridAxesRepresentation'
    slice2Display.PolarAxes = 'PolarAxesRepresentation'
    slice2Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.009999999776482582, 1.0, 0.5, 0.0]
    slice2Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.009999999776482582, 1.0, 0.5, 0.0]
    # turn off scalar coloring
    ColorBy(slice2Display, None)

    # change solid color
    slice2Display.AmbientColor = [1.0, 1.0, 1.0]
    slice2Display.DiffuseColor = [1.0, 1.0, 1.0]
    
    ss_setAndUpdate('y', sv)
    
    return slice2
    
def ss_yarrows(slice2, sv):
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
    glyph1Display.Representation = 'Surface'
    glyph1Display.ColorArrayName = [None, '']
    glyph1Display.OSPRayScaleArray = 'U'
    glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    glyph1Display.SelectOrientationVectors = 'None'
    glyph1Display.ScaleFactor = 0.03
    glyph1Display.SelectScaleArray = 'None'
    glyph1Display.GlyphType = 'Arrow'
    glyph1Display.GlyphTableIndexArray = 'None'
    glyph1Display.GaussianRadius = 3.018014947883785e-05
    glyph1Display.SetScaleArray = ['POINTS', 'U']
    glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
    glyph1Display.OpacityArray = ['POINTS', 'U']
    glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
    glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
    glyph1Display.PolarAxes = 'PolarAxesRepresentation'
    glyph1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.009999999776482582, 1.0, 0.5, 0.0]
    glyph1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.009999999776482582, 1.0, 0.5, 0.0]   
    ss_uColorBar(sv.renderView1, glyph1Display, 0.02)
    ss_setAndUpdate('y', sv)
    
def ss_sliceandyarrows(sv):
    ss_hideAll()
    slice2 =  ss_plainslice(sv)
    ss_yarrows(slice2, sv)
    
def ss_tube(sv):
    streamTracer1 = StreamTracer(Input=sv.casevtmseries, SeedType='High Resolution Line Source')
    streamTracer1.Vectors = ['POINTS', 'U']
    streamTracer1.MaximumStreamlineLength = 0.006
    streamTracer1.SeedType.Point1 = [-0.003014999907463789, -0.0021104998886585236, tubeh]
    streamTracer1.SeedType.Point2 = [0.003014999907463789, 0.0021104998886585236, tubeh]
    streamTracer1.SeedType.Resolution = 50

    # toggle 3D widget visibility (only when running from the GUI)
    Hide3DWidgets(proxy=streamTracer1.SeedType)

    # show data in view
    streamTracer1Display = Show(streamTracer1, sv.renderView1, 'GeometryRepresentation')


    # trace defaults for the display properties.
    streamTracer1Display.Representation = 'Surface'
    streamTracer1Display.ColorArrayName = [None, '']
    streamTracer1Display.OSPRayScaleArray = 'AngularVelocity'
    streamTracer1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    streamTracer1Display.SelectOrientationVectors = 'Normals'
    streamTracer1Display.ScaleFactor = 0.001738997269421816
    streamTracer1Display.SelectScaleArray = 'AngularVelocity'
    streamTracer1Display.GlyphType = 'Arrow'
    streamTracer1Display.GlyphTableIndexArray = 'AngularVelocity'
    streamTracer1Display.GaussianRadius = 8.69498634710908e-05
    streamTracer1Display.SetScaleArray = ['POINTS', 'AngularVelocity']
    streamTracer1Display.ScaleTransferFunction = 'PiecewiseFunction'
    streamTracer1Display.OpacityArray = ['POINTS', 'AngularVelocity']
    streamTracer1Display.OpacityTransferFunction = 'PiecewiseFunction'
    streamTracer1Display.DataAxesGrid = 'GridAxesRepresentation'
    streamTracer1Display.PolarAxes = 'PolarAxesRepresentation'

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

    # trace defaults for the display properties.
    tube1Display.Representation = 'Surface'
    tube1Display.ColorArrayName = [None, '']
    tube1Display.OSPRayScaleArray = 'AngularVelocity'
    tube1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    tube1Display.SelectOrientationVectors = 'Normals'
    tube1Display.ScaleFactor = 0.001094136107712984
    tube1Display.SelectScaleArray = 'AngularVelocity'
    tube1Display.GlyphType = 'Arrow'
    tube1Display.GlyphTableIndexArray = 'AngularVelocity'
    tube1Display.GaussianRadius = 5.4706805385649205e-05
    tube1Display.SetScaleArray = ['POINTS', 'AngularVelocity']
    tube1Display.ScaleTransferFunction = 'PiecewiseFunction'
    tube1Display.OpacityArray = ['POINTS', 'AngularVelocity']
    tube1Display.OpacityTransferFunction = 'PiecewiseFunction'
    tube1Display.DataAxesGrid = 'GridAxesRepresentation'
    tube1Display.PolarAxes = 'PolarAxesRepresentation'
    tube1Display.ScaleTransferFunction.Points = [-3.4285312251856137, 0.0, 0.5, 0.0, 3.581233990897969, 1.0, 0.5, 0.0]
    tube1Display.OpacityTransferFunction.Points = [-3.4285312251856137, 0.0, 0.5, 0.0, 3.581233990897969, 1.0, 0.5, 0.0]

    Hide(streamTracer1, sv.renderView1)
    sv.renderView1.Update()

    ss_uColorBar(sv.renderView1, tube1Display, 0.01)

    # update the view to ensure updated data information
    sv.renderView1.Update()
    
def ss_surfandtube(sv):
    ss_hideAll()
    ss_plainalphasurface(sv)
    ss_tube(sv)
    
    
################ scripting #########################
class ss_stateVars():
    def __init__(self, folder):
        ss_mksubdirs(folder)
        self.folder = folder        
        self.casevtmseries = ""
        self.renderView1 = ""
        self.animationScene1 = ""
        self.timeKeeper1 = ""
        self.timeStamp = ""
        self.times = ""
        self.initialized = False
        

def ss_runthrough(tlist, viewlistin, coloring, subfolder, sv, f, mode):
    # tlist is the list of times to go through
    # viewlist is the list of views to go through
    # coloring is the coloring label
    # subfolder is 'frames' or 'videos'
    # sv is a stateVars object
    # f is a function to create the rendering
    tlist2 = []
    fnlist = []
    viewlist = []

    for view in viewlistin: # for each view point
        # if mode==1:
        #     # if we're creating video frames, check if the video already exists, and quit if it does
        #     vidfn = vidFN(view, coloring, folder)
        #     if os.path.exists(vidfn):
        #         return 
        for t in tlist:
            # go through each time step
            fn = ss_fullFN(t, view, coloring, folder, subfolder)
                # find the full file name for this set of circumstances
            if (not os.path.exists(fn)) and t in sv.times:
                # if the file does not exist, add this set to the list
                tlist2.append(t)
                fnlist.append(fn)
                viewlist.append(view)

    
    
    # if we already have all the files we need, don't run this
    if len(tlist2)>0:
        if not sv.initialized:
            try:
                ss_initializeSV(sv)
            except:
                raise
        f(sv)
        Show(sv.timeStamp, sv.renderView1)     # add time stamp

        # iterate through times and take snapshots of the surfaces
        for [i,t] in enumerate(tlist2):
            if t-round(t)==0:
                print(f't = {t}, view = {viewlist[i]}, filename = {fnlist[i]}\n' )
                ss_setTime(t/10, sv)
                ss_setView(viewlist[i], sv)
                SaveScreenshot(fnlist[i], sv.renderView1)
    return sv
 
        


def ss_cleansession():
    Disconnect()
    Connect()
    
def ss_parseVTKSeries(folder):
    seriesfile = ss_series(folder)
    times = []
    if os.path.exists(seriesfile):
        with open(seriesfile, 'r') as f:
            for line in f:
                if 'name' in line:
                    times.append(int(round(10*float(re.split('time\" : | }', line)[1]))))
        return times
    else:
        return []

def ss_folderscript(folder, modes, volumes, meshes, vectors, tubes, tlist):
    print(folder)
    try:
        if not os.path.exists(folder):
            return
        sv = ss_stateVars(folder)
        sv.times = ss_parseVTKSeries(folder) # these times are ints
        for mode in modes:
            if mode==0:
                subfolder = 'frames'
                trange = tlist
            else:
                subfolder = 'videos'  
                trange = sv.times

            coloringlist = []
            flist = []
            vollist = []

            # ink-support interfaces
            if volumes:
                coloringlist.append('umag')
                flist.append(vol)
                vollist.append(volviewlist)

            # mesh from the y view
            if meshes:
                coloringlist.append('mesh')
                flist.append(mesh)
                vollist.append(['y'])

            # flow vectors from the y view
            if vectors:
                coloringlist.append('vecs')
                flist.append(sliceandyarrows)
                vollist.append(['y'])

            # streamlines with interface
            if tubes:
                coloringlist.append('stre')
                flist.append(surfandtube)
                vollist.append(volviewlist)

            for i in range(len(coloringlist)):
                try:
                    sv = ss_runthrough(trange, vollist[i], coloringlist[i], subfolder, sv, flist[i], mode)
                except:
                    return
        ss_cleansession()
        return
    except Exception as e:
        print(e)
        ss_cleansession()
        return


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                             SCRIPT
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for folder in folders:
    ss_folderscript(folder, modes, volumes, meshes, vectors, tubes, tlist)
print('---------------\n--Done exporting images--\n---------------\n')