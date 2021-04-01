import os
import numpy as np

# t0 = 5
# tf = 55
# dt = 5
t0 = 1
tf = 51
dt = 1
y = False
z = False
a = True
b = False
tubes = True
tubeh = 0.001

folder = 'C:/Users/lmf1/Documents/OpenFOAM/nozzlebath3d/nb16'

for folder in os.listdir("C:/Users/lmf1/Documents/OpenFOAM/nozzlebath3d/"):
    if os.path.isdir(folder) and folder.startswith(nb):

seriesfile = os.path.join(folder, 'case/VTK/case.vtm.series')

yims = []
zims = []
aims = []
bims = []

#### import the simple module from the paraview
from paraview.simple import *
LoadPalette(paletteName='WhiteBackground')
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

ss = GetSources()
for s in ss:
    Hide(ss[s])

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
renderView1.ViewSize = [1216, 1216]
# get layout
layout1 = GetLayout()
renderView1.CameraParallelProjection = 1
# get the material library
materialLibrary1 = GetMaterialLibrary()

# get active source.
#casevtmseries = GetActiveSource()
casevtmseries = XMLMultiBlockDataReader(FileName=seriesfile)

# # Properties modified on casevtmseries
casevtmseries.CellArrayStatus = ['U', 'alpha.ink', 'p_rgh']
casevtmseries.PointArrayStatus = ['U', 'alpha.ink', 'p_rgh']

# show data in view
casevtmseriesDisplay = Show(casevtmseries, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
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

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
casevtmseriesDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.009999999776482582, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
casevtmseriesDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.009999999776482582, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(casevtmseriesDisplay, ('FIELD', 'vtkBlockColors'))

# show color bar/color legend
casevtmseriesDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'vtkBlockColors'
vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')

# get opacity transfer function/opacity map for 'vtkBlockColors'
vtkBlockColorsPWF = GetOpacityTransferFunction('vtkBlockColors')


# create a new 'Clip'
clip1 = Clip(Input=casevtmseries)
clip1.ClipType = 'Plane'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['POINTS', 'alpha.ink']
clip1.Value = 0.5

# Properties modified on clip1
clip1.ClipType = 'Scalar'
clip1.Invert = 0

# show data in view
clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'
clip1Display.ColorArrayName = [None, '']
clip1Display.OSPRayScaleArray = 'U'
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.SelectOrientationVectors = 'None'
clip1Display.ScaleFactor = 0.00018132181139662862
clip1Display.SelectScaleArray = 'None'
clip1Display.GlyphType = 'Arrow'
clip1Display.GlyphTableIndexArray = 'None'
clip1Display.GaussianRadius = 9.066090569831431e-06
clip1Display.SetScaleArray = ['POINTS', 'U']
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityArray = ['POINTS', 'U']
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display.DataAxesGrid = 'GridAxesRepresentation'
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityUnitDistance = 5.84412949766512e-05
clip1Display.ExtractedBlockIndex = 1

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# hide data in view
Hide(casevtmseries, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(clip1Display, None)
HideScalarBarIfNotNeeded(vtkBlockColorsLUT, renderView1)

renderView1.OrientationAxesVisibility = 0





# create a new 'Annotate Time Filter'
annotateTimeFilter1 = AnnotateTimeFilter(Input=clip1)
# Properties modified on annotateTimeFilter1
annotateTimeFilter1.Format = '%2.1f s'
# show data in view
annotateTimeFilter1Display = Show(annotateTimeFilter1, renderView1, 'TextSourceRepresentation')
# Properties modified on annotateTimeFilter1Display
annotateTimeFilter1Display.FontSize = 12

# Properties modified on animationScene1
animationScene1.AnimationTime = 5.0

# Properties modified on timeKeeper1
timeKeeper1.Time = 5.0


# reset view to fit data
renderView1.ResetCamera()

if tubes:
    # create a new 'Stream Tracer'
    streamTracer1 = StreamTracer(Input=casevtmseries, SeedType='High Resolution Line Source')
    streamTracer1.Vectors = ['POINTS', 'U']
    streamTracer1.MaximumStreamlineLength = 0.004


    # init the 'High Resolution Line Source' selected for 'SeedType'
    streamTracer1.SeedType.Point1 = [-0.003014999907463789, -0.0021104998886585236, tubeh]
    streamTracer1.SeedType.Point2 = [0.003014999907463789, 0.0021104998886585236, tubeh]
    streamTracer1.SeedType.Resolution = 50

    # toggle 3D widget visibility (only when running from the GUI)
    Hide3DWidgets(proxy=streamTracer1.SeedType)

    # show data in view
    streamTracer1Display = Show(streamTracer1, renderView1, 'GeometryRepresentation')


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
    renderView1.Update()
    
    # create a new 'Tube'
    tube1 = Tube(Input=streamTracer1)
    tube1.Scalars = ['POINTS', 'AngularVelocity']
    tube1.Vectors = ['POINTS', 'Normals']
    tube1.Radius = 1e-05

    # show data in view
    tube1Display = Show(tube1, renderView1, 'GeometryRepresentation')

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

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    tube1Display.ScaleTransferFunction.Points = [-3.4285312251856137, 0.0, 0.5, 0.0, 3.581233990897969, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    tube1Display.OpacityTransferFunction.Points = [-3.4285312251856137, 0.0, 0.5, 0.0, 3.581233990897969, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(streamTracer1, renderView1)

    # update the view to ensure updated data information
    renderView1.Update()

    # rescale color and/or opacity maps used to include current data range
    tube1Display.RescaleTransferFunctionToDataRange(True, False)

    # set scalar coloring
    ColorBy(tube1Display, ('POINTS', 'U', 'Magnitude'))

    # show color bar/color legend
    tube1Display.SetScalarBarVisibility(renderView1, True)

    # get color transfer function/color map for 'U'
    uLUT = GetColorTransferFunction('U')

    # get opacity transfer function/opacity map for 'U'
    uPWF = GetOpacityTransferFunction('U')

    # Rescale transfer function
    uLUT.RescaleTransferFunction(0.0, 0.01)
    uLUTColorBar = GetScalarBar(uLUT, renderView1)
    uLUTColorBar.AutoOrient = 0
    uLUTColorBar.Orientation = 'Horizontal'
    uLUTColorBar.WindowLocation = 'LowerCenter'
    uLUTColorBar.ScalarBarLength = 0.7
    uLUTColorBar.UseCustomLabels = 1
    uLUTColorBar.CustomLabels = np.arange(0, 0.012, 0.002)
    uLUTColorBar.TitleFontSize = 12
    uLUTColorBar.LabelFontSize = 12
    uLUTColorBar.ComponentTitle = 'Mag (mm/s)'
    uLUTColorBar.RangeLabelFormat = '%-#1.3f'

    # Rescale transfer function
    uPWF.RescaleTransferFunction(0.0, 0.01)
    
    # update the view to ensure updated data information
    renderView1.Update()

# t = 5

# surfaces

for t in range(t0, tf, dt):
    tstr = str(t)
    if t<10:
        tstr = "00" + tstr
    elif t<100:
        tstr = "0" + tstr

    # Properties modified on animationScene1
    animationScene1.AnimationTime = t/10

    # Properties modified on timeKeeper1
    timeKeeper1.Time = t/10

    if tubes:
        s2 = "stream"
    else:
        s2 = "surface"
    
    if z:
        # z view
        renderView1.CameraFocalPoint = [0.001, 0,0]
        renderView1.CameraPosition = [0.001, 0, 10]
        renderView1.CameraViewUp = [0,1,0]
        SaveScreenshot(folder+'/stills/z_'+s2+'_'+tstr+'.png', renderView1)

    if y:
        # y view
        renderView1.CameraPosition = [0.001, -10, 0]
        renderView1.CameraViewUp = [0,0,1]
        SaveScreenshot(folder+'/stills/y_'+s2+'_'+tstr+'.png', renderView1)

    if a:
        # view a
        renderView1.CameraFocalPoint = [0.0006, -0.0002, 0.0005]
        renderView1.CameraPosition = [2, 1, 2]
        renderView1.CameraViewUp = [0,0,1]
        SaveScreenshot(folder+'/stills/a_'+s2+'_'+tstr+'.png', renderView1)

    if b:
        # view b
        renderView1.CameraFocalPoint = [0.001, 0, 0.001]
        renderView1.CameraPosition = [-1,-1,1]
        renderView1.CameraViewUp = [0,0,1]
        SaveScreenshot(folder+'/stills/b_'+s2+'_'+tstr+'.png', renderView1)
