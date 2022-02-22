#!/usr/bin/env pvpython
'''Collecting interface points into csvs from vtk files'''

import os
import numpy as np
import csv
import re
from paraview.simple import * # import the simple module from the paraview
import time
from datetime import datetime
from paraview_general import *
from typing import List, Dict, Tuple, Union, Any, TextIO
import logging
logger = logging.getLogger(__name__)

__author__ = "Leanne Friedrich"
__copyright__ = "This data is publicly available according to the NIST statements of copyright, fair use and licensing; see https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software"
__credits__ = ["Leanne Friedrich"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Leanne Friedrich"
__email__ = "Leanne.Friedrich@nist.gov"
__status__ = "Production"

############################
## CSV ##############
##################

def initializeAll(folder, mode:str):
    '''initialize all of paraview'''
    logging.info(f'CSVs: Initializing paraview for {folder}')
    sv = stateVars(folder)        # make subdirectories
    sv =  initializeP(sv)        # initialize Paraview
    if mode=='interface': # RG
        sv =  initSeriesInter(sv) # import the vtk series files
    elif mode=='nozzle':
        sv =  initSeriesNoz(sv, mode='volume') # import the vtk series files
    elif mode=='nozzleSlice':
        sv =  initSeriesNoz(sv, mode='slice') # import the vtk series files
    return sv


def initSeriesInter(sv:stateVars):
    '''Import the case.vtm.series file and just get the interface'''
    caseVTMSeries = initSeries0(sv)

    clip1 = Clip(Input=caseVTMSeries)
    clip1.Scalars = ['POINTS', 'alpha.ink']
    clip1.Value = 0.4
    clip1.ClipType = 'Scalar'
    clip1.Invert = 0
    clip2 = Clip(Input=clip1)
    clip2.Scalars = ['POINTS', 'alpha.ink']
    clip2.Value = 0.6
    clip2.ClipType = 'Scalar'
    clip2.Invert = 1
    slice1 = Slice(Input=clip2)
    slice1.SliceType = 'Plane'
    slice1.SliceOffsetValues = [0.0]
    slice1.SliceType.Origin = [0.00334, 0, 0]
    slice1.SliceType.Normal = [1, 0, 0] # look down the y axis
    Hide(caseVTMSeries, sv.renderView1)
    Hide(clip1, sv.renderView1)
    Hide(clip2, sv.renderView1)
    sv.caseVTMSeries = caseVTMSeries
    sv.slice = slice1
    slice1Display = Show(slice1, sv.renderView1, 'GeometryRepresentation')
    return sv

def initSeriesNoz(sv:stateVars, mode:str='slice'): # RG
    '''Import the case.vtm.series file and just get the ink in the nozzle. edges mode uses planar boundaries around the '''
    caseVTMSeries = initSeries0(sv)
    renderView1 = GetActiveViewOrCreate('RenderView')

    le = fp.legendUnique(sv.folder)
    ink_rho = float(le['ink_rho'])
    sup_rho = float(le['sup_rho'])
    di = float(le['nozzle_inner_width'])/1000 # inner diameter
    z0 = float(le['nozzle_bottom_coord'])/1000 # tip of nozzle
    zf = (float(le['bath_top_coord'])-float(le['nozzle_bottom_coord']))*0.97/1000 + z0
    x0 = float(le['nozzle_center_x_coord'])/1000
    y0 = float(le['nozzle_center_y_coord'])/1000
    xr = x0+di/2
    xl = x0-di/2
    na = np.deg2rad(float(le['nozzle_angle'])) # nozzle angle in radians

    # clips ink
    clip1 = Clip(Input=caseVTMSeries)
    clip1.Scalars = ['POINTS', 'alpha.ink'] 
    clip1.Value = 0.01
    clip1.ClipType = 'Scalar'
    clip1.Invert = 0

    if mode=='slice':

        # # clip bottom of nozzle
        # clip2 = Clip(Input=clip1)
        # clip2.ClipType = 'Plane'
        # clip2.HyperTreeGridClipper = 'Plane'
        # clip2.ClipType.Origin = [0, 0, z0]
        # clip2.HyperTreeGridClipper.Origin = [0, 0, z0]
        # clip2.Invert = 0
        # clip2.ClipType.Normal = [0.0, 0.0, 1.0]
        # Hide3DWidgets(proxy=clip2.ClipType)

        # clip right edge
        clip3 = Clip(Input=clip1)
        clip3.ClipType = 'Plane'
        clip3.HyperTreeGridClipper = 'Plane'
        clip3.ClipType.Origin = [xr, 0, z0]
        clip3.HyperTreeGridClipper.Origin = [xr, 0, z0]
        clip3.Invert = 0
        clip3.ClipType.Normal = [-np.cos(na), 0.0, np.sin(na)]
        Hide3DWidgets(proxy=clip3.ClipType)

        # clip left edge
        clip4 = Clip(Input=clip3)
        clip4.ClipType = 'Plane'
        clip4.HyperTreeGridClipper = 'Plane'
        clip4.ClipType.Origin = [xl, 0, z0]
        clip4.HyperTreeGridClipper.Origin = [xl, 0, z0]
        clip4.Invert = 0
        clip4.ClipType.Normal = [np.cos(na), 0.0, np.sin(na)]
        Hide3DWidgets(proxy=clip4.ClipType)

        


        # below here is specific to shear stress and should be changed/removed if not looking at shear stress
        calculator1 = Calculator(Input=clip4)

    else:
        # take whole volume
        clip2 = Clip(Input=clip1)
        clip2.Scalars = ['POINTS', 'p'] # clips ink more than just in nozzle, but the excess is ignored in processing
        clip2.Value = 70.0 # may need to change value if max value of pressure in nozzle is low or in bath is high
        clip2.ClipType = 'Scalar'
        clip2.Invert = 0
        # below here is specific to shear stress and should be changed/removed if not looking at shear stress
        calculator1 = Calculator(Input=clip2)

    calculator1.ResultArrayName = 'magU'
    calculator1.Function = 'mag(U)'
    computeDerivatives1 = ComputeDerivatives(Input=calculator1)
    computeDerivatives1.Scalars = ['POINTS', 'magU']
    computeDerivatives1.Vectors = ['POINTS', 'U']
    cellDatatoPointData1 = CellDatatoPointData(Input=computeDerivatives1)
    #     cellDatatoPointData1.CellDataArraytoprocess = ['ScalarGradient', 'U', 'VectorGradient', 'alpha.ink', 'cellID', 'nu1', 'nu2', 'p', 'p_rgh', 'rAU']
    cellDatatoPointData1.CellDataArraytoprocess = ['ScalarGradient', 'U', 'alpha.ink', 'nu1', 'nu2', 'p', 'p_rgh', 'rAU']
    calculator2 = Calculator(Input=cellDatatoPointData1)
    calculator2.ResultArrayName = 'ShearStress'
    calculator2.Function = stressFunc(le)
    calculator3 = Calculator(Input=calculator2)
    calculator3.ResultArrayName = 'ShearStressMag'
    calculator3.Function = 'mag(ShearStress)'


    if mode=='slice':
        # slice through y center of nozzle
        slice1 = Slice(registrationName='Slice1', Input=calculator3)
        slice1.SliceType = 'Plane'
        slice1.HyperTreeGridSlicer = 'Plane'
        slice1.SliceOffsetValues = [0.0]
        slice1.SliceType.Origin = [0,0,0]
        slice1.HyperTreeGridSlicer.Origin = [0,0,0]
        slice1.SliceType.Normal = [0.0, 1.0, 0.0]
        Hide3DWidgets(proxy=slice1.SliceType)

        # below here is needed
        slice2 = Slice(Input=slice1)

    else:
        slice2 = Slice(Input=calculator3)
    slice2.SliceType = 'Plane'
    slice2.SliceOffsetValues = [0.0]
    slice2.SliceType.Origin = [0, 0, zf] # 3% the length of the nozzle below the nozzle top to avoid innacurate shear stress at the ink-atmosphere interface
    slice2.SliceType.Normal = [0, 0, 1] # look down the z axis

    Hide(caseVTMSeries, sv.renderView1)
    # Hide(clip1, sv.renderView1)
    # Hide(clip2, sv.renderView1)
    # Hide(clip3, sv.renderView1)
    # Hide(clip4, sv.renderView1)
    # Hide(slice1, sv.renderView1)
    Hide(calculator1, sv.renderView1)
    Hide(computeDerivatives1, sv.renderView1)
    Hide(cellDatatoPointData1, sv.renderView1)
    Hide(calculator2, sv.renderView1)
    Hide(calculator3, sv.renderView1)
    sv.caseVTMSeries = caseVTMSeries
    sv.slice = slice2
    slice1Display = Show(slice2, sv.renderView1, 'GeometryRepresentation')
    return sv



################
# creating csvs

def addToFile(x:float, tempFile:str, w:csv.writer, sv:stateVars, skipHeader:bool, mode:str, xdisp:float=0) -> bool:
    '''Add more points to the file. x is the x position that we take the slice at. tempFile is the full path of the temp file that we store points in temporarily before we put them into the actual csv file. We already opened a csvwriter w in csvFolder. sv is a stateVars object that holds the actual Paraview objects that we need to analyze. skipHeader is true if we don't want to include the header and false if we want to write the header to w. mode is either interface or nozzle. xdisp is a displacement that we can add to the saved x values, if we want to make x relative to a certain point. It is best to leave this at 0.'''

    ret = skipHeader
    if mode=='interface':
        sv.slice.SliceType.Origin = [x, 0, 0]
        SaveData(tempFile, FieldAssociation="Point Data", ChooseArraysToWrite=1, AddTime=1) # saves the slice to the temp file
    else:
        sv.slice.SliceType.Origin = [0, 0, x] # moving down the z axis
        SaveData(tempFile, FieldAssociation="Point Data", ChooseArraysToWrite=1, AddTime=1, PointDataArrays=['alpha.ink', 'ShearStressMag', 'magU', 'nu1', 'nu2', 'p']) # RG
    with open(tempFile, mode='r') as f2: # scrape the data out of the temp file and move it to the permanent file
        ftemp = csv.reader(f2)
        header=next(ftemp)
        if 'Points:0' not in header:
            return False
        xpos = header.index('Points:0')
        if not skipHeader:
            if len(header)>4:
                headerdict = {'alpha.ink':'alpha', 'Points:0':'x', 'Points:1':'y', 'Points:2':'z', 'U:0':'vx', 'U:1':'vy', 'U:2':'vz', 'Time':'time', 'nu1':'nu_ink', 'nu2':'nu_sup'}
                if not xdisp==0:
                    if xdisp>0:
                        headerdict['Points:0'] = f'x+{xdisp}'
                    else:
                        headerdict['Points:0'] = f'x-{abs(xdisp)}'
                header2 = [headerdict.get(h,h) for h in header]
                if mode=='interface':
                    unitdict = {'time':'s', 'x':'m', 'y':'m', 'z':'m', 'vx':'m/s', 'vy':'m/s', 'vz':'m/s', 'alpha':'', 'nu1':'m^2/s', 'nu2':'m^2/s', 'nu_ink':'m^2/s', 'nu_sup':'m^2/s', 'p':'kg/(m*s^2)', 'p_rgh':'kg/(m*s^2)', 'rAU':'m^3*s/kg'}
                else:
                    unitdict = {'time':'s', 'x':'m', 'y':'m', 'z':'m', 'alpha':'', 'nu_ink':'m^2/s', 'nu_sup':'m^2/s', 'p':'kg/(m*s^2)', 'magU':'m/s', 'ShearStressMag':'Pa'} # RG
                unitline = [unitdict.get(s, s) for s in header2]
                w.writerow(header2)
                w.writerow(unitline)
                ret = True
            else:
                ret = False
        for row in ftemp:
            r = row
            if len(r)>1:
                r[xpos] = str(float(r[xpos])+xdisp)
                w.writerow(r)
    return ret

####
# script

def ipFN(f1:str, mode:str, time:float) -> str:
    '''get the filename for an interfacepoints file. f1 is the folder name, mode is interface or nozzle, time is actual time in s'''
    tstr = int(time*10)
    return os.path.join(f1, f"{mode}Points_t_{tstr}.csv")


def csvFolderMode(folder:str, mode:str, forceOverwrite:bool, times0:List[float]=[2.5]):
    '''goes through all the timesteps and generates either interfacePoints or nozzlePoints'''
    try:
        if not os.path.exists(folder):
            return

        le = fp.legendUnique(folder)
        if mode=='interface':
            # xmin = -0.004824 # left edge of bath
            # xmax = -xmin # right edge of bath
            xmin = float(le['bath_left_coord'])/1000
            xmax = float(le['bath_right_coord'])/1000
            dx = (xmax-xmin)/50
            # dx = 0.0002
            # times = fp.times(folder)
        elif mode=='nozzleSlice' or mode=='nozzle':
            
            xmin = float(le['nozzle_bottom_coord'])/1000 # tip of nozzle
            # xmin = 0.0003015 # tip of nozzle RG
            xmax = (float(le['bath_top_coord'])-float(le['nozzle_bottom_coord']))*0.97/1000 + xmin
            # xmax = 0.00205623 # 3% below top of nozzle RG
            dx = (xmax-xmin)/35
            times = times0 # very storage intensive to do all 25 timesteps
        else:
            raise Exception("Valid modes are 'interface' and 'nozzleSlice'")

        f1 = os.path.join(folder, f'{mode}Points')

        # only collect files that don't already exist, or overwrite if we're forced to
        times = []
        for time in times0:
            ipfile =  ipFN(f1, mode, time) # this is the file that all points for this time will be saved in RG
            if not os.path.exists(ipfile) or forceOverwrite:
                times.append(time)
        if len(times)==0:
            return

        xchunk = 0.0015
        initialized = False
        tempFile = os.path.join(f1, "temp.csv")

        for time in times:
            ipfile =  ipFN(f1, mode, time)
            if not initialized: # if paraview hasn't already been initialized, initialize it
                sv = initializeAll(folder, mode)
                times = [i for i in times if i in sv.times] # only select time steps that are in the file
                logging.info(f'-----{times}, {sv.times}')
                # sv.times = times
                initialized = True
            if time in times:
                # if this time has now been eliminated from the list, skip it
                setTime(time, sv) 
                # error checking: sometimes we lose the ability to set the time, for inexplicable reasons.
                logging.info(f'--t {time}, file: {ipfile}' )
                if not sv.timeKeeper1.Time==time:
                    raise Exception(f'{time}, {sv.timeKeeper1.Time}, Timekeeper broken. Start again.')
                with open(ipfile, mode='w', newline='') as f:
                    w = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                    skipHeader = False # addToFile generates a list of points at a given x position, and each of those tables has a header. We only want to include the header once, so skipHeader tells us if we already have a header in w
                    for x in (np.arange(xmin, xmax, dx)):
                        skipHeader = addToFile(x, tempFile, w, sv, skipHeader, mode, xdisp=0)
                os.remove(tempFile)
        ResetSession()
        cleanSession()
        try:
            del sv
        except:
            return
    except Exception as e:
        logging.error(folder+': '+str(e))
        ResetSession()
        cleanSession()
        return

def csvFolder(folder:str, modes:List[str], forceOverwrite:bool, times0:List[float]=[2.5]):
    '''Goes through all of the timesteps in the simulation and generates an interfacePoints/nozzlePoints table, containing all the points at the interface.'''
    for mode in modes: # run through for each mode RG
        csvFolderMode(folder, mode, forceOverwrite, times0)