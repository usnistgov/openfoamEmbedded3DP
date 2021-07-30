#!/usr/bin/env python
'''Functions for generating legends for OpenFOAM simulations of embedded 3D printing of single filaments. Written for OpenFOAM v1912 and OpenFOAM 8. Scrapes input files for input variables.'''

# external packages
import os
import re
import csv
import shutil 
import errno
from typing import List, Dict, Tuple, Union, Any, TextIO
from datetime import datetime
import time
import logging, platform, socket, sys
from backwardsRead import fileReadBackwards

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from folderparser import *

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
for s in ['matplotlib']:
    logging.getLogger(s).setLevel(logging.WARNING)

# info
__author__ = "Leanne Friedrich"
__copyright__ = "This data is publicly available according to the NIST statements of copyright, fair use and licensing; see https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software"
__credits__ = ["Leanne Friedrich"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Leanne Friedrich"
__email__ = "Leanne.Friedrich@nist.gov"
__status__ = "Production"

#-------------------------------------------------------------------------------------------------  


class scrape:
    '''the scrape class is used to store variables scraped from the log files
        scrape contains a function scrape.table(), which converts all the data in the class into a 2 column table'''
        
    def __init__(self, folder:str):
        '''folder is a folder name (string) which should point either to the case folder or one level above it'''
        
        self.fold = folder # full path
        self.meshfold = meshFolder(folder)
        self.casefold = caseFolder(folder)
        self.folder = ['folder', os.path.basename(folder)] # just the folder name
        self.compareto = ['compare to', '']

        self.initRates()
        self.initGeo()
                        
        # these lists store values ripped from the dictionaries
        self.initMesh()
        self.initDMD()                
        self.initTransport()
        self.initControl()                
        self.initFV()
    
    #-------------------------
    def initRates(self):
        '''initial (empty) values for simulation rates'''
        
        # these values are ripped from log files
        self.shmtimes = ['snappyHexMesh time (s)', ''] 
            # how long it takes snappyHexMesh to run
        self.shmtimem = ['snappyHexMesh time (min)', '']
        self.iftimes = ['interFoam time (s)', ''] 
            # how long it takes interfoam to run
        self.iftimehr = ['interFoam time (hr)', '']
        self.simTime = ['simulation time (s)', ''] 
            # how many seconds within the simulation we ran
        self.simrate = ['simulation rate (hr/s)', ''] 
            # how many real hours it takes to run each second within the simulation
        
        
    #-------------------------
    def initGeo(self):
        '''these values measure the geometry of the system and are stored in geometry.csv, which is created by ncreate3d.py/noz3dscript.ipynb'''
        self.GEOniw = ['nozzle inner width (mm)', '0.603']
        self.GEOnt = ['nozzle thickness (mm)', '0.152']
        self.GEObw = ['bath width (mm)', '']
        self.GEObd = ['bath depth (mm)', '']
        self.GEOnl = ['nozzle length (mm)', '']
        self.GEOblc = ['bath left coord (mm)', ''] # x coord
        self.GEObrc = ['bath right coord (mm)', ''] # x coord
        self.GEObfc = ['bath front coord (mm)', ''] # y
        self.GEObbackc = ['bath back coord (mm)', ''] # y
        self.GEObbotc = ['bath bottom coord (mm)', ''] # z
        self.GEObtc = ['bath top coord (mm)', ''] # z
        self.GEOnbc = ['nozzle bottom coord (mm)', ''] # z
        self.GEOncxc = ['nozzle center x coord (mm)', '']
        self.GEOncyc = ['nozzle center y coord (mm)', '']
        self.GEOna = ['nozzle angle (degrees)', ''] # RG
        self.GEObathv = ['bath velocity (m/s)', '']
        self.GEOinkv = ['ink velocity (m/s)', '']
    
    #-------------------------
    def initMesh(self):
        '''initial values for mesh variables'''
        
        self.SHMlist = [['castellatedMesh', ''],\
                        ['snap', ''],\
                        ['addLayers', '']]
        self.CMClist = [['maxLocalCells', ''], \
                        ['maxGlobalCells', ''], \
                        ['minRefinementCells', ''], \
                        ['nCellsBetweenLevels', ''],\
                        ['resolveFeatureAngle', ''],\
                        ['locationInMesh', ''], \
                        ['allowFreeStandingZoneFaces', '']]
        self.CMCfixedWallsLevel = ['fixedWalls level', '']
        self.SClist = [['nSmoothPatch', ''], \
                       ['tolerance', ''], \
                       ['nSolveIter', ''], \
                       ['nRelaxIter', ''],\
                       ['nFeatureSnapIter', ''], \
                       ['implicitFeatureSnap', ''],\
                       ['explicitFeatureSnap', ''], \
                       ['multiRegionFeatureSnap', '']]
        self.ALClist = [['relativeSizes', ''], \
                        ['expansionRatio', ''], \
                        ['finalLayerThickness', ''],\
                        ['minThickness', ''], \
                        ['nGrow', ''], \
                        ['featureAngle', ''], \
                        ['nRelaxIter', ''], \
                        ['nSmoothSurfaceNormals', ''], \
                        ['nSmoothNormals', ''], \
                        ['nSmoothThickness', ''],\
                        ['maxFaceThicknessRatio', ''],\
                        ['maxThicknessToMedialRatio', ''],\
                        ['minMedialAxisAngle', ''],\
                        ['nBufferCellsNoExtrude', ''], \
                        ['nLayerIter', ''], \
                        ['nRelaxedIter', '']]
        self.ALCfixedWallsLayers = ['fixedWalls nSurfaceLayers', '']
        self.MQClist = [['nSmoothScale', ''], \
                        ['errorReduction', ''], \
                        ['maxNonOrtho', '']]
        self.SHMmergeTolerance = ['mergeTolerance', '']
                        
        self.blocksdims = ['blocks dims', '']
    
    #-------------------------
    def initDMD(self):
        '''initial values for dynamic mesh dictionary'''
        self.DMDlist = [['refineInterval', ''], \
                        ['lowerRefineLevel', ''], \
                        ['upperRefineLevel', ''], \
                        ['unrefineLevel', ''],\
                        ['nBufferLayers', ''], \
                        ['maxRefinement', ''], \
                        ['maxCells', '']]
    
    #-------------------------
    def initTransport(self):
        '''initial values for transport model'''
        self.TPinklist = [['ink', ''],\
                        ['transportModel', ''],\
                          ['nu', ''], \
                          ['nu0', ''], \
                          ['tau0', ''], \
                          ['k', ''],\
                          ['n', ''], \
                          ['rho', '']]
        self.TPsuplist = [['sup', ''],\
                        ['transportModel', ''],\
                          ['nu', ''], \
                          ['nu0', ''], \
                          ['tau0', ''], \
                          ['k', ''], \
                          ['n', ''], \
                          ['rho', '']]  
        self.inkLabel = ''
        self.supLabel = ''
        self.TPsigma = ['sigma', '']
    
    #-------------------------
    def initControl(self):
        '''initial values for controlDict'''
        self.controlDictList =  [['application', ''],\
                                 ['endTime', ''], \
                                 ['deltaT', ''], \
                                 ['adjustTimeStep', ''], \
                                 ['maxCo', ''],\
                                 ['maxAlphaCo', ''],\
                                 ['maxDeltaT', ''], \
                                 ['writeControl', ''], \
                                 ['writeInterval', ''], \
                                 ['purgeWrite', ''],\
                                 ['runTimeModifiable', '']]
    
    #-------------------------
    def initFV(self):
        '''initial values for fvSolution, fvSchemes'''
        self.fvSailist = [['nAlphaSubCycles', ''],\
                          ['cAlpha', ''], \
                          ['nAlphaCorr', ''], \
                          ['MULESCorr', ''], \
                          ['nLimiterIter', ''], \
                          ['solver', ''], \
                          ['smoother', ''],\
                          ['tolerance', ''],\
                          ['relTol', ''], \
                          ['isofaceTol', ''], \
                          ['surfCellTol', ''], \
                          ['nAlphaBounds', ''], \
                          ['snapTol', ''],\
                          [ 'clip', '']]                
        self.fvSpcorrlist = [['solver', ''], \
                             ['preconditioner', ''], \
                             ['tolerance', ''], \
                             ['relTol', '']]
        self.fvSprghlist = [['solver', ''], \
                            ['preconditioner', ''],\
                            ['tolerance', ''],\
                            ['relTol', '']]
        self.fvSprghfinallist = [['relTol', '']]
        self.fvSUlist = [['solver', ''],\
                         ['smoother', ''], \
                         ['tolerance', ''], \
                         ['relTol', '']]
        self.fvSPIMPLElist = [['momentumPredictor', ''], \
                              ['nOuterCorrectors', ''], \
                              ['nCorrectors', ''], \
                              ['nNonOrthogonalCorrectors', '']]
        
        
    #-------------------------

    def table(self) -> List[List[str]]:
        '''this function takes all the values we ripped and stored as variables in the scrape object and turns them into a table with 2 columns: the variable name and the variable value'''
        col = [] # start with an empty table, then add to it sequentially
        for i in [self.folder, self.compareto, self.shmtimes, \
                  self.shmtimem, self.iftimes, self.iftimehr, \
                  self.simTime, self.simrate]:
            col.append(i) # RG
        ca(col, ['', '', 'mesh','GEOMETRY'],\
           [self.GEOniw, self.GEOnt, self.GEObw, self.GEObd, \
            self.GEOnl, self.GEOblc, self.GEObrc, self.GEObfc,\
            self.GEObbackc, self.GEObbotc, self.GEObtc, self.GEOnbc,\
            self.GEOncxc, self.GEOncyc, self.GEOna, self.GEObathv, self.GEOinkv])
        ca(col, ['', 'SYSTEM'], [])    
        ca(col, ['snappyHexMeshDict'], self.SHMlist)
        ca(col, ['castellatedMeshControls'], self.CMClist)
        col.append(self.CMCfixedWallsLevel)
        ca(col, ['snapControls'], self.SClist)
        ca(col, ['addLayersControls'], self.ALClist)
        col.append(self.ALCfixedWallsLayers)
        ca(col, ['meshQualityControls'], self.MQClist)
        col.append(self.SHMmergeTolerance)
        ca(col, ['', 'blockMeshDict'], [self.blocksdims])
        ca(col, ['', 'case', 'CONSTANT', 'dynamicMeshDict'], self.DMDlist)
        ca(col, ['', 'transportProperties'], self.TPinklist)
        ca(col, [], self.TPsuplist)
        col.append(self.TPsigma)
        ca(col, ['', 'SYSTEM', 'controlDict'], self.controlDictList)
        ca(col, ['', 'fvSolution', 'alpha.ink'], self.fvSailist)
        ca(col, ['', 'pcorr'], self.fvSpcorrlist)
        ca(col, ['', 'p_rgh'], self.fvSprghlist)
        ca(col, ['', 'p_rghFinal'], self.fvSprghfinallist)
        ca(col, ['', 'U'], self.fvSUlist)
        ca(col, ['', 'PIMPLE'], self.fvSPIMPLElist)
        return col




def ca(col:List[List[str]], hlist:List[str], li:List[List[str]]) -> None:
    '''ca is used by scrape.table() to compile all of the data in the object into a table
    col is a table with 2 columns. 
    hlist is a list of headers that describe this chunk of data. they are added as rows with an empty value 
    li is a list of [variable name, value] to add to col'''
    for i in hlist:
        col.append([i, ''])
    for i in li:
        col.append(i)
        
 
 ###################### SCRAPING TOOLS ##################

        

def placeInList(l:List[List[str]], s:str, v:str) -> None: 
    '''placeInList puts a value v into a list with 2 columns, where s is the name of the variable
        this function will only place the value into the list if the variable name s is in the list and there isn't already a value there
        l is a list
        s is a string
        v is a value'''
    try:
        i = l.index([s, '']) # i is the position in the list
    except:
        return
    else:
        l[i][1] = v # this puts the value in the list
    return


def cancelUnits(s:str) -> str:
    '''cancelUnits removes the units list (e.g. [ 0 2 -1 0 0 0 0]) from a value
    s is a string'''
    strs = re.split(' ', s)
    return strs[-1]


def listLevel(startString:str, endString:str, line:str, f:TextIO, l:List[List[str]]) -> str:
    '''listLevel is a tool that looks for sections of interest within files and scrapes values out of them
    startString is a string that tells us that we've reached the section of interest. It should be at the beginning of the line
    endString tells us that we've reached the end of the section of interest. It should be at the beginning of the line.
    line is the starting line
    f is a file stream created by open()
    l is a list of variable names and values that we're going to scrape values into
    returns the line we just read'''
    while not line.startswith(startString):
        line = f.readline()
    while not line.startswith(endString):
        strs = re.split(';|\t', line) # split the line at tabs and semicolons
        ii = 0
        si = 0
        while ii<len(strs) and len(strs[ii])==0:
            ii+=1 # find the first non-empty string in the list
        if ii+1<=len(strs)-1:
            s0 = strs[ii] # if there are enough entries in the split line to contain a name and value, place this entry into l
            s1 = cancelUnits(strs[ii+1])
            placeInList(l, s0, s1)
        line = f.readline()
    return line  


def readLevel0(s:str, line:str, f:TextIO, obj:List[List[str]]) -> str:
    '''readLevel0 is a simpler version of listLevel, where instead of placing many values in a list,
    we're looking for a single value. This function targets lines in files that have no tabs at the beginning, just "name\tvalue"
    s is a trigger string that tells us we've found the value we're looking for
    line is a starting line
    f is a file stream created by open()
    obj is a [1x2] list into which we'll store our value
    returns the line we just read'''
    while not line.startswith(s):
        line = f.readline()
    strs = re.split(';|\t', line) # split the line at ; and tabs
    obj[1] = strs[1] # the value will always be the second item
    return line


###################### SCRAPING FILES ##################


def scrapeSHMLog(s:scrape) -> None:
    '''scrape the snappyHexMesh log file
    s is a scrape object'''
    shmlog = os.path.join(s.meshfold, 'log_snappyHexMesh')
    if os.path.exists(shmlog):
        # it's useful to go backwards here because we're just looking for the time reported at the end of the file
        for line in fileReadBackwards(shmlog):
            if line.startswith('Finished meshing'):
                strs = re.split('Finished meshing in = | s', line)
                shmtime = float(strs[1])
                s.shmtimes[1] = "%.2f" % shmtime
                s.shmtimem[1] = "%.2f" % (shmtime/60)
                return 
    

              
def scrapeIFLog(s:scrape) -> None:
    '''scrape the interFoam log file
        Because interFoam can take hours to days, sometimes runs get split into pieces. 
        Each interFoam run adds onto the existing log file. 
        s is a scrape object  '''
    ifLog = os.path.join(s.casefold, 'log_interFoam')
    if not os.path.exists(ifLog):
        ifLog = os.path.join(s.fold, 'log_interFoam')
    if not os.path.exists(ifLog):
        return
    iftime = 0 # this variable adds up all the times for separate interFoam runs
    waitfortop = False # this variable tells us what to look for, so we only collect the last reported time from each run
    simtime = 0
    # it's useful to go backwards here because we're just looking for the time reported at the end of the run
    for line in fileReadBackwards(ifLog):
        if simtime==0 and line.startswith('Time = '):
            strs = re.split('Time = ', line)
            simtime = float(strs[1])
        if (waitfortop and line.startswith('fileModificationChecking')): 
            # when we hit fileModificationChecking, it's the end of the run, so now we should look for the next reported time
            waitfortop = False
        if (not waitfortop and line.startswith('ExecutionTime')):
            # read the last reported ExecutionTime
            strs = re.split('ExecutionTime = | s', line)
            iftime+=float(strs[1])
            waitfortop = True # now that we've read the time, look for the next end of run
    # store the time in the scrape object in seconds and hr
    s.iftimes[1] = "%.2f" % iftime 
    s.iftimehr[1] = "%.2f" % (iftime/60/60)
    try:
        if simtime==0:
            simtime = "%.2f" % float(s.simTime[1])
        if simtime>0:
            # if we've already measured a simulation time, calculate the simulation speed
            s.simrate[1] = "%.2f" % (iftime/60/60/simtime)
    except:
        pass
    return

def scrapeRunTime(s:scrape) -> None:
    '''Get the simulation time from the folder'''
    ti = times(s.fold)
    if len(ti)>0:
        s.simTime[1] = str(max(ti))
    else:
        s.simTime[1] = '0'
    return


def scrapeLogs(s:scrape) -> None:
    '''scrape all of the times (run time, simulation time, etc.)
    s is a scrape object'''
    scrapeRunTime(s)
    scrapeSHMLog(s)
    scrapeIFLog(s)
    

def scrapeBlockMeshDict(s:scrape) -> None:
    '''scrape blockMeshDict
    s is a scrape object'''
    bm = os.path.join(s.meshfold, 'system', 'blockMeshDict')
    if os.path.exists(bm):
        with open (bm, "r") as f:
            line = f.readline()
            while not line.startswith('vertices'):
                line = f.readline()
            # now we have reached the list of vertices
            # because we establish a basic mesh with blockmeshDict and refine with snappyHexMesh, this list only contains 8 vertices
            for i in range(2):
                line = f.readline() # read coords from the first list of points
            strs = re.split('\(|\)| ', line)
            s.GEOblc[1] = strs[1] # bath left coord
            s.GEObbackc[1] = strs[2] # bath back coord
            s.GEObbotc[1] = strs[3] # bath bottom coord
            for i in range(7):
                line = f.readline() # read coords from the last list of points
            strs = re.split('\(|\)| ', line)
            s.GEObrc[1] = strs[1] # bath right coord
            s.GEObfc[1] = strs[2] # bath front coord
            s.GEObtc[1] = strs[3] # bath top coord
            s.GEObw[1] = str(float(s.GEObrc[1]) - float(s.GEOblc[1])) # bath width
            s.GEObd[1] = str(float(s.GEObfc[1]) - float(s.GEObbackc[1])) # bath depth
            while not line.startswith('blocks'):
                line = f.readline()    
            # now we have reached the list of blocks
            # again, there is only one block because we're using snappyHexMesh
            for i in range(2):
                line = f.readline() # read coords from the first list of blocks
            strs = re.split('\) \(|\) s', line)
            s.blocksdims[1] = strs[1] # number of cells in the blocks: this should look like (# # #)
            return


def scrapeSetFieldsDict(s:scrape) -> None:
    '''scrape setFieldsDict
    s is a scrape object'''
    bm = os.path.join(s.casefold, 'system', 'setFieldsDict')
    if os.path.exists(bm):
        with open (bm, "r") as f:
            line = f.readline()
            while not line.startswith('\t\tp1'):
                line = f.readline()
            # now we have reached the bottom point of the nozzle
            strs = re.split('\( | \)| ', line)
#             strs = re.split('\(|\)| ', line) # RG
            s.GEOncxc[1] = str(1000*float(strs[2])) # nozzle center x
            s.GEOncyc[1] = str(1000*float(strs[3])) # nozzle center y
            s.GEOnbc[1] = (1000*float(strs[4])) # nozzle bottom 
            try:
                btc = float(s.GEObtc[1]) # bath top coord
            except:
                s.GEOnl[1] = ""
            else:
                s.GEOnl[1] = str(btc - s.GEOnbc[1]) # nozzle length
            s.GEOnbc[1] = str(s.GEOnbc[1])
            return
        

def scrapeU(s:scrape) -> None:
    '''scrape 0/U
    s is a scrape object'''
    bm = os.path.join(s.casefold, '0', 'U')
    if os.path.exists(bm):
        with open (bm, "r") as f:
            line = f.readline()
            while not line.startswith('\tbathFlow'):
                line = f.readline()
            # now we have reached the bath flow section
            for i in range(3):
                line = f.readline()
            strs = re.split('\(|\)| ', line)
            s.GEObathv[1] = strs[2] # bath velocity: this will be reported in m/s
            while not line.startswith('\tinkFlow'):
                line = f.readline()
            # now we have reached the ink flow section
            for i in range(3):
                line = f.readline()
            strs = re.split('\(|\)| |-', line)
            s.GEOinkv[1] = strs[5] # ink velocity: this will be reported in m/s
            return
    else:
        logging.info(f'path {bm} does not exist')


def scrapeSHM(s:scrape) -> None:
    '''scrape snappyHexMeshDict
    here we use listLevel because snappyHexMeshDict has a lot of sections and we're scraping all of the data
    to change which fields we collect, go back to the scrape class definition
    we're collecting data for SHMlist, CMClist, CMCfixedWallsLevel, SClist, ALClist, ALCfixedWallsLayers, MQClist, and SHMmergeTolerance
    s is a scrape object'''
    bm = os.path.join(s.meshfold, 'system', 'snappyHexMeshDict')
    if os.path.exists(bm):
        with open (bm, "r") as f:
            line = f.readline()
            # first determine if castellatedMesh, snapping, and layers are active
            line = listLevel('castellatedMesh', 'geometry', line, f, s.SHMlist)
            # then read in castellated mesh controls
            if s.SHMlist[0][1]=='true':
                # only collect castellatedMesh variables if we're using a castellated mesh
                line = listLevel('castellatedMeshControls', '\tfeatures', line, f, s.CMClist)
                while not line.startswith('\t\t\tfile\t\"fixed'):
                    line = f.readline()
                line = f.readline()
                strs = re.split(';|\t', line)
                s.CMCfixedWallsLevel[1] = strs[4]
            if s.SHMlist[1][1]=='true':
                # only collect snap variables if we're using snapping
                line = listLevel('snapControls', 'addLayersControls', line, f, s.SClist)
            if s.SHMlist[2][1]=='true':
                # only collect layers variables if we're using layers
                listLevel('addLayersControls', '\tlayers', line, f, s.ALClist)
                while not line.startswith('\t\tfixed'):
                    line = f.readline()
                for i in range(2):
                    line = f.readline()
                strs = re.split(';|\t', line)
                s.ALCfixedWallsLayers[1] = strs[4]
            line = listLevel('meshQualityControls', 'writeFlags', line, f, s.MQClist)
            line = readLevel0('mergeTolerance', line, f, s.SHMmergeTolerance)
            return
    else:
        logging.warning(f'path {bm} does not exist')


def scrapeDMD(s:scrape) -> None:
    '''scrape dynamicMeshDict
    s is a scrape object'''
    bm = os.path.join(s.casefold, 'constant', 'dynamicMeshDict')
    if os.path.exists(bm):
        with open (bm, "r") as f:
            line = f.readline()
            # just collect the variables in the dynamicRefineFvMeshCoeffs section
            line = listLevel('dynamicRefineFvMeshCoeffs', '\tcorrectFluxes', line, f, s.DMDlist)
            return
    else:
        logging.warning(f'path {bm} does not exist')


def scrapeTP(s:scrape) -> None:
    '''scrape transportProperties
    s is a scrape object'''
    bm = os.path.join(s.casefold, 'constant', 'transportProperties')
    if os.path.exists(bm):
        with open (bm, "r") as f:
            line = f.readline()
            line = listLevel('ink', '}', line, f, s.TPinklist)
            line = listLevel('sup', '}', line, f, s.TPsuplist)
            line = readLevel0('sigma', line, f, s.TPsigma)
            return
    else:
        logging.warning(f'path {bm} does not exist')
        
def scrapeLabels(s:scrape) -> None:
    '''scrape the labels.csv document'''
    bm = os.path.join(s.fold, 'labels.csv')
    if os.path.exists(bm):
        with open(bm, "r") as f:
            data = list(csv.reader(f))
            s.TPinklist[0][1]=data[0][1]
            s.TPsuplist[0][1]=data[1][1]
        

def scrapeCD(s:scrape) -> None:
    '''scrape controlDict
    s is a scrape object'''
    bm = os.path.join(s.casefold, 'system', 'controlDict')
    if os.path.exists(bm):
        with open (bm, "r") as f:
            line = f.readline()
            line = listLevel('application', '//', line, f, s.controlDictList)
            return
    else:
        logging.warning(f'path {bm} does not exist')
        


def scrapeFV(s:scrape) -> None:
    '''scrape fvSolution
        fvSolution files are broken into sections by the variable we're solving for, e.g. alpha, pcorr, p_rgh. 
        scrape each section into a different list stored in s
        s is a scrape object'''
    bm = os.path.join(s.casefold, 'system', 'fvSolution')
    if os.path.exists(bm):
        with open (bm, "r") as f:
            line = f.readline()
            line = listLevel('\t\"alpha', '\t}', line, f, s.fvSailist)
            line = listLevel('\t\"pcorr', '\t}', line, f, s.fvSpcorrlist)
            line = listLevel('\tp_rgh', '\t}', line, f, s.fvSprghlist)
            line = listLevel('\tp_rghFinal', '\t}', line, f, s.fvSprghfinallist)
            line = listLevel('\tU', '\t}', line, f, s.fvSUlist)
            line = listLevel('PIMPLE', '\t}', line, f, s.fvSPIMPLElist)
            return
    else:
        logging.warning(f'path {bm} does not exist')


def scrapeGeo(s:scrape) -> None:
    '''scrape geometry.csv
        geometry.csv was created by ncreate3d.py/noz3dscript.ipynb when we generated the whole folder
        s is a scrape object'''
    bm = os.path.join(s.fold, 'geometry.csv')
    if os.path.exists(bm):
        with open(bm, "r") as f:
            data = list(csv.reader(f))
            s.GEOniw = data[0] # nozzle inner width
            s.GEOnt = data[1] # nozzle thickness
            s.GEObw = data[2] # bath width
            s.GEObd = data[3] # bath depth
            s.GEOnl = data[4] # nozzle length
            s.GEOblc = data[5] # bath left coordinate
            s.GEObrc = data[6] # bath right coordinate
            s.GEObfc = data[7] # bath front coordinate
            s.GEObbackc = data[8] # bath back coordinate
            s.GEObbotc = data[9] # bath bottom coordinate
            s.GEObtc = data[10] # bath top coordinate
            s.GEOnbc = data[11] # nozzle bottom coordinate
            s.GEOncxc = data[12] # nozzle center x coordinate
            s.GEOncyc = data[13] # nozzle center y coordinate
            s.GEOna = data[14] # nozzle angle
            s.GEObathv = data[15] # bath velocity
            s.GEOinkv = data[16] # ink velocity
            return
    else:
        return


#----------------------------------------------------------------------------

def populate(folder:str, *varargin) -> List[List[str]]:
    '''populate scrapes all the data from the folder and exports it to a table called legend.csv
    folder is a full path name
    can also add strings that evaluate specific functions, e.g. scrapeTP, so that you can just scrape the transportproperties if there is already a legend'''
    if not isSimFolder(folder):
        raise Exception("Not a simulation folder")
    s = scrape(folder)   # create an object to store variables
    fn = os.path.join(folder, 'legend.csv')     # export file name
    scrapeLogs(s)   # scrape the logs
    if os.path.exists(fn):
        for func in varargin:
            try:
                eval(func+'(s)')
            except Exception as e:
                print(e)
                pass
        # if a legend file already exists, keep that 
        # legend file and just replace the processing times
        # all of the other variables are set at the folder generation and don't change
        with open(fn, 'r') as f:
            t = list(csv.reader(f))
            t2 = s.table()
            t[2:8] = t2[2:8]
            
            # keep any collected variables
            if len(varargin)>0:
                for i in range(len(t)):
                    if not t[i][1]==t2[i][1] and not t2[i][1]=='' :
                        t[i]=t2[i]
    else:
        # if there is no legend file, we have to go through all the files and scrape data
        scrapeGeo(s)
        scrapeBlockMeshDict(s)
        scrapeSetFieldsDict(s)
        scrapeU(s)
        scrapeSHM(s)
        scrapeDMD(s)
        scrapeTP(s)
        scrapeLabels(s)
        scrapeCD(s)
        scrapeFV(s)
        t = s.table() # generate a table from all the data we scraped
    exportCSV(fn, t) # export the updated table
    return t



def populateToTable(folder:str, repopulate:bool = False) -> List[List[str]]:
    '''take one folder and either import the existing legend or create a new one
    folder is a full path name
    repopulate is true if you want to overwrite existing legend.csv files'''
    if repopulate:
        # overwrite times in legend.csv, or if there is no file, create a new legend file
        t2 = populate(folder)
    else:
        # don't overwrite anything, just import the existing file
        t2 = importIf(folder)
    return t2


def populateList(liInit:List[str], exportFilename:str, repopulate:bool = False) -> None:
    '''populateList scrapes all data for all files in a list and creates a combined table
    liInit is a list of folders to scrape
    exportFilename is the destination to export the summary file
    repopulate is true if you want to overwrite existing legend.csv files'''
    li = []
    for folder in liInit:
        if os.path.exists(folder):
            li.append(folder)
    t1 = populateToTable(li[0], repopulate)   
        # first collect one table so you have the format
    tbig = [['' for i in range(len(li)+1)] for j in range(len(t1))] 
        # tbig combines all the legend files into one table
    for i,row in enumerate(t1):
        tbig[i][0] = row[0] 
            # import the variable names and the values 
            # for the first file into the big table
        tbig[i][1] = row[1]
    for j,l in enumerate(li[1:]): 
        # go through the rest of the files and put the values in the big table
        t2 = populateToTable(l, repopulate)
        for i,row in enumerate(t2):
            tbig[i][j+2] = row[1]
    exportCSV(exportFilename, tbig) 
        # export the combined table
    return
