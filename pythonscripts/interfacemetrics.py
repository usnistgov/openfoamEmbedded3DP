#!/usr/bin/env python
'''Analyzing simulated single filaments'''

# external packages
import sys
import os
import csv
import numpy as np
import pandas as pd
from shapely.geometry import Polygon
import re
from typing import List, Dict, Tuple, Union, Any
import logging
import traceback

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import folderparser as fp
from pvCleanup import addUnits

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

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

# BATHRIGHT = 4.834
# NOZZLEBOTTOM = 0.3015
# NOZZLEDIAMETER = 0.603
# NOZZLETHICKNESS = 0.152
# NOZZLEXCENTER = -2.4120
# NOZZLERIGHTEDGE = NOZZLEXCENTER + NOZZLEDIAMETER/2+NOZZLETHICKNESS
# BEHINDNOZZLE = NOZZLERIGHTEDGE + NOZZLEDIAMETER

#################################################################

class folderStats:
    '''the folderStats class is used to store critical information about the nozzle geometry'''
    
    def __init__(self, folder:str):
        self.folder = folder
        try:
            self.geo = importLegend(folder)
        except:
            raise Exception('No legend file')
        
        i0 = self.geo[self.geo['title']=='nozzle inner width (mm)'].index[0] # find the index where the nozzle inner width is
        self.niw = float(self.geo.loc[i0, 'val']) # nozzle inner width
        self.nt = float(self.geo.loc[i0+1, 'val']) # nozzle thickness
        self.bd = float(self.geo.loc[i0+3, 'val']) # bath depth RG
        self.nl = float(self.geo.loc[i0+4, 'val']) # nozzle length RG
        self.brx = float(self.geo.loc[i0+6, 'val']) # bath right x
        self.nbz = float(self.geo.loc[i0+11, 'val']) # nozzle bottom z
        self.ncx = float(self.geo.loc[i0+12, 'val']) # nozzle center x
        self.ncy = float(self.geo.loc[i0+13, 'val']) # nozzle center y RG
        self.na = float(self.geo.loc[i0+14, 'val']) # nozzle angle RG
        self.bv = float(self.geo.loc[i0+15, 'val'])*1000 # bath velocity (mm/s)
        self.iv = float(self.geo.loc[i0+16, 'val'])*1000 # ink velocity (mm/s)
        self.nre = self.ncx + self.niw/2 + self.nt # nozzle right edge
        self.behind = self.nre + self.niw # 1 inner diameter behind nozzle
        self.intentzcenter = self.nbz - self.niw/2
        self.intentzbot = self.nbz - self.niw
        self.xlist = []
       

    
########### FILE HANDLING ##############


def importLegend(folder:str) -> pd.DataFrame:
    '''import the legend file'''
    file = os.path.join(folder, 'legend.csv')
    if os.path.exists(file):
        return pd.read_csv(os.path.join(folder, 'legend.csv'), names=['title', 'val'])
    else:
        raise Exception('No legend file')
    

def plainIm(file:str, ic:Union[int, bool]) -> Tuple[Union[pd.DataFrame, List[Any]], Dict]:
    '''import a csv to a pandas dataframe. ic is the index column. Int if there is an index column, False if there is none'''
    if os.path.exists(file):
        try:
            d = pd.read_csv(file, index_col=ic, skiprows=[1], dtype=np.float64) # skiprows skips units, dtype defines as floats to save memory RG
            d.columns = map(str.lower, d.columns)
            u = pd.read_csv(file, nrows=1) # read the units seperately since they are not numbers RG
            row1 = list(u.iloc[0])
            if type(row1[0]) is str and ('m' in row1 or 's' in row1):
                unitdict = dict(u.iloc[0])
                skiprows = [1]
            else:
                unitdict = dict([[s,'undefined'] for s in toprows])
                skiprows = []
            d = pd.read_csv(file, index_col=ic, dtype=float, skiprows=skiprows)
            d.columns = map(str.lower, d.columns)
        except:
            traceback.print_exc()
            return [],{}
        return d, unitdict
    else:
        return [], {}


def importFilemm(file:str, slist:List[str]) -> Tuple[Union[pd.DataFrame, List[Any]], Dict]:
    '''Import a file and convert m to mm. File is a full path name. slist is the list of column names for the columns that need mm conversion'''
    d,units = plainIm(file, False)
    mdict = {'m':'mm', 'm/s':'mm/s'}
    if len(d)==0:
        return d,units
    else:
        for s in slist:
            d[s.lower()]*=1000
            units[s]=mdict.get(units[s], units[s]+'*10^3')
        return d,units
    

def importLine(folder:str, time:float, x:float=1.5) -> pd.DataFrame:
    '''import a csv of a line trace pulled from ParaView. time and x are the time and position the line trace are taken at'''
    file = os.path.join(folder, f'line_t_{int(round(time*10))}_x_{x}.csv')
    try:
        data, units = importPointsFile(file) 
    except:
        addUnits(file)
        data, units = importPointsFile(file) 
    return  data, units

def importPointsFile(file:str) -> Tuple[Union[pd.DataFrame, List[Any]], Dict]:
    '''This is useful for importing interfacePoints.csv files. '''
    d,units = plainIm(file, False)
    if len(d)==0:
        return d,units
    try:
        d = d[pd.to_numeric(d.time, errors='coerce').notnull()] # remove non-numeric times
    except:
        pass
    mdict = {'m':'mm', 'm/s':'mm/s'}
    for s in ['x', 'y', 'z', 'vx', 'vy', 'vz', 'arc_length']:
        if s in d:
            try:
                d[s] = d[s].astype(float)
            except:
                raise Exception('Non-numeric values for '+s+' in '+file)
            else:
                d[s]*=1000
            units[s]=mdict.get(units[s], units[s]+'*10^3')
    d = d.sort_values(by='x')
    return d,units
    

def importPoints(folder:str, time:float) -> Tuple[Union[pd.DataFrame, List[Any]], Dict]:
    '''import all points in the interface at a given time'''
    file = os.path.join(folder, 'interfacePoints', 'interfacePoints_t_'+str(int(round(time*10)))+'.csv')
    return importPointsFile(file)

def importPtsNoz(folder:str, time:float) -> Tuple[Union[pd.DataFrame, List[Any]], Dict]: # RG
    '''import all points in the nozzle at a given time'''
    file = os.path.join(folder, 'nozzlePoints', 'nozzlePoints_t_'+str(int(round(time*10)))+'.csv')
    return importPointsFile(file)


def importPtsSlice(folder:str, time:float, x:float) -> Union[pd.DataFrame, List[Any]]:
    '''import points just from one slice. folder is full path name, time is in s, x is absolute position in mm. Finds the closest position to the requested position and gives up if there is no x value within 0.2 mm'''
    pts,units = importPoints(folder, time)
    if len(pts)>0:
        xlist = xpts(pts)
        xreal = closest(xlist, x)
        if abs(xreal-x)>0.2:
            return []
        ptsx = pts[pts['x']==xreal]
        return ptsx
    else:
        return []
    
def importSS(folder:str) -> pd.DataFrame:
    '''import slice summaries. folder is full path name'''
    file = os.path.join(folder, 'sliceSummaries.csv')
    d, units = plainIm(file, 0)
    if len(d)==0:
        return [], []
    try:
        for s in d:
            d[s] = pd.to_numeric(d[s], errors='coerce') # remove non-numeric times
    except Exception as e:
        pass
    d = d.dropna() # drop NA values
    return d, units

def imFn(exportfolder:str, labels:str, topfolder:str, **kwargs) -> str:
    '''Construct an image file name with no extension. Exportfolder is the folder to export to. Label is any given label. Topfolder is the folder this image refers to, e.g. HBHBsweep. Insert any extra values in kwargs as keywords'''
    bn = os.path.basename(topfolder)
    s = ''
    for k in kwargs:
        if not k in ['adjustBounds', 'svg', 'png', 'overwrite', 'split'] and 'list' not in k:
            s = s + k + '_'+str(kwargs[k])+'_'
    s = s[0:-1]
    s = s.replace('*', 'x')
    s = s.replace('/', 'div')
    
    # RG
    label = ''
    labels = [labels] if isinstance(labels, str) else labels
    for i in labels:
        label += i+'_'
    
    return os.path.join(exportfolder, bn, label+bn+'_'+s) # RG


def exportIm(fn:str, fig, svg:bool=True, png:bool=True, **kwargs) -> None:
    '''export an image. fn is a full path name, without the extension. fig is a matplotlib figure'''
    if svg:
        slist = ['.svg']
    else:
        slist = []
    if png:
        slist.append('.png')
    for s in slist:
        fig.savefig(fn+s, bbox_inches='tight', dpi=300)
    logging.info(f'Exported {fn}')



    ########## DATA HANDLING ############
def xpts(data:pd.DataFrame) -> pd.Series:
    '''List of all of the x positions in the dataframe'''
    xl = np.sort(data.x.unique())
    return xl


def splitCentroid(p:Polygon) -> Tuple[float, float, float]:
    '''get the coordinates of a polygon centroid
    p is a Polygon from shapely.geometry'''
    c = p.centroid.wkt
    slist = re.split('\(|\)| ', c)
    x = float(slist[2])
    y = float(slist[3])
    return x,y

def xyFromDf(sli: Union[np.ndarray, pd.DataFrame]) -> np.ndarray:
    '''Get a list of y,z positions from points in a dataframe'''
    if type(sli) is np.ndarray:
        sli2 = sli
    else:
        # sli is a Pandas DataFrame
        sli2 = sli[['y','z']].values.tolist()
    return sli2


def centroid(sli: Union[np.ndarray, pd.DataFrame]) -> Tuple[float, float]:
    '''find centroid and other slice measurements
    sli is a list of points in a slice or a pandas dataframe that contains points'''
    sli2 = xyFromDf(sli)
    p = Polygon(sli2).convex_hull
    return splitCentroid(p)


def centroidAndArea(sli: Union[np.ndarray, pd.DataFrame]) -> Tuple[float, float, float]:
    '''get the centroid and area of a slice, given as a list of points or a dataframe'''
    sli2 = xyFromDf(sli)
    p = Polygon(sli2).convex_hull
    x,y=splitCentroid(p)
    a = p.area
    return x,y,a


def pdCentroid(xs:pd.DataFrame) -> Tuple[float, float, float]:
    '''find the centroid of a slice represented as a pandas dataframe'''
    return centroid(np.array(xs[['y', 'z']]))

 
def closest(lst:List[float], K:float) -> float: # RG
    '''find the closest value in a list to the float K'''
    lst = np.asarray(lst)
    idx = (np.abs(lst - K)).argmin()
    return lst[idx]

############ CREATE FILES ############


#### SUMMARIZING SLICES


def sliceSummary(fs:folderStats, sli:pd.DataFrame) -> Dict[str, float]:
    '''sliceSummary collects important stats from a slice of a filament at a certain x and time and returns as a dictionary
    sli is a subset of points as a pandas dataframe
    fs is a folderStats object'''
    ev = -100 # error value
    rv = {'x':ev, 'xbehind':ev, 'time':ev, \
          'centery':ev, 'centerz':ev, 'area':ev, 'maxheight':ev, 'maxwidth':ev, \
           'centeryn':ev, 'centerzn':ev, 'arean':ev, 'maxheightn':ev, 'maxwidthn':ev,\
          'vertdisp':ev, 'vertdispn':ev, 'aspectratio':ev, 'speed':ev, 'speeddecay':ev} # return value
        # x is in mm
        # time is in s
        # centery, centerz, minz, vertdisp, maxz, maxheight, maxwidth are in mm
        # speed is in mm/s
        # centeryn, centerzn, vertdispn, maxheightn, maxwidthn are normalized by nozzle inner diameter
        # aspectratio, topbotratio are dimensionless
        # speeddecay is normalized by the bath speed
    rv['x'] = sli.iloc[0]['x']
    rv['xbehind'] = rv['x']-fs.ncx
    rv['time'] = np.float64(sli.iloc[0]['time'])
    
    if len(sli)<10:
        #logging.error('Not enough points')
        raise Exception('Not enough points')

    try:
        rv['centery'], rv['centerz'], rv['area'] = centroidAndArea(sli)
    except:
        #logging.error('centroid error')
        raise Exception('centroid error')
    rv['maxheight'] = sli['z'].max()-sli['z'].min()
    rv['maxwidth'] = sli['y'].max()-sli['y'].min()
    if rv['maxheight']==0 or rv['maxwidth']==0:
        #logging.error('Cross-section is too small')
        raise Exception('Cross-section is too small')
    
    rv['centerzn'] = rv['centerz']/fs.niw
    rv['centeryn'] = rv['centery']/fs.niw
    rv['arean'] = rv['area']/(np.pi*(fs.niw/2)**2) # normalize area by nozzle area
    rv['maxheightn'] = rv['maxheight']/fs.niw # normalize height by nozzle diameter
    rv['maxwidthn'] = rv['maxwidth']/fs.niw # normalized width by intended width
    
    rv['vertdisp'] = (sli['z'].min() - fs.intentzbot)
    rv['vertdispn'] = rv['vertdisp']/fs.niw
    rv['aspectratio'] = rv['maxheight']/rv['maxwidth']
    rv['speed'] = sli['vx'].mean() # speed of interface points in x
    rv['speeddecay'] = rv['speed']/fs.bv # relative speed of interface relative to the bath speed

    return rv

def sliceUnits(ipUnits:Dict) -> Dict:
    '''Find the units for a slice summary dictionary based on the units of the interfacePoints csv'''
    xu = ipUnits['x']
    return {'x':xu, 'xbehind':xu, 'time':ipUnits['time'], \
          'centery':xu, 'centerz':xu, 'area':xu+'^2', 'maxheight':xu, 'maxwidth':xu, \
           'centeryn':'', 'centerzn':'', 'arean':'', 'maxheightn':'', 'maxwidthn':'',\
          'vertdisp':xu, 'vertdispn':'', 'aspectratio':'', 'speed':ipUnits['vx'], 'speeddecay':''}
    

def summarizeSlices(fs: folderStats) -> pd.DataFrame:
    '''go through all the interface points files and summarize each x and time slice
    fs should be a folderStats object
    outputs a pandas DataFrame'''
    ipfolder = os.path.join(fs.folder, 'interfacePoints')
    if not os.path.exists(ipfolder):
        raise Exception('No interface points')
    ipfiles = os.listdir(ipfolder)
    if len(ipfiles)==0:
        return [], {}
    s = []
    for f in ipfiles:
        data, units = importPointsFile(os.path.join(ipfolder, f))
        if len(data)>0:
            xlist = xpts(data)
            try:
                xlist = xlist[xlist>fs.behind]
            except Exception as e:
                logging.error(f+': '+e+', '+str(xlist))
                raise e
            for x in xlist:
                sli = data[data['x']==x]
                if len(sli)>9:
                    try:
                        ss1 = sliceSummary(fs, sli)
                    except Exception as e:
                        pass
                    else:
                        s.append(ss1)
    sliceSummaries = pd.DataFrame(s, dtype=np.float64)
    sliceSummaries = sliceSummaries.dropna()
    return sliceSummaries, units


def ssFile(folder:str) -> str:
    '''slice Summaries file name'''
    return os.path.join(folder, 'sliceSummaries.csv')


def summarize(folder:str, overwrite:bool) -> int:
    '''given a simulation, get critical statistics on the filament shape
    folder is the full path name to the folder holding all the files for the simulation
        there must be an interfacePoints folder holding csvs of the interface points
    overwrite true to overwrite existing files, false to only write new files
    a return value of 0 indicates success
    a return value of 1 indicates failure'''
    cf = fp.caseFolder(folder)
    if not os.path.exists(cf):
        return
    if os.path.exists(folder):
        sspath = ssFile(folder)
        if not os.path.exists(sspath) or overwrite:
            
            # this trycatch makes sure that there is a legend file
            try:
                fs = folderStats(folder)
            except:
                return 1
            
            # summarizeSlices will raise an exception if there are no interface
            # points files
            try:
                sliceSummaries, ipunits = summarizeSlices(fs)
            except Exception as e:
                logging.error(str(e))
                return 1
            
            # if there are points in the summary, save them
            if len(sliceSummaries)>0:
                sliceSummaries = sliceSummaries.sort_values(by=['time', 'x'])
                if len(sliceSummaries)>0:
                    # add a units header
                    units = sliceUnits(ipunits)
                    unitrow = pd.DataFrame(dict([[s,units[s]] for s in sliceSummaries]), index=[0])
                    sliceSummaries = pd.concat([unitrow, sliceSummaries]).reset_index(drop=True)
                sliceSummaries.to_csv(sspath)
                logging.info(f'    Exported {sspath}')
            else:
                logging.info(f'No slices recorded in {folder}')
                header = ['x', 'xbehind', 'time', 'centery', 'centerz', 'area', 'maxheight', 'maxwidth', 'centeryn', 'centerzn', 'arean', 'maxheightn', 'maxwidthn','vertdisp', 'vertdispn', 'aspectratio', 'speed', 'speeddecay']
                sliceSummaries = pd.DataFrame([], columns=header)
                sliceSummaries.to_csv(sspath)
                return 1
    else:
        return 1
    return 0



####### DETERMINING STEADY STATE FILAMENT SHAPE


def steadyList(folder:str, dother:float, vdcrit:float, col:str, mode:str) -> pd.DataFrame:
    '''steadyList determines a list of times and positions which have reached steady state, for either steady in time or steady in position.
    folder is a full path name
    dt is the size of the chunk of time in s over which we want to evaluate if the chunk is "steady", e.g. 1
    vdcrit is the maximum range of the variable of interest to be considered "steady", e.g. 0.01
    col is the column name of the variable we're watching to see if it's steady, e.g. 'vertdispn'
    mode is 'xbehind' or 'time': 
        'xbehind' means that for each x, we're finding a range of steady times. 
        'time' means that for each time, we're finding a range of steady xbehinds'''
    try:
        ss, ssunits = importSS(folder)
    except:
        raise Exception('Problem with summaries')
    
    if mode=='xbehind': # mode is the variable that we use to split into groups
        other='time' # other is the variable that we scan across
        outlabel = ['x', 't0', 'tf']
        flatlist = [[ssunits['xbehind'], ssunits['time'], ssunits['time']]]
    else:
        other='xbehind'
        outlabel = ['t', 'x0', 'xf']
        flatlist = [[ssunits['time'], ssunits['xbehind'], ssunits['xbehind']]]
        
    if len(ss)<2:
        return pd.DataFrame([], columns=outlabel)
    
    list1 = ss[mode].unique() # this gets the list of unique values for the mode variable
    for xtmode in list1:
        l1 = ss[ss[mode]==xtmode] # get this slice in x if mode is x, time if mode is time
        l1 = l1.sort_values(by=other) # sort the slice by time if mode is x, x if mode is time
        list2 = l1[other].unique()
        vdrange = 100
        i = -1
        modevar = xtmode
        while vdrange>vdcrit and i+1<len(list2):
            i+=1
            xtother = list2[i]
            l2 = l1[(l1[other]>=xtother-dother/2) & (l1[other]<=xtother+dother/2)] # get a chunk of size dt centered around this time
            vdrange = l2[col].max()-l2[col].min()
        if i+1<len(list2):
            other0 = xtother
            while vdrange<=vdcrit and i+1<len(list2):
                i+=1
                xtother = list2[i]
                l2 = l1[(l1[other]>=xtother-dother/2) & (l1[other]<=xtother+dother/2)] # get a chunk of size dt centered around this time
                vdrange = l2[col].max()-l2[col].min()
            if i+1<len(list2):
                otherf = xtother
            else:
                otherf = 1000
            if mode=='xbehind':
                modevar=round(modevar,3)
            else:
                other0 = round(other0,3)
                otherf = round(otherf,3)
            flatlist.append([modevar, other0, otherf]) # stable within this margin
    return pd.DataFrame(flatlist, columns=outlabel)


def steadyTime(folder:str, dt:float, vdcrit:float, col:str) -> pd.DataFrame:
    '''get steady in time stats'''
    return steadyList(folder, dt, vdcrit, col, 'xbehind')


def steadyPos(folder:str, dx:float, vdcrit:float, col:str) -> pd.DataFrame:
    '''get steady in position stats'''
    return steadyList(folder, dx, vdcrit, col, 'time')



def stFile(folder:str) -> str:
    '''steady time file name'''
    return os.path.join(folder, 'steadyTimes.csv')

def spFile(folder:str) -> str:
    '''steady positions file name'''
    return os.path.join(folder, 'steadyPositions.csv')


def steadyMetric(folder:str, mode:int, overwrite:bool) -> int:
    '''exports steadyTimes or steadyPositions
    folder is full path name
    mode 0 for steady times, 1 for steady positions
    overwrite true to overwrite existing files'''
    if mode==0:
        fn = stFile(folder)
        f = steadyTime
    else:
        fn = spFile(folder)
        f = steadyPos
    if not os.path.exists(fn) or overwrite:
        try:
            tab = f(folder, 1, 0.01, 'vertdispn') # returns a table
        except Exception as e:
            logging.error(str(e))
            raise e
        tab.to_csv(fn)
        logging.info(f'    Exported {fn}')
    return 0


def steadyMetrics(folder:str, overwrite:bool) -> int:
    '''exports steadyTimes and steadyPositions
    folder is full path name
    overwrite true to overwrite existing files'''
    cf = fp.caseFolder(folder)
    if not os.path.exists(cf):
        return
    for mode in [0,1]:
        try:
            steadyMetric(folder, mode, overwrite)
        except:
            return 1
    return 0


def sumAndSteady(folder:str, overwrite:bool) -> None:
    '''summarize slices and find steady state for a simulation
    folder is full path name
    overwrite true to overwrite existing files, false to only write new ones'''
    
    if not overwrite:
        ssfn = ssFile(folder)
        stfn = stFile(folder)
        spfn = spFile(folder)

        if os.path.exists(ssfn) and os.path.exists(stfn) and os.path.exists(spfn):
            return
    
    # if there are no interface points in the folder, 
    # we won't be able to analyze this folder
    cf = os.path.join(folder, 'interfacePoints')
    if not os.path.exists(cf):
        return 

    logging.debug(folder)
    
    # get initial summaries of the slices in time and position
    sret = summarize(folder, overwrite)
    
    # if sret is 1, summarizing failed, and we won't be able to gather steady metrics
    if sret==0:
        steadyMetrics(folder, overwrite)
    return 
    
    
    
##### ARCHIVE

    
# def correctPointsFiles(folder:str):
#     ipfolder = os.path.join(folder, 'interfacePoints')
#     if not os.path.exists(ipfolder):
#         return
#     for file in os.listdir(ipfolder):
#         correctPointsFile(os.path.join(ipfolder, file))
    
# def correctPointsFile(csvfile:str):
#     blist = []
#     with open(csvfile, 'r') as b:
#         csv_reader = csv.reader(b)
#         line1 = next(csv_reader)
#         line2 = next(csv_reader)
#         if len(line1)==len(line2):
#             return
#         for i in range(len(line2)-len(line1)):
#             line1.append('col'+str(i))
#         blist.append(line1)
#         blist.append(line2)
#         blist = blist+list(csv_reader)
#     with open(csvfile, 'w', newline='', encoding='utf-8') as b:
#         writer = csv.writer(b)
#         for row in blist:
#             #print(row)
#             writer.writerow(row)
            
# def deleteBadPointsFile(csvfile:str):
#     if os.path.exists(csvfile):
#         if 'temp' in csvfile:
#             os.remove(csvfile)
#             return 1
#         else:
#             with open(csvfile, 'r') as b:
#                 csv_reader = csv.reader(b)
#                 try:
#                     line1 = next(csv_reader)
#                     line2 = next(csv_reader)
#                 except:
#                     os.remove(csvfile)
#                     return 1
#                 if 'col' in line1:
#                     os.remove(csvfile)
#                     return 1
#                 if len(line1)==len(line2):
#                     return 0
#             os.remove(csvfile)
#             return 1
#     return 0


# def exportImage(exportfolder:str, label:str, topfolder:str, cp) -> None:
#     '''export an image from a comboplot or gridofplots object, given an export folder, a label, a folder this came from, and a comboPlot object'''
    
#     if len(cp.xlistreal)==0: # don't export an empty plot
#         return
#     fn = imFn(exportfolder, label, topfolder)
#     exportIm(fn, cp.fig)



    
    
# def addUnitsAll(serverfolder):
#     for topfolder in [os.path.join(serverfolder, s) for s in ['HBHBsweep', 'HBnewtsweep', 'newtHBsweep', 'newtnewtsweep']]:
#         for cf in fp.caseFolders(topfolder):
#             print(cf)
#             for ipfile in os.listdir(cf):
#                 if 'line' in ipfile and 'csv' in ipfile:
#                     f = os.path.join(cf, ipfile)
#                     if os.path.exists(f):
#                         try:
#                             addUnits(f)
#                         except Exception as e:
#                             logging.error(f'ERROR in {f}: {str(e)}')
#                             return
#     return


# def correctName(csvfile:str) -> None:
#     d,units = importPointsFile(csvfile)
#     if len(d)==0:
#         return
#     xabs = (d['x'].unique())[0]
#     xrel = xabs+2.412
#     folder = os.path.dirname(csvfile)
#     bn = os.path.basename(csvfile)
#     csvsplit = re.split('x_', bn)
#     newname = os.path.join(folder, csvsplit[0]+'x_{:.1f}.csv'.format(xrel))
#     if newname==csvfile:
#         return
#     idx = 2
#     while os.path.exists(newname):
#         d2, units2 = importPointsFile(newname)
#         if len(d2)==0:
#             os.remove(newname)
#         elif d2.equals(d) and units2==units:
#             os.remove(csvfile)
#             return
#         else:
#             newname = newname.replace('.csv', f' ({idx}).csv')
#             idx+=1
#     os.rename(csvfile, newname)

