import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import seaborn as sns
import statistics as st
from PIL import Image
import scipy as sp
from shapely.geometry import Polygon
import re
from folderparser import caseFolder, caseFolders
from typing import List, Dict, Tuple, Union, Any
import interfacemetricsplots as intmp

# BATHRIGHT = 4.834
# NOZZLEBOTTOM = 0.3015
# NOZZLEDIAMETER = 0.603
# NOZZLETHICKNESS = 0.152
# NOZZLEXCENTER = -2.4120
# NOZZLERIGHTEDGE = NOZZLEXCENTER + NOZZLEDIAMETER/2+NOZZLETHICKNESS
# BEHINDNOZZLE = NOZZLERIGHTEDGE + NOZZLEDIAMETER

# the folderStats class is used to store critical information about the nozzle geometry
class folderStats:
    def __init__(self, folder:str):
        self.folder = folder
        try:
            self.geo = importLegend(folder)
        except:
            raise Exception('No legend file')
        
        i0 = self.geo[self.geo['title']=='nozzle inner width (mm)'].index[0] # find the index where the nozzle inner width is
        self.niw = float(self.geo.loc[i0, 'val']) # nozzle inner width
        self.nt = float(self.geo.loc[i0+1, 'val']) # nozzle thickness
        self.brx = float(self.geo.loc[i0+6, 'val']) # bath right x
        self.nbz = float(self.geo.loc[i0+11, 'val']) # nozzle bottom z
        self.ncx = float(self.geo.loc[i0+12, 'val']) # nozzle center x
        self.bv = float(self.geo.loc[i0+14, 'val'])*1000 # bath velocity (mm/s)
        self.iv = float(self.geo.loc[i0+15, 'val'])*1000 # ink velocity (mm/s)
        self.nre = self.ncx + self.niw/2 + self.nt # nozzle right edge
        self.behind = self.nre + self.niw # 1 inner diameter behind nozzle
        self.intentzcenter = self.nbz - self.niw/2
        self.intentzbot = self.nbz - self.niw
        self.xlist = []
       

    
########### FILE HANDLING ##############

# import the legend file
def importLegend(folder:str) -> pd.DataFrame:
    file = os.path.join(folder, 'legend.csv')
    if os.path.exists(file):
        return pd.read_csv(os.path.join(folder, 'legend.csv'), names=['title', 'val'])
    else:
        raise Exception('No legend file')
    
# import a csv to a pandas dataframe
# ic is the index column. Int if there is an index column, False if there is none
def plainIm(file:str, ic:Union[int, bool]) -> Union[pd.DataFrame, List[Any]]:
    if os.path.exists(file):
        try:
            d = pd.read_csv(file, index_col=ic)
            d.columns = map(str.lower, d.columns)
        except:
            return []
        return d
    else:
        return []

# file is a full path name
# slist is the list of column names for the columns that need mm conversion
def importFilemm(file:str, slist:List[str]) -> Union[pd.DataFrame, List[Any]]:
    d = plainIm(file, False)
    if len(d)==0:
        return d
    else:
        for s in slist:
            d[s.lower()]*=1000
        return d
    
def correctPointsFiles(folder:str):
    ipfolder = os.path.join(folder, 'interfacePoints')
    if not os.path.exists(ipfolder):
        return
    for file in os.listdir(ipfolder):
        correctPointsFile(os.path.join(ipfolder, file))
    
def correctPointsFile(csvfile:str):
    blist = []
    with open(csvfile, 'r') as b:
        csv_reader = csv.reader(b)
        line1 = next(csv_reader)
        line2 = next(csv_reader)
        if len(line1)==len(line2):
            return
        for i in range(len(line2)-len(line1)):
            line1.append('col'+str(i))
        blist.append(line1)
        blist.append(line2)
        blist = blist+list(csv_reader)
    with open(csvfile, 'w', newline='', encoding='utf-8') as b:
        writer = csv.writer(b)
        for row in blist:
            #print(row)
            writer.writerow(row)
    
def importPointsFile(file:str) -> Union[pd.DataFrame, List[Any]]:
    d = plainIm(file, False)
    if len(d)==0:
        return []
    try:
        d = d[pd.to_numeric(d.time, errors='coerce').notnull()] # remove non-numeric times
    except:
        pass
    for s in ['x', 'y', 'z', 'vx', 'vy', 'vz']:
        d[s] = d[s].astype(float)
        if type(d[s][0]) is str:
            raise Exception('Non-numeric values for '+s)
        else:
            d[s]*=1000
    d = d.sort_values(by='x')
    return d
    
# import all points in the interface at a given time
def importPoints(folder:str, time:float) -> Union[pd.DataFrame, List[Any]]:
    file = os.path.join(folder, 'interfacePoints', 'interfacePoints_t_'+str(int(round(time*10)))+'.csv')
    return importPointsFile(file)

# import points just from one slice
# folder is full path name
def importPtsSlice(folder:str, time:float, x:float) -> Union[pd.DataFrame, List[Any]]:
    pts = importPoints(folder, time)
    if len(pts)>0:
        xlist = xpts(pts)
        xreal = closest(xlist, x)
        #print(xreal, x, abs(xreal-x))
        if abs(xreal-x)>0.2:
            return []
        ptsx = pts[pts['x']==xreal]
        return ptsx
    else:
        return []
    


# import slice summaries
# folder is full path name
def importSS(folder:str) -> pd.DataFrame:
    file = os.path.join(folder, 'sliceSummaries.csv')
    return plainIm(file, 0)

def imFn(exportfolder:str, label:str, topfolder:str, **kwargs) -> str:
    bn = os.path.basename(topfolder)
    s = ''
    for k in kwargs:
        s = s + k + '_'+str(kwargs[k])+'_'
    s = s[0:-1]
    s = s.replace('*', 'x')
    s = s.replace('/', 'div')
    return os.path.join(exportfolder, bn, label+'_'+bn+'_'+s)

# export an image
# fn is a full path name, without the extension
# fig is a matplotlib figures
def exportIm(fn:str, fig:plt.Figure) -> None:
    for s in ['.svg', '.png']:
        fig.savefig(fn+s, bbox_inches='tight', dpi=300)
    print('Exported ', fn)

# export an image from a comboplot or gridofplots object
def exportImage(exportfolder:str, label:str, topfolder:str, mode:int, cp:Union[intmp.comboPlot, intmp.gridOfPlots]) -> None:
    # don't export an empty plot
    if len(cp.xlistreal)==0:
        return
    fn = imFn(exportfolder, label, topfolder, mode)
    exportIm(fn, cp.fig)
 
    
    # import a csv of a line trace pulled from ParaView. 
# These are not automatically created for all folders
def importLine(folder:str, time:float) -> pd.DataFrame:
    file = os.path.join(folder, 'line_t_'+str(int(round(time*10)))+'_x_1.5.csv')
    return importFilemm(file, ['U:0', 'U:1', 'U:2', 'arc_length', 'Points:0', 'Points:1', 'Points:2'])  


    ########## DATA HANDLING ############
def xpts(data:pd.DataFrame) -> pd.Series:
    xl = np.sort(data.x.unique())
    return xl

# get the coordinates of a polygon centroid
    # p is a Polygon from shapely.geometry
def splitCentroid(p:Polygon) -> Tuple[float, float, float]:
    c = p.centroid.wkt
    slist = re.split('\(|\)| ', c)
    x = float(slist[2])
    y = float(slist[3])
    return x,y

def xyFromDf(sli: Union[np.ndarray, pd.DataFrame]) -> np.ndarray:
    if type(sli) is np.ndarray:
        sli2 = sli
    else:
        # sli is a Pandas DataFrame
        sli2 = sli[['y','z']].values.tolist()
    return sli2

# find centroid and other slice measurements
    # sli is a list of points in a slice
def centroid(sli: Union[np.ndarray, pd.DataFrame]) -> Tuple[float, float]:
    sli2 = xyFromDf(sli)
    p = Polygon(sli2).convex_hull
    return splitCentroid(p)

# get the centroid and area
def centroidAndArea(sli: Union[np.ndarray, pd.DataFrame]) -> Tuple[float, float, float]:
    sli2 = xyFromDf(sli)
    p = Polygon(sli2).convex_hull
    x,y=splitCentroid(p)
    a = p.area
    return x,y,a

# find the centroid of a slice represented as a pandas dataframe
def pdCentroid(xs:pd.DataFrame) -> Tuple[float, float, float]:
    return centroid(np.array(xs[['y', 'z']]))

# find the closest value in a list to the float K    
def closest(lst:List[float], K:float) -> float:       
    return lst[min(range(len(lst)), key = lambda i: abs(lst[i]-K))]   

############ CREATE FILES ############


#### SUMMARIZING SLICES

# sliceSummary collects important stats from a slice of a filament at a certain x and time
    # sli is a subset of points as a pandas dataframe
    # fs is a folderStats object
def sliceSummary(fs:folderStats, sli:pd.DataFrame) -> Dict[str, float]:
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
    rv['time'] = sli.iloc[0]['time']
    
    if len(sli)<10:
        #print('Not enough points')
        raise Exception('Not enough points')

    try:
        rv['centery'], rv['centerz'], rv['area'] = centroidAndArea(sli)
    except:
        #print('centroid error')
        raise Exception('centroid error')
    rv['maxheight'] = sli['z'].max()-sli['z'].min()
    rv['maxwidth'] = sli['y'].max()-sli['y'].min()
    if rv['maxheight']==0 or rv['maxwidth']==0:
        #print('Cross-section is too small')
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


# go through all the interface points files and summarize each x and time slice
# fs should be a folderStats object
# outputs a pandas DataFrame
def summarizeSlices(fs: folderStats) -> pd.DataFrame:
    ipfolder = os.path.join(fs.folder, 'interfacePoints')
    if not os.path.exists(ipfolder):
        raise Exception('No interface points')
    ipfiles = os.listdir(ipfolder)
    s = []
    for f in ipfiles:
        data = importPointsFile(os.path.join(ipfolder, f))
        if len(data)>0:
            xlist = xpts(data)
            try:
                xlist = xlist[xlist>fs.behind]
            except Exception as e:
                print(f, e, xlist)
                raise e
            for x in xlist:
                sli = data[data['x']==x]
                if len(sli)>9:
                    try:
                        ss1 = sliceSummary(fs, sli)
                    except Exception as e:
                        #print(e)
                        pass
                    else:
                        s.append(ss1)
    slicesummaries = pd.DataFrame(s)
    return slicesummaries

# slice Summaries file name
def ssFile(folder:str) -> str:
    return os.path.join(folder, 'sliceSummaries.csv')

# given a simulation, get critical statistics on the filament shape
# folder is the full path name to the folder holding all the files for the simulation
# there must be an interfacePoints folder holding csvs of the interface points
# overwrite true to overwrite existing files, false to only write new files
# a return value of 0 indicates success
# a return value of 1 indicates failure
def summarize(folder:str, overwrite:bool) -> int:
    cf = caseFolder(folder)
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
                slicesummaries = summarizeSlices(fs)
            except Exception as e:
                print(e)
                return 1
            
            # if there are points in the summary, save them
            if len(slicesummaries)>0:
                slicesummaries = slicesummaries.sort_values(by=['time', 'x'])
                slicesummaries.to_csv(sspath)
                print('    Exported', sspath)
            else:
                print('No slices recorded in ', folder)
                header = ['x', 'xbehind', 'time', 'centery', 'centerz', 'area', 'maxheight', 'maxwidth', 'centeryn', 'centerzn', 'arean', 'maxheightn', 'maxwidthn','vertdisp', 'vertdispn', 'aspectratio', 'speed', 'speeddecay']
                slicesummaries = pd.DataFrame([], columns=header)
                slicesummaries.to_csv(sspath)
                return 1
    else:
        return 1
    return 0



####### DETERMINING STEADY STATE FILAMENT SHAPE

# steadyTime determines a list of times and positions which have reached steady state
    # folder is a full path name
    # dt is the size of the chunk of time in s over which we want to evaluate if the chunk is "steady", e.g. 1
    # vdcrit is the maximum range of the variable of interest to be considered "steady", e.g. 0.01
    # col is the column name, e.g. 'vertdispn'
    # mode is 'x' or 'time': 
        # 'x' means that for each x, we're finding a steady time. 
        # 'time' means that for each time, we're finding a steady x
def steadyList(folder:str, dother:float, vdcrit:float, col:str, mode:str) -> pd.DataFrame:
    if mode=='xbehind': # mode is the variable that we use to split into groups
        other='time' # other is the variable that we scan across
        outlabel = ['x', 't0', 'tf']
    else:
        other='xbehind'
        outlabel = ['t', 'x0', 'xf']
    try:
        ss = importSS(folder)
    except:
        raise Exception('Problem with summaries')
    if len(ss)<2:
        return pd.DataFrame([], columns=outlabel)
    
    list1 = ss[mode].unique()
    flatlist = []
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

# steady in time
def steadyTime(folder:str, dt:float, vdcrit:float, col:str) -> pd.DataFrame:
    return steadyList(folder, dt, vdcrit, col, 'xbehind')

# steady in position
def steadyPos(folder:str, dx:float, vdcrit:float, col:str) -> pd.DataFrame:
    return steadyList(folder, dx, vdcrit, col, 'time')



def stFile(folder:str) -> str:
    return os.path.join(folder, 'steadyTimes.csv')

def spFile(folder:str) -> str:
    return os.path.join(folder, 'steadyPositions.csv')

# exports steadyTimes or steadyPositions
# folder is full path name
# mode 0 for steady times, 1 for steady positions
# overwrite true to overwrite existing files
def steadyMetric(folder:str, mode:int, overwrite:bool) -> int:
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
            print(e)
            raise e
        tab.to_csv(fn)
        print('    Exported', fn)
    return 0

# exports steadyTimes and steadyPositions
# folder is full path name
# overwrite true to overwrite existing files
def steadyMetrics(folder:str, overwrite:bool) -> int:
    cf = caseFolder(folder)
    if not os.path.exists(cf):
        return
    for mode in [0,1]:
        try:
            steadyMetric(folder, mode, overwrite)
        except:
            return 1
    return 0


# summarize and find steady state for a simulation
# folder is full path name
# overwrite true to overwrite existing files, false to only write new ones
def sumAndSteady(folder:str, overwrite:bool) -> None:
    
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

    print(folder)
    
    # get initial summaries of the slices in time and position
    sret = summarize(folder, overwrite)
    
    # if sret is 1, summarizing failed, and we won't be able to gather steady metrics
    if sret==0:
        steadyMetrics(folder, overwrite)
    return 
    

