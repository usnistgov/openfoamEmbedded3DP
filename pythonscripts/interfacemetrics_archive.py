
###### ARCHIVE

 ########## subfunctions from folderstats
    
class folderstats:
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
        self.data = []
        self.xlist = []
        
    def xslice(self, x):
        return self.data[self.data['x']==x]
    
    def xslices(self, xlist):
        return self.data[self.data['x'].isin(xlist)].sort_values(by='x')
    
    def xslices2(self, xmin, xmax, numslices):
        xlist2 = self.xlist[(self.xlist>xmin)*(self.xlist<xmax)]
        if len(xlist2)<numslices:
            spacing = 1
        else:
            spacing = int(np.floor(len(xlist2)/numslices))
        return self.xslices(xlist2[0:-1:spacing])

def importCSV(folder:str) -> Union[pd.DataFrame, List[Any]]:
    file = os.path.join(folder, 'interfacePoints.csv')
    return importpointsfile(file)

def summarize(sli):
        # locations in the sli list
    il = int(sli[['y']].idxmin()) # left
    ir = int(sli[['y']].idxmax()) # right
    it = int(sli[['z']].idxmax()) # top
    ib = int(sli[['z']].idxmin()) # bottom
    
    # extreme points
    ptl = sli.loc[il] # left
    ptr = sli.loc[ir] # right
    ptt = sli.loc[it] # top
    ptb = sli.loc[ib] # bottom  

    
        
    toppts = quadtr+quadtl
    botpts = quadbl+quadbr

    quadtr = sli[(sli['y']>ptt['y']) & (sli['z']>ptr['z'])] # points in top right quadrant of filament
    quadtl = sli[(sli['y']<ptt['y']) & (sli['z']>ptr['z'])] # points in top left quadrant of filament
    quadbl = sli[(sli['y']<ptt['y']) & (sli['z']<ptr['z'])] # points in bottom left quadrant of filament
    quadbr = sli[(sli['y']>ptt['y']) & (sli['z']<ptr['z'])] # points in bottom right quadrant of filament

#     intt = intzdy(quadtr) + intzdy(quadtl) # integral of the top segment of the filament
#     intb = intzdy(quadbr) + intzdy(quadbl) # integral of the bottom segment of the filament
#     rv['centerz'] = (intt-intb)/(intdy(quadtr)+intdy(quadtl)+intdy(quadbr)+intdy(quadbl)) # z coord of centroid
        
#     intl = intydz(quadtl) + intydz(quadbl) # integral of the left segment of the filament
#     intr = intydz(quadbr) + intydz(quadtr) # integral of the left segment of the filament
#     rv['centery'] = (intr+intl)/(intdz(quadtr)+intdz(quadtl)+intdz(quadbr)+intdz(quadbl)) # z coord of centroid
    
#     intmidhor = np.trapz([ptl['z'], ptr['z']], x=[ptl['y'], ptr['y']]) # integral of the midline between the left and right points
#     at = intt - intmidhor # area of the top segment
#     ab = intmidhor - intb # area of the bottom segment
#     if ab>0:
#         rv['topbotratio'] = at/ab

  def sortcircle(data):
    d = data.loc[:, ['x', 'y','z']]
    dleft = (d[d['y']<0]).sort_values(by='z')
    dright = ((d[d['y']>=0]).sort_values(by='z')).iloc[::-1]
    d = pd.concat([dleft, dright])
    return d

def ptangle(c, yc, zc):
    v1 = [1,0]
    z = c['z']-zc
    v2 = [c['y'] - yc, z]
    u1 = v1/np.linalg.norm(v1)
    u2 = v2/np.linalg.norm(v2)
    dp = np.dot(u1, u2)
    angle = np.arccos(dp)*360/(2*np.pi)
    if z>0:
        angle = 360-angle
    return angle

def closest(lst, K):       
    return lst[min(range(len(lst)), key = lambda i: abs(lst[i]-K))]      
    
    
    def intzdy(pts):
    p2 = pts.sort_values(by='y')
    return np.trapz(p2['z'], x=p2['y'])

def intdy(pts):
    p2 = pts.sort_values(by='y')
    return np.trapz(np.ones(len(pts)), x=p2['y'])

def intydz(pts):
    p2 = pts.sort_values(by='z')
    return np.trapz(p2['y'], x=p2['z'])

def intdz(pts):
    p2 = pts.sort_values(by='z')
    return np.trapz(np.ones(len(pts)), x=p2['z'])
    
    ###### archive
def plotSummaries(slicesummaries, t, columns):
    grps = [['centery', 'centerz', 'minz', 'maxz'], ['vertdisp', 'maxwidth', 'maxheight'], ['maxwidthn', 'maxheightn', 'aspectratio', 'topbotratio', 'speeddecay']] # plot 1,2,3
    axlist = []
    for c in columns:
        for i in range(len(grps)):
            if c in grps[i]:
                axlist.append(i)  # axlist: for each requested column, give the plot number
    plots = np.unique(axlist) # find which plots are requested (out of 1,2,3)
    plotcols = [[] for i in range(len(plots))]
    for i, p in enumerate(plots): # i is actual plot number, p is initial plot number
        for j, a in enumerate(axlist): # j is column number, a is the initial plot number
            if a==p:
                 plotcols[i].append(columns[j])
    fig, axes = plt.subplots(nrows=len(plots), ncols=1, sharex='col', figsize=(4, 4*len(plots)))
    l = slicesummaries
    l = l[l['time']==t]
    l = l[l['topbotratio']>0]    
    for i,p in enumerate(plotcols): # i is plotnum, p is col number list
        colors = sns.color_palette("Paired", len(p))
        for j, col in enumerate(p):
            l.plot.scatter(x='x', y=col, color=colors[j], ax=axes[i], label=col)
        axes[i].legend(bbox_to_anchor=(1, 1))
        axes[i].set_xlabel('x (mm)')
        axes[i].set_ylabel(None)
        axes[i].set_aspect(1.0/axes[i].get_data_ratio(), adjustable='box')
    
        
def plotST(steadytime, columns):
    grps = [['centery', 'centerz', 'minz', 'maxz', 'centerx', 'x0', 'xf'], ['vertdisp', 'maxwidth', 'maxheight'], ['maxwidthn', 'maxheightn', 'aspectratio', 'topbotratio', 'speeddecay']] # plot 1,2,3
    axlist = []
    for c in columns:
        for i in range(len(grps)):
            if c in grps[i]:
                axlist.append(i)  # axlist: for each requested column, give the plot number
    plots = np.unique(axlist) # find which plots are requested (out of 1,2,3)
    plotcols = [[] for i in range(len(plots))]
    for i, p in enumerate(plots): # i is actual plot number, p is initial plot number
        for j, a in enumerate(axlist): # j is column number, a is the initial plot number
            if a==p:
                 plotcols[i].append(columns[j])
    fig, axes = plt.subplots(nrows=len(plots), ncols=1, sharex='col', figsize=(4, 4*len(plots)))
    l = steadytime  
    for i,p in enumerate(plotcols): # i is plotnum, p is col number list
        colors = sns.color_palette("Paired", len(p))
        for j, col in enumerate(p):
            ster = col+'_sterr'
            if ster in l.columns:
                ye = ster
            else:
                ye = None
            l.plot.scatter(x='time', y=col, yerr=ye, color=colors[j], ax=axes[i], label=col)
        axes[i].legend(bbox_to_anchor=(1, 1))
        axes[i].set_xlabel('time (s)')
        axes[i].set_ylabel(None)
        axes[i].set_aspect(1.0/axes[i].get_data_ratio(), adjustable='box')
        



def steadytimeplot(t1, col):
    p = t1.plot.scatter(x='time', y=col, yerr=(col+'_sterr'))
    p.set_aspect(1.0/p.get_data_ratio(), adjustable='box')
    
    
#### cross-section plot defined by plateau position and arbitrary plateau time
def XSplott(folder, t, xpv):
    # t is a time
    # xpv is a XSplotvars object
    if not os.path.exists(folder):
        return
    xs = importXS(folder, t)
    if len(xs)>0:
        XSplot(xs, folder, xpv)
    
def XSplotst(folderlist, t, xpv):
    for folder in folderlist:
        XSplott(folder, t, xpv)
    viscplotclean(xpv)
    xpv.ax.set_title('Steady state cross sections at t='+str(t))
    
#### cross-section plot defined by arbitrary x position and t time  \

def XSplotxt(folder, x, t, xpv):
    # x is a distance behind the back of the nozzle
    # t is a time
    # xpv is a XSplotvars object
    xs = XSselection(folder, x, t)
    if len(xs)>0:
        XSplot(xs, folder, xpv)


def XSplotsxt(folderlist, x, t, xpv):
    for folder in folderlist:
        XSplotxt(folder, x, t, xpv)
    viscplotclean(xpv)
    xpv.ax.set_title('Cross sections at x='+str(x)+' and t='+str(t))

    
def XSselection(folder, x, t):
    if not os.path.exists(folder):
        return []
    xs = importCSV(folder)
    if len(xs)==0:
        return []
    xi = BEHINDNOZZLE+x
    xlist = xpts(xs)
    xf = closest(xlist, xi)
    if abs(xf-xi)>0.2:
        return []
    xs = xs[(xs['x']==xf)&(xs['time']==t)]
    return xs

def exportslice(folder, distance, time):
    if os.path.exists(folder):
        sspath = os.path.join(folder, 'slice_d'+str(distance)+'_t'+str(int(round(time*10)))+'.csv')
        if not os.path.exists(sspath):
            xs = XSselection(folder, (distance-1)*NOZZLEDIAMETER, time)
            if len(xs)>0:
                xs.to_csv(sspath, index_label=False)
                print('Exported ', sspath)
    return
            
def stdevchunk(ss, xmin, xmax, scrit):
    subset = ss[(ss['x']>=xmin) & (ss['x']<xmax)]
    stdevs = [0,0,0]
    for i,a in enumerate(['vertdisp', 'maxheight', 'maxwidth']):
        stdevs[i] = st.stdev(subset[a])
    return stdevs[0]<scrit and stdevs[1]<scrit and stdevs[2]<scrit
            
def findsteady(ssfull, t, scrit, dx):
    # ss is a sliceSummaries pandas dataFrame
    # scrit is the critical standard deviation to qualify as steady state
    # dx is the chunk size to use when calculating standard deviations
    xmeets = []
    ss = ssfull[ssfull['time']==t]
    xlist = xpts(ss)
    # get a list of x locations where the chunk standard deviations are below scrit
    for x in xlist[xlist<xlist[-1]-dx]:
        if stdevchunk(ss, x, x+dx, scrit):
            xmeets.append(x+dx/2)
    # if there are no chunks, return an error value
    if len(xmeets)<1:
        return -1
    else:
        # take the median chunk location that worked
        centerx = st.median(xmeets)
        s = 0
        dxf = dx
        # expand the chunk size until we go above scrit
        while stdevchunk(ss, centerx-dxf/2, centerx+dxf/2, scrit):
            dxf = dxf+dx/10
        dxf = dxf-dx/10
        dict1 = {'time':t, 'centerx':centerx, 'thickness':dxf, 'x0':centerx-t/2, 'xf':centerx+t/2}
        subset = ss[(ss['x']>=centerx-t/2) & (ss['x']<centerx+t/2)]
        subset = subset.drop(columns=['x', 'time'])
        dict2 = (subset.mean(axis=0)).to_dict()
        dict3 = ((subset.sem(axis=0)).add_suffix('_sterr')).to_dict()
        return {**dict1, **dict2, **dict3}
    
def findsteadytime(folder, scrit, dx):
    ppath = os.path.join(folder, 'steadyTime.csv')
    if os.path.exists(ppath):
        return
    ssfull = importSS(folder)
    if len(ssfull)==0:
        return -1
    t1 = []
    tlist = ssfull.time.unique()
    for t in tlist:
        out = findsteady(ssfull, t, scrit, dx)
        if type(out)==dict:
            t1.append(out)
    if len(t1)>0:
        p = pd.DataFrame(t1)
        p.to_csv(ppath, index_label=False)
        print('Exported ', ppath)
        return 
    else:
        return


def averagechunk(chunk):
    center = chunk.mean(axis=0)
    if not 'angle' in chunk.columns:
        # chunk.loc[:, 'angle'] = [0 for i in range(len(chunk))]
        chunk = chunk.assign(angle=pd.Series([0 for i in range(len(chunk))]).values)
        yc = center['y']
        zc = center['z']
        for i in chunk.index:
            chunk.at[i, 'angle'] = ptangle(chunk.loc[i], yc, zc)
    finalpts = pd.DataFrame()
    for a in range(360):
        pts = chunk[(chunk['angle']>=a)&(chunk['angle']<a+1)]
        if len(pts)>0:
            finalpts = finalpts.append(pts.mean(axis=0), ignore_index=True)
    return finalpts
    
def tracksteadyregion(folder, t):
    # t is a time in s that we use to define the steady state time
    ppath = os.path.join(folder, 'trackRegion_'+str(int(round(t*1000)))+'ms.csv')
    XSpath = os.path.join(folder, 'steadyXSpts_'+str(int(round(t*1000)))+'ms.csv')
    if os.path.exists(ppath) and os.path.exists(XSpath):
        return
    st1 = importST(folder)
    if len(st1)<1:
        return
    steady = st1[st1['time']==t] # get the row in steadyTime.csv that matches time t
    if len(steady)==0:
        return
    
    print(folder)
    x0 = steady.iloc[0]['x0']
    xf = steady.iloc[0]['xf']
    allpts = importCSV(folder)
    allptsreg = allpts[(allpts['x']>=x0)&(allpts['x']<=xf)]   
    fs = folderstats(folder)
    s = []
    for ti in np.sort(allptsreg.time.unique()):
        print('t=',ti)
        thesepts = allptsreg[allptsreg['time']==ti]
        chunkf = averagechunk(thesepts)
        if ti==t:
            chunkf.to_csv(XSpath)
            print('Exported ', XSpath)
        s.append(slicesummary(fs, chunkf))
    steadytimetrack = pd.DataFrame(s)
    steadytimetrack.to_csv(ppath)
    print('Exported ', ppath)
    return
        
def processfile(folder):
    summarize(folder)
    findsteadytime(folder, 0.01, 1)
    tracksteadyregion(folder, 2.5)
    exportslice(folder, 3, 2.5)

def importST(folder):
    file = os.path.join(folder, 'steadyTime.csv')
    return plainim(file, False)

def importXS(folder, t):
    # t is a time that the slices are taken at
    file = os.path.join(folder, 'steadyXSpts_'+str(int(round(1000*t)))+'ms.csv')
    return plainim(file, False)

def importTR(folder, t):
    # t is a time that the slices are taken at
    file = os.path.join(folder, 'trackRegion_'+str(int(round(1000*t)))+'ms.csv')
    return plainim(file, False)

# from an excel table containing all of the legend entries, transposed, get the list of ink viscosities, support viscosities, and surface tensions
    # t1 is a 
def tpfromxl(t1):
    ivlist = kintodyn(list(t1['nuink'].unique()))
    svlist = kintodyn(list(t1['nusup'].unique()))
    sigmalist = list(t1['sigma'].unique())
    ivlist.sort()
    svlist.sort()
    sigmalist.sort()
    return ivlist,svlist,sigmalist


###################################
## FOLDERPARSER ARCHIVE

##### if we run interFoam then foamToVTK on a folder that only contains the last time step, but we previously generated vtk files for previous time steps, the new .vtk.series file will not reference the old time steps we generated. However, when we run foamToVTK, we should be careful to always add the output onto the log_foamToVTK file so we know what timesteps correspond to which vtk files. By scraping that log, we can reconstruct the .vtk.series file with all of the time steps, not just the ones in the folder at the time
    

    
# redoVTKSeries goes into a folder and regenerates the .vtk.series file with all time steps included
    # folder is a full path name to either the case folder or its parent
def redoVTKSeries(folder:str) -> None:
    cf = caseFolder(folder) 
        # sometimes the case file is below the folder and named 'case', 
        # and sometimes it is the folder and named something like 'nb64'. 
        # This allows flexibility for those two situations
    cfbasename = os.path.basename(cf) # e.g. 'case' or 'nb64'
    vtklog = os.path.join(cf, 'log_foamToVTK')
    
    flist = [] # folder numbers
    tlist = [] # times
    if os.path.exists(vtklog):
        with open(vtklog, 'r') as f:
            line = f.readline()
            while len(line)>0:
                # the foamToVTK log will contain chunks that start with 
                # time: and include a line that contains Internal and the folder name
                while not line.startswith('Time:') and len(line)>0:
                    line = f.readline()
                if len(line)>0:
                    strs = re.split('Time: |\n', line)
                    time = strs[1]
                    while not 'Internal' in line:
                        line = f.readline()
                    strs = re.split(cfbasename+'_|/internal', line)
                    folderlabel = strs[1]
                    if not folderlabel in flist:
                        flist.append(folderlabel)
                        tlist.append(time)
                    line = f.readline()
        generateVTKSeries(tlist, flist, cf)
    return 
