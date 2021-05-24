#!/usr/bin/env python
'''Functions for moving folders between computers, servers, for OpenFOAM simulations of embedded 3D printing of single filaments. '''

import os
import re
import shutil 
from typing import List, Dict, Tuple, Union, Any, TextIO
import time
import logging, platform, socket, sys

from folderscraper import populate
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from folderparser import *

__author__ = "Leanne Friedrich"
__copyright__ = "This data is publicly available according to the NIST statements of copyright, fair use and licensing; see https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software"
__credits__ = ["Leanne Friedrich"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Leanne Friedrich"
__email__ = "Leanne.Friedrich@nist.gov"
__status__ = "Production"

#-------------------------------------------------------------------------------------------------  


#########################################################################
########## COPY FUNCTIONS ###############################################
#########################################################################


def copyFolder(fromFolder:str, targetFolder:str, f:str, st:str) -> None:
    '''fromFolder is the parent directory to copy from
    targetFolder is the parent directory to copy to
    f is the folder within fromFolder to copy
    st is folder within f to copy'''
    ipfrom = os.path.join(fromFolder, f, st)
    if not os.path.exists(ipfrom):
        return
    ipfolder = os.path.join(targetFolder, f, st)
    if not os.path.exists(ipfolder):
        mkdirif(ipfolder)
    for file in os.listdir(ipfrom):
        targetFile = os.path.join(ipfolder, file)
        fromFile = os.path.join(ipfrom, file)
        if os.path.isdir(fromFile):
            copyFilesInFolder(fromFile, targetFile)
        else:
            if not os.path.exists(targetFile):
                shutil.copyfile(fromFile, targetFile)
            else:
                sizediff = os.stat(targetFile).st_size - os.stat(fromFile).st_size
                if not sizediff==0:
                    shutil.copyfile(fromFile, targetFile)
    return
                

def copyFilesInFolder(fromFolder:str, targetFolder:str) -> None:
    '''directly copy the files in one folder to another'''
    mkdirif(targetFolder)
    printOut = False
    checkSizeDiff = False
    for file in os.listdir(fromFolder):
        if 'geometry' not in file and not 'mesh'==file:
            fromFile = os.path.join(fromFolder, file)
            targetFile = os.path.join(targetFolder, file)
            # if this is a directory, make the directory in the 
            # target folder if it doesn't exist, and copy all files
            if os.path.isdir(fromFile):
                if not os.path.exists(targetFile):
                    mkdirif(targetFile)
                copyFilesInFolder(fromFile, targetFile)
            else:
                # if this is a file, copy the file if it doesn't exist
                if not os.path.exists(targetFile):
                    shutil.copyfile(fromFile, targetFile)
                elif 'legend.csv' in fromFile:
                    copyLegend(fromFolder, targetFolder) 
                elif '.series' in fromFile:
                    copySeries(fromFile, targetFile) 
                elif checkSizeDiff:
                    # if the file already exists, determine if the size changed
                    sizediff = os.stat(targetFile).st_size - os.stat(fromFile).st_size
                    if not sizediff==0:
                        shutil.copyfile(fromFile, targetFile)
                        printOut = True                       
    if printOut:
        print([fromFolder, targetFolder])
    return

#--------------------------------------------------
#--------------------------------------------------
#### COPY TYPES OF FILES


def toDoList(topfolder:str) -> List[str]:
    '''Simulations still to run'''
    with open (os.path.join(topfolder, 'runallfiles.sh'), "r") as f:
        line = f.readline()
        while not line.startswith('dirlist'):
            line = f.readline()
        lspl = re.split('\'| ', line)[1:-1]
    lspl = [x for x in lspl if len(x)>0]
    lspl.sort()    
    return lspl

def doneFolder(topfolder:str, tfinal:float, loopTime:float=0) -> bool:
    '''Go through files in the folder and determine which are done, in progress, and aborted, and abort any simulations that are too slow. Returns true if all simulations are done, false if not done.'''
    abort = False
    while not abort:
        done = []
        aborted = []
        aborted1 = []
        progress = []
        notstarted = []
        tdl = toDoList(topfolder)
        for f in caseFolders(topfolder):
            ct = currentTime(f)
            f0 = os.path.basename(f)
            if ct==0:
                if os.path.exists(os.path.join(f, 'log_interFoam')) or os.path.exists(os.path.join(f, 'case', 'log_interFoam')):
                    if f0 in tdl:
                        progress.append(f0)
                    else:
                        aborted.append(f0)
                else:
                    notstarted.append(f0)
            elif ct>=tfinal:
                done.append(f0)
            else:
                if f0 in tdl:
                    progress.append(f0)
                else:
                    if ct>=1:
                        aborted1.append(f0)
                    else:
                        aborted.append(f0)
        plog = []

        print('Populating legends')

        for p in progress:
            o = 0
            endtime = 2.5
            abort = '     '
            f = os.path.join(topfolder, p)
            t = populate(f)
            t = importIf(f)
            runtime = t[5][1]
            ct = t[6][1]
            rate = t[7][1]
            if len(rate)>0:
                try:
                    rate = float(rate)
                    ct = float(ct) # ct = current time in simulation
                    runtime = float(runtime) # actual time in real life (hr)
                except:
                    pass
                else:
                    if rate>120 and runtime>1:
                        if ct<0.8:
                            #logging.info('Modifying controldict to 0 for '+f)
                            endtime = 0.0
                            abort = 'ABORT'
                        else:
                            #logging.info('Modifying controldict to 1 for '+f)
                            endtime = 1.0
                    elif rate>60:
                        #logging.info('Modifying controldict to 1 for '+f)
                        endtime = 1.0
                        if 2.1>=ct>=1:
                            abort = 'ABORT'
                    elif rate<60:
                        endtime = 2.5
                    o = modifyControlDict(f, endtime)
            if o==1 and endtime<=ct:
                # if the controlDict was already modified
                if ct>=1:
                    aborted1.append(os.path.basename(f))
                else:
                    aborted.append(os.path.basename(f))
            else:
                if type(ct) is float:
                    ct = '{:.1f}'.format(ct)
                if type(rate) is float:
                    rate = '{:.2f}'.format(rate)
                plog.append('{!s}\t{!s}\t{!s}\t\t{:1.1f}\t{!s}'.format(p, ct, rate, endtime, abort))

        aborted.sort()
        logging.info('###### done: '+', '.join(done))
        logging.info('###### not started: '+', '.join(notstarted))
        logging.info('###### aborted at 1: '+', '.join(aborted1))
        logging.info('###### aborted before 1: '+', '.join(aborted))
        logging.info('###### in progress:')
        logging.info('folder\ttime\trate (hr/s)\tendtime(s)')
        for p in plog:
            logging.info(p)
        logging.info('')

        if len(notstarted)==0 and len(plog)==0:
            finished=True
            abort=True
        else:
            finished=False
            if loopTime>0:
                abort=False
                printCurrentTime()
                logging.info(f'Waiting {loopTime} hours\n\n\n')
                time.sleep(60*60*loopTime)
            else:
                abort=True
    return finished
    


#---------
# legend

def copyAndPrint(fromfile:str, tofile:str) -> None:
    '''Copy the file and print a status update.'''
    shutil.copyfile(fromfile, tofile)
    logging.info('Copied '+ tofile)
    return

def copyLegend(fromFolder:str, toFolder:str) -> None:
    '''Copy most up to date legend from one folder to the other'''
    try:
        fromL = os.path.join(fromFolder, 'legend.csv')
        toL = os.path.join(toFolder, 'legend.csv')
        forward = True
        backward = True
        if not os.path.exists(fromL):
            # if there's no original file, return
            return
        if not os.path.exists(toL):
            # if there's no existing file, copy and return
            return copyAndPrint(fromL, toL)
    
        # if there are legends in both folders
        fromim = importIf(fromFolder)
        fromrate = fromim[7][1]
        fromtime = fromim[6][1]
        fromruntime = fromim[4][1]
        toim = importIf(toFolder)
        torate = toim[7][1]
        totime = toim[6][1]
        toruntime = toim[4][1]
        if not fromrate.isnumeric:
            # if the from legend doesn't have a valid rate, don't copy forward
            forward = False
        if not torate.isnumeric:
            # if the to legend doesn't have a valid rate, don't copy backward
            backward = False
        if not fromtime.isnumeric:
            # if the from legend doesn't have a valid time, don't copy forward
            forward = False 
        if not totime.isnumeric:
            # if the new legend doesn't have a valid time, don't copy backward
            backward = False 
        if float(fromtime)>float(totime) or float(fromruntime)>float(toruntime):
            # if the from legend is more up to date, don't copy backward
            backward = False
        elif float(totime)>float(fromtime) or float(toruntime)>float(fromruntime):
            # if the to legend is more up to date, don't copy forward
            forward = False
        else:
            # if the legends have the same times, don't copy either
            backward = False
            forward = False
            
        # if we've flagged this to copy the legend in both directions
        # or to copy neither legend, return
        if forward==backward:
            return 
        if forward:
            return copyAndPrint(fromL, toL)
        elif backward:
            return copyAndPrint(toL, fromL)
    except Exception as e:
        return
    else:
        return
            
#----------------------------
# logs


def moveLog(folder:str, cf:str, s:str, moveSource:bool) -> str:
    '''move the log files to the parent folder (folder), out of the case folder (cf)
    s is the name of the log file
    returns the name of the log
    moveSource true to move around files in folder, false to leave it alone'''
    # we only want to keep the version in the parent folder, not the case folder
    filefromcf = os.path.join(cf, s)
    filefrom = os.path.join(folder, s)

    # if the log doesn't exist anywhere, return
    if not os.path.exists(filefrom) and not os.path.exists(filefromcf):
        return ''

    # if there is only a version in the case folder, move it to the main folder
    if os.path.exists(filefromcf) and not os.path.exists(filefrom):
        if moveSource:
            os.rename(filefromcf, filefrom)
            return filefrom
        else:
            return filefromcf

    # if there are two versions of the log, 
    # keep the larger one and move it to the parent folder
    if os.path.exists(filefrom) and os.path.exists(filefromcf):
        sizecf = os.stat(filefromcf).st_size
        size = os.stat(filefrom).st_size

        # if the one in the main folder is bigger, delete the one in the case folder
        if size>=sizecf:
            if moveSource:
                os.remove(filefromcf)
                return filefrom
            else:
                return filefrom
        else:
            # of the version in the case folder is bigger,
            # move it to the main folder and overwrite the one in the main folder
            if moveSource:
                os.remove(filefrom)
                os.rename(filefromcf, filefrom)
            else:
                return filefromcf
            
    return filefrom


def moveLogs(folder:str, moveSource:bool) -> List[str]:
    '''move the log files to the parent folder, out of the case folder
        returns the names of the logs
        moveSource true to move around files in fromFolder, false to leave it alone'''
    cf = caseFolder(folder)
    files = []
    for s in ['log_foamToVTK', 'log_interFoam']:
        files.append(moveLog(folder, cf, s, moveSource))
    return files

def moveLogsTop(topFolder:str) -> None:
    '''move all of the logs from the case folder to the parent folder for all simulations in topFolder'''
    for f in os.listdir(topFolder):
        folder = os.path.join(topFolder, f)
        if isSimFolder(folder):
            moveLogs(folder, True)
    return
            

def copyLogs(fromFolder:str, toFolder:str, f:str, moveSource:bool) -> None:
    '''copy just the logs for one folder
        f is subfolder name
        moveSource true to move around files in fromFolder, false to leave it alone'''
    fromf = os.path.join(fromFolder, f)
    if not os.path.isdir(fromf):
        return
    logs = moveLogs(fromf, moveSource)
    for log in logs:
        if os.path.exists(log):
            s = os.path.basename(log)
            fileto = os.path.join(toFolder, f, s)
            if os.path.exists(fileto):
                sizefrom = os.stat(log).st_size
                sizeto = os.stat(fileto).st_size
                if sizeto>sizefrom:
                    shutil.copyfile(fileto, log)
                    logging.info('Copied back '+log)
                    return
                elif sizeto==sizefrom:
                    return
            # if there is no existing file or it's smaller, copy
            shutil.copyfile(log, fileto)
            logging.info('Copied '+fileto)


def copyLogsTop(fromFolder:str, toFolder:str, moveSource:bool) -> None:
    '''copy just the logs for many folders
        moveSource true to move around files in fromFolder, false to leave it alone'''
    for f in os.listdir(fromFolder):
        if isSimFolder(os.path.join(fromFolder, f)):
            copyLogs(fromFolder, toFolder, f, moveSource)
    return


##############

def vtkCheck(file):
    '''check if the .vtk.series file is up to date, and redo it if it's not'''
    folder = os.path.dirname(os.path.dirname(file))
    times = parseVTKSeries(folder)
    files = vtkfiles(folder)
    if len(times)<files:
        redoVTKSeriesNoLog(folder)
        times = parseVTKSeries(folder)
    return files

def copySeries(fromFile:str, toFile:str) -> None:
    '''update and copy vtk.series file from one place to another'''
    fromFiles = vtkCheck(fromFile)
    toFiles = vtkCheck(toFile)
    if fromFiles>toFiles:
        shutil.copyfile(fromFile, toFile)
        print('Copied '+toFile)
    elif toFiles>fromFiles:
        shutil.copyfile(toFile, fromFile)
        print('Copied back '+fromFile)

#############################

def copySubFolder(fromFolder:str, toFolder:str, s:str) -> None:
    '''copy subfolder within fromFolder to subfolder in toFolder'''
    cf = caseFolder(fromFolder) # folder that the case files are in
    vtkfolder = os.path.join(cf, s)
    if not os.path.exists(vtkfolder):
        vtkfolder = os.path.join(fromFolder, s)
    if os.path.exists(vtkfolder):
        copyFilesInFolder(vtkfolder, os.path.join(toFolder, s))  

            

        
########## COPY WHOLE FOLDERS



def copyEtoServer(EFolder:str, serverFolder:str) -> None:
    '''copy the legend, interface points, vtk, and images
    this is most useful for moving files from E to server'''
    logging.info(' ---------- Copying files from E to server')
    for f in os.listdir(EFolder):
        fromf = os.path.join(EFolder, f)
        tof = os.path.join(serverFolder, f)
        mkdirif(tof)
        copyLegend(fromf, tof)
        for s in ['images', 'interfacePoints']: # copy images, points, VTK files
            copySubFolder(fromf, tof, s)
        for f1 in os.listdir(fromf):
            if f1.endswith('.csv') and not 'legend' in f1:
                # this is for slice summaries, other csvs created by interfaceMetrics
                shutil.copyfile(os.path.join(fromf, f1), os.path.join(tof, f1))
    logging.info(' ---------- Done copying files from E to server')
    return


def copyToE(fromFolder:str, Efolder:str, fromServer:bool) -> None:
    '''copy folders to E. fromServer is true if we're copying from the server'''
    if fromServer:
        logging.info(' ---------- Copying files from server to E')
    else:
        logging.info(' ---------- Copying files from C to E')
    for f in os.listdir(fromFolder):
        fromf = os.path.join(fromFolder, f)
        if os.path.isdir(fromf):            
            tof = os.path.join(Efolder, f) # to folder
            mkdirif(tof) # create E folder
            for s in ['images', 'interfacePoints', 'VTK']: # copy images, points, VTK files
                copySubFolder(fromf, tof, s)
            copyLogs(fromFolder, Efolder, f, fromServer) 
                # copy log files, only overwrite if from server
            copyLegend(fromf, tof)                  
                # copy legend
    logging.info(' ---------- Done copying files to E')
    

def moveCtoServer(CFolder:str, serverFolder:str, Efolder:str='') -> None:
    '''move files from C to E and server. If you only want C to server, leave Efolder empty'''
    if not os.path.exists(CFolder) and os.path.exists(Efolder):
        copyToE(serverFolder, Efolder, True)
        logging.info('------ all files copied in '+os.path.basename(Efolder))
    else:
        s = os.path.basename(CFolder)
        try:
            print(CFolder, serverFolder)
            doneFolder(CFolder, 2.5)
            if os.path.exists(Efolder):
                copyToE(CFolder, Efolder, False)
                copyToE(serverFolder, Efolder, True)
            logging.info(' ---------- Copying files from C to server')
            copyFilesInFolder(CFolder, serverFolder)
            logging.info(' ---------- Done copying files to server, moving logs')
            moveLogsTop(serverFolder)
            logging.info(' ---------- Done moving logs')
        except Exception as e:
            logging.error('----- ERROR in '+s+', '+str(e))
        else:
            logging.info('------ all files copied in '+s)
        


def copyCtoServer(CFolder:str, serverFolder:str, timex:float, Efolder:str='') -> None:
    '''this loops infinitely, copying files to the server
        time is the time in hours between loops'''
    if timex>1:
        waittime = str(timex)+' hours'
    else:
        waittime = str(timex)+' hour'
    while True:
        moveCtoServer(CFolder, serverFolder, Efolder)
           
        if timex==0:
            break
        else:
            printCurrentTime()
            logging.info('Waiting '+waittime+'\n\n\n')
            time.sleep(60*60*timex)
            
        
def copyCtoServerFolders(SERVERFOLDER:str, CFOLDER:str, slist:List[str], timex:float, EFOLDER:str='') -> None:
    ''' loop indefinitely, for multiple folders. Leave EFOLDER empty to only copy from cfolder to serverfolder. slist is a list of subfolders, e.g. ['HBHBsweep']    '''
    if timex>1:
        waittime = str(timex)+' hours'
    else:
        waittime = str(timex)+' hour'
    while True:
        for s in slist:
            serverfolder = os.path.join(SERVERFOLDER,s)
            cfolder = os.path.join(CFOLDER, s)
            Ef = os.path.join(EFOLDER,s)
            if os.path.exists(cfolder):
                moveCtoServer(cfolder,serverfolder, Efolder=Ef)
            
        if timex==0:
            break
        else:
            printCurrentTime()
            logging.info('Waiting '+waittime+'\n\n\n')
            time.sleep(60*60*timex)


    
####################################

# ---------------------------------
# ---------------------------------
# ---------------------------------
# ---------------------------------
# ---------------------------------

### ARCHIVE #####



# def consolidateImages(folder:str) -> None:
#     imf = os.path.join(folder, 'images')
#     if not os.path.exists(imf):
#         return
#     for fo in os.listdir(imf):
#         fofull = os.path.join(imf, fo)
#         if os.path.isdir(fofull):
#             for f in os.listdir(fofull):
#                 fold = os.path.join(fofull, f)
#                 fnew = fold.replace(r'\videos', '')
#                 fnew = fnew.replace(r'\frames', '')
#                 if os.path.exists(fnew):
#                     os.remove(fold)
#                 else:
#                     os.rename(fold, fnew)
#             os.rmdir(fofull)