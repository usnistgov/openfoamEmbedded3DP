#!/usr/bin/env python
'''Script for combining images into videos. Loops continuously every 6 hours.'''

import os
import glob
from IPython.display import Image
import imageio
import subprocess, json
import re
import time

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import folderparser as fp
from config import cfg

__author__ = "Leanne Friedrich"
__copyright__ = "This data is publicly available according to the NIST statements of copyright, fair use and licensing; see https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software"
__credits__ = ["Leanne Friedrich"]
__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Leanne Friedrich"
__email__ = "Leanne.Friedrich@nist.gov"
__status__ = "Production"

#----------------------------------------------


def get_length(filename:str) -> float:
    '''Get the length of the video'''
    rl = ["ffprobe", "-v", "error", "-show_entries", "format=duration", "-of", "default=noprint_wrappers=1:nokey=1", filename]
    result = subprocess.run(rl, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    duration =float(result.stdout)
    return duration

def get_time(filename:str) -> float:
    '''Get the time of the simulation timepoint'''
    return float(re.split('t|_', os.path.basename(filename))[1])/10

def saveVid(folder:str, s:str, p:str) -> None:
    '''compile images into a video, where images contain the string {s}_{p}.png'''
    zfiles = glob.glob(os.path.join(folder,'images', '*'+s+'_'+p+'.png'))
    if len(zfiles)>0:
        fn = os.path.join(folder,'images',  s+'_'+p+'.mp4')
        if os.path.exists(fn):
            #print(fn)
            fnlength = get_length(fn)
            times = [get_time(f) for f in zfiles]
            #print(fn, fnlength, max(times))
            if fnlength < max(times)-0.1:
                print('remove', fn)
                os.remove(fn)
        if not os.path.exists(fn):
            print('write', fn)
            mp4writer = imageio.get_writer(fn, fps=10)
            for filename in zfiles:
                image = imageio.imread(filename)
                mp4writer.append_data(image)
            mp4writer.close()

            
#-----------------------------------------------            
            
            
if __name__ == "__main__":
    serverfolder = os.path.join(cfg.server, 'yieldingsweep', 'HBHByielded')
    topfolders = [os.path.join(serverfolder,s) for s in ['k', 'n', 'tau0']]
    while True:
        try:
            for topfolder in topfolders:
                for folder in fp.caseFolders(topfolder):
                    for p in ['umag']:
                        for s in ['y']:
                            try:
                                saveVid(folder, s, p)      
                            except:
                                pass
        except:
            pass
        print('waiting 6 hrs for more files...')
        fp.printCurrentTime()
        time.sleep(60*60*6)