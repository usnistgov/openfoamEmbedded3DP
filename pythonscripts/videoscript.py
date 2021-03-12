import os
import glob
from IPython.display import Image
import imageio
import folderparser as fp
import subprocess, json
import re
import time
from config import cfg


def get_length(filename):
    rl = ["ffprobe", "-v", "error", "-show_entries", "format=duration", "-of", "default=noprint_wrappers=1:nokey=1", filename]
    result = subprocess.run(rl, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    duration =float(result.stdout)
    return duration

def get_time(filename):
    return float(re.split('t|_', os.path.basename(filename))[1])/10

def saveVid(folder, s:str, p:str):
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