#!/usr/bin/env python
'''Script for combining images into videos. Loops continuously every 6 hours.'''

# external packages
import os, sys
import time
import logging

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(currentdir)
sys.path.append(parentdir)
import folderparser as fp
from videofuncs import saveVid
from config import cfg

# logging
LOGGERDEFINED = fp.openLog(os.path.realpath(__file__), False)

# info
__author__ = "Leanne Friedrich"
__copyright__ = "This data is publicly available according to the NIST statements of copyright, fair use and licensing; see https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software"
__credits__ = ["Leanne Friedrich"]
__license__ = "NIST"
__version__ = "1.0.0"
__maintainer__ = "Leanne Friedrich"
__email__ = "Leanne.Friedrich@nist.gov"
__status__ = "Production"

            
#-----------------------------------------------            
      
            
if __name__ == "__main__":
    logging.info('Compiling frames into videos...')
    
#     serverfolder = os.path.join(cfg.path.server, 'yieldingsweep', 'HBHByielded')
#     topfolders = [os.path.join(serverfolder,s) for s in ['k', 'n', 'tau0']]

    serverfolder = os.path.join(cfg.path.server, 'viscositysweep')
    topfolders = [os.path.join(serverfolder,s) for s in ['HBHBsweep', 'HBnewtsweep', 'newtHBsweep', 'newtnewtsweep']]
    
    while True:
        try:
            for topfolder in topfolders:
                for folder in fp.caseFolders(topfolder):
                    logging.info(f'Checking {fp.shortName(folder)}')
                    for p in ['umag']:
                        for s in ['y']:
                            try:
                                saveVid(folder, s, p)      
                            except:
                                pass
        except:
            pass
        
        if len(sys.argv)>1:
            waittime = float(sys.argv[1])
        else:
            waittime = 0
        if waittime>0:
            logging.info(f'waiting {waittime} hrs for more files...')
            fp.printCurrentTime()
            time.sleep(60*60*waittime)
        else:
            break