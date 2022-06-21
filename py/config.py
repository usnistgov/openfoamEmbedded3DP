#!/usr/bin/env python
'''Tools for loading settings'''

# external packages
import yaml
import sys
import os
try:
    from box import Box
except ModuleNotFoundError:
    nobox = True
else:
    nobox = False
import shutil


#----------------------------------------------------

class Struct:
    def __init__(self, **entries):
        for key,val in entries.items():
            if type(val) is dict:
                setattr(self, key, Struct(**val))
            else:
                setattr(self, key, val)


def getConfigDir() -> str:
    '''find the configs directory'''
    currentdir = os.path.dirname(os.path.realpath(__file__))
    configdir = os.path.join(currentdir, 'configs')
    if not os.path.exists(configdir):
        parentdir = os.path.dirname(currentdir)
        configdir = os.path.join(parentdir, 'configs')
        if not os.path.exists(configdir):
            raise FileNotFoundError(f"No configs directory found")
    return configdir

def dumpConfigs(cfg, path:str) -> int:
    '''Saves config file. cfg could be a Box or a dict'''
    with open(path, "w") as ymlout:
        if type(cfg) is Box:
            cout = cfg.to_dict()
        elif type(cfg) is dict:
            cout = cfg
        else:
            return 1
        yaml.safe_dump(cout, ymlout)
        return 0
    
    
def findConfigFile() -> str:
    '''find the config file and return the path'''
    configdir = getConfigDir()
    path = os.path.join(configdir,"config.yml")
    if os.path.exists(path):
        return path
    
    # config.yml does not exist: find default or template
    path2 = os.path.join(configdir, 'config_default.yml')
    if os.path.exists(path2):
        shutil.copy2(path2, path)
        return path
    else:
        path2 = os.path.join(configdir, 'config_template.yml')
        if os.path.exists(path2):
            shutil.copy2(path2, path)
            return path
            
    # no config_default or config_template either: find any file that looks like a config file
    llist = os.listdir(configdir)
    while not os.path.exists(path):
        l = llist.pop(0)
        if (l.endswith('yml') or l.endswith('yaml')) and ('config' in l):
            path = os.path.join(configdir, l)
            return path
        
    return path
    
        
def loadConfigFile(path:str) -> Box:
    '''open the config file and turn it into a Box'''
    with open(path, "r") as ymlfile:
        y = yaml.safe_load(ymlfile)
        if nobox:
            cfg = Struct(**y)
        else:
            cfg = Box(y)
        return cfg

def loadConfig() -> Box:
    path = findConfigFile()
    cfg = loadConfigFile(path)
    return cfg
        
#----------------------------------------------------

cfg = loadConfig()
    
    
