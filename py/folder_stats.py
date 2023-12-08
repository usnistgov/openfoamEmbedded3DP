#!/usr/bin/env python
'''Collecting metadata about files'''

# external packages
import sys
import os
import csv
import numpy as np
import re
from typing import List, Dict, Tuple, Union, Any
import logging
import traceback

# local packages
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import tools.val_tools as vt
import tools.strings as st
import file.file_handling as fh

# logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)



#################################################################

class fluidStats:
    '''holds stats about a single fluid'''
    
    def __init__(self, fluid:str):
        self.fluid = fluid          # name of the fluid: ink or bath
        self.kinematic = {}         # values in OpenFOAM's preferred units, in m and kinematic units
        self.kinematicU = {}        # units
        self.dynamic = {}           # values in practical units, in mm and dynamic units
        self.dynamicU = {}          # units
        self.convertLater = []
        self.transportModel = ''
        
        
    def addVal(self, key:str, val:Any, units:str):
        '''add the value to the dictionary'''
        if val=='':
            return
        varname = key.replace(f'{self.fluid}_', '') 
        if varname=='velocity':
            self.addSpeed(val, units)
        elif varname=='rho':
            self.addDensity(val, units)
        elif varname=='transportModel':
            self.transportModel = val
        elif varname=='n':
            self.addScaledUnit('n', val, '', '', 1, 1)
        elif varname in ['nu', 'nu0', 'tau0', 'k', 'n']:
            self.addKinUnit(varname, val, units)
            
    def metaRow(self) -> Tuple[dict,dict]:
        '''get a dict describing the useful metadata about this fluid'''
        d = {f'{self.fluid}_transportModel':self.transportModel}
        u = {f'{self.fluid}_transportModel':''}
        self.localVisc()  # calculate the reynolds number, etc.
        for key,val in self.dynamic.items():
            d[f'{self.fluid}_{key}'] = val
            u[f'{self.fluid}_{key}'] = self.dynamicU[key]
        return d,u
            
    def addScaledUnit(self, name:str, val:Any, munit:str, mmunit:str, kscale:float, dscale:float) -> None:
        '''add a unit that's in meters to the variable dictionaries, and convert to mm'''
        if val=='':
            return
        val = float(val)
        self.kinematic[name] = kscale*val
        self.kinematicU[name] = munit
        self.dynamic[name] = dscale*val
        self.dynamicU[name] = mmunit

    def addSpeed(self, val:Any, units:str):
        '''convert speeds to mm/s'''
        if units=='m/s':
            self.addScaledUnit('v', val, 'm/s', 'mm/s', 1, 1000)
        elif units=='mm/s':
            self.addScaledUnit('v', val, 'm/s', 'mm/s', 1/1000, 1)
        else:
            raise ValueError(f'Unexpected speed units {units}')
                
    def addDensity(self, val:Any, units:str):
        '''convert densities to g/mL'''
        if units=='kg/m^3':
            self.addScaledUnit('rho', val, 'kg/m^3', 'g/mL', 1, 1/1000)
        elif units=='g/mL':
            self.addScaledUnit('rho', val, 'kg/m^3', 'g/mL', 1000, 1)

    def addKinUnit(self, name:str, val:Any, units) -> None:
        '''add kinematic rheology units'''
        if val=='':
            return
        self.kinematic[name] = float(val)
        self.kinematicU[name] = units
        if not 'rho' in self.kinematic:
            self.convertLater.append(name)
        else:
            self.convertKin(name)
            
    def convertKin(self, name:str) -> None:
        '''convert kinematic units to dynamic units'''
        val = self.kinematic[name]
        units = self.kinematicU[name] 
        self.dynamic[name] = self.kinToDyn(val)
        self.dynamicU[name] = {'m^2/s':'Pa*s', 'm^2/s^2':'Pa', 'm^2*s^(n-2)':'Pa*s^n'}[units]
        
    def convertAllKin(self) -> None:
        '''now that we have a density, convert all kinematic units to dynamic units'''
        for name in self.convertLater:
            self.convertKin(name)

    def kinToDyn(self, val:float) -> float:
        '''converts kinematic viscosity to dynamic viscosity
        val is a float in units of m and s.
        multiply by density in kg/m^3'''
        val = float(val)
        density = self.kinematic['rho']
        v2 = val*density
        if abs(round(v2)-v2)<10**-20:
            v2 = int(v2)            
        return v2
    
    def localVisc(self) -> dict:
        '''calculate local viscosity, reynolds number'''
        if hasattr(self, 'localViscDict'):
            return self.localViscDict
        out = {}
        tm = self.transportModel

        if not 'rho' in self.dynamic:
            return {}
        rho = self.dynamic['rho']     # g/mL
        out[f'{self.fluid}_rho'] = rho
        d = self.dynamic['d']   # mm
        v = self.dynamic['v']  # mm/s
        

        if tm=='Newtonian':
            out[f'{self.fluid}_visc'] = self.dynamic['nu'] # Pa*s
            self.dynamic['visc0'] = self.dynamic['nu']
            self.dynamicU['visc0'] = self.dynamicU['nu']
        else:
            # Herschel-Bulkley
            nu0 = self.dynamic['nu0']   # dynamic plateau viscosity
            tau0 = self.dynamic['tau0'] # yield stress
            k = self.dynamic['k']       # consistency index
            n = self.dynamic['n']           # power law index

            gdot = v/d # shear rate in 1/s
                                                                                                
            out[f'{self.fluid}_gdot'] = gdot
            self.dynamic['gdot'] = gdot
            self.dynamicU['gdot'] = '1/s'
            if gdot>0:
                visc0 = min([nu0, tau0/gdot+k*gdot**(n-1)]) # HB equation
            else:
                visc0 = nu0
            self.dynamic['visc0'] = visc0
            self.dynamicU['visc0'] = self.dynamicU['nu0']
            out[f'{self.fluid}_visc'] = visc0

        # reynolds number
        Re = (self.kinematic['rho']*self.kinematic['v']*self.kinematic['d'])/out[f'{self.fluid}_visc']
        self.dynamic['Re'] = Re
        self.dynamicU['Re'] = ''
        out[f'{self.fluid}_Re'] = Re

        self.localViscDict = out
        return out
        
class geoStats:
    '''for holding the nozzle geometry'''
    
    def __init__(self, folder:str):
        self.folder = folder
        self.vn = st.varNicknames()
        self.u = {}  # units
        self.orig = {}
        
    def metaRow(self) -> Tuple[dict, dict]:
        '''get a dictionary of the values that are useful to describe the inputs into this system'''
        d = {}
        u = {}
        for s in ['niw', 'nozzle_angle'
                  , 'adjacent_filament_orientation', 'adjacent_filament_offset'
                  , 'spacing', 'dEst']:
            if hasattr(self, s):
                d[s] = getattr(self, s)
                u[s] = self.u[s]
        return d,u
        
    def addVal(self, name:str, val:Any, units:str) -> None:
        val = vt.tryfloat(val)
        if units=='m':
            val = val*1000
            units = 'mm'
        setattr(self, name, val)
        self.orig[name] = val    # keep a list of original values
        self.u[name] = units
        shortname = self.vn.shorten(name)
        if shortname==name:
            return
        setattr(self, shortname, val)
        self.u[shortname] = units
        
    def complete(self, inkv:float, bathv:float):
        self.addCoords(inkv, bathv)
        self.addSpacing(inkv)
        
    def addDefaultUnits(self, s:str) -> None:
        '''add this variable to the units list, using the default length units (mm)'''
        self.u[s] = self.u['niw']

    def addCoords(self, inkv:float, bathv:float):
        '''calculate extra coordinates'''
        self.now = self.niw + self.nt*2
        self.nre = self.ncx + self.niw/2 + self.nt # nozzle right edge
        self.behind = self.nre + self.niw # 1 inner diameter behind nozzle
        if inkv>0:
            self.dEst = self.niw*np.sqrt(inkv/bathv)
        else:
            self.dEst = self.niw  # diameter of a single filament
        self.intentztop = self.nbz
        self.intentzbot = self.nbz - self.dEst
        self.intentzcenter = (self.intentztop+self.intentzbot)/2
        self.intentyleft = self.ncy - self.dEst/2
        self.intentyright = self.ncy + self.dEst/2
        self.intentycenter = (self.intentyleft+ self.intentyright)/2
        self.intentarea = np.pi*(self.niw/2)**2
        
        for s in ['now', 'nre', 'behind', 'intentzcenter', 'intentzbot', 'intentyleft', 'intentyright', 'intentycenter', 'dEst']:
            self.addDefaultUnits(s)
        
    def addSpacing(self, inkv:float):
        '''calculate relative spacing between filaments'''
        if hasattr(self, 'adjacent_filament_offset'):
            self.spacing = self.adjacent_filament_offset/self.dEst
            self.u['spacing'] = 'nid'
            slist = [0.5, 0.625, 0.75, 0.875, 1, 1.25]
            for s in slist:
                if abs(self.spacing-s)<0.003:
                    self.spacing = s
            if self.adjacent_filament_orientation=='z':
                self.intentzbot = self.intentzbot-self.adjacent_filament_offset
                if inkv==0:
                    self.intentztop = self.intentztop-self.adjacent_filament_offset
                self.intentzcenter = (self.intentzbot+self.intentztop)/2
            elif self.adjacent_filament_orientation=='y':
                self.intentyleft = self.intentyleft - self.adjacent_filament_offset
                if inkv==0:
                    self.intentyright = self.intentyright - self.adjacent_filament_offset
                self.intentycenter = (self.intentyleft + self.intentyright)/2
            self.intentw = self.intentyright-self.intentyleft
            self.intenth = self.intentztop-self.intentzbot
            if inkv>0:
                self.intentarea = self.niw*self.adjacent_filament_offset + self.intentarea  # rectangle between + 2 semicircles
    


class folderStats:
    '''the folderStats class is used to store critical information about the nozzle geometry'''
    
    def __init__(self, folder:str):
        self.vn = st.varNicknames()
        self.folder = folder
        self.bn = os.path.basename(self.folder)
        self.ink = fluidStats('ink')
        self.sup = fluidStats('sup')
        self.geo = geoStats(self.folder)
        self.times = {}
        self.u = {}
        self.readLegend() 
        self.ink.convertAllKin()
        self.sup.convertAllKin()
        self.geo.complete(self.ink.dynamic['v'], self.sup.dynamic['v'])
        self.addDiameters()
        self.fh1 = fh.folderHandler(self.folder)
        self.xlist = []
        
    def addDiameters(self):
        if self.geo.u['niw']=='mm':
            self.ink.addScaledUnit('d', self.geo.niw, 'm', 'mm', 1/1000, 1)
        elif self.geo.u['niw']=='m':
            self.ink.addScaledUnit('d', self.geo.niw, 'm', 'mm', 1, 1000)
        if self.geo.u['now']=='mm':
            self.sup.addScaledUnit('d', self.geo.now, 'm', 'mm', 1/1000, 1)
        elif self.geo.u['now']=='m':
            self.sup.addScaledUnit('d', self.geo.now, 'm', 'mm', 1, 1000)
                                   
    def addTiming(self, name:str, val:Any, units:str) -> None:
        '''add simulation time data'''
        if name=='compare_to':
            self.compare_to = val.strip()
        elif name=='folder' or name=='openfoam_version':
            return
        else:
            setattr(self, name, val)
            self.u[name] = units
        
    def addSigma(self, val:Any, units:str) -> None:
        '''add surface tension data'''
        val = float(val)
        if units=='J/m^2':
            self.sigma = val*1000
            self.u['sigma'] = 'mJ/m^2'
        elif units=='mJ/m^2':
            self.sigma = val
            self.u['sigma'] = units
        else:
            raise ValueError(f'Unexpected sigma units {units}')
            
    def readLegendRow(self, row:list) -> None:
        '''read a single row from the legend.csv file'''
        row[0] = row[0].replace(' ', '_')
        if row[0] in ['folder', 'openfoam_version', '']:
            return
        elif row[0]=='GEOMETRY' or row[0]=='mesh':
            self.section = 'geo'
            return
        elif row[0]=='SYSTEM_SYSTEM':
            self.section = 'system'
            return
        elif row[0]=='bath_velocity' or row[0]=='bath velocity':
            self.section = 'system'
            
        if len(row)<3:
            row.append('')
        if self.section=='time':
            self.addTiming(row[0], row[1], row[2])
        elif self.section=='geo':
            self.geo.addVal(row[0], row[1], row[2])
        elif row[0].startswith('ink'):
            self.ink.addVal(row[0], row[1], row[2])
        elif row[0].startswith('sup') or row[0].startswith('bath'):  
            self.sup.addVal(row[0].replace('bath', 'sup'), row[1], row[2])
        elif row[0]=='sigma':
            self.addSigma(row[1], row[2])

    def readLegend(self):
        '''read all values from the legend'''
        legendFile = os.path.join(self.folder, 'legend.csv')
        if not os.path.exists(legendFile):
            raise FileNotFoundError(f'No legend file in {self.folder}')
        self.section = 'time'
        with open(legendFile, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='|')
            for row in reader:
                # save all rows as class attributes
                if row[0]=='controlDict':
                    return
                self.readLegendRow(row)    
       
    def getVal(self, s:str) -> Any:
        '''find a value'''
        if ',' in s:
            # this is a list of values
            return ', '.join([str(self.getVal(i)) for i in re.split(',', s)])
        if '*' in s or '/' in s or '+' in s or '-' in s:
            # this is an expression
            spl = re.split('\*|/|\+|-', s)
            for sub in spl:
                s = s.replace(sub, f'self.getVal(\'{sub}\')')
            return eval(s)
        
        
        if s=='const':
            return 0
        elif s=='iv' or s=='vink' or s=='ink_velocity':
            out = self.ink.dynamic['v']
        elif s=='bv' or s=='vsup' or s=='sup_velocity':
            out = self.sup.dynamic['v']
        elif s=='ink_transportModel':
            out = self.ink.transportModel
        elif s=='sup_transportModel':
            out = self.sup.transportModel
        elif len(s)>3 and s[-3:]=='ink':
            out = self.ink.dynamic[s[:-3]]
        elif len(s)>3 and s[-3:]=='sup':
            out = self.sup.dynamic[s[:-3]] 
        elif s=='sigma':
            out = getattr(self, s)
        elif hasattr(self.geo, s):
            out = getattr(self.geo, s)
        elif hasattr(self, s):
            out = getattr(self, s)
        else:
            return np.nan
        if type(out) is float and abs(round(out)-out)<10**-10:
            return int(out)
        else:
            return out
        
        
    def getUnits(self, s:str) -> str:
        '''get the units for a variable'''
        if ',' in s:
            # this is a list of values
            l = [self.getUnits(i) for i in re.split(',', s)]
            j = ', '.join(l)
            if len(j)>len(l)*2:
                return j
            else:
                return ''
        if '*' in s or '/' in s or '+' in s or '-' in s:
            # this is an expression
            spl = re.split('\*|/|\+|-', s)
            for sub in spl:
                s = s.replace(sub, self.getUnits(sub))
            return s
        
        if s=='const':
            return ''
        elif s=='iv' or s=='vink' or s=='ink_velocity':
            return self.ink.dynamicU['v']
        elif s=='bv' or s=='vsup' or s=='sup_velocity':
            return self.sup.dynamicU['v']
        elif len(s)>3 and s[-3:]=='ink':
            return self.ink.dynamicU[s[:-3]]
        elif len(s)>3 and s[-3:]=='sup':
            return self.sup.dynamicU[s[:-3]] 
        elif s=='sigma':
            return self.u[s]
        elif hasattr(self.geo, s):
            return self.geo.u[s]
        else:
            return ''
            
    def meta(self, xvar:str, yvar:str, cvar:str, splitxvar:str, splityvar:str, mvar:str, restrictions:dict) -> dict:
        '''get just the values that matter'''
        self.xval = self.getVal(xvar)
        self.yval = self.getVal(yvar)
        self.cval = self.getVal(cvar)
        self.splitxval = self.getVal(splitxvar)
        self.splityval = self.getVal(splityvar)
        self.mval = self.getVal(mvar)
        out = {'folder':self.folder, 'bn':self.bn, 'xvar':self.xval, 'yvar':self.yval
               , 'cvar':self.cval, 'splitxvar':self.splitxval, 'splityvar':self.splityval
              , 'mvar':self.mval}
        for key in restrictions:
            out[key] = self.getVal(key)
        return out
    
    def getLabel(self, var:str, short:bool=True, name:bool=True, units:bool=True, **kwargs) -> str:
        '''get a label for a variable name. If name, include the variable name. If units, include the units.'''
        if name:
            if units:
                u = self.getUnits(var)
                if len(u)>0:
                    if 'val' in kwargs:
                        v = kwargs['val']
                        label = f'{self.vn.toEnglish(var)}={v} {u}'
                    else:
                        label = f'{self.vn.toEnglish(var)} ({u})'
                else:
                    if 'val' in kwargs:
                        v = kwargs['val']
                        label = f'{self.vn.toEnglish(var)}={v}'
                    else:
                        label = self.vn.toEnglish(var)
            else:
                if 'val' in kwargs:
                    v = kwargs['val']
                    label = f'{self.vn.toEnglish(var)}={v}'
                else:
                    label = self.vn.toEnglish(var)
        else:
            if units:
                label = self.getUnits(var)
            else:
                return ''
        if short:
            label = self.vn.symbolic(self.vn.shorten(label))
        return label
    
    def viscRatio(self) -> float:
        nuink = self.ink.localVisc()
        nusup = self.sup.localVisc()
        return {'folder':self.folder
                 , 'viscRatio':nuink['ink_visc']/nusup['sup_visc']
                , 'ReRatio':nuink['ink_Re']/nusup['sup_Re']}
    
    
    def currentTime(self) -> float:
        '''Get the latest simulated time
        Input folder should be a simulation folder, e.g. "C:\\...\\nb30"
        If there is no vtk file or simulation folders, it takes the value from the legend'''
        rdict = self.currentRate()
        endTime = rdict['end_time']
        
        t2 = self.fh1.times()
        if len(t2)>0:
            simtime = max(t2)
        else:
            simtime = rdict['simulation_time']

        try:
            simtime = float(simtime)
        except:
            pass
        try:
            endTime = float(endTime)
        except:
            pass
        return simtime, endTime

    def currentRate(self) -> dict:
        '''get the current time, end time, and simulation rate from the folder'''
        out =  {'simulation_time':'', 'end_time':'', 'rate':'', 'run_time':''}
        for s in ['endTime', 'endTime_s', 'endTime_(s)']:
            if hasattr(self, s):
                out['end_time'] = getattr(self, s)
        for s in ['simulation_time', 'simulation_time_(s)', 'simulation_time_s']:
            if hasattr(self, s):
                out['simulation_time'] = getattr(self, s)
        for s in ['_interFoam_time', 'interFoam_time_(s)', 'interFoam_time_s']:
            if hasattr(self, s):
                out['run_time'] = getattr(self, s)
        for s in ['simulation_rate', 'simulation_rate_(hr/s)', 'simulation_rate_hr/s']:
            if hasattr(self, s):
                out['rate'] = getattr(self, s)
        return out
        
    def metaRow(self) -> Tuple[dict, dict]:
        '''get a row describing the inputs into this system'''
        d = {'folder':self.folder, 'bn':self.bn, 'sigma':self.sigma}
        u = {'folder':'', 'bn':'', 'sigma':self.u['sigma']}
        for o in self.ink, self.sup, self.geo:
            di, ui = o.metaRow()
            d = {**d, **di}
            u = {**u, **ui}
        for fluid in ['ink', 'sup']:
            fo = getattr(self, fluid)
            if 'visc0' in fo.dynamic:
                d[f'{fluid}_CaInv'] = self.sigma/(fo.dynamic['visc0']*fo.dynamic['d'])
                u[f'{fluid}_CaInv'] = ''
        if 'visc0' in self.ink.dynamic and 'visc0' in self.sup.dynamic:
            d['viscRatio'] = self.ink.dynamic['visc0']/self.sup.dynamic['visc0']
            u['viscRatio'] = ''
            d['ReRatio'] = self.ink.dynamic['Re']/self.sup.dynamic['Re']
            u['ReRatio'] = ''
        return d,u
