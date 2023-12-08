#!/usr/bin/env python
'''Functions to generate OpenFOAM input files for a nozzle in a 3D bath'''


# external packages
from typing import List, Dict, Tuple, Union, Any, TextIO, Callable
import logging

# local packages

# logging
logging.basicConfig(level=logging.INFO)

#------------------------------------------------------------------------------------------------- 


class CDVars:
    '''controlDict variables '''
        
    def __init__(self, startTime:float, endTime:float, dt:float, writeDt:float):
        '''Inputs: startTime, endTime, dt, writeDt
        startTime is in seconds
        endTime is in seconds
        dt is the initial time step of the simulation in seconds
        writeDt is the time step at which to write results to file'''
        # solver
        self.application = "interFoam"
        
        # time control
        self.startFrom = "latestTime"
        self.startTime = startTime
        self.stopAt = "endTime"
        self.endTime = endTime
        self.deltaT = dt # Time step of the simulation
        
        # time step control
        self.adjustTimeStep = "yes" # yes/no† to adjust time step according to maximum Courant number in transient simulation.
        self.maxCo = 1 # Maximum Courant number allowed.
        self.maxAlphaCo = 1
        self.maxDeltaT = 1 # Maximum time step allowed in transient simulation
        
        # data writing
        #self.writeControl = "adjustable" # Writes data every writeInterval seconds of simulated time.
        self.writeControl = "adjustableRunTime"
        self.writeInterval = writeDt # Scalar used in conjunction with writeControl described above.
        self.purgeWrite = 0 # Integer representing a limit on the number of time directories that are stored by overwriting time directories on a cyclic basis. Example of t0 = 5s, Δt = 1s and purgeWrite 2;: data written into 2 directories, 6 and 7, before returning to write the data at 8 s in 6, data at 9 s into 7, etc. To disable the time directory limit, specify purgeWrite 0;† For steady-state solutions, results from previous iterations can be continuously overwritten by specifying purgeWrite 1;
        self.writeFormat = "ascii" # ASCII format, written to writePrecision significant figures.
        self.writePrecision = 6 # Integer used in conjunction with writeFormat described above, 6† by default
        self.writeCompression = "uncompressed" # Specifies the compression of the data files.
        self.timeFormat = "general" # Specifies scientific format if the exponent is less than -4 or greater than or equal to that specified by timePrecision.
        self.timePrecision = 6 # Integer used in conjunction with timeFormat described above, 6† by default
        
        # data reading
        self.runTimeModifiable = "yes" # yes†/no switch for whether dictionaries, e.g.controlDict, are re-read by OpenFOAM at the beginning of each time step.

    
    def varList(self) -> List[List[str]]:
        '''Constructs a table with all of the variables held in this object. The table has two columns: the variable name and the value'''
        l = []        
        for attr, value in self.__dict__.items():
            l.append([attr, value])
        return l