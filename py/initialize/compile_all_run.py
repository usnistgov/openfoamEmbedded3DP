#!/usr/bin/env python
'''Functions to generate OpenFOAM input files for a nozzle in a 3D bath'''

# external packages
from typing import List, Dict, Tuple, Union, Any, TextIO, Callable
import logging
import os
import numpy as np

# local packages


# logging
logging.basicConfig(level=logging.INFO)

#-------------------------------------------------------------------------------------------------  

class bashCompiler:
    '''class for bash files'''
    
    def __init__(self):
        return
    
    def opening(self):
        s = self.binbash()
        s = s + 'cd \"$(dirname \"$0\")\"\n'
        return s
    
    def binbash(self):
        return "#!/bin/bash\n\n"
    
    def fListLoop(self, s:str, functionlist:List[str], folder:str, ifstarted:bool=False) -> str:
        '''Write a bash script to go through multiple functions. ifstarted True if you do not want to run this function if the sim has started'''
        for f in functionlist:
            s1 = f'echo \"running {f} in {os.path.basename(folder)}\";\n{f}>>log_{f}'
            if ifstarted:
                s = f'{s}[ ! d \"0.1\"] && ({s1}); '
            else:
                s = f'{s}{s1};\n'
        return s


class compileAllClean(bashCompiler):
    '''this removes residual data from previous runs'''
    
    def __init__(self, endTime:float, writeDt:float):
        super().__init__()
        s = self.opening()
        s = s + 'rm -f -r 0.* '
        for t in range(1, int(np.ceil(endTime+writeDt)+1), 1):
            s = f'{s}{t}* '
        s = s + "slurm*.out log_*"
        self.s = s

class compileAllAllClean(bashCompiler):
    '''for folders that contain both a mesh and case folder, create a function that runs both allrun functions'''
    
    def __init__(self):
        super().__init__()
        s = self.opening()
        s = s + 'rm -f -r slurm*.out log_*\n'
        s = s + ("./case/Allclean.sh")
        self.s = s


class compileAllAllRun(bashCompiler):
    '''for folders that contain both a mesh and case folder, create a function that runs both allrun functions'''
    
    def __init__(self):
        super().__init__()
        s = self.opening()
        s = s + ("./mesh/Allrun.sh; ./case/Allrun.sh")
        self.s = s

class compileAllRun(bashCompiler):
    '''this is the allrun bash script for the case folder'''
    
    def __init__(self, folder:str, solver:str): # RG
        super().__init__()
        f = os.path.basename(folder)
        s = self.opening()
        s = s + '. $WM_PROJECT_DIR/bin/tools/RunFunctions;\n'
        s = s + f'echo {f}\n'
        s = s + 'if [ ! -d "0.1" ]; then\n'
        s = s + '\tcp -r ../mesh/constant/polyMesh constant;\n'
        s = s + '\tcp 0/alpha.ink.orig 0/alpha.ink;\n'
        s = s + f'\techo \"running setFields in {f}\";\n'
        s = s + '\tsetFields>>log_setFields;\n'
        s = s + 'fi \n'
        s = self.fListLoop(s, [solver, "foamToVTK"], f, ifstarted=False)
        self.s = s

class compileSlurm(bashCompiler):
    '''this is the slurm script for the case folder'''
    
    def __init__(self, folder:str, parentdir:str):
        super().__init__()
        # workdir = (os.path.join(parentdir, os.path.basename(folder))).replace("\\","/")
        workdir = parentdir.replace("\\","/") # RG
        s = f'#!/bin/bash\n#SBATCH -p local\n'
        s = s + '#SBATCH --time=42-00:00:00\n'
        s = s + '#SBATCH --nodes=1\n'
        s = s + '#SBATCH --cpus-per-task=1\n'
        s = s + f'#SBATCH --job-name={os.path.basename(folder)}\n'
        s = s + f'#SBATCH --chdir={workdir}\n'
        s = s + f'#SBATCH --partition=local\n\n'
        s = s + 'export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}\n\n'
        s = s + 'srun bash mesh/Allrun.sh\n'
        s = s + 'srun bash case/Allrun.sh'
        self.s = s

class compileAllRunMesh(bashCompiler):
    '''this script runs the meshing functions in the mesh folder'''
    
    def __init__(self, folder:str):
        super().__init__()
        s = self.binbash()
        s = s + 'cd \"$(dirname \"$0\")\" || exit; \n'
        s = s + '. $WM_PROJECT_DIR/bin/tools/RunFunctions;\n'
        s = s + 'if [ ! -d "VTK" ]; then\n'
        functionlist = ["surfaceFeatures", "blockMesh", "snappyHexMesh -overwrite", "foamToVTK"]
            # use surfaceFeatures for OpenFOAM 8 and later
            # use surfaceFeatureExtract otherwise
        s = self.fListLoop(s, functionlist, folder)
        s = s + 'fi \n'
        self.s = s