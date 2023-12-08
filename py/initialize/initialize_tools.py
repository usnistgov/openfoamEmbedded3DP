#!/usr/bin/env python
'''Functions to generate OpenFOAM input files for a nozzle in a 3D bath'''


# external packages
from typing import List, Dict, Tuple, Union, Any, TextIO, Callable
import logging
# local packages


# logging
logging.basicConfig(level=logging.INFO)



#-------------------------------------------------------------------------------------------------  


def scale() -> float:
    '''scale all units by this amount * m'''
    return 0.001 

        
class OpenFOAMFile:
    '''tools for creating openfoam files'''
    
    def __init__(self):
        return
    

    def header(self, cl:str, obj:str) -> str:
        '''header for openfoam files
        cl is the class name (string)
        obj is the object name (string)
        '''
        s = ("/*--------------------------------*- C++ -*----------------------------------*\n"
            +"| =========                 |                                                 |\n"
            +"| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
            +"|  \\\\    /   O peration     | Version:  8                               |\n"
            +"|   \\\\  /    A nd           | Website:  www.openfoam.com                      |\n"
            +"|    \\\\/     M anipulation  |                                                 |\n"
            +"*---------------------------------------------------------------------------*/\n"
            +"FoamFile\n"
            +"{\n"
            +"    version     2.0;\n"
            +"    format      ascii;\n"
            +f"    class       {cl};\n"
            +f"    object      {obj};\n"
            +"}\n"
            +"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //;\n\n"
            )
        return s

    def closeLine(self) -> str:
        '''put this at the end of files'''
        return "// ************************************************************************* //";