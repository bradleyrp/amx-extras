#!/usr/bin/env python

from .mesh import *
from .structural_biology import *
from .geometry_tools import  *
from .Multimer import *
#! trying to use magic_importer from amx here
from amx import state
mesh.state = state
structural_biology.state = state
#! etc

#from ortho.imports import magic_importer
#import os
#magic_importer(globals(),os.path.dirname(__file__))
