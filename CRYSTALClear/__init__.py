from CRYSTALClear import adsorb
from CRYSTALClear import base
from CRYSTALClear import calculate
from CRYSTALClear import convert
from CRYSTALClear import execute
from CRYSTALClear import geometry
from CRYSTALClear import crystal_io
from CRYSTALClear import plot
from CRYSTALClear import thermodynamics
from CRYSTALClear import units
#from CRYSTALClear import utils

import sys

if sys.version_info >= (3, 8):
    from importlib import metadata
else:
    from importlib_metadata import metadata

__author__ = "CRYSTALClear Development Team"
__email__ = "crystalunito@gmail.com"
__version__ = metadata.version('CRYSTALClear')
