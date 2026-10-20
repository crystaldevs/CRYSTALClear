from CRYSTALClear import base
from CRYSTALClear import calculate
from CRYSTALClear import convert
from CRYSTALClear import execute
from CRYSTALClear import geometry
from CRYSTALClear import crystal_io
from CRYSTALClear import plot
from CRYSTALClear import units

import sys

if sys.version_info >= (3, 8):
    from importlib import metadata
else:
    from importlib_metadata import metadata

__author__ = "CRYSTAL devs"
__version__ = metadata.version('CRYSTALClear')
