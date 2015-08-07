# -*- coding: utf-8 -*-
"""
MultiBlock
======

**MultiBlock** is a python package that solves (AB)n like multiblock copolymer melts and analyzes the loop and
bridge chain conformation. 

References
----------

* K. Rasmussen et al Journal of Polymer Science: Part B: Polymer Physics, Vol. 41, 104 –111 (2003)

"""

__author__ = "Jiuzhou Tang <tangjiuzhou@iccas.ac.cn>"
__license__ = "BSD License"
__version__ = "Alpha"

from .common import *
from .unitcell import *
from .sim_box import *
#from .symmetry import *
#from .space_group import *
#from .group import *
from .grid import *
from .chem import *
from .chain import *
from .mde import *
from .utility import *
from .density_FE import *
from .scft import *
from .scft_io import *
from .field_init import *
from .iterator import *

