# -*- coding: utf-8 -*-
"""
gyroid.common
=============

This moudle define all global constants.

It includes::

    EPS, SMALL, LARGE
    BRAVAIS, CARTESIAN
    LAMELLAR
    SQUARE, RECTANGULAR, HEXAGONAL, OBLIQUE
    CUBIC, TETRAGONAL, ORTHORHOMBIC, TRIGONAL, MONOCLINIC, TRICLINIC
    DEFAULT
    CRYSTAL_SYSTEM1, CRYSTAL_SYSTEM2, CRYSTAL_SYSTEM3

"""

import numpy as np
from numpy.linalg import inv

# To let the G2 comparison pass, EPS should not be too small
EPS = 1e-6

SMALL = 1e-10
LARGE = 1e+10

BRAVAIS = "Bravais"
CARTESIAN = "Cartesian"

LAMELLAR = "Lamellar"
SQUARE = "Square"
RECTANGULAR = "Rectangular"
HEXAGONAL = "Hexagonal"
OBLIQUE = "Oblique"
CUBIC = "Cubic"
TETRAGONAL = "Tetragonal"
ORTHORHOMBIC = "Orthorhombic"
TRIGONAL = "Trigonal"
MONOCLINIC = "Monoclinic"
TRICLINIC = "Triclinic"
DEFAULT = "Default"
CRYSTAL_SYSTEM1 = [LAMELLAR,DEFAULT]
CRYSTAL_SYSTEM2 = [SQUARE,RECTANGULAR,HEXAGONAL,OBLIQUE,DEFAULT]
CRYSTAL_SYSTEM3 = [CUBIC,TETRAGONAL,ORTHORHOMBIC,HEXAGONAL,
                   TRIGONAL,MONOCLINIC,TRICLINIC,DEFAULT]





