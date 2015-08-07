# -*- coding: utf-8 -*-
"""
MultiBlock
=========

:copyright: (c) 2015 by Jiuzhou Tang
:license: BSD, see LICENSE for more details.

"""
import pylab as pl

import numpy as np
#import scipy.io
#import matplotlib.pyplot as plt
#from mayavi import mlab
from scft import Scft_Obj
from sim_box import SimBox, BoxShape
from grid import Grid
from chem import Chem
from chain import Chain
from field_init import Fields_Init
from mde import pesudo_spectrum_mde
from utility import simposon_int_pbc
from density_FE import *
from iterator import *
from scft_io import Scft_IO
def run_scft():
    abn_melt=Scft_Obj("ABn_melt")
    cell=[1.70]
    box=SimBox(1,cell)
    grid=Grid(box,box.shape,[128])
    XN=80.0 
    AB=Chem(2,grid,XN)
    ABn_multiblk=Chain(10,500,"ABn",grid)   
# prepare the initial lamellar field
    field_param=[0,2] # on X axis, one period     
    Fields_Init(ABn_multiblk,AB,grid,"lamellar",[0,1])
    for i in range(100):      
        pesudo_spectrum_mde(ABn_multiblk,AB,grid)
        density(ABn_multiblk,AB,grid,abn_melt)
        free_energy_abmelt(ABn_multiblk,AB,grid,abn_melt)
        SimpleMixing_AB(AB,ABn_multiblk,grid)  
    
    #io=Scft_IO(["Dump_Field","Dump_Result"],chem=AB)
    print "FE:",abn_melt.f_tot, abn_melt.f_FH
    #pl.plot(grid.x,AB.R[0,:,0,0])
    #pl.show()
    qshow=np.reshape(ABn_multiblk.qf,(501,128))
    pl.imshow(qshow)
    pl.show 

if __name__ == '__main__':
    run_scft()
