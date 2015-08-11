# -*- coding: utf-8 -*-
"""
MultiBlock
=========

:copyright: (c) 2015 by Jiuzhou Tang
:license: BSD, see LICENSE for more details.

"""
import pylab as pl

import numpy as np
from scipy import optimize
#import scipy.io
#import matplotlib.pyplot as plt
#from mayavi import mlab
from scft import Scft_Obj
from sim_box import SimBox, BoxShape
from grid import Grid
from chem import Chem
from chain import Chain
from field_init import Fields_Init
from mde import *
from utility import simposon_int_pbc
from density_FE import *
from iterator import *
from scft_io import Scft_IO
def run_scft(D):
    abn_melt=Scft_Obj("ABn_melt")
    cell=[D]
    #cell=[1.70]
    print "D=",D
    box=SimBox(1,cell)
    Nx=128
    grid=Grid(box,box.shape,[Nx])
    XN=100.0 
    AB=Chem(2,grid,XN)
    N_intpo=1
    Ns=500
    N_blk=10
    ABn_multiblk=Chain(N_blk,Ns,"ABn",N_intpo,grid)   
# prepare the initial lamellar field
    field_param=[0,2] # on X axis, one period     
    Fields_Init(ABn_multiblk,AB,grid,"lamellar",field_param)
    #pl.plot(grid.x,AB.W[0,:,0,0],'ro')
    for i in range(300):      
        pesudo_spectrum_mde(ABn_multiblk,AB,grid)
        density(ABn_multiblk,AB,grid,abn_melt)
        SimpleMixing_AB(AB,ABn_multiblk,grid)  
    
    free_energy_abmelt(ABn_multiblk,AB,grid,abn_melt)
    #io=Scft_IO(["Dump_Field","Dump_Result"],chem=AB)
    print "FE:",abn_melt.f_tot, abn_melt.f_FH
    pl.plot(grid.x,AB.R[0,:,0,0],'ro')
    pl.plot(grid.x,AB.R[1,:,0,0],'b^')
    pl.show()
    #crank_nickson_mde(ABn_multiblk,AB,grid)
    #t=np.arange(Ns+1)
    #pl.plot(grid.x,AB.R[0,:,0,0]  ,'bo',label="RA")
    #pl.plot(grid.x,AB.R[1,:,0,0]  ,'r^',label="RB")
    #pl.plot(t[0:Ns+1],np.log(ABn_multiblk.qf[0:Ns+1,0,0,0])  ,'ro',label="qf")
    #pl.plot(t[0:501],ABn_multiblk.qf[0:501,120,0,0],'go',label="qf")
    #pl.plot(t[0:Ns+1],np.log(ABn_multiblk.qf1[0:Ns+1,0,0,0]),'b^',label="qf1")
    #data=ABn_multiblk.qf[0:Ns,:,0,0]-ABn_multiblk.qf1[0:Ns,:,0,0]
    #print "shape",np.shape(data)
    #data=np.reshape(data,(Ns,len(grid.x)))  
    #pl.imshow(data)
    #pl.plot(grid.x,data[100,:],'r^',label="qf")
    #pl.plot(grid.x,data[300,:],'b^',label="qf")
    #pl.show()
    return abn_melt.f_tot



  

    pl.plot(grid.x,AB.R[0,:,0,0])
    t=np.arange(501)
    t1=np.arange(500*N_intpo+1)
    t1=t1/N_intpo
    #pl.plot(t[0:501],np.log(ABn_multiblk.qf[0:501,0,0,0]),'bo',label="qf")
    #pl.show()
    #pl.plot(t[0:501],np.log(ABn_multiblk.qf[0:501,0,0,0]),'bo',label="qf")
    crank_nickson_mde(ABn_multiblk,AB,grid)
    ABn_multiblk.qf[:,:,:,:]=ABn_multiblk.qf1[:,:,:,:] 
    ABn_multiblk.qb[:,:,:,:]=ABn_multiblk.qb1[:,:,:,:] 
    #density(ABn_multiblk,AB,grid,abn_melt)
    pl.plot(grid.x,AB.R[0,:,0,0],'ro')
    pl.plot(grid.x,AB.R[1,:,0,0],'b^')
    #pl.show()
    #pl.plot(t,np.log(ABn_multiblk.qloop[:,0,0,0]))
    #pl.show()
    #crank_nickson_mde(ABn_multiblk,AB,grid)
    #print "q_loop",ABn_multiblk.qloop[:,0,0,0]
    #pl.plot(t,ABn_multiblk.qloop[:,0,0,0]-ABn_multiblk.qf[:,0,0,0])
    #pl.plot(t1[0:1000],ABn_multiblk.qloop[0:1000,0,0,0],label="qloop")
    #pl.plot(t[0:501],ABn_multiblk.qf[0:501,0,0,0],'bo',label="qf")
    #pl.plot(t[0:501],np.log(ABn_multiblk.qf1[0:501,0,0,0]),'r^',label="qf1")
    #pl.plot(grid.x,AB.R[0,:,0,0],'ro',label="RA")
    #pl.plot(grid.x,AB.R[1,:,0,0],'bo',label="RB")
    pl.showbar()
    pl.show()
    print "RA",AB.R[0,:,0,0]
    





if __name__ == '__main__':
    run_scft(1.7)
    #optimize.brent(run_scft, brack=(1.7,1.9,2.1), tol=1.0e-5,full_output=0, maxiter=500)
