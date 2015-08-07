# -*- coding: utf-8 -*-

import numpy as np
from common import EPS
from scipy.integrate import simps
from utility import simposon_int_pbc

__all__ = ["density","free_energy_abmelt"]


def density(chain,chem,grid,scft_obj):
    nx=len(grid.x) 
    ny=len(grid.y) 
    nz=len(grid.z)
    ds=1.0/chain.Ns 
    if grid.dim==1:
        Mv=nx*grid.x[1]
    if grid.dim==2:
        Mv=nx*grid.x[1]*ny*grid.y[1]
    if grid.dim==3:
        Mv=nx*grid.x[1]*ny*grid.y[1]*nz*grid.z[1]
    if np.abs(Mv)<EPS:
        raise ValueError('Unkonwn Mv.')
      
    ar=np.zeros((nx,ny,nz))
    ar[:,:,:]=chain.qf[chain.Ns,:,:,:]
    chain.Q=simposon_int_pbc(ar,grid)/Mv
    scft_obj.logQ=np.log(chain.Q)
 
    R_s=np.zeros(chain.Ns+1)
    for x in np.arange(nx) :
        for y in np.arange(ny) :  
            for z in np.arange(nz) :
                ar[x,y,z]=0.0  
                chem.R[:,x,y,z]=0.0
                for blk in np.arange(chain.n_blk):
                    sum_blk=0.0
                    sta=chain.blk_sta[blk]
                    end=chain.blk_end[blk]
                    for s in np.arange(sta,end+1): # this is a trick to make sure we have odd number when doing simps
                        R_s[s]=chain.qf[s,x,y,z]*chain.qb[s,x,y,z]
                    if blk!=chain.n_blk-1:
                        end=end+1
                        R_s[end]=chain.qf[end,x,y,z]*chain.qb[end,x,y,z]
                    sum_blk=simps(R_s[sta:end+1:1],dx=ds)
                    ar[x,y,z]=ar[x,y,z]+sum_blk 
                    nch=chain.blk_sp[blk]
                    chem.R[nch,x,y,z]=chem.R[nch,x,y,z]+sum_blk


    totden=simposon_int_pbc(ar,grid)/Mv
    chem.R[:,:,:,:]=chem.R[:,:,:,:]/totden  

                        
def free_energy_abmelt(chain,chem,grid,scft_obj):
    nx=len(grid.x) 
    ny=len(grid.y) 
    nz=len(grid.z)
    if grid.dim==1:
        Mv=nx*grid.x[1]
    if grid.dim==2:
        Mv=nx*grid.x[1]*ny*grid.y[1]
    if grid.dim==3:
        Mv=nx*grid.x[1]*ny*grid.y[1]*nz*grid.z[1]
    if np.abs(Mv)<EPS:
        raise ValueError('Unkonwn Mv.')
                 
    ar=np.zeros((nx,ny,nz))
    ar[:,:,:]=chem.XN*chem.R[0,:,:,:]*chem.R[1,:,:,:]
    ar[:,:,:]=ar[:,:,:]-chem.W[0,:,:,:]*chem.R[0,:,:,:]-chem.W[1,:,:,:]*chem.R[1,:,:,:]
    ar[:,:,:]=ar[:,:,:]+0.5*(chem.W[0,:,:,:]+chem.W[1,:,:,:])*(chem.R[0,:,:,:]+chem.R[1,:,:,:]-1.0)
    
    scft_obj.f_FH=simposon_int_pbc(ar,grid)/Mv
    scft_obj.f_tot=scft_obj.f_FH-scft_obj.logQ
    #print "FE in time",scft_obj.f_FH
    
















