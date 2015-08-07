# -*- coding: utf-8 -*-

import numpy as np

iter_eps=1e-4
lambda_t=0.1

__all__ = ["SimpleMixing_AB"]

def SimpleMixing_AB(chem,chain,grid):
                 
    nx=len(grid.x) 
    ny=len(grid.y) 
    nz=len(grid.z)
    field_err=np.zeros(2)
    wA_tmp=np.zeros((nx,ny,nz))
    wB_tmp=np.zeros((nx,ny,nz))
    yita=np.zeros((nx,ny,nz))

    yita[:,:,:]=0.5*(chem.W[0,:,:,:]+chem.W[1,:,:,:])
    wA_tmp[:,:,:]=chem.XN*(chem.R[1,:,:,:]-chain.f_sp[1])+yita[:,:,:]
    wB_tmp[:,:,:]=chem.XN*(chem.R[0,:,:,:]-chain.f_sp[0])+yita[:,:,:]

    field_err[0]=np.max(np.abs(wA_tmp[:,:,:]-chem.W[0,:,:,:]))
    field_err[1]=np.max(np.abs(wB_tmp[:,:,:]-chem.W[1,:,:,:]))

    chem.W[0,:,:,:]=chem.W[0,:,:,:]+lambda_t*(wA_tmp[:,:,:]-chem.W[0,:,:,:])
    chem.W[1,:,:,:]=chem.W[1,:,:,:]+lambda_t*(wB_tmp[:,:,:]-chem.W[1,:,:,:])
    print "field_err",field_err 
    return field_err

    
  
    
















