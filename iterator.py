# -*- coding: utf-8 -*-

import numpy as np

iter_eps=1e-4
lambda_t=0.18

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

    
class AndersonMixing(object):
    ''' A :class:`Grid` object represents an evenly Discretized sim_box cell in real and reciprocal space.

    '''

    def __init__(self,chem,grid,anderson_dim,maxium,lambda_t,tol):

        self.num=chem.num
        self.an_dim=anderson_dim
        self.i=0 # counter  
        self.maxium=maxium # maxium counter  steps
        self.lambda_t=lambda_t
        self.tol=tol
        nx=len(grid.x) 
        ny=len(grid.y) 
        nz=len(grid.z)
        W_an=np.zeros((self.num,an_dim,nx,ny,nz))
        dW_an=np.zeros((self.num,an_dim,nx,ny,nz))
        
    def __iterate(self,chem,chain,grid,i):
                 
        nx=len(grid.x) 
        ny=len(grid.y) 
        nz=len(grid.z)
        field_err=np.zeros(2)
        wA_tmp=np.zeros((nx,ny,nz))
        wB_tmp=np.zeros((nx,ny,nz))
        yita=np.zeros((nx,ny,nz))

        # doing simple mixing first
        total_err=1.0
        if i < self.nsimple:
            total_err=SimpleMixing_AB(chem,chain,grid)
        elif i==self.nsimple:
            total_err=SimpleMixing_AB(chem,chain,grid)
            yita[:,:,:]=0.5*(chem.W[0,:,:,:]+chem.W[1,:,:,:]-chem.XN)
            self.W_an[0,0,:,:,:]=chem.XN*chem.R[1,:,:,:]+yita[:,:,:] 
            self.W_an[1,0,:,:,:]=chem.XN*chem.R[0,:,:,:]+yita[:,:,:] 
            self.dW_an[0,0,:,:,:]=self.W_an[0,0,:,:,:]-chem.W[0,:,:,:]
            self.dW_an[1,0,:,:,:]=self.W_an[1,0,:,:,:]-chem.W[1,:,:,:]
        else:
            nr_tmp=min((i-self.nsimple),self.an_dim)
        pass     
                   





    
  
    
















    
















