# -*- coding: utf-8 -*-

import numpy as np

import pylab as pl
__all__ = ["Fields_Init"]

def Fields_Init(chain,chem,grid,field_type,field_param=None):
    nx=len(grid.x) 
    ny=len(grid.y) 
    nz=len(grid.z)
    
    if field_type == "lamellar" :
        if len(field_param) !=2:
            raise ValueError('lamellar field takes two params.')
        else:  
            orient  = field_param[0]
            periods = field_param[1]
    else:
        raise ValueError('currently only lamellar field can be initialized.')
         
  
    for x in np.arange(nx) :
        for y in np.arange(ny) :  
            for z in np.arange(nz) :
                chem.R[0,x,y,z]=chain.f_sp[0]*(1.0+0.8*np.cos(2*np.pi*periods*x/nx))  
                chem.R[1,x,y,z]=1.0-chem.R[0,x,y,z]
                #chem.R[0,x,y,z]=np.random.rand(1)
                #chem.R[1,x,y,z]=1.0-chem.R[0,x,y,z]

    chem.W[0,:,:,:]=chem.XN*chem.R[1,:,:,:]        
    chem.W[1,:,:,:]=chem.XN*chem.R[0,:,:,:]        

 #   pl.plot(grid.x,chem.R[0,:,0,0])
 #   pl.show()                    
