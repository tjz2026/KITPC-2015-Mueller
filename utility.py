# -*- coding: utf-8 -*-

from scipy.integrate import simps
import numpy as np

def simposon_int_pbc(f,grid):
    if grid.dim==1 :
        integral=simps(np.append(f[:,0,0],f[0,0,0]),dx=grid.x[1])
    elif grid.dim==2 :
        integral=0.0
        pass 
    elif grid.dim==3 :
        integral=0.0
        pass
    else:
        raise ValueError('Unkonwn dimension for simposon integral.')
    
    return integral     
