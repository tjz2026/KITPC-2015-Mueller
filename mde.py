# -*- coding: utf-8 -*-
import numpy as np
from numpy.linalg import inv
import scipy as sc

import scipy.sparse as sparse
import scipy.sparse.linalg
import pylab as pl

def pesudo_spectrum_mde(chain,chem,grid):

    ds=1.0/chain.Ns
    exp_w=np.zeros((chem.num,len(grid.x),len(grid.y),len(grid.z)))
    exp_k2=np.zeros((len(grid.x),len(grid.y),len(grid.z)))
    exp_w[:,:,:,:]=np.exp(chem.W[:,:,:,:]*(-0.5*ds))
    local_data=np.zeros((len(grid.x),len(grid.y),len(grid.z)),dtype=complex)

    for x in np.arange(len(grid.x)) :
        for y in np.arange(len(grid.y)) :  
            for z in np.arange(len(grid.z)) :  
                if len(grid.x)==1:
                    raise ValueError('Unkonwn dimension for mde.')
                    
                if len(grid.y)==1 and len(grid.z)==1:
                    k2=grid.kx[x]**2
                if len(grid.y)!=1 and len(grid.z)==1:
                    k2=grid.kx[x]**2+grid.ky[y]**2
                if len(grid.y)!=1 and len(grid.z)!=1:
                    k2=grid.kx[x]**2+grid.ky[y]**2+grid.kz[z]**2
                exp_k2[x][y][z]=np.exp(-0.5*ds*k2)
    chain.qf[0,:,:,:]=1.0
     
    for s in np.arange(1,chain.Ns+1):
        sp=chain.blk_sp[chain.blk_s[s]] 
        local_data[:,:,:]=chain.qf[s-1,:,:,:]*exp_w[sp,:,:,:]+0.0j
        local_data=np.fft.fftn(local_data)
        local_data[:,:,:]=local_data[:,:,:]*exp_k2[:,:,:] 
        local_data=np.fft.ifftn(local_data)
        chain.qf[s,:,:,:]=local_data[:,:,:].real*exp_w[sp,:,:,:]
    
    chain.qb[chain.Ns,:,:,:]=1.0
     
    for s in np.arange(0,chain.Ns)[::-1]:
        sp=chain.blk_sp[chain.blk_s[s]] 
        local_data[:,:,:]=chain.qb[s+1,:,:,:]*exp_w[sp,:,:,:]+0.0j
        local_data=np.fft.fftn(local_data)
        local_data[:,:,:]=local_data[:,:,:]*exp_k2[:,:,:] 
        local_data=np.fft.ifftn(local_data)
        chain.qb[s,:,:,:]=local_data[:,:,:].real*exp_w[sp,:,:,:]




def crank_nickson_mde(chain,chem,grid,q_init):

    ds=1.0/chain.Ns
    q_init=np.append(chain.qf[s,:,0,0],chain.qf[s,0,0,0])
    for i in range(1,chain.n_blk+1,2):
        start = chain.blk_sta[i]
        end   = chain.blk_end[i]  
    for s in np.arrange(1:chain.Ns+1):
       
        q_init=np.append(chain.qf[s,:,0,0],chain.qf[s,0,0,0])
        blk=chain.blk_s[s]
        nch=chain.blk_sp[blk]
        Ws[:]=np.append(chem.W[nch,:,0,0],chem.W[nch,0,0,0])
        size=grid.x[1]*len(grid.x)
        crank_nicolson_onestep(len(grid.x),size,chain.Ns,q_init)
       
    


def crank_nicolson_onestep(Nx,Lx,Ns,Ws,q_init,q_loop):
"""
	This program solves the 1D modified diffusion equation
		u_t = u_xx-w*u
	with reflective boundary condition
		u(-1,t) = u(1,t) = 0
		u(Nx+1,t) = u(Nx-1,t) = 0
                  
	with the Initial Conditions
		u(x,0) = 1.0
	over the domain x = [0, Lx] Nx+1 points starts at 0 and t= [0,1] Ns+1 point
 
	The program solves the heat equation using a finite difference method
	where we use a center difference method in space and Crank-Nicolson 
	in time.
"""
 
# Number of internal points
 
# Calculate Spatial Step-Size
dx = Lx/Nx
 
# Create Temporal Step-Size, TFinal, Number of Time-Steps
dt=1/Ns

x=np.linspace
k = h/2
TFinal = 1
NumOfTimeSteps = int(TFinal/k)
 
# Create grid-points on x axis
x = np.linspace(0,Lx,Nx+1)
t = np.linspace(0,1,Ns+1)

 
# Initial Conditions
u = np.transpose(np.mat(q_init))
 
# Second-Derivative Matrix
data = np.ones((3, Nx+1))
data[1] = -2*data[1]
data[1,:]=data[1,:]-0.5*dt*Ws[:]

# Reflective boundary 
data[2,0]=2*data[2,0]
data[0,Nx]=2*data[0,Nx]
diags = [-1,0,1]
D2 = sparse.spdiags(data,diags,Nx+1,Nx+1)/(dx**2)
 
# Identity Matrix
I = sparse.identity(Nx+1)
 
# Data for each time-step
data = []
 
# Solve the System: (I - k/2*D2) u_new = (I + k/2*D2)*u_old
A = (I -dt/2*D2)
b = ( I + dt/2*D2 )*u
u = np.transpose(np.mat( sparse.linalg.spsolve( A,  b ) ))
 
# Save Data
chem.qf[]
data.append(u)



 
             








