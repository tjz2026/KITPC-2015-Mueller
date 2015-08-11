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




def crank_nickson_mde(chain,chem,grid):
    q_init=np.zeros((len(grid.x)+1,len(grid.y),len(grid.z)))
    ds=1.0/chain.Ns
    start = chain.blk_sta[1]
    end   = chain.blk_end[1]
    #end   = chain.blk_sta[1]+1
    q_loop=np.zeros(((chain.Ns+1)*10,len(grid.x)+1,len(grid.y),len(grid.z)))
    Ws=np.append(chem.W[1,:,0,0],chem.W[1,0,0,0])
    #for i in range(1,chain.n_blk+1,2):
    for i in range(0,1,1):
        print "the ith block",i
        start = chain.blk_sta[i]
        end   = chain.blk_end[i]
        start=0
        end=chain.Ns 
        #end   = chain.blk_sta[i]+1
        print "start,end",start,end
        Nt=end-start+1
        if i== chain.n_blk-1:
            Nt=Nt-1
        Nt=chain.Ns   
        #q_init=np.append(chain.qf[start-1,:,0,0],chain.qf[start-1,0,0,0])
        q_init=np.append(chain.qf[0,:,0,0],chain.qf[0,0,0,0])
        print "qinit[0],[128]",q_init[0],q_init[128]
        size=grid.x[1]*len(grid.x)
        #crank_nicolson(len(grid.x),size,Nt,ds,Ws,chem,chain,q_init,q_loop)
        Euler(len(grid.x),size,Nt*10,ds*0.1,Ws,chem,chain,q_init,q_loop)
        #print "q_loop(0,0,0,0)",q_loop[0,0,0,0]
        #chain.qloop[start:end+1,:]=q_loop[:,0:len(grid.x)]
        chain.qloop[0:Nt*10+1:]=q_loop[0:Nt*10+1,0:len(grid.x)]
        #if i== chain.n_blk-1:
        #    chain.qloop[start:end,:]=q_loop[:,0:len(grid.x)]
        #else: 
        #    chain.qloop[start:end+1,:]=q_loop[:,0:len(grid.x)]
    
    err=np.sum(np.abs(chain.qloop[1,:,0,0]-chain.qf[1,:,0,0]))   
    print "total err1 in qf and q_bar:",err 
    err=np.sum(np.abs(chain.qloop[5,:,0,0]-chain.qf[5,:,0,0]))   
    print "total err5 in qf and q_bar:",err     
    err=np.sum(np.abs(chain.qloop[20,:,0,0]-chain.qf[20,:,0,0]))   
    print "total err20 in qf and q_bar:",err     
    err=np.sum(np.abs(chain.qloop[30,:,0,0]-chain.qf[30,:,0,0]))   
    print "total err30 in qf and q_bar:",err     
    err=np.sum(np.abs(chain.qloop[200,:,0,0]-chain.qf[200,:,0,0]))   
    print "total err200 in qf and q_bar:",err     


def crank_nicolson(Nx,Lx,Nt,dt,Ws,chem,chain,q_init,q_loop):
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
     
    # Create grid-points on x axis
    x = np.linspace(0,Lx,Nx+1)
     
    # Initial Conditions
    u = np.transpose(np.mat(q_init))
     
    # Second-Derivative Matrix
    data = np.ones((3, Nx+1))
    data[1] = -2*data[1]
    #data[1,:]=data[1,:]-0.5*dt*Ws[:]*(2*dx**2/dt)
    
    # Reflective boundary 
    #data[0,0]=0.0
    data[2,0]=2*data[2,0]
    data[0,0]=2*data[0,0]
    data[2,Nx]=2*data[2,Nx]
    data[0,Nx]=2*data[0,Nx]
    #data[2,Nx]=0.0
    diags = [-1,0,1]
    #D2 = sparse.spdiags(data,diags,Nx+1,Nx+1)/(dx**2)
     
    # Identity Matrix
    I = sparse.identity(Nx+1)
     
    # Data for each time-step
    #data = []
     
    # Solve the System: (I - k/2*D2) u_new = (I + k/2*D2)*u_old
    q_loop[0,0:128,0,0]=chain.qf[0,:,0,0]   
    q_loop[0,128,0,0]=chain.qf[0,0,0,0]   
    for i in range(1,Nt+1,1):
        nch=chain.blk_sp[chain.blk_s[i]]
        Ws=np.append(chem.W[nch,:,0,0],chem.W[nch,0,0,0])
        data[1,:]=data[1,:]-0.5*dt*Ws[:]*(2*dx**2/dt)
        D2 = sparse.spdiags(data,diags,Nx+1,Nx+1)/(dx**2)
        A = (I -dt/2*D2)
        b = ( I + dt/2*D2 )*u
        u = np.transpose(np.mat( sparse.linalg.spsolve( A,  b ) ))
        q_loop[i,:,0,0]=np.reshape(u[:,0],(Nx+1))
    print " q_loop[0,0,0,0] in crank",q_loop[0,0,0,0]



 
def Euler(Nx,Lx,Nt,dt,Ws,chem,chain,q_init,q_loop):
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
     
    # Create grid-points on x axis
    x = np.linspace(0,Lx,Nx+1)
     
    # Initial Conditions
    u = np.transpose(np.mat(q_init))
     
    # Second-Derivative Matrix
    data = np.ones((3, Nx+1))
    data[1] = -2*data[1]
    #data[1,:]=data[1,:]-0.5*dt*Ws[:]*(2*dx**2/dt)
    
    # Reflective boundary 
    #data[0,0]=0.0
    data[2,0]=2*data[2,0]
    data[0,0]=2*data[0,0]
    data[2,Nx]=2*data[2,Nx]
    data[0,Nx]=2*data[0,Nx]
    #data[2,Nx]=0.0
    diags = [-1,0,1]
    #D2 = sparse.spdiags(data,diags,Nx+1,Nx+1)/(dx**2)
     
    # Identity Matrix
    I = sparse.identity(Nx+1)
     
    # Data for each time-step
    #data = []
     
    # Solve the System: (I - k/2*D2) u_new = (I + k/2*D2)*u_old
    q_loop[0,0:128,0,0]=chain.qf[0,:,0,0]   
    q_loop[0,128,0,0]=chain.qf[0,0,0,0]   
    for i in range(1,Nt+1,1):
        nch=chain.blk_sp[chain.blk_s[i/10]]
        Ws=np.append(chem.W[nch,:,0,0],chem.W[nch,0,0,0])
        q_loop[i,0,0,0]=q_loop[i-1,0,0,0]+(dt/dx**2)*(q_loop[i-1,1,0,0]+q_loop[i-1,1,0,0]- \
        2*q_loop[i-1,0,0,0])-Ws[0]*dt*q_loop[i-1,0,0,0]
        for xi in range(1,128):
            q_loop[i,xi,0,0]=q_loop[i-1,xi,0,0]+(dt/dx**2)*(q_loop[i-1,xi-1,0,0]+q_loop[i-1,xi+1,0,0]- \
            2*q_loop[i-1,xi,0,0])-dt*Ws[xi]*q_loop[i-1,xi,0,0]
        q_loop[i,128,0,0]=q_loop[i-1,128,0,0]+(dt/dx**2)*(q_loop[i-1,127,0,0]+q_loop[i-1,127,0,0]- \
        2*q_loop[i-1,128,0,0])-dt*Ws[0]*q_loop[i-1,128,0,0]

    print " q_loop[2,0,0,0]",q_loop[2,0,0,0]



 
             








             








