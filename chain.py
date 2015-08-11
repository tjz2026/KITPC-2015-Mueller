# -*- coding: utf-8 -*-
"""
chain.py
===========

"""

import numpy as np


__all__ = ["Chain"]

class Chain(object):
    ''' A :class:`Grid` object represents an evenly Discretized sim_box cell in real and reciprocal space.

    '''

    def __init__(self,n_blk,Ns,chain_squence,N_interpo,grid):

        self.n_blk=n_blk
        self.Ns=Ns
        self.blk_ns=Ns/n_blk
        self.topo=chain_squence
        if chain_squence == "ABn":
            self.blk_f=np.arange(n_blk,dtype=float)
            self.blk_f[:]=1.0/n_blk
            self.f_sp=np.zeros(2)
            self.f_sp[0]=0.5
            self.f_sp[1]=0.5
            self.blk_sp=np.arange(n_blk)
            self.blk_s=np.arange(Ns+1)
            self.blk_sp[0:n_blk:2]=0
            self.blk_sp[1:n_blk+1:2]=1
            self.blk_sta=np.arange(n_blk)
            self.blk_end=np.arange(n_blk)
            self.intpo=N_interpo
             
        else:
            raise ValueError(' only ABn are considered currently')

        for i in np.arange(n_blk):
            self.blk_sta[i]=i*self.blk_ns
            self.blk_end[i]=(i+1)*self.blk_ns-1
            print "block start:end",self.blk_sta[i],self.blk_end[i]
        self.blk_end[n_blk-1]=self.Ns
        print "correction for last block end",self.blk_end[n_blk-1]
    
        for s in np.arange(0,Ns+1,1):
            for i in np.arange(n_blk):
                if s>=self.blk_sta[i] and s<=self.blk_end[i]:
                    self.blk_s[s]=i
 
        self.__create_propagator(grid) 


    def __create_propagator(self,grid):
        '''

        '''
        self.qf=np.zeros((self.Ns+1,len(grid.x),len(grid.y),len(grid.z)))
        self.qb=np.zeros((self.Ns+1,len(grid.x),len(grid.y),len(grid.z)))
        self.qf1=np.zeros((self.Ns+1,len(grid.x),len(grid.y),len(grid.z)))
        self.qb1=np.zeros((self.Ns+1,len(grid.x),len(grid.y),len(grid.z)))
        if self.topo == "ABn":
            self.qloop=np.zeros(((self.Ns+1)*self.intpo,len(grid.x),len(grid.y),len(grid.z)))
        self.Q=0.0


