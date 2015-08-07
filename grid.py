# -*- coding: utf-8 -*-
"""
gyroid.grid
===========

"""

import numpy as np
from numpy.linalg import inv


__all__ = ["Grid"]

class Grid(object):
    ''' A :class:`Grid` object represents an evenly Discretized sim_box cell in real and reciprocal space.
        :param dim: space dimension, can be 1, 2 or 3
        :type dim: integer
        :param box
        :type box: class :'SimBox'
        :param shape: the shape of a unit cell
        :type boxshape: :class:`BoxShape`
        :grid_parameter (the grid number on each axis of the simulation box).
        Note: Grid.mesh is always constructed as 3D numpy array using meshgrid method

    '''

    def __init__(self,box,boxshape,grid_param=None):

        self.dim = box.dim
        if self.dim == 1 and np.size(grid_param)==1:
            self.x=np.linspace(0.0,boxshape.ll[0],grid_param[0],endpoint=False)
            self.y=np.array([0.0])
            self.z=np.array([0.0])
            self.n_grid=len(self.x) 
        elif self.dim == 2 and np.size(grid_param)==2:
            self.x=np.linspace(0.0,boxshape.ll[0],grid_param[0],endpoint=False)
            self.y=np.linspace(0.0,boxshape.ll[1],grid_param[1], endpoint=False)
            self.z=np.array([0.0])
            self.n_grid=len(self.x)*len(self.y) 
        elif self.dim == 3 and np.size(grid_param)==3:
            self.x=np.linspace(0.0,boxshape.ll[0],grid_param[0],endpoint=False)
            self.y=np.linspace(0.0,boxshape.ll[1],grid_param[1],endpoint=False)
            self.z=np.linspace(0.0,boxshape.ll[2],grid_param[2],endpoint=False)
            self.n_grid=len(self.x)*len(self.y)*len(self.z) 
        else:
            raise ValueError(' dimension mismatch for SimBox and grid_param.')
        #print "xgrid info",self.x
        #print "ygrid info",self.y
        #print "zgrid info",self.z

        self.__create_kspace_grid(box) 

    def __create_kspace_grid(self,box):
        '''
        :var kx,ky,kz: mesh grid on the fourier space on kx/ky/kz axis.note that the last half frequency of each grid is negtive, 
         as is addressed in FFTW.

        '''
        d_kx=2*np.pi/box.shape.ll[0]
        half_kx=len(self.x)/2 
        self.kx=np.arange(len(self.x),dtype=float)
        self.kx[:half_kx+1]=self.kx[:half_kx+1]*d_kx
        self.kx[half_kx+1:]=d_kx*(len(self.x)-self.kx[half_kx+1:])
        if self.dim==2 :
            d_ky=2*np.pi/box.shape.ll[1]
            half_ky=len(self.y)/2 
            self.ky=np.arange(len(self.y),dtype=float)
            self.ky[:half_kx+1]=self.ky[:half_ky+1]*d_ky
            self.ky[half_ky+1:]=d_ky*(len(self.y)-self.ky[half_ky+1:])
         
        if self.dim==3 :
            d_ky=2*np.pi/box.shape.ll[1]
            half_ky=len(self.y)/2 
            self.ky=np.arange(len(self.y),dtype=float)
            self.ky[:half_kx+1]=self.ky[:half_ky+1]*d_ky
            self.ky[half_ky+1:]=d_ky*(len(self.y)-self.ky[half_ky+1:])

            d_kz=2*np.pi/box.shape.ll[2]
            half_kz=len(self.z)/2 
            self.kz=np.arange(len(self.z),dtype=float)
            self.kz[:half_kz+1]=self.kz[:half_kz+1]*d_kz
            self.kz[half_kz+1:]=d_kz*(len(self.z)-self.kz[half_kz+1:])

















