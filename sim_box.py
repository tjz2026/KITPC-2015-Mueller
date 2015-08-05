# -*- coding: utf-8 -*-
"""
multiBlock.sim_box
===============

Define a the simulation box infomation when doing a non unitcell calculation.
Currently, only a standard XYZ cartesian system is supported, ie. X, Y,Z axises are all perpendicular to 
each other.
Note that for the last dimension (in numpy, the last dimension varies fastest!!!) is always picked as X .

"""

import numpy as np


__all__ = ["SimBox","BoxShape"]

class SimBox(object):
    ''' A :class:`SimBox` object contains the shape information of a nonunit.
    cell.

    1D, 2D and 3D spaces are supported.

    To construct a SimBox of, one should provide
    the space dimension, and the cell
    parameters (the lengths on each axis of the simulation box.

    '''

    def __init__(self,dim,cell_param=None):

        self.dim = dim
        if dim == 1:
            self.__standard_cell_1D(cell_param)
        elif dim == 2:
            self.__standard_cell_2D(cell_param)
        elif dim == 3:
            self.__standard_cell_3D(cell_param)
        else:
            raise ValueError('Unkonwn dimension for SimBox.')

        self.shape = self.__create_shape()

    def __create_shape(self):
        '''
        :var x: length of SimBox on X axis
        :var y: length of SimBox on Y axis
        :var z: length of SimBox on Z axis
        :retrun: the shape of the unit cell
        :rtype: :class:`Shape`

        '''

        if self.dim == 1:
            return BoxShape(1,np.array([self.x]))

        if self.dim == 2:
            x = np.array([self.x, 0.0])
            y = np.array([0.0,self.y])
            return BoxShape(2,x,y)

        if self.dim == 3:
            x = np.array([self.x, 0.0, 0.0])
            y = np.array([0.0, self.y, 0.0])
            z = np.array([0.0, 0.0, self.z])
            return BoxShape(3,x,y,z)

    def __standard_cell_1D(self,cp):
        if np.size(cp) != 1:
            raise ValueError('1D SimBox requires 1 parameters.')
        self.x = cp[0]

    def __standard_cell_2D(self,cp):
        if np.size(cp) != 2:
            raise ValueError('2D SimBox requires 2 parameters.')
        self.x = cp[0]
        self.y = cp[1]

    def __standard_cell_1D(self,cp):
        if np.size(cp) != 3:
            raise ValueError('3D SimBox requires 3 parameters.')
        self.x = cp[0]
        self.y = cp[1]
        self.z = cp[2]


class BoxShape(object):
    ''' A :class:`BoxShape` object fully describes the shape of a SimBox.

    The heart of a :class:`BoxShape` object is a shape matrix which can be
    constructed from unit vectors of a unit cell in Cartesian Coordinate.

    The Morse convention is used to construct the shape matrix. That is each row
    in the shape matrix represents a unit vector, e.g.::
       h = [X 0 0
            0 Y 0
            0 0 Z]

    '''

    def __init__(self,dim,a1,a2=None,a3=None):
        self.dim = dim
        if dim==3 and np.size(a1)==3 and np.size(a2)==3 and np.size(a3)==3:
            self.m = np.array([a1,a2,a3])
        elif dim==2 and np.size(a1)==2 and np.size(a2)==2:
            self.m = np.array([a1,a2])
        elif dim==1 and np.size(a1)==1:
            self.m = np.array([a1])
        else:
            raise ValueError('Dimension and the number'
                             'of Box cell vector not match')
        self.gg = 2.0 * np.pi * inv(self.m).T
        self.ll = np.sqrt(np.sum(self.m**2,axis=1))


    @property
    def l(self):
        ''' (a,b,c), the length vectors of each dimension. '''
        return self.ll

    @property
    def h(self):
        ''' The shape matrix in real space. '''

        return self.m

    @property
    def g(self):
        ''' The shape matrix in reciprocal space. '''

        return self.gg


