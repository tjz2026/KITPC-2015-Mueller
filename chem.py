# -*- coding: utf-8 -*-
"""
chemical.py
===========

"""
import numpy as np

__all__ = ["Chem"]

class Chem(object):

    def __init__(self,n_sp,grid,XN):

        self.num = n_sp
        if self.num!=2 and len[XN]!=1:
            raise ValueError(' currently,only two species are supported')
        else:
            self.XN=XN
            print "Flory-Huggins param:",XN
        self.__create_operator(grid) 


    def __create_operator(self,grid):
        '''
        ''' 
        
        self.W=np.zeros((self.num,len(grid.x),len(grid.y),len(grid.z)))
        self.R=np.zeros((self.num,len(grid.x),len(grid.y),len(grid.z)))
        
class Loop(object):

    def __init__(self,n_sp,grid,XN):

        self.num = n_sp
        if self.num!=2 and len[XN]!=1:
            raise ValueError(' currently,only two species are supported')
        else:
            self.XN=XN
            print "Flory-Huggins param:",XN
        self.__create_operator(grid) 


    def __create_operator(self,grid):
        '''

        ''' 
        
        self.W=np.zeros((self.num,len(grid.x),len(grid.y),len(grid.z)))
        self.R=np.zeros((self.num,len(grid.x),len(grid.y),len(grid.z)))
        






















