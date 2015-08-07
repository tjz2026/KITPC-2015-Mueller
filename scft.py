# -*- coding: utf-8 -*-
import numpy as np
__all__ = ["scft_obj"]

class Scft_Obj(object):
    ''' A :class:`Grid` object represents an evenly Discretized sim_box cell in real and reciprocal space.
        :param dim: space dimension, can be 1, 2 or 3
        :type boxshape: :class:`BoxShape`
    '''

    def __init__(self,sim_system):

        self.f_tot=0.0
        self.f_FH=0.0
        self.f_logQ=0.0
        if sim_system!="ABn_melt":
            raise ValueError(' currently,only ABn melts scft driver is supported')
        else:
            print "set up scft obj named ABn_melt"   
        self.res=np.zeros(4)
         
          
