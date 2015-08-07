# -*- coding: utf-8 -*-

import numpy as np


__all__ = ["Scft_IO"]

IO_TYPES = ["Read_Param","Read_Field","Dump_Field","Dump_Result"]

class Scft_IO(object):
    """ Representation of a whole SABF set.

    """
    def __init__(self,io_type,chem,filename=None):
           
        self.__IO(io_type,chem=chem)
    
    def __IO(self,io_type,chem,io_param=None):
        if len(io_type) == 0:
            print "empty io_type"
            return
        else:
            for x in io_type:
                if x in IO_TYPES:
                    if x == "Read_Param":
                        scft_obj.param=np.loadtxt('param.dat').T
                    elif x == "Read_Field":
                       chem.W=np.loadtxt('fields.dat').T
                    elif x== "Dump_Field":
                       #np.savetxt('fields.dat',chem.W.T,fmt='%3.6e')    
                       #np.savetxt('densities.dat',chem.R.T,fmt='%3.6e')    
                       np.savetxt('fields.dat',chem.W.T)    
                       np.savetxt('densities.dat',chem.R.T)    
                    elif x== "Dump_Result":
                       np.savetxt('result.dat',scft_obj.res)    
                else:
                    raise ValueError('Unkonwn IO type in io_type list.')
                
           
  
