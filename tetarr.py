#!/usr/bin/python
"""
Created on Thu May 23 00:00:25 2013
@author: joe
"""

import numpy as np

class Tetarr(np.ndarray):
    """"""
    def __new__(cls,L=100,dtype=int,data=None):
        """"""
        if data == None:
            input_array = np.zeros((47,L))
        elif type(data) == np.ndarray:
            input_array = data
        obj = np.asarray(input_array).view(cls)
        obj.germline = obj[:2]
        obj.somatic = obj[2:]
        obj.L = obj.size[1]
        return obj
        
    def __array_finalize__(self, obj):
        """"""
        if obj is None: return
        self.germline = getattr(obj, 'germline', None)
        self.somatic = getattr(obj, 'somatic', None)
        self.L = getattr(obj, 'L', None)
        
    def mutate(self,mu,value = -1):
        """"""
        mutants = (np.random.binomial(1,mu,self.shape) == 1) & (self == 0)
        self[mutants] = value
        
    def reproduce(self):
        """"""
        new_tet = Tetarr(data=self)
        dup_somatic = np.repeat(self.somatic,2,axis=0)
        for i in range(self.L):
            new_tet.somatic[:,i] = np.random.choice(dup_somatic[:,i],45,False)
        return new_tet
        