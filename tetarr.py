#!/usr/bin/python
"""
Created on Thu May 23 00:00:25 2013
@author: joe
"""

import numpy as np

class Tetarr(np.ndarray):
    """"""
    def __new__(cls,L=100,data=None):
        """"""
        if data == None:
            input_array = np.zeros((47,L))
        elif hasattr(data,'shape'):
            input_array = data.copy()
        obj = np.asarray(input_array).view(cls)
        obj.germline = obj[:2]
        obj.somatic = obj[2:]
        obj.L = obj.shape[1]
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
        
    def fitness(self,s=0.2,h=0.5):
        """"""
        deleterious = (self.somatic == -1)
        beneficial = (self.somatic == 1)
        d = sum((np.mean(deleterious,axis=0)**(np.log(h)/np.log(0.5)))*s)
        b = sum((np.mean(beneficial,axis=0)**(np.log(h)/np.log(0.5)))*s)
        w = 1 - d + b
        return w
        
    def genomic_exclusion(self):
        """"""
        which_allele = np.random.binomial(1,0.5,self.L)
        #which_allele = (np.random.binomial(1,0.5,self.L) == 0)
        #hap_germ = (np.where(which_allele,self.germline[0,:],self.germline[1,:])).reshape((1,self.L))
        hap_germ = which_allele.choose([self[0,:],self[1,:]]).reshape((1,self.L))        
        homozygous = np.repeat(hap_germ,47,axis=0)
        homozygote = Tetarr(data=homozygous)
        return homozygote
        