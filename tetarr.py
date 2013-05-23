#!/usr/bin/python
"""
Created on Thu May 23 00:00:25 2013
@author: joe
"""

import numpy as np

class Tetarr(np.ndarray):
    """numpy.ndarray subclass for representing individual organisms. A 47xL
    array is used to represent the genomes. The first two rows of the array
    correspond to the germline genome and the last 45 rows correspond to the
    somatic genome. Each column corresponds to a locus, and each element in the
    array specifies an allele (wild-type=0, deleterious mutation= -1,
    beneficial mutation=1)"""
    def __new__(cls,L=100,data=None):
        """Uses view casting of np.ndarray to make a Tetarr object.
        L - number of loci
        data - (should be) an array with shape=(47,L)
        If data is None, makes an array of zeros with shape=(47,L)
            """
        if data == None:
            input_array = np.zeros((47,L),dtype=np.int8)
        elif hasattr(data,'shape'):
            input_array = data.copy()
        obj = np.asarray(input_array,dtype=np.int8).view(cls)
        obj.germline = obj[:2]
        obj.somatic = obj[2:]
        obj.L = obj.shape[1]
        return obj
        
    def __array_finalize__(self, obj):
        """Attributes:
            L - number of loci (int)
            germline - two rows of genomes array, corresponds to germline
            somatic - 45 rows of genomes array, corresponds to somatic"""
        if obj is None: return
        self.germline = getattr(obj, 'germline', None)
        self.somatic = getattr(obj, 'somatic', None)
        self.L = getattr(obj, 'L', None)
        
    def mutate(self,mu,value = -1):
        """Mutates each allele at each locus in each genome to the given value
        with probability mu. Uses numpy.random.binomial to determine which
        alleles mutate, then uses a boolean array to assign the value those
        alleles. Right now, mutations are irreversible...
        (currently using value = 0 for wild-type, value = -1 for deleterious
        mutations, and value = 1 for beneficial mutations)."""
        mutants = (np.random.binomial(1,mu,self.shape) == 1) & (self == 0)
        self[mutants] = value
        
    def reproduce(self):
        """Returns a new Tetarr object. The germline genome of the new
        individual is copied from the progenater. To form the new somatic
        genome, the somatic genome of the progenater is duplicated and each
        locus is sampled without replacement 45 times via np.random.choice."""
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
        