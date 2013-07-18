#!/usr/bin/python
"""
Created on Fri Jun 21 11:55:54 2013
@author: joe
"""

import numpy as np
from numpy.random import gamma
from numpy import inf

def mutational_effects(n, s, beta):
    """
    Generate n mutations from a gamma distribution with mean effect s and shape parameter beta.
    Negative (positive) values of s indicate deleterious (beneficial) mutations. 
    
    Arguments:
    n -- number of mutations
    s -- mean effect of a mutation
    beta -- shape parameter of the gamma distribution (inf indicates equal effects)
    """
    if sign(s) == 1:
        beneficial = True
    elif sign(s) == -1:
        beneficial = False
    else:
        return "Invalid s: must be nonzero."
    if beta > 0:
        if beta == inf:
            mutations = np.repeat(s, n)            
        else:
            alpha = beta / abs(s)
            if beneficial:
                mutations = gamma(shape=beta, scale=1/alpha, size=n)
            else:
                mutations = - gamma(shape=beta, scale=1/alpha, size=n)
        return mutations
    else:
        return "Invalid beta: must be 0 < beta < inf."

class Tet:
    
    """
    Represent individual organism.
    
    Attributes:
    L - Number of loci.
    germline - 2xL array representing the diploid germline genome.
    somatic - 45XL array representing the somatic genome.
    
    """
    
    def __init__(self, L=100, germline=None, somatic=None):
        """
        Initialize a Tet object.
        
        Arguments:
        L - Number of loci (default is 100).
        germline - An array with shape (2,L) or None (default is None).
            If germline argument is None, then a 2xL array of zeros.
        somatic - An array with shape (45,L) or None (default is None).
            If somatic argument is None, then a 45xL array of zeros.
        
        NOTE: Doesn't check germline/somatic type or shape. FIX.
        """
        self.L = L
        if germline is None:
            self.germline = np.zeros((2,self.L))
        else:
            self.germline = germline
        if somatic is None:
            self.somatic = np.zeros((45,self.L))
        else:
            self.somatic = somatic
            
    def __repr__(self):
        """Return the string representation of each genome."""
        s = "Germline:\n" + str(self.germline) + "\nSomatic:\n" + str(self.somatic)
        return s
        
    def mutate_germline(self,mu,value):
        """Mutates each allele at each locus in the germline to the given value
        with probability mu. Uses numpy.random.binomial to determine which
        alleles mutate, then uses a boolean array to assign the value those
        alleles. Right now, mutations are irreversible...
        (currently using value = 0 for wild-type, value = -1 for deleterious
        mutations, and value = 1 for beneficial mutations)."""
        mutants = (np.random.binomial(1,mu,(47,self.L)) == 1) & (self.germline == 0)
        self.germline[mutants] = value
        
    def mutate_somatic(self,mu,value):
        """Mutates each allele at each locus in the somatic to the given value
        with probability mu. Uses numpy.random.binomial to determine which
        alleles mutate, then uses a boolean array to assign the value those
        alleles. Right now, mutations are irreversible...
        (currently using value = 0 for wild-type, value = -1 for deleterious
        mutations, and value = 1 for beneficial mutations)."""
        mutants = (np.random.binomial(1,mu,(47,self.L)) == 1) & (self.somatic == 0)
        self.somatic[mutants] = value
        
    def reproduce(self):
        """Returns a new Tet object. The germline genome of the new individual
        is copied from the progenater. To form the new somatic genome, the
        somatic genome of the progenater is duplicated and each locus is
        sampled without replacement 45 times using numpy.random.choice."""
        dup_somatic = np.zeros((90,self.L))
        dup_somatic[:45] = self.somatic.copy()
        dup_somatic[45:] = self.somatic.copy()
        new_tet = Tet(self.L)
        for i in range(self.L):
            new_tet.somatic[:,i] = np.random.choice(dup_somatic[:,i],45,False)
        new_tet.germline = self.germline
        return new_tet
        
    def fitness(self,s,h):
        """Calculates additive fitness..."""
        deleterious = (self.somatic == -1)
        beneficial = (self.somatic == 1)
        d = sum((np.mean(deleterious,axis=0)**(np.log(h)/np.log(0.5)))*s)
        b = sum((np.mean(beneficial,axis=0)**(np.log(h)/np.log(0.5)))*s)
        w = 1 - d + b
        return w
        
    def genomic_exclusion(self):
        """"""
        
    
        
        
        
        
        