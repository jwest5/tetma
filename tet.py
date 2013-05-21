#!/usr/bin/python
"""
Created on Tue May 21 12:54:26 2013
@author: joe
"""

import numpy as np

class Tet:
    """Class for representing individual organisms. A 47xL array is used to
    represent the genomes. The first two rows of the array correspond to the
    germline genome and the last 45 rows correspond to the somatic genome.
    Each column corresponds to a locus, and each element in the array specifies
    an allele (wild-type=0, deleterious mutation= -1, beneficial mutation=1)"""
    def __init__(self,L=100,genomes=None):
        """Initializes a Tet object.
        attributes:
            L - number of loci (int)
            genomes - numpy array (if none is provided, fills with zeros
            germline - two rows of genomes array, corresponds to germline
            somatic - 45 rows of genomes array, corresponds to somatic"""
        self.L = L
        if genomes == None:
            self.genomes = np.zeros((47,self.L))
        else:
            self.genomes = genomes
        self.germline = self.get_germline()
        self.somatic = self.get_somatic()
        
    def get_germline(self):
        """Returns the diploid germline genome"""
        return self.genomes[:2]
        
    def set_germline(self,germline):
        """Sets the diploid germline genome to the array provided."""
        self.genomes[:2] = germline
        
    def get_somatic(self):
        """Returns the 45-ploid somatic genome."""
        return self.genomes[2:]
        
    def set_somatic(self,somatic):
        """Sets the 45-ploid somatic genome to the array provided."""
        self.genomes[2:] = somatic
        
    def mutate(self,mu,value= -1):
        """Mutates each allele at each locus in each genome to the given value
        with probability mu. Uses numpy.random.binomial to determine which
        alleles mutate, then uses a boolean array to assign the value those
        alleles. Right now, mutations are irreversible...
        (currently using value = 0 for wild-type, value = -1 for deleterious
        mutations, and value = 1 for beneficial mutations)."""
        mutants = (np.random.binomial(1,mu,(47,self.L)) == 1) & (self.genomes == 0)
        self.genomes[mutants] = value
        
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
        new_tet.set_germline(self.germline)
        return new_tet

