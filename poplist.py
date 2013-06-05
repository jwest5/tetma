#!/usr/bin/python
"""
Created on Wed May 22 23:13:56 2013
@author: joe
"""

from tet import *
import numpy as np
import pickle

class Population(list):
    """Class used to contain the individuals in a population"""
    def __init__(self):
        """Initializes an empty Population list. N is population size,
        and generation keeps track of the generation."""
        list.__init__(self)
        self.N = len(self)
        self.generation = 0
        self.somatic_fitness = None
        self.germline_fitness = None
        
    def append(self,value):
        """Modified append so population size (N) updates when individuals
        are added."""
        list.append(self,value)
        self.N = len(self)

    def populate(self,N=10,L=100):
        """Creates the initial population of N individuals each with L loci."""
        while len(self) < N:
            self.append(Tet(L))
        self.N = len(self)
            
    def reproduction(self):
        """Applies Tet.reproduce() to each member of the population."""
        for i in range(self.N):
            self[i] = self[i].reproduce()
        self.generation += 1
            
    def mutation(self,mu,value=-1):
        """Applies Tet.mutate(mu,value) to each member of the population."""
        for i in range(self.N):
            self[i].mutate(mu,value)
            
    def fitness(self,s=0.2,h=0.5):
        """Returns a list containing the result of Tet.fitness() for each
        member of the population."""
        w = []
        for i in range(self.N):
            w.append(self[i].fitness(s,h))
        return w
        
    def ge(self):
        """Returns a new Population object containing the result of
        Tet.genomic_exclusion() for each member of the population."""
        gep = Population()
        for i in range(self.N):
            gep.append(self[i].genomic_exclusion())
        return gep
        
    def next_generation(self,mu_d=0.000025,mu_b=0.0,s=0.2,h=0.5):
        """"""
        self.reproduction()
        self.mutation(mu_d)
        self.mutation(mu_b,1)
        w = self.fitness(s,h)
        print w
        pw = []
        for i in range(self.N):
            pw.append(w[i]/sum(w))
        new_gen = Population()
        new = np.random.choice(range(self.N),self.N,replace=True,p=pw)
        for i in range(self.N):
            new_gen.append(self[new[i]])
        new_gen.generation = self.generation
        return new_gen
        
    def pickle_save(self,outfile):
        """"""
        stream = open(outfile,'a')
        pickle.dump(self,stream)
        stream.close()
        
    def save_fitness(self,outfile,s=0.2,h=0.5):
        """"""
        fitness = self.fitness(s,h)
        stream = open(outfile,'a')
        pickle.dump(fitness,stream)
        stream.close()
        
        
        
