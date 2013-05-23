#!/usr/bin/python
"""
Created on Wed May 22 23:13:56 2013
@author: joe
"""

from tet import *
import numpy as np

class Population(list):
    """"""
    def __init__(self):
        """"""
        list.__init__(self)
        self.N = len(self)
        self.generation = 0
        
    def append(self,value):
        """"""
        list.append(self,value)
        self.N = len(self)

    def populate(self,N=10,L=100):
        """"""
        while len(self) < N:
            self.append(Tet(L))
        self.N = len(self)
            
    def reproduction(self):
        """"""
        for i in range(self.N):
            self[i] = self[i].reproduce()
        self.generation += 1
            
    def mutation(self,mu):
        """"""
        for i in range(self.N):
            self[i].mutate(mu)
            
    def fitness(self,s=0.2,h=0.5):
        """"""
        w = []
        for i in range(self.N):
            w.append(self[i].fitness(s,h))
        return w
        
    def next_generation(self,mu,s=0.2,h=0.5):
        """"""
        self.reproduction()
        self.mutation(mu)
        w = self.fitness(s,h)
        pw = []
        for i in range(self.N):
            pw.append(w[i]/sum(w))
        new_gen = Population()
        new = np.random.choice(self,self.N,replace=True,p=pw)
        for i in range(self.N):
            new_gen.append(new[i])
        new_gen.generation = self.generation
        return new_gen
        
    def ge(self):
        """"""
        gep = Population()
        for i in range(self.N):
            gep.append(self[i].genomic_exclusion())
        return gep
        
    def germline_fitness(self,s=0.2,h=0.5):
        """"""
        gep = self.ge()
        return gep.fitness(s,h)
        
        
