#!/usr/bin/python
"""
Created on Sat May 25 13:41:10 2013
@author: joe
"""

from poplist import *
from pandas import Series,DataFrame
import pandas as pd

class Evolution:
    """"""
    def __init__(self,N=10,L=100,gens=1000,mu_d=0.000025,mu_b=0.0,s=0.2,h=0.5):
        """"""
        self.population = Population()
        self.population.populate(N,L)
        self.generations = gens
        self.mu_d = mu_d
        self.mu_b = mu_b
        self.s = s
        self.h = h
        self.m_somatic_fitness = dict()
        self.v_somatic_fitness = dict()
        self.m_germline_fitness = dict()
        self.v_germline_fitness = dict()
        
    def evolve(self):
        """"""
        while self.population.generation < self.generations:
            print self.population.generation
            w = self.population.fitness(self.s,self.h)
            self.m_somatic_fitness[self.population.generation] = np.mean(w)
            self.v_somatic_fitness[self.population.generation] = np.var(w)
            gep = self.population.ge()
            w_gep = gep.fitness(self.s,self.h)
            self.m_germline_fitness[self.population.generation] = np.mean(w_gep)
            self.v_germline_fitness[self.population.generation] = np.var(w_gep)
            self.population = self.population.next_generation(self.mu_d,self.mu_b,self.s,self.h)
        self.m_somatic_fitness = Series(self.m_somatic_fitness)
        self.v_somatic_fitness = Series(self.v_somatic_fitness)
        self.m_germline_fitness = Series(self.m_germline_fitness)
        self.v_germline_fitness = Series(self.v_germline_fitness)
                    
        