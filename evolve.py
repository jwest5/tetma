#!/usr/bin/python
"""
Created on Tue Jun  4 18:12:38 2013
@author: joe
"""

from poplist import *
import pandas
from pandas import Series,DataFrame

# add runabunch function, readabunch function, and DA

class Evo:
    """"""
    def __init__(self,popfile,somafile,germfile,L=100,mu_d=0.000025,mu_b=0.0,mu_effect=0.02,N=10,gens=1000):
        """"""
        self.popfile = popfile
        self.somafile = somafile
        self.germfile = germfile
        self.L = L
        self.mu_d = mu_d
        self.mu_b = mu_b
        self.mu_effect = mu_effect
        self.N = N
        self.generations = gens
        self.population = Population()
        self.population.populate(self.N,self.L)
        self.somatic_fitness = []
        self.germline_fitness = []
        
    def evolve(self):
        """"""
        #self.population.pickle_save(self.popfile)
        self.somatic_fitness.append(self.population.fitness(self.mu_effect,h=0.5))
        self.germline_fitness.append(self.population.fitness(self.mu_effect,h=0.5))
        while self.population.generation <= self.generations:
            next_gen = self.population.next_generation(mu_d=self.mu_d,mu_b=self.mu_b,s=self.mu_effect,h=0.5)
            self.somatic_fitness.append(next_gen.fitness(self.mu_effect,h=0.5))
            next_gen_ge = next_gen.ge()
            self.germline_fitness.append(next_gen_ge.fitness(self.mu_effect,h=0.5))
            #next_gen.pickle_save(self.popfile)
            self.population = next_gen
            print self.population.generation
        sfitness = DataFrame(self.somatic_fitness)  #change to series
        gfitness = DataFrame(self.germline_fitness)
        sfitness.to_csv(self.somafile)
        gfitness.to_csv(self.germfile)
        
        

