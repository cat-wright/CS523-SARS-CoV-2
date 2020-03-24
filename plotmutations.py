# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 13:34:27 2020

@author: cathe
"""

## To calculate the number of genomes each mutation away
from scipy import special
import matplotlib.pyplot as plt

length_of_seq = 15

total_combinations = 4**length_of_seq

def calc_mutations(k):  # kth mutation
    choose = special.comb(length_of_seq, k)
    return choose * 3**k

def plot_mutations():
    mutations = []
    for muts in range(length_of_seq):
        mutations.append(calc_mutations(muts+1))
    
    plt.plot(list(range(1, length_of_seq+1)), mutations, 'bo', markersize=12)
    plt.yscale('log')
    plt.xlabel('mutations away')
    plt.ylabel('number of genomes (log-scale)')
    plt.xticks(list(range(1, length_of_seq+1)))
    plt.rc('font', size=25)
    plt.grid()
    plt.show()
    
plot_mutations()