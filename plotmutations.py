# -*- coding: utf-8 -*-
"""
CS523 Project 2
Part 1b: "Calculate the number of genomes 1, 2, 3, all the way up to
15 mutations away from a length 15 genome.
Show these in some sort of figure. A log scale may be helpful."
"""

# special.comb performs binomial n choose k
from scipy import special
import matplotlib.pyplot as plt

length_of_seq = 15
total_combinations = 4**length_of_seq


# Calculate the number of genomes at the kth mutation away
def calc_mutations(k):  # kth mutation
    choose = special.comb(length_of_seq, k)
    return choose * 3**k


# Creates an array storing the genomes at each mutation away. (y-values)
def mutation_list():
    mutations = []
    for muts in range(length_of_seq):
        mutations.append(calc_mutations(muts+1))
    return mutations


# Plots the number of mutations in a y-log scale
def plot_mutations():
    mutations = mutation_list()
    plt.plot(list(range(1, length_of_seq+1)), mutations, 'bo', markersize=12)
    plt.yscale('log')
    plt.xlabel('mutations away from root genome')
    plt.ylabel('number of genomes (log-scale)')
    plt.xticks(list(range(1, length_of_seq+1)))
    plt.rc('font', size=25)
    plt.grid()
    plt.show()


# Checks to ensure the sum of all mutations is equal to the total possible
# mutations.
def check_sum():
    mutations = mutation_list()
    assert sum(mutations) == total_combinations - 1  # minus 1 for the root


plot_mutations()
# check_sum()
