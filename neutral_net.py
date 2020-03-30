# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 17:00:50 2020

@author: cathe
"""

import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import cm
import random as rd
import copy
import csv

# Graph to hold neutral network
G = nx.Graph()
# List to hold all current strands in neutral network
ALL = []
ALL_base = set()
covid_base = 'LFQQN'
sars_base = 'YLNYT'
civet_base = 'YLKYS'
bat_base = 'SFNYN'


# Function authored by Amartya Ranjan Saikia taken from GeeksForGeeks
# https://www.geeksforgeeks.org/dna-protein-python-3/
def translate(seq):

    table = {
            'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
            'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
            'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
            'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
            'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
            'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
            'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
            'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
            'UAC': 'Y', 'UAU': 'Y', 'UAA': '_', 'UAG': '_',
            'UGC': 'C', 'UGU': 'C', 'UGA': '_', 'UGG': 'W',
            }
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += table[codon]
    return protein


# Modification of the above function.  Translates a string of amino acids
# to a possible sequence of codons
def translate_AA(seq):
    table = {
            'I': ['AUA', 'AUC', 'AUU'],
            'M': ['AUG'],
            'T': ['ACA', 'ACC', 'ACG', 'ACU'],
            'N': ['AAC', 'AAU'],
            'K': ['AAA', 'AAG'],
            'S': ['AGC', 'AGU', 'UCA', 'UCC', 'UCG', 'UCU'],
            'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGU'],
            'L': ['CUA', 'CUC', 'CUG', 'CUU', 'UUA', 'UUG'],
            'P': ['CCA', 'CCC', 'CCG', 'CCU'],
            'H': ['CAC', 'CAU'],
            'Q': ['CAA', 'CAG'],
            'V': ['GUA', 'GUC', 'GUG', 'GUU'],
            'A': ['GCA', 'GCC', 'GCG', 'GCU'],
            'D': ['GAC', 'GAU'],
            'E': ['GAA', 'GAG'],
            'G': ['GGA', 'GGC', 'GGG', 'GGU'],
            'F': ['UUC', 'UUU'],
            'Y': ['UAC', 'UAU'],
            'C': ['UGC', 'UGU'],
            '_': ['UAA', 'UAG', 'UGA'],
            'W': ['UGG'],
            }
    codon = ""
    rd.seed(2)
    for i in range(len(seq)):
        codon += rd.choice(table[seq[i]])
    return codon


def make_mutations(seq):
    seq_list = list(seq)
    mutations = []
    nucleotides = ['A', 'C', 'G', 'U']
    for i in range(len(seq_list)):
        for n in nucleotides:
            if not n == seq_list[i]:
                strand_copy = copy.copy(seq_list)
                strand_copy[i] = n
                strand = "".join(strand_copy)
                mutations.append(strand)
    return mutations


def add_mutations(seq, k=None, target_base=None):
    queue = []
    b = 0.15
    queue.append(seq)
    ALL.append(seq)
    ALL_base.add(translate(seq))
    mutation = 0
    target_count = 0
    while queue:
        mutation += 1
        strand = queue.pop(0)
        mutations = make_mutations(strand)
        for m in mutations:
            silent_mutation = (translate(m) == translate(strand))
            functional_variation = (rd.uniform(0, 1) < b)
            flag = silent_mutation
            if (not silent_mutation) and functional_variation:
                if '_' not in translate(m):
                    count = 0
                    for base in range(len(translate(m))):
                        if (translate(m))[base] == target_base[base]:
                            count += 1
                    if count >= target_count:
                        if count > target_count:
                            print(translate(m), count, mutation)
                            target_count = count
                        flag = True
            if flag and m not in ALL:
                queue.append(m)
                ALL.append(m)
                ALL_base.add(translate(m))
                G.add_edge(strand, m)
                if target_count == len(translate(m)):
                    distance = nx.shortest_path_length(G, source=seq, target=m)
                    print('Found target variation: ', distance)
                    queue.clear()
        if mutation == k:
            break
        if mutation % 5000 == 0:
#            cw = csv.writer(open('all_strands.csv', 'w'))
#            cw.writerow(['Mutations: '+str(mutation)])
#            cw.writerow(['Length of ALL: '+str(len(ALL))])
#            cw.writerow(['Length of queue: '+str(len(queue))])
#            cw.writerow([list(ALL_base)])
            print(mutation, len(queue))


def make_colors2():
    names = list(mcolors.CSS4_COLORS)
    possible_colors = rd.sample(names, len(ALL_base))
    colors = {}
    for index, strand in enumerate(ALL_base):
        colors[strand] = possible_colors[index]
    if len(colors.values()) != len(set(colors.values())):
        extra_colors = len(colors.values()) - len(set(colors.values()))
        print('Repeated colors: ', extra_colors)
    return colors

def make_colors():
    colormap = list((cm.get_cmap('tab20c', 80)).colors)
    possible_colors = rd.sample(colormap, len(ALL_base))
    colors = {}
    for index, strand in enumerate(ALL_base):
        colors[strand] = tuple(possible_colors[index][0:3])
    if len(colors.values()) != len(set(colors.values())):
        extra_colors = len(colors.values()) - len(set(colors.values()))
        print('Repeated colors: ', extra_colors)
    print(type(colors))
    return colors


def plot_graph(covid_seq):
    color_map = []
    colors = make_colors2()
    for node in G:
        if node == covid_seq:
            color_map.append('black')
        elif translate(node) == covid_base:
            color_map.append('white')
        elif translate(node) == sars_base:
            print('SARS BASE')
            color_map.append('black')
        elif translate(node) == civet_base:
            print('CIVET BASE')
            color_map.append('black')
        elif translate(node) == bat_base:
            print('BAT BASE')
            color_map.append('black')
        else:
            color = colors[translate(node)]
            color_map.append(color)
    nx.draw(G, node_color=color_map)
    ax = plt.gca()
    ax.collections[0].set_edgecolor('black')
    plt.show()


def main():
    covid_seq = translate_AA(covid_base)
    add_mutations(covid_seq, target_base=bat_base)
    print(len(ALL_base))
    #plot_graph(covid_seq)


main()
