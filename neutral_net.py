# -*- coding: utf-8 -*-
"""
CS523 Project 2
Neutral Network implementation used in part 1.
main() function can be used to generate all values/plots seen in the report.

def main(task, b_val=0.0, k=None, target_base=None, iterations=40):
Arguments:
    task -
        'neutral clusters': all other parameters can be default.
                    Creates and plots the clusters of neutral networks that
                    form the COVID-19 5 critical residues base (LFQSN)
        'neutral clusters 1e': all other parameters can be default.
                    Used for generating results to part 1e.
                    Creates and plots the clusters of neutral networks that
                    form 5 chosen spike protein amino acids in pangolin
                    SARS-CoV and an initial COVID-19 case in Wuhan.
        'mutations': Specify b_val and the number of iterations k.
                    Used to create Figure 3.  Generates a random possible base
                    sequence for covid_base.  Allows functional variations
                    with percentage b and plots the network after k iterations.
                    For figure 3 we used b = 0.015 (1.5%) and k=40.
        'target base': Specify b_val and the target_base [bat_base, civet_base
                    , sars_base (human)].  Generates a random possible base
                    sequence for covid_base and performs "productive searching"
                    to find the target base.  Returns the distance from the
                    root node to the target base.
        'target base with iterations': Specify b_val, target_base, and
                    iterations. Used to calculate data for table 3.
                    Runs the same as 'target base' for specified
                    iterations, each time appending the minimum distance to
                    a csv file.
"""

import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import random as rd
import copy
import math
from part1e import return_diffs

# List to hold all current strands in neutral network
ALL = []
ALL_base = set()

covid_base = 'LFQSN'
sars_base = 'YLNDT'  # Human SARS-CoV
civet_base = 'YLKDS'
bat_base = 'SFNDN'

# PART 1e
pangolin_base = 'STESI'
wuhan_base = 'AADNV'


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
def translate_AA(seq, all_possibilities=False, get_table=False):
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
    if all_possibilities:  # Only works for 5 AA seqences
        list1 = table[seq[0]]
        list2 = table[seq[1]]
        list3 = table[seq[2]]
        list4 = table[seq[3]]
        list5 = table[seq[4]]
        codons = []
        [codons.append(i+j+k+l+m) for i in list1
         for j in list2
         for k in list3
         for l in list4
         for m in list5]
        return codons
    elif get_table:
        return table
    else:
        codon = ""
        # rd.seed(2)
        for i in range(len(seq)):
            codon += rd.choice(table[seq[i]])
        return codon


def get_size_of_network(nt_seq):
    aa_seq = translate(nt_seq)
    table = translate_AA(aa_seq, get_table=True)
    network_size = 0  # 1
    for i in range(len(aa_seq)):
        if aa_seq[i] == 'S':
            nt = nt_seq[3*i:3*i+3]
            if nt in ['AGC', 'AGU']:
                network_size += math.log10(2)
            else:
                network_size += math.log10(4)
        else:
            network_size += math.log10(len(table[aa_seq[i]]))
    return network_size  # Returned as the log base 10 of the network size


# Returns all possible genomes 1 mutation away from input genome "seq"
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


# Creates and returns the NetworkX graph for our neutral network.
def add_mutations(seq, b_val, k=None):
    G = nx.Graph()
    # Used similar to the queue in BFS.  Holds "unvisited" nodes
    queue = []
    queue.append(seq)
    ALL.append(seq)
    ALL_base.add(translate(seq))
    mutation = 0
    # Iterates until the queue is empty (no more nodes to search) or the kth
    # iteration.
    while queue:
        mutation += 1
        # Pulls random node from list to mutate
        strand = queue.pop(rd.randrange(len(queue)))
        # Pulls first node from list to mutate (BFS)
        # strand = queue.pop(0)
        mutations = make_mutations(strand)
        for m in mutations:
            silent_mutation = (translate(m) == translate(strand))
            functional_variation = (rd.uniform(0, 1) < b_val)
            flag = silent_mutation
            if (not silent_mutation) and functional_variation:
                if '_' not in translate(m):
                    flag = True
            if flag and m not in ALL:
                queue.append(m)
                ALL.append(m)
                ALL_base.add(translate(m))
                G.add_edge(strand, m)
        if mutation == k:
            print('Break')
            return G
    return G


# Very similar to add_mutations function.  Searches for a target base, using
# our "productive searching" approach.
def add_mutations_target(seq, target_base, b_val):
    G = nx.Graph()
    queue = []
    queue.append(seq)
    ALL.append(seq)
    ALL_base.add(translate(seq))
    mutation = 0
    # Fitness function: number of amino acids correctly mutated to target base
    target_count = 0
    while queue:
        mutation += 1
        # Pulls random node from list to mutate
        strand = queue.pop(rd.randrange(len(queue)))
        # Pulls first node from list to mutate (BFS)
        # strand = queue.pop(0)
        mutations = make_mutations(strand)
        for m in mutations:
            silent_mutation = (translate(m) == translate(strand))
            functional_variation = (rd.uniform(0, 1) < b_val)
            flag = silent_mutation
            if (not silent_mutation) and functional_variation:
                if '_' not in translate(m):
                    # Fitness function for mutation m
                    count = 0
                    for base in range(len(translate(m))):
                        if (translate(m))[base] == target_base[base]:
                            count += 1
                    if count >= target_count:
                        # New maximum fitness
                        if count > target_count:
                            # print(translate(m), count, mutation)
                            target_count = count
                        # Will only append to queue if mutation is of equal
                        # or better fitness than maximum fitness seen
                        flag = True
            if flag and m not in ALL:
                queue.append(m)
                ALL.append(m)
                ALL_base.add(translate(m))
                G.add_edge(strand, m)
                # If fitness function = 5, target base has been found.
                # Returns the min-hop distance from the root to the target
                if target_count == len(translate(m)):
                    distance = nx.shortest_path_length(G, source=seq, target=m)
                    return distance
    return 0


# Generates the colors for plotting networks.  Based on how many different
# amino acid sequences are in the network.
def make_colors():
    names = list(mcolors.CSS4_COLORS)
    possible_colors = rd.sample(names, len(ALL_base))
    colors = {}
    for index, strand in enumerate(ALL_base):
        colors[strand] = possible_colors[index]
    if len(colors.values()) != len(set(colors.values())):
        extra_colors = len(colors.values()) - len(set(colors.values()))
        print('Repeated colors: ', extra_colors)
    return colors


# Plots our neutral network.  Uses make_colors to generate colors needed for
# graph G and creates a color map, making the root black and neutral genomes
# white.  All other genomes that are neutral to eachother are given the same
# color in the color map.
def plot_graph(G, covid_seq):
    color_map = []
    colors = make_colors()
    print(len(G.nodes))
    for node in G:
        if node == covid_seq:
            color_map.append('black')
        elif translate(node) == covid_base:
            color_map.append('white')
        elif translate(node) == pangolin_base:
            color_map.append('springgreen')
        elif translate(node) == wuhan_base:
            color_map.append('plum')
        else:
            color = colors[translate(node)]
            color_map.append(color)
    plt.figure()
    plt.axis('off')
    pos = nx.kamada_kawai_layout(G)
    nx.draw(G, pos, node_color=color_map, alpha=0.75)
    ax = plt.gca()
    ax.collections[0].set_edgecolor('black')
    plt.show()


# Main function, defined in initial comment.
def main(task, b_val=0.0, k=None, target_base=None, iterations=40):
    if task == 'neutral clusters':
        # All sequences of covid19
        covid_seqs = translate_AA(covid_base, all_possibilities=True)
        for covid_seq in covid_seqs:
            if covid_seq not in ALL:
                G = add_mutations(covid_seq, b_val)
                plot_graph(G, covid_seq)
    if task == 'neutral clusters 1e':
        p_spike, wu_spike, eng_spike, wa_spike = return_diffs(get_spikes=True)
        assert len(p_spike) == len(wu_spike) == len(eng_spike) == len(wa_spike)
        print(get_size_of_network(p_spike))
        print(get_size_of_network(wu_spike))
        print(get_size_of_network(eng_spike))
        print(get_size_of_network(wa_spike))
        return
        rna_seqs = return_diffs()  # Function from part1e.py
        pang_seq = rna_seqs['Pangolin']
        wuhan_seq = rna_seqs['Wuhan']
        eng_seq = rna_seqs['England']
        print(get_size_of_network(pang_seq))
        print(get_size_of_network(wuhan_seq))
        print(get_size_of_network(eng_seq))
    elif task == 'mutations':
        # Random seqence of covid19
        covid_seq = translate_AA(covid_base)
        G = add_mutations(covid_seq, b_val, k=k)
        plot_graph(G, covid_seq)
    elif task == 'target base':
        covid_seq = translate_AA(covid_base)
        hops, G = add_mutations_target(covid_seq, target_base, b_val)
        print(hops)
        print(len(G.nodes))
    elif task == 'target base with iterations':
        if target_base == bat_base:
            filename = 'bat_base_' + str(iterations) + '.csv'
        elif target_base == civet_base:
            filename = 'civet_base_' + str(iterations) + '.csv'
        elif target_base == sars_base:
            filename = 'sars_base_' + str(iterations) + '.csv'
        covid_seqs = translate_AA(covid_base, all_possibilities=True)
        for it in range(iterations):
            dist = 0
            search = 0
            while dist == 0:
                search += 1
                tried_genomes = []
                # Random seqence of covid19
                covid_seq = translate_AA(covid_base)
                if covid_seq not in tried_genomes:
                    dist = add_mutations_target(covid_seq, target_base, b_val)
                    tried_genomes.append(covid_seq)
                    print('Search', search)
                if search == 50:
                    break
            with open(filename, 'a') as file:
                file.write(str(dist))
                file.write('\n')
            print('Iteration', it)
    else:
        print('Not a valid task')


# main('neutral clusters')
# main('mutations', b_val=0.001, k=200)
# main('target base', b_val=0.15, target_base=bat_base)
# main('target base with iterations', b_val=0.1, target_base=sars_base)

# PART 1e
main('neutral clusters 1e')
