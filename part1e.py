# -*- coding: utf-8 -*-
"""
CS523 Project 2
Used for part 1e.  Gets the spike protein sequences from 4 cases of COVID-19:
    Pangolin (2019) - GISAID Accession ID: EPI_ISL_410721
    Wuhan (12/24/2019) - GISAID Accession ID: EPI_ISL_402120
    England (3/28/2020) - GISAID Accession ID: EPI_ISL_420250
    Washington (3/24/2019) - GISAID Accenssion ID: EPI_ISL_418915
"""
d = 'fasta_files\\'


# Returns the spike protein genome.
def get_spike(filename, s_min, s_max, c_min, c_max):
    with open(filename, 'r') as infile:
        filename = (infile.readline()).rstrip().split('|')[1] + '.txt'
        # print(filename)
        spike_protein = ""
        for i, line in enumerate(infile):
            if i == s_min:
                spike_protein += line[c_min:-1]
            if i > s_min and i < s_max:
                spike_protein += line.rstrip()
            if i == s_max:
                spike_protein += line[:c_max]
    spike_protein = list(spike_protein)
    for i in range(len(spike_protein)):
        if spike_protein[i] == 'T':
            spike_protein[i] = 'U'
    return "".join(spike_protein)


# Same function as in neutral_net.py, translates a sequence of nucleotides
# into amino acids
def translate(seq):

    table = {
            'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
            'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
            'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
            'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
            'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
            'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
            }
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            if codon in table:
                protein += table[codon]
            else:
                protein += 'X'  # Way to track N's
    else:
        print('Error in modulo: ', len(seq) % 3)
    return protein


# Used in error checking, prints differences in two spikes
def print_diffs(spike1, spike2):
    source = translate(spike1)
    target = translate(spike2)
    indxs = [i for i in range(len(source)) if source[i] != target[i]]
    print('Spike 1', [source[i] for i in indxs])
    print('Spike 2', [target[i] for i in indxs])


def return_strand(spike, indxs):
    strand = ""
    for i in indxs:
        strand += spike[3*i:3*i+3]
    return strand


# Only used for returning spike proteins in final version of the report.
# Called in neutral_net.py
def return_diffs(get_spikes=False):
    val = -1998  # 1692
    pan_spike = get_spike(d+'pangolin1.fasta', 268, 316, -5, 58)
    wuhan_spike = get_spike(d+'12_24_Wuhan.fasta', 269, 317, 42, 24)
    eng_spike = get_spike(d+'3_25_England.fasta', 269, 317, 42, 24)
    wa_spike = get_spike(d+'3_24_WA.fasta', 267, 315, 78, 60)
    if get_spikes:
        return pan_spike, wuhan_spike, eng_spike, wa_spike
    tran_pan = translate(pan_spike[-val:])
    tran_wuhan = translate(wuhan_spike[-val:])
    # Index 613 is mutated from wuhan to england/wa
    pang_613 = return_strand(pan_spike, [613])
    wuhan_613 = return_strand(wuhan_spike, [613])
    eng_613 = return_strand(eng_spike, [613])
    wa_613 = return_strand(wa_spike, [613])

    # Indexes of different values in amino acids.  There are only differences
    # in the pangolin and wuhan cases, the other human cases are neutral.
    indices = [i for i in range(len(tran_pan)) if tran_pan[i] != tran_wuhan[i]]
    return {'Pangolin': pang_613 + return_strand(pan_spike[-val:], indices),
            'Wuhan': wuhan_613 + return_strand(wuhan_spike[-val:], indices),
            'England': eng_613 + return_strand(eng_spike[-val:], indices),
            'WA': wa_613 + return_strand(wa_spike[-val:], indices)
            }
