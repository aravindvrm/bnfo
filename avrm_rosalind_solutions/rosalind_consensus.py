"""
ROSALIND - Translating RNA into Protein

Given: A collection of at most 10 DNA strings of equal length (at most 1 kbp) in FASTA format.

Return: A consensus string and profile matrix for the collection.
(If several possible consensus strings exist, then you may return any one of them.)
"""
from __future__ import print_function
from collections import defaultdict


def profile_matrix(input_file):
    name = ''
    seq_dict = defaultdict(str)
    s = open(input_file, 'r')
    for line in s:
        line = line.strip()
        if line.startswith('>'):
            name = line[1:]
            continue
        else:
            seq_dict[name] += line

    seq_len = len(seq_dict.values()[0].strip())
    matrix = {'A': [0]*seq_len, 'C': [0]*seq_len, 'G': [0]*seq_len, 'T': [0]*seq_len}

    for k in seq_dict:
        pos = 0
        for nt in seq_dict[k]:
            if nt == 'A':
                matrix['A'][pos] += 1
            elif nt == 'C':
                matrix['C'][pos] += 1
            elif nt == 'G':
                matrix['G'][pos] += 1
            elif nt == 'T':
                matrix['T'][pos] += 1
            pos += 1

    return matrix


def find_consensus(profile):

    consensus_seq = []
    keys = profile.keys()

    for i in range(len(profile[keys[0]])):
        max_v = 0
        max_k = None

        for k in keys:
            v = profile[k][i]

            if v > max_v:
                max_v = v
                max_k = k
        consensus_seq.append(max_k)

    return ''.join(consensus_seq)


m = profile_matrix('rosalind_cons.txt')
print(find_consensus(m))
for k in m.keys():
    print(k, ':', ' '.join(str(c) for c in m[k]))
