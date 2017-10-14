"""
ROSALIND - Counting Point Mutations

Given: Two DNA strings  and  of equal length (not exceeding 1 kbp).

Return: The Hamming distance
"""
from __future__ import print_function
from __future__ import division


def hamming_dist(input_file):
    h_dist = 0
    sequences = open(input_file, 'r').readlines()
    seq_len = len(sequences[0])
    for i in range(seq_len):
        if sequences[0][i] != sequences[1][i]:
            h_dist += 1
    return h_dist


seq_1 = 'GAGCCTACTAACGGGAT'
seq_2 = 'CATCGTAATGACGGCCT'
print(hamming_dist('rosalind_hamm.txt'))