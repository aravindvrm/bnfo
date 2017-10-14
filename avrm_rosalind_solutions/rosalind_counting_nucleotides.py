"""
ROSALIND - Counting DNA Nucleotides

Given: A DNA string  of length at most 1000 nt.

Return: Four integers (separated by spaces) counting the respective number of times
that the symbols 'A', 'C', 'G', and 'T' occur in.
"""
from __future__ import print_function
from __future__ import division


def nucleotide_counter(input_file):
    for line in input_file:
        line = line.strip()
        total_counts = [0]*4

        for nt in line:
            if nt == 'A':
                total_counts[0] += 1
            elif nt == 'C':
                total_counts[1] += 1
            elif nt == 'G':
                total_counts[2] += 1
            else:
                total_counts[3] += 1

    return ' '.join(str(c) for c in total_counts)


infile = open('rosalind_dna.txt', 'r')
print('A\tC\tG\tT')
print(nucleotide_counter(infile))