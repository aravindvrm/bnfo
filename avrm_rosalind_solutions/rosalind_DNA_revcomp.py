"""
ROSALIND - Complementing a Strand of DNA

Given: A DNA string  of length at most 1000 bp.

Return: The reverse complement
"""
from __future__ import print_function
from __future__ import division


def revcomp(input_file):
    revcomped = []

    for line in input_file:
        line = line.strip()

        for nt in line:
            if nt == 'A':
                revcomped.append('T')
            elif nt == 'C':
                revcomped.append('G')
            elif nt == 'G':
                revcomped.append('C')
            elif nt == 'T':
                revcomped.append('A')
    revcomped.reverse()
    return ''.join(revcomped)


infile = open('rosalind_revc.txt', 'r')
print(revcomp(infile))