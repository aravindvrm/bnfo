"""
ROSALIND - Finding a Motif in DNA

Given: Two DNA strings, target sequence and motif

Return: Start positions of all locations of motif
"""
from __future__ import print_function
from __future__ import division
import re


def simple_search(sequence, query):
    output = []
    matches = re.compile(query)

    for m in matches.finditer(sequence):
        output.append(m.start() + 1)

    return ' '.join(str(c) for c in output)


# Overlapping motifs are counted in this version
def simple_search_2(sequence, query):
    output = []
    query_len = len(query)
    sequence_len = len(sequence)

    for i in range(0, sequence_len - query_len + 1):
        if sequence[i:i+query_len] == query:
            output.append(i + 1)

    return ' '.join(str(c) for c in output)

seq = 'AGTGTTATCCCAATGTATAGAGACCCAATGCCCAATGTGAGCTCCCAATGCCCCAATGTCCC'
query = 'CCCAATGCC'
print(simple_search_2(seq, query))