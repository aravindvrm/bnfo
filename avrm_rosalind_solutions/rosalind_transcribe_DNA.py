"""
ROSALIND - Transcribing DNA to RNA

Given: A DNA string  having length at most 1000 nt.

Return: The transcribed RNA string
"""
from __future__ import print_function
from __future__ import division


def transcribe_sequence(input_file):
    rna_sequence = []
    for line in input_file:
        line = line.strip()

        for nt in line:
            if nt == 'T':
                rna_sequence.append('U')
            else:
                rna_sequence.append(nt)

    return ''.join(rna_sequence)


infile = open('rosalind_rna.txt', 'r')
print(transcribe_sequence(infile))

