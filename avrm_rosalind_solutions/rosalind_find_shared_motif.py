"""
ROSALIND - Finding a Shared Motif

Given: A collection of DNA strings of length at most 1 kbp each in FASTA format.

Return: A longest common substring of the collection.
(If multiple solutions exist, you may return any single solution.)
"""
from __future__ import print_function
from collections import defaultdict


def shared_motif(input_file):
    name = ''
    seq_dict = defaultdict(str)

    for line in input_file:
        line = line.strip()
        if line.startswith('>'):
            name = line[1:]
            continue
        else:
            seq_dict[name] += line
    # for m in seq_dict.keys():
    #     print(m, seq_dict[m])
    srt_seq = sorted(seq_dict.values(), key=len)
    short_seq = srt_seq[0]
    comp_seq = srt_seq[1:]
    motif = ''

    for i in range(len(short_seq)):

        for j in range(i, len(short_seq)):
            m = short_seq[i:j + 1]
            found = False

            for sequ in comp_seq:

                if m in sequ:
                    found = True

                else:
                    found = False
                    break

            if found and len(m) > len(motif):
                motif = m

    return motif


infile = open('rosalind_lcsm.txt', 'r')
print(shared_motif(infile))

