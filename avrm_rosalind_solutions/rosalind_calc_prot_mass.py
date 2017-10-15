"""
ROSALIND - Calculating Protein Mass

Given: A protein string  of length at most 1000 aa.

Return: The total weight . Consult the monoisotopic mass table.
"""
from __future__ import print_function


def protein_mass(aa_seq):

    tot_weight = 0
    # monoisotopic mass table
    mass = {'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
            'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
            'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
            'M': 131.04049, 'F': 147.06841, 'P': 97.05276, 'S': 87.03203,
            'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841}

    for aa in aa_seq:
        tot_weight += mass[aa]

    return tot_weight


sample_seq = 'SKADYEK'
full_seq = open('rosalind_prtm.txt').read().strip()

print(protein_mass(sample_seq))
# print(protein_mass(full_seq))

