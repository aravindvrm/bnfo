"""
ROSALIND - Translating RNA into Protein

Given: An RNA string corresponding to a strand of mRNA (of length at most 10 kbp).

Return: The RNA sequence translated into a string of amino acids
"""
from __future__ import print_function


def translate_rna(input_file):

    rna_seq = open(input_file, 'r').read()

    codon_table = {'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
                   'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
                   'UAU': 'Y', 'UAC': 'Y', 'UAA': '-', 'UAG': '-',
                   'UGU': 'C', 'UGC': 'C', 'UGA': '-', 'UGG': 'W',
                   'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
                   'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                   'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                   'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                   'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
                   'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                   'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
                   'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                   'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
                   'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                   'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
                   'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

    triplets = []
    aa_seq = ''
    for i in range(0, len(rna_seq), 3):
        triplets.append(rna_seq[i:i + 3])
    for codon in triplets:
        if codon in codon_table:
            aa_seq += (codon_table[codon][0])

    return aa_seq


print(translate_rna('rosalind_prot.txt'))

