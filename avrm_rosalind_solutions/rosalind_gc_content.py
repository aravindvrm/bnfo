"""
ROSALIND - Computing GC Content

Given: At most 10 DNA strings in FASTA format (of length at most 1 kbp each).

Return: The ID of the string having the highest GC-content, followed by the GC-content of that string.
Rosalind allows for a default error of 0.001 in all decimal answers unless otherwise stated;
"""
from __future__ import print_function
from __future__ import division

def GC_content(input_file):

    output = {}

    for line in input_file:
        line = line.strip()

        if line.startswith('>'):
            header = line[1:]
            output[header] = ''
        else:
            output[header] += line

    for sequence in output.keys():

        total_counts = [0] * 4

        for nt in output[sequence]:
            if nt == 'A':
                total_counts[0] += 1
            elif nt == 'C':
                total_counts[1] += 1
            elif nt == 'G':
                total_counts[2] += 1
            else:
                total_counts[3] += 1
        gc_ratio = ((total_counts[1] + total_counts[2]) /
                    (total_counts[0] + total_counts[1] +
                     total_counts[2] + total_counts[3]))
        output[sequence] = round(gc_ratio, 8) * 100

    return output


infile = open('rosalind_gc.txt', 'r')
gc_content = GC_content(infile)
for seq in gc_content:
    print(seq,  gc_content[seq])
