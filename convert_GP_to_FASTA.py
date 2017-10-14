"""
BNFO 601 - Exam 1 - 09/30/17
converts GenPept format file to FASTA format file
"""
from __future__ import print_function
import re
from string import digits


class GenPeptToFASTA(object):

    def __init__(self, input_file):
        self.input_file = input_file
        self.output = []

    def read_input(self):

        # regex to capture required fields
        LOCUS_RE = re.compile("LOCUS\s*(\S*)\s")
        DEFINITION_RE = re.compile("DEFINITION(.*|\n)(.|\s{8}.+)")
        SEQUENCE_RE = re.compile("^\s+\d{1,3}\s[a-z].*")
        GI_RE = re.compile("GI:(\d+)")
        elements = []  # list to collect fields required for FASTA format
        complete_seq = ''

        # begin reading input file
        for line in self.input_file:

            seq = SEQUENCE_RE.search(line)
            if seq:
                # convert sequence to FASTA format
                # remove whitespace and digits, convert to uppercase
                partial_seq = seq.group().translate(None, digits).replace(" ", "").upper()
                complete_seq += partial_seq

            line = line.strip()
            gb = LOCUS_RE.search(line)
            desc = DEFINITION_RE.search(line)
            gi = GI_RE.search(line)

            if gb:
                gb_out = '|gb|' + gb.group(1) + '|'
            elif gi:
                gi_out = '>gi|' + gi.group(1) + gb_out  # build output string
                elements.insert(0, gi_out)
            elif desc:
                s = desc.group(1) + desc.group(2)
                str.join(" ", s.splitlines())
                elements.insert(1, s)  # insert description in element list
            # '//' indicates end of entry
            elif line.startswith("//"):
                elements.insert(2, complete_seq)  # insert sequence in element list
                complete_seq = ''  # reset sequence string
                self.output.append(elements)  # add list of elements to the output list
                elements = []  # reset element list

    def write_output(self, output_file):
        # open output file and print the elements of the output list
        with open(output_file, 'w') as outfile:
            for q in self.output:
                print(q[0] + q[1], file=outfile, end='\n')
                print(q[2], file=outfile, end='\n\n')


def main():

    # set input file
    filename = "O104_H4_GP.txt"
    outfile = filename.replace("GP.txt", "FASTA.txt")

    with open(filename) as infile:
        convert = GenPeptToFASTA(infile)
        convert.read_input()
        convert.write_output(outfile)

if __name__ == "__main__":
    main()
