"""Counts the number of CDS and creates a fasta file that will be used 
for building a Kallisto index. 

Author: Jimmy Capecci
"""

import sys
import argparse

def main():

    # Create parser object
    parser = argparse.ArgumentParser(description="Creates fasta files of CDS and counts number of CDS")
    parser.add_argument("-i", "--input", help="input file", required=True)
    parser.add_argument("-o", "--output", help="output file",required=True)
    parser.add_argument("-f", "--file", help="output file",required=True)

    # Allow parser to see command line
    arguments = parser.parse_args(sys.argv[1:])

    # Set infile and outfile
    infile = arguments.input
    outfile = arguments.output
    file = arguments.file

    # Open the infile and collect records in a dictionary
    with open(infile, 'r') as f:
        d = {}

        for line in f:

            # If line is a header collect protein_id and add it to dictionary
            if line.startswith('>'):
                protein_index = line.index('protein_id=') + 11 # add 11 because we want everything after protein_id
                protein_id = line[protein_index:].split(']')[0] # this remvoes everythng after ] so that we only collect protein id
                d[protein_id] = '' # add to dictionary

            # If line is fasta sequence add to dictionary for the appropriat protein_id
            else:
                d[protein_id] += line.strip()

    # Open outfile and pipeline report
    with open(outfile, 'w') as f, \
         open(file, 'w') as f1:
        
        # For each protein id write its correspondent fasta sequence
        for key in d:
            f.write(f'>{key}\n')
            f.write(f'{d[key]}\n')
        
        # Write the number of CDS into the pipeline report
        f1.write(f'The HCMV genome (GCF_000845245.1) has {len(d)} CDS\n')
        

if '__main__' == __name__:
    main()
                   