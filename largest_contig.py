import argparse
import sys


def main():

    parser = argparse.ArgumentParser(description="Creates fasta files of CDS and counts number of CDS")
    parser.add_argument("-i", "--input", help="input file", required=True)
    parser.add_argument("-o", "--output", help="output file",required=True)

    arguments = parser.parse_args(sys.argv[1:])
    infile = arguments.input
    outfile = arguments.output

    with open(infile, 'r') as f:
        header = f.readline().strip()
        longest_contig = ''
        for line in f:
            if line.startswith('>'):
                break  
            else:
                longest_contig += line.strip() + '\n'
    
    with open(outfile, 'w') as f:
        f.write(f'{header}\n{longest_contig}\n')

if '__main__' == __name__:
    main()






