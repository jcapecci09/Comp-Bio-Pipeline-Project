import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description="Creates fasta files of CDS and counts number of CDS")
    parser.add_argument("-s", "--samples", help="input file", required=True)
    parser.add_argument("-c", "--condition", help="input file", required=True)
    parser.add_argument("-o", "--output", help="output file",required=True)


    # Allow parser to see command line
    arguments = parser.parse_args(sys.argv[1:])

    samples = arguments.samples
    condition = arguments.condition
    outfile = arguments.output

    with open(outfile, 'w') as f:
        for sample in samples:
            f.write(f'{sample}\t{condition[sample]}\tresults\{sample}\n')

if '__main__' == __name__:
    main()
