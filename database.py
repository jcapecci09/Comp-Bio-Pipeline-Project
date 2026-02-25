from Bio import Entrez
import argparse
import sys

def main():

    # Supply email to Entrez
    Entrez.email = 'jcapecci@luc.edu'

    # create parser object and set arguments
    parser = argparse.ArgumentParser(description="Creates fasta file of Betaherpesvirinae complete genomes")
    parser.add_argument("-o", "--output", help="output file",required=True)

    # grab the output
    arguments = parser.parse_args(sys.argv[1:])
    outfile = arguments.output


    # Search NCBI for Betaherpersvirinae 
    handle = Entrez.esearch(db='nucleotide',
                            term='Betaherpesvirinae[Organism] AND complete genome',
                            retmax=1000)
    record = Entrez.read(handle)
    handle.close()

    # Find ids
    ids = record['IdList']

    # fetch in fasta format the organisms
    fetch_handle = Entrez.efetch(db='nucleotide',
                                id=ids, 
                                rettype='fasta',
                                retmode='text')

    # Write to outfile
    with open(outfile, 'w') as f:
        f.write(fetch_handle.read())

    fetch_handle.close()

if '__main__' == __name__:
    main()
