import argparse
import sys

def main():

    # create argparse object and set arguments
    parser = argparse.ArgumentParser(description="turns tsv into text and adds appropriate headers")
    parser.add_argument("-o", "--output", help="output file",required=True)
    parser.add_argument("-i", "--input", help="output file",required=True)

    # Grab the arguments
    arguments = parser.parse_args(sys.argv[1:])

    # Set outfile and infile
    outfile = arguments.output
    infile = arguments.input

    # determine samples in pipeline
    samples = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']

    # Collect each line from tsv
    with open(infile, 'r') as f:
        collect = []
        for line in f:
            collect.append(line)

    # Write to the outfile
    with open(outfile, 'w') as f:

        # Intiate an index for samples
        # and a count for line count
        count = 4
        index = 0
        for line in collect:

            # If count == 4 we've reached a header write the appropriate lines
            if count == 4:
                f.write(f'{samples[index]}:\n')
                f.write(f'sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n')
                f.write(line)
                index += 1
                count = 0

            # Else write the rest
            else:
                f.write(line)
                count += 1

if '__main__' == __name__:
    main()
            


        