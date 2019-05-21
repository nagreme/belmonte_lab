import argparse
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

# Tutorial/Cookbook
# http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc100

# Setup pgm args
parser = argparse.ArgumentParser()

parser.add_argument('-q', '--query_file', dest='query_file', required=True, metavar='query', help='File to blastn')
parser.add_argument('-n', '--num_threads', dest='threads', default=1, metavar='threads', help="Number of threads to use")
parser.add_argument('-t', '--taxa_id', dest='taxa_id', metavar='taxaID', help="Filter to inly search whithin specified taxa (NCBI IDs)")


def main():
    args = parser.parse_args()
    # help(NcbiblastnCommandline)

    if args.taxa_id is not None:
        blastn_cline = NcbiblastnCommandline(db="refseq_rna", query=args.query_file, window_masker_taxid=args.taxa_id, num_threads=args.threads, task="blastn", out="test_blastn_results.xml", outfmt=5) # has to be in .xml for the parser
    else:
        blastn_cline = NcbiblastnCommandline(db="refseq_rna", query=args.query_file, num_threads=args.threads, task="blastn", out="test_blastn_results.xml", outfmt=5) # has to be in .xml for the parser

    print("\nRunning:\n", blastn_cline)

    # Run the blastn
    stdout, stderr = blastn_cline()

    result_handle = open("test_blastn_results.xml")

    # parse records
    blast_records = NCBIXML.parse(result_handle)

    q_num = 1

    print("\n\nTop hits:\n")

    # print title, e value, and identity for top hit
    for record in blast_records:
        # check if we have at least one hit
        print("Query:", q_num) # to be friendly to humans
        if record.alignments:
            alignment = record.alignments[0] # assuming first one is top hit (seems right from limited testing)
            for hsp in alignment.hsps: # there seems to only ever be one? At least in our case? Unclear
                print(alignment.title)
                print("e value:", hsp.expect)
                print("identity:", hsp.identities/record.query_length)
        else:
            print("No hits")

        q_num += 1
        print()

    result_handle.close()


if __name__ == '__main__':
    main()
