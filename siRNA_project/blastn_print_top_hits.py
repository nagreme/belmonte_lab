import argparse
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

# Note: I had to modify the NcbiblastnCommandline object ot add the -taxids option
# (see ~/.local/lib/python3.6/site-packages/Bio/Blast/Applications.py)

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
        blastn_cline = NcbiblastnCommandline(db="refseq_rna_v5", query=args.query_file, num_threads=args.threads, task="blastn-short", out="blastn_results.xml", outfmt=5, num_alignments=3, taxids=args.taxa_id) # has to be in .xml for the parser

    else:
        # It's an object (not a str) so I'm not sure how to just append taxid. We'll deal with the redundancy for now
        blastn_cline = NcbiblastnCommandline(db="refseq_rna_v5", query=args.query_file, num_threads=args.threads, task="blastn-short", out="blastn_results.xml", outfmt=5, num_alignments=3) # has to be in .xml for the parser

    print("\nRunning:\n", blastn_cline)

    # Run the blastn
    stdout, stderr = blastn_cline()

    result_handle = open("blastn_results.xml")

    # parse records
    blast_records = NCBIXML.parse(result_handle)

    print("\n\nTop hits:\n")

    # print title, e value, and identity for top hit
    for record in blast_records:
        # check if we have at least one hit
        print("Query:", record.query)
        if record.alignments:
            alignment = record.alignments[0] # assuming first one is top hit
            hsp = alignment.hsps[0] # there should be one if there's an alignment. there's sometimes more than one hsp but the first one has higher or equal score to the others
            print(alignment.title)
            print("e value:", hsp.expect)
            print("identity:", hsp.identities/record.query_length)
        else:
            print("No hits")

        print()

    result_handle.close()


if __name__ == '__main__':
    main()
