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

    query_file = args.query_file
    threads = args.threads

    print("query file:", query_file)
    print("threads: ", threads)

    # you can check if is not None
    if args.taxa_id is not None:
        print("taxa id:", args.taxa_id)

    # if args.taxa_id:
    #     taxa_id = args.taxa_id
    # else:
    #     taxa_id = False


    # blastn_cline = NcbiblastnCommandline(db="refseq_rna", query="test_query.fa", window_masker_taxid="3708", num_threads=10, task="blastn", out="test_blastn_results.xml", outfmt=5) # has to be in .xml for the parser
    #
    # print(blastn_cline)
    #
    # blastn_cline()
    #
    # result_handle = open("test_blastn_results.xml")

    # parse records

    # print title, e value, and identity for top hit




if __name__ == '__main__':
    main()
