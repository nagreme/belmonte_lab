import argparse
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

# Setup pgm args
parser = argparse.ArgumentParser()

parser.add_argument('-q', '--query_file', required=True, metavr='query', help='File to blastn')
parser.add_argument('-p', '--num_threads', default=1, metavar='threads', help="Number of threads to use")
parser.add_argument('-t', '--taxa_id', metavar='taxaID', help="Filter to inly search whithin specified taxa (NCBI IDs)")


def main():
    args = parser.parse_args()
    # help(NcbiblastnCommandline)

    blastn_cline = NcbiblastnCommandline(db="refseq_rna", query="test_query.fa", window_masker_taxid="3708", num_threads=10, task="blastn", out="test_blastn_results.xml", outfmt=5) # has to be in .xml for the parser

    print(blastn_cline)

    blastn_cline()

    result_handle = open("test_blastn_results.xml")
