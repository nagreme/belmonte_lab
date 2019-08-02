# --------------------------------------------------------------
# Date: 2019-07-30
# Nadège Pulgar-Vidal
#
# Make bed4 files from the blast result data of ShortStack output
# (siRNAs) so that it can be compared with Deirdre's methylation
# data using bedIntersect
#
# --------------------------------------------------------------

import sys
from Bio.Blast import NCBIXML

def main():
    if len(sys.argv) < 2:
        print("Usage $python3 shortstack_blastn_to_bed4.py blastn_filename.xml")
        sys.exit()

    infile = sys.argv[1]

    result_handle = open(infile)

    # parse records
    blast_records = NCBIXML.parse(result_handle)

    for record in blast_records:
        query = record.query
        col4 = "nohit"

        if record.alignments:
            alignment = record.alignments[0] # look only at top hit
            col4 = alignment.accession

        print("{}\t{}\t{}\t{}".format(query[:query.index(':')], query[query.index(':')+1:query.index('-')], query[query.index('-')+1:], col4))


    result_handle.close()



if __name__ == '__main__':
    main()
