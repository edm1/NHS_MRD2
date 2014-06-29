#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Contains simple parsers for fasta and fastq files, taken directly from the
# biopython source code. Also, includes a homemade BLAST xml output parser that
# is written specifically with NHS MRD needs in mind.
#

import sys
from bio_file_parsers import write_fasta, fastq_parser


def main():

    convert_fastq_to_fasta(args['in_fastq'], args['out_fasta'])
    
    return 0

def convert_fastq_to_fasta(fastq_in, fasta_out):
    
    with open(fastq_in, 'r') as in_h:
        with open(fasta_out, 'w') as out_h:
            for title, seq, _ in fastq_parser(in_h):
                write_fasta(out_h, title, seq)

if __name__ == '__main__':
    
    global args
    args = {}
    args['in_fastq'] = sys.argv[1]
    args['out_fasta'] = sys.argv[2]
    
    main()

 
