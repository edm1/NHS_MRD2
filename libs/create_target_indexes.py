#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import os
import sys
import subprocess
from bio_file_parsers import fasta_parser, write_fasta
import re
import pandas as pd

def main(args):
    
    make_indexes(args['in_clusters'], args['out_prefix'])
    
    return 0

def make_indexes(cluster_fasta, out_prefix):
    
    bowtie_build = './libs/bowtie2-2.2.1/bowtie2-build'
    
    # Get targets
    target_df = get_target_seqs(cluster_fasta, 0.01)
    
    # Write fasta
    fasta_name = '{0}.fasta'.format(out_prefix)
    with open(fasta_name, 'w') as out_handle:
        for idx in target_df.index:
            write_fasta(out_handle, target_df.loc[idx, 'name'],
                        target_df.loc[idx, 'seq'])
    
    # Run bowtie2-build
    cmd = '{0} {1} {2}'.format(bowtie_build, fasta_name, out_prefix)
    subprocess.call(cmd, shell=True)

def get_target_seqs(clusters_file, threshold):
    """ Takes the fasta file containing clusters from the first stage
        and makes a bowtie2 index for clusters above specified proportion.
    """
    records = []
    pattern = re.compile(r'.*;size=([0-9]+)')
    with open(clusters_file, 'r') as in_handle:
        for header, seq in fasta_parser(in_handle):
            size = long(pattern.match(header).group(1))
            records.append([header, size, seq])
    
    df = pd.DataFrame(records, columns=['name', 'size', 'seq'])
    total_reads = df['size'].sum()
    df['prop'] = df['size'].apply(lambda x: float(x) / total_reads)
    target_df = df.loc[df['prop'] >= threshold, :]
    
    return target_df

if __name__ == '__main__':
    
    args = {}
    args['in_clusters'] = sys.argv[1]
    args['out_prefix'] = sys.argv[2]
    
    main(args)
