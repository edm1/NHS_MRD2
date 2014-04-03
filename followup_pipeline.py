#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import os
import sys
import subprocess

def main(args):
    
    bowtie_script = './libs/bowtie2-2.2.1/bowtie2'
    
    # Make bowtie indexes
    from libs.create_target_indexes import make_indexes
    targets_prefix = os.path.join(args['out_dir'], 'targets')
    target_names = make_indexes(args['in_clusters'], targets_prefix, 0.05)
    
    # Derep fastq
    from libs.derep_fastq import derep_fastq
    print 'Dereplicating fastq...'
    fastq_derep = os.path.join(args['out_dir'], 'follow_up_reads.fastq')
    derep_fastq(args['in_fastq'], fastq_derep)
    
    # Run bowtie
    print 'Running bowtie2...'
    sam_out = os.path.join(args['out_dir'], 'followup_mapped.sam')
    cmd = '{0} -p 4 --very-sensitive-local --reorder -x {1} -U {2} -S {3}'
    cmd = cmd.format(bowtie_script, targets_prefix, fastq_derep, sam_out)
    subprocess.call(cmd, shell=True)
    
    return 0



if __name__ == '__main__':
    
    args = {}
    args['in_clusters'] = sys.argv[1]
    args['in_fastq'] = sys.argv[2]
    args['out_dir'] = sys.argv[3]
    
    main(args)
