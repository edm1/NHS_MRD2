#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import os
import sys
from shutil import rmtree
import subprocess
from multiprocessing import cpu_count

def main(args):
    
    # Set initial variables
    bowtie_exec = 'libs/bowtie2-2.2.1/bowtie2'
    cdhit_exec = 'libs/cd-hit-v4.6.1-2012-08-27/cd-hit-est'
    out_dir = args['out_dir']
    num_cpu = min(4, cpu_count())
    
    # Creat out directory
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    else:
        resp = raw_input('Out directory already exists. Overwrite? [y/n]')
        if resp == 'y':
            rmtree(out_dir)
            os.mkdir(out_dir)
        else:
            print 'Exiting.'
            return 1
    
    # Derep fastq
    from libs.derep_fastq import derep_fastq
    print 'Dereplicating raw fastq...'
    fastq_file = os.path.join(out_dir, 'derep_reads.fastq')
    derep_fastq(args['in_reads'], fastq_file)
    
    # Run bowtie
    print 'Running bowtie2 alignment...'
    j_sam = os.path.join(out_dir, 'J_align.sam')
    v_sam = os.path.join(out_dir, 'V_align.sam')
    j_bowtie_cmd = '{0} --very-sensitive-local --reorder -x bowtie_indexes/J_w_phix_indexes -U {1} -S {2} -p {3}'
    v_bowtie_cmd = '{0} --very-sensitive-local --reorder -x bowtie_indexes/V_indexes -U {1} -S {2} -p {3}'
    j_bowtie_cmd = j_bowtie_cmd.format(bowtie_exec, fastq_file, j_sam, num_cpu)
    v_bowtie_cmd = v_bowtie_cmd.format(bowtie_exec, fastq_file, v_sam, num_cpu)
    subprocess.call(j_bowtie_cmd, shell=True)
    subprocess.call(v_bowtie_cmd, shell=True)
    
    # Process sams
    from libs.process_target_sam import parse_sams
    print 'Processing SAM files...'
    ref_names, metrics, ndn_fastq = parse_sams(j_sam, v_sam, out_dir, fastq_file)
    
    # Derep fastq
    print 'Dereplicating N-D-N sequences...'
    ndn_derep_fastq = os.path.join(out_dir, 'derep_ndn.fastq')
    derep_fastq(ndn_fastq, ndn_derep_fastq)
    
    # Convert fastq to fasta for cd-hit
    from libs.convert_fastq_to_fasta import convert_fastq_to_fasta
    ndn_derep_fasta = os.path.join(out_dir, 'NDN_reads.fasta')
    convert_fastq_to_fasta(ndn_derep_fastq, ndn_derep_fasta)
    
    # Make CD-HIT command template
    cdhit_templ = ['{0} -i {1} -o {2}',        # Script, input, output
                   '-c 0.90',                  # Identity threshold
                   '-G 1',                     # Use global alignment
                   '-d 0',                     # Report whole seq name
                   '-s 0.9',                   # Length difference cutof
                   '-r 0',                     # Only search the +/+ strand
                   '-T {3}',                   # Num threads
                   '-M 2000',                  # RAM limit
                   '-p 1',                     # Print alignment in output
                   '-l 8 -n 8']                # Length throw-away, word length
    cdhit_templ = ' '.join(cdhit_templ)
    
    # Run multiple iterations of clustering
    from libs.make_cluster_consensus import make_consensus
    from libs.update_cdhit_sizes import update_sizes
    in_fasta = ndn_derep_fasta
    i = 1
    while i <= 3:
        # Cluster
        print 'Clustering N-D-N sequences step{0}...'.format(i)
        clstr_out = os.path.join(out_dir, 'NDN_clusters_step{0}.fasta'.format(i))
        clstr_meta = clstr_out + '.clstr'
        cdhit_cmd = cdhit_templ.format(cdhit_exec, in_fasta, clstr_out, num_cpu)
        subprocess.call(cdhit_cmd, shell=True)
        # Make consensus
        print 'Making consensus sequences step{0}...'.format(i)
        cons_out = os.path.join(out_dir, 'NDN_clusters_step{0}.consensus'.format(i))
        num_of_clusters = make_consensus(in_fasta, clstr_meta, cons_out)
        # Update sizes
        print 'Updating cluster sizes step{0}...'.format(i)
        total_clusters_size = update_sizes(cons_out, clstr_meta)
        # Set input for next round
        in_fasta = cons_out
        i += 1
    
    # Process clusters
    print 'Processing clusters...'
    from libs.process_clusters import process_clusters
    targets_out = os.path.join(out_dir, 'target_results.txt')
    process_clusters(cons_out, ref_names, metrics, targets_out)
    
    # Write a log file
    log_file = os.path.join(out_dir, 'log_file.txt')
    with open(log_file, 'w') as out_h:
        out_h.write('Total number of reads:\t{0}\n'.format(metrics["total_count"]))
        out_h.write('Num of mapped reads:\t{0}\n'.format(metrics["mapped_count"]))
        out_h.write('Num of unmapped reads:\t{0}\n'.format(metrics["unmapped_count"]))
        out_h.write('phiX mapped reads:\t{0}\n'.format(metrics["phiX_count"]))
        out_h.write('CD19 CAR mapped reads:\t{0}\n'.format(metrics["pUPATrap_count"]))
        out_h.write('\n')
        out_h.write('Number of clusters:\t{0}\n'.format(num_of_clusters))
        out_h.write('Total reads in clusters:\t{0}\n'.format(total_clusters_size))
        
    
    return 0


if __name__ == '__main__':
    
    args = {}
    args['in_reads'] = sys.argv[1]
    args['out_dir'] = sys.argv[2]
    
    main(args)
