#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import os
import sys
from shutil import rmtree
import subprocess

def main(args):
    
    # Set initial variables
    bowtie_exec = 'libs/bowtie2-2.2.1/bowtie2'
    cdhit_exec = 'libs/cd-hit-v4.6.1-2012-08-27/cd-hit-est'
    out_dir = args['out_dir']
    
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
    j_bowtie_cmd = '{0} -p 4 --very-sensitive-local --reorder -x bowtie_indexes/J_w_phix_indexes -U {1} -S {2}'
    v_bowtie_cmd = '{0} -p 4 --very-sensitive-local --reorder -x bowtie_indexes/V_indexes -U {1} -S {2}'
    j_bowtie_cmd = j_bowtie_cmd.format(bowtie_exec, fastq_file, j_sam)
    v_bowtie_cmd = v_bowtie_cmd.format(bowtie_exec, fastq_file, v_sam)
    subprocess.call(j_bowtie_cmd, shell=True)
    subprocess.call(v_bowtie_cmd, shell=True)
    
    # Process sams
    from libs.process_target_sam import parse_sams
    print 'Processing SAM files...'
    ref_names, metrics, ndn_fasta = parse_sams(j_sam, v_sam, out_dir)
    
    # Derep fasta
    from libs.derep_fasta import derep_fasta
    print 'Dereplicating N-D-N sequences...'
    ndn_derep = os.path.join(out_dir, 'derep_ndn.fasta')
    derep_fasta(ndn_fasta, ndn_derep)
    
    # Make CD-HIT command template
    cdhit_templ = ['{0} -i {1} -o {2}', # Script, input, output
                   '-c 0.90', # Identity threshold
                   '-G 1',    # Use global alignment
                   '-d 0',    # Report whole seq name
                   '-s 0.9',  # Length difference cutof
                   '-r 0',    # Only search the +/+ strand
                   '-T 0',    # Num threads
                   '-M 2000', # RAM limit
                   '-p 1',    # Print alignment in output
                   '-l 8 -n 8'] # Length throw-away, word length
    cdhit_templ = ' '.join(cdhit_templ)
    
    # Run multiple iterations of clustering
    from libs.make_cluster_consensus import make_consensus
    from libs.update_cdhit_sizes import update_sizes
    in_fasta = ndn_derep
    for i in range(1, 4):
        # Cluster
        print 'Clustering N-D-N sequences step{0}...'.format(i)
        clstr_out = os.path.join(out_dir, 'NDN_clusters_step{0}.fasta'.format(i))
        clstr_meta = clstr_out + '.clstr'
        cdhit_cmd = cdhit_templ.format(cdhit_exec, in_fasta, clstr_out)
        subprocess.call(cdhit_cmd, shell=True)
        # Make consensus
        print 'Making consensus sequences step{0}...'.format(i)
        cons_out = os.path.join(out_dir, 'NDN_clusters_step{0}.consensus'.format(i))
        num_of_clusters = make_consensus(in_fasta, clstr_meta, cons_out)
        # Update sizes
        print 'Updating cluster sizes step{0}...'.format(i)
        total_clusters_size = update_sizes(cons_out, clstr_meta)
        
        in_fasta = cons_out
    
    #~ # Run cd-hit
    #~ print 'Clustering N-D-N sequences step1...'
    #~ clstr_out = os.path.join(out_dir, 'NDN_clusters_step1.fasta')
    #~ clstr_meta = clstr_out + '.clstr'
    #~ cdhit_cmd = cdhit_templ.format(cdhit_exec, ndn_derep, clstr_out)
    #~ subprocess.call(cdhit_cmd, shell=True)
    #~ 
    #~ # Make consensus sequences (and return number of clusters)
    #~ print 'Making consensus sequences step1...'
    #~ from libs.make_cluster_consensus import make_consensus
    #~ cons_out1 = os.path.join(out_dir, 'NDN_clusters_step1.consensus')
    #~ num_of_clusters = make_consensus(ndn_derep, clstr_meta, cons_out1)
    #~ 
    #~ # Update cd-hit cluster sizes
    #~ from libs.update_cdhit_sizes import update_sizes
    #~ print 'Updating cluster sizes step1...'
    #~ total_clusters_size = update_sizes(cons_out1, clstr_meta)
    #~ 
    #~ # Cluster again on consensus seqs
    #~ print 'Clustering N-D-N sequences step2...'
    #~ clstr_out = os.path.join(out_dir, 'NDN_clusters_step2.fasta')
    #~ clstr_meta = clstr_out + '.clstr'
    #~ cdhit_cmd = cdhit_templ.format(cdhit_exec, cons_out1, clstr_out)
    #~ subprocess.call(cdhit_cmd, shell=True)
    #~ 
    #~ # Make consensus step 2
    #~ print 'Making consensus sequences step2...'
    #~ cons_out2 = os.path.join(out_dir, 'NDN_clusters_step2.consensus')
    #~ num_of_clusters = make_consensus(cons_out1, clstr_meta, cons_out2)   
    #~ 
    #~ # Update cd-hit cluster sizes step2
    #~ print 'Updating cluster sizes step2...'
    #~ total_clusters_size = update_sizes(cons_out2, clstr_meta)
    
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
