#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
import pickle
import libs.stage1_funcs as s1f
import libs.stage1_classes as s1c
import libs.motif_filter as mtf
import libs.clustering as clf
import subprocess

def main(args):
    
    start_time = time.time()
    
    # Print start time
    if not args['--quiet']:
        print "\n{1} Job name {0} ---------------".format(args['<id>'],
                time.strftime('%H:%M:%S', time.localtime()))
    
    # Prepare folders for analysis
    s1f.prepare_folders(args)
    
    # Write parameters to log file
    s1f.log_parameters('./results/{0}/{0}.log'.format(args['<id>']), args)

    #
    # Filter forward and reverse reads down to those containing motif, and
    # splice clone-specific subsequence
    #
    
    # Print progress to stdout
    if not args['--quiet']:
        s1f.print_timestamp('Searching reads for motif...')
    
    # Execute motif filter
    reads_path, total_input_seqs = mtf.filter_reads(
            args['<fastq>'],
            args['<id>'],
            './tmp/{0}'.format(args['<id>']),
            args['--r'],
            args['--germ'],
            args['--ab'],
            args['--bb'],
            args['--m'],
            './results/{0}/{0}.log'.format(args['<id>']))

    #
    # Clustering: dereplication, clustering and consensus calculation
    #

    if not args['--quiet']:
        s1f.print_timestamp('Clustering reads...')
    
    clus_start_time = time.time()
    
    # Dereplicate the reads
    record_list, num_unique_seqs, total_filtered_reads = \
            clf.dereplication(reads_path, args['--min_dup'])
    
    record_list_len = len(record_list)
    
    # DEBUG **********************
    record_tot_size = 0
    for rec in record_list:
        record_tot_size += rec['size']
    print 'Total size check 1:', record_tot_size
    
    # Do clustering
    clusters = clf.cluster(record_list[:],
                           args['--identity'])
    del record_list
    num_clusters = len(clusters)
    
    # DEBUG **********************
    record_tot_size = 0
    for key in clusters:
        for rec in clusters[key]:
            record_tot_size += rec['size']
    print 'Total size check 2:', record_tot_size
    
    # Convert cluster groups into an ordered list of consensus records
    consensus_list = clf.clusters_to_consensus(clusters)
    del clusters
    
    # DEBUG **********************
    record_tot_size = 0
    for consen in consensus_list:
        record_tot_size += consen['size']
    print 'Total size check 3:', record_tot_size
    
    # Write clustering stats to log file
    clf.log_stats('./results/{0}/{0}.log'.format(args['<id>']),
                  total_filtered_reads,
                  num_unique_seqs,
                  record_list_len,
                  num_clusters,
                  clus_start_time)
    
    # Check that more than 1 cluster exists (else BLAST will crash)
    if not num_clusters > 1:
        print 'Warning: Less than 2 clusters were produced. Exiting!'
        return 1
    
    #
    # Initiate records class, load clusters and write target fastas
    #

    if not args['--quiet']:
        s1f.print_timestamp('Writing top clusters...')
    
    cluster_records = s1c.Records()
    cluster_records.add_records_from_clustering(consensus_list)
    del consensus_list

    # Write total number of reads in input fastq to Records class for use in
    # sample comparison
    
    cluster_records.total_input_seqs = total_input_seqs
    
    # DEBUG ****************************
    print 'Total size check 4:', cluster_records.total_clusters_size
    
    # Write the top N clusters to fasta for BLASTing
    cluster_records.write_top_to_fasta(
            args['--top'],
            './tmp/{0}/{0}.top'.format(args['<id>']),
            size_out=False)
            
    # Write all clusters to fasta in results folder
    cluster_records.write_top_to_fasta(
            None,
            './results/{0}/{0}.cons'.format(args['<id>']),
            size_out=True)

    #
    # BLAST consensus sequences for V, D and J hits.
    #
    
    if not args['--quiet']:
        s1f.print_timestamp('BLASTing segments...')
    
    # BLAST D regions
    s1f.run_D_blast('./tmp/{0}/{0}.top'.format(args['<id>']), args['--t'])
    
    # BLAST J regions
    s1f.run_J_blast('./tmp/{0}/{0}.top'.format(args['<id>']), args['--t'])
    
    # BLAST V regions
    s1f.run_V_blast('./tmp/{0}/{0}.top'.format(args['<id>']), args['--t'])
    
    if not args['--quiet']:
        s1f.print_timestamp('Parsing BLAST results...')
                
    # Parse the D and J segment hits from XML file. J must be done before D.
    # NB that J and V must be done before D
    for segment in ['J', 'V', 'D']:
        in_file = './tmp/{0}/{0}.top_{1}.blast'.format(args['<id>'], segment)
        cluster_records.add_VDJ_BLAST_xml(in_file,
                                          args['--e'],
                                          args['--t'],
                                          segment)

    #
    # Write results files
    #
    
    if not args['--quiet']:
        s1f.print_timestamp('Writing outputs...')
        
    # Write the .dat summary file
    with open('./results/{0}/{0}.dat'.format(args['<id>']), 'w') as out_handle:
        pickle.dump(cluster_records, out_handle)
    
    # Write the .tab summary file
    tab_file = './results/{0}/{0}.tab'.format(args['<id>'])
    cluster_records.write_tabulated(tab_file, args['--top'])
    
    # Write the detailed summary report
    cluster_records.write_detail_output(
            './results/{0}/{0}.detail'.format(args['<id>']),
            args['--top'])

    # Try to make bubble plot using Rscript
    cmd = ["Rscript", "extras/make-bubbleplot.R", tab_file]
    ret = subprocess.call(cmd)
    if ret == 0:
        print("Bubble plot done.")
    else:
        print("Failed to make bubble plot.")



    
    return 0

if __name__ == '__main__':
    pass
