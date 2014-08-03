#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import os
import sys
from shutil import rmtree
import subprocess
from multiprocessing import cpu_count
import argparse

def main():

    # Parse args
    global args
    args = parse_arguments()

    # Set initial variables
    num_cpu = min(args.numCPU, cpu_count())

    # Creat out directory
    if not os.path.exists(args.OutDir):
        os.mkdir(args.OutDir)
    else:
        resp = raw_input('Out directory already exists. Overwrite? [y/n]')
        if resp == 'y':
            rmtree(args.OutDir)
            os.mkdir(args.OutDir)
        else:
            print 'Exiting.'
            return 1

    # Derep fastq
    from libs.derep_fastq import derep_fastq
    print 'Dereplicating raw fastq...'
    fastq_file = os.path.join(args.OutDir, 'derep_reads.fastq')
    derep_fastq(args.InFastq, fastq_file)

    # Run bowtie
    print 'Running bowtie2 alignment...'
    j_sam = os.path.join(args.OutDir, 'J_align.sam')
    v_sam = os.path.join(args.OutDir, 'V_align.sam')
    j_bowtie_cmd = '{0} --very-sensitive-local --reorder -x bowtie_indexes/J_w_phix_indexes -U {1} -S {2} -p {3}'
    v_bowtie_cmd = '{0} --very-sensitive-local --reorder -x bowtie_indexes/V_indexes -U {1} -S {2} -p {3}'
    j_bowtie_cmd = j_bowtie_cmd.format(args.Bowtie2exe, fastq_file, j_sam, num_cpu)
    v_bowtie_cmd = v_bowtie_cmd.format(args.Bowtie2exe, fastq_file, v_sam, num_cpu)
    subprocess.call(j_bowtie_cmd, shell=True)
    subprocess.call(v_bowtie_cmd, shell=True)

    # Process sams
    from libs.process_target_sam import parse_sams
    print 'Processing SAM files...'
    ref_names, metrics, ndn_fastq = parse_sams(j_sam, v_sam, args.OutDir, fastq_file)

    # Derep fastq
    print 'Dereplicating N-D-N sequences...'
    ndn_derep_fastq = os.path.join(args.OutDir, 'NDN_reads_derep.fastq')
    derep_fastq(ndn_fastq, ndn_derep_fastq)

    # Convert fastq to fasta for cd-hit
    from libs.convert_fastq_to_fasta import convert_fastq_to_fasta
    ndn_derep_fasta = os.path.join(args.OutDir, 'NDN_reads_derep.fasta')
    convert_fastq_to_fasta(ndn_derep_fastq, ndn_derep_fasta)

    # Make CD-HIT command template
    cdhit_templ = ['{0} -i {1} -o {2}',        # Script, input, output
                   '-c {4}',                  # Identity threshold
                   '-G 1',                     # Use global alignment
                   '-d 0',                     # Report whole seq name
                   '-s {5}',                   # Length difference cutof
                   '-r 0',                     # Only search the +/+ strand
                   '-T {3}',                   # Num threads
                   '-M 2000',                  # RAM limit
                   '-p 1',                     # Print alignment in output
                   '-l 8 -n 8']                # Length throw-away, word length
    cdhit_templ = ' '.join(cdhit_templ)

    # Run multiple iterations of clustering
    from libs.make_cluster_consensus_futures import make_consensus
    in_fasta = ndn_derep_fasta
    in_fastq = ndn_derep_fastq
    i = 1
    while i <= args.ClusIter:
        # Cluster
        print 'Clustering N-D-N sequences step{0}...'.format(i)
        clstr_out = os.path.join(args.OutDir, 'NDN_clusters_step{0}.fasta'.format(i))
        clstr_meta = clstr_out + '.clstr'
        cdhit_cmd = cdhit_templ.format(args.CDhitexe, in_fasta, clstr_out,
                                       num_cpu, args.ClusIdentity,
                                       args.ClusLen)
        subprocess.call(cdhit_cmd, shell=True)
        # Make consensus
        print 'Making consensus sequences step{0}...'.format(i)
        cons_fastq_out = os.path.join(args.OutDir, 'NDN_clusters_step{0}.fastq.consensus'.format(i))
        num_of_clusters, total_clusters_size = make_consensus(in_fastq, clstr_meta, cons_fastq_out, num_cpu)
        # Convert fastq consensus to fasta
        cons_fasta_out = os.path.join(args.OutDir, 'NDN_clusters_step{0}.fasta.consensus'.format(i))
        convert_fastq_to_fasta(cons_fastq_out, cons_fasta_out)
        # Set input for next round
        in_fasta = cons_fasta_out
        in_fastq = cons_fastq_out
        i += 1

    # Process clusters
    print 'Processing clusters...'
    from libs.process_clusters import process_clusters
    targets_out = os.path.join(args.OutDir, 'target_results.txt')
    process_clusters(cons_fasta_out, ref_names, metrics, targets_out)

    # Write a log file
    log_file = os.path.join(args.OutDir, 'log_file.txt')
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

def parse_arguments():
    """ Load command line args.
    """

    parser = argparse.ArgumentParser(description='Minimal resdiual disease target detection pipeline.')

    # Required args
    parser.add_argument('InFastq',
                        metavar='<fastq>',
                        type=str,
                        # nargs='1',
                        help='Input fastq reads.')
    parser.add_argument('OutDir',
                        metavar='<dir>',
                        type=str,
                        # nargs='1',
                        help='Output directory.')

    # Optional arguments
    parser.add_argument('--ClusIter',
                        metavar='<int>',
                        type=int,
                        required=False,
                        default=2,
                        help="Number of clustering iterations to perform. (2)")
    parser.add_argument('--ClusIdentity',
                        metavar='<float>',
                        type=float,
                        required=False,
                        default=0.9,
                        help="CD-hit clustering identity threshold. (0.9)")
    parser.add_argument('--ClusLen',
                        metavar='<float>',
                        type=float,
                        required=False,
                        default=0.9,
                        help="CD-hit clustering length similarity threshold. (0.9)")
    parser.add_argument('--Bowtie2exe',
                        metavar='<str>',
                        type=str,
                        required=False,
                        default='libs/bowtie2-2.2.1/bowtie2',
                        help="Location of Bowtie2 executable. (libs/bowtie2-2.2.1/bowtie2)")
    parser.add_argument('--CDhitexe',
                        metavar='<str>',
                        type=str,
                        required=False,
                        default='libs/cd-hit-v4.6.1-2012-08-27/cd-hit-est',
                        help="Location of CD-hit-est executable. (libs/cd-hit-v4.6.1-2012-08-27/cd-hit-est)")
    parser.add_argument('--numCPU',
                        metavar='<int>',
                        type=int,
                        required=False,
                        default=4,
                        help='Number of CPUs to use. (4)')


    return parser.parse_args()

if __name__ == '__main__':

    main()
