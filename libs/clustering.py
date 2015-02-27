#!/usr/bin/env python
# -*- coding: utf-8 -*-

import shelve
import bio_file_parsers as parser
import sys
import time

def dereplication(in_file, min_size, derep_ref_file=None):
    """Will take a fasta file and return a list of records with each containing
       only unique sequences. It will also return a dictionary with pointers
       from each read to its centroid.
    """
    # Use seq as dict key and list of ids as values
    total_filtered_reads = 0
    unique_dict = {}
    with open(in_file, 'r') as in_handle:
        for record in parser.fasta_parser(in_handle):
            total_filtered_reads += 1
            header = str(record[0])
            seq = str(record[1])
            try:
                unique_dict[seq].append(header)
            except KeyError:
                unique_dict[seq] = [header]
    
    # Make a list of records, in size decending order
    record_list = []
    for entry in sorted(unique_dict.items(), key=lambda x: len(x[1]),
                        reverse=True):
        size = len(entry[1])
        if size >= min_size:
            record = {}
            record['seq'] = entry[0]
            record['id'] = entry[1][0]
            record['size'] = size
            record_list.append(record)
        else:
            break
    
    if derep_ref_file:
        # Make a reference dictionary with pointers from reads to centroid
        derep_refs = shelve.open(derep_ref_file)
        for seq in unique_dict:
            headers = unique_dict[seq]
            centroid = headers[0]
            for identical_copy in headers:
                derep_refs[identical_copy] = centroid
        derep_refs.close()
    
    return record_list, len(unique_dict), total_filtered_reads

def cluster(record_list, min_identity, cluster_refs_file=None):
    """Takes a list of records in decesending group size order and clusters 
       them according to similarity (caluclated using edit distance). 
       
       If the distance between two sequences is < auto_merge_dist the 
       sequences will be merged irrespective of the size differences 
       between them.
       
       However if the distance is > auto_merge_dist, then merging will only 
       occur if merge_function retruns True.
       
       Return: A dictionary entry for each cluster, containing a list of 
       records that make up the cluster.
    """
    # Initiate variables
    total_records_len = len(record_list)
    clusters = {}
    
    # All sequences should be the same size because of motif_filter.
    seq_len = len(record_list[0]['seq'])
    seq_len_range = range(seq_len)
    max_dist = int(seq_len - (seq_len * (float(min_identity) / 100)))
    
    not_merged = record_list
    
    while len(not_merged) > 0:
        centroid = not_merged[0]
        # Print progress to stdout
        if len(not_merged)%20 == 0:
            print_progress(len(not_merged), total_records_len)          
        # Make new entry in clusters dict
        clusters[centroid['id']] = [centroid]
        # Create list of targets
        target_list = not_merged[1:]
        not_merged = []
        # Iterate through all other records
        for target in target_list:
            # Calculate distance between centroid and target
            dist = distance(str(centroid['seq']),
                            str(target['seq']),
                            seq_len_range,
                            max_dist)
            # Merge clusters if dist is <= to max_dist
            if dist <= max_dist:
                clusters[centroid['id']].append(target)
            else:
                not_merged.append(target)
    
    # Update progress readout to 100%
    if total_records_len:
        print_progress(len(not_merged), total_records_len)
        print
    
    if cluster_refs_file:
        # Make reference dict with pointers from reads to cluster centroid
        cluster_refs = shelve.open(cluster_refs_file)
        for centroid_id in clusters:
            for target in clusters[centroid_id]:
                cluster_refs[target['id']] = centroid_id
        cluster_refs.close()
    
    return clusters

def print_progress(not_merged, total):
    """Prints progress to the stdout to keep user informed.
    """
    progress = 100 - ((100 * not_merged) / total)
    sys.stdout.write('\rProgress: {0}%'.format(progress))
    sys.stdout.flush()

def distance(seq1 , seq2, seq_len_range, max_dist=None):
    """ Very simple distance function. I kept it simple to hopefully speed
        things up.
    """
    
    dist = 0
    for i in seq_len_range:
        if seq1[i] != seq2[i]:
            dist += 1
            # Exit early if max_dist has been reached
            if max_dist and dist > max_dist:
                return dist
    
    return dist

def merge_function(c, dist, size1, size2):
    """Function will return True if clusters should be merged and False
       if they should be separate.
    """
    if ((float(size2)/size1) < (c**dist)):
        return True
    else:
        return False

def consensus(record_list):
    """Will take a list of records and return a record containing the consensus
       sequence.
    """
    consensus_id = record_list[0]['id']
    seq_len = len(record_list[0]['seq'])
    
    # Build the dictionary
    position_matches = {}
    for i in range(seq_len):
        position_matches[i] = {'A': 0, 'C': 0, 'T': 0, 'G': 0, 'N': 0}
    # Record how many of each base there are at each position
    total_cluster_size = 0
    for record in record_list:
        num_identical_reads = record['size']
        total_cluster_size += num_identical_reads
        sequence = str(record['seq']).upper()
        for i in range(seq_len):
            position_matches[i][sequence[i]] += num_identical_reads
    # Calculate consensus sequence
    consensus_seq = ''
    for i in range(seq_len):
        consensus_seq += highest_key(position_matches[i])
    # Make a consensus record
    consensus_record = {}
    consensus_record['id'] = consensus_id
    consensus_record['seq'] = consensus_seq
    consensus_record['size'] = total_cluster_size
    
    return consensus_record

def highest_key(dic):
    """Given a dictionary with str(key) and int(values) it will return
       the key with the highest value.
    """
    top_key = None
    top_value = None
    for key in dic:
        if dic[key] > top_value:
            top_value = dic[key]
            top_key = key
    return top_key

def clusters_to_consensus(cluster_dict):
    """Will take the dictionary of clusters produced by clustering() and 
       turn it into an ordered list of consensus records
    """
    consensus_record_list = []
    for entry in sorted(cluster_dict.items(), key=lambda x: len(x[1]),
                              reverse=True):
        consensus_record_list.append(consensus(entry[1]))
    return consensus_record_list

def write_consensus_fasta(consensus_record_list, out_file):
    """Will write all the consensus records to a fasta file with size appended
       to the header
    """
    with open(out_file, 'w') as out_handle:
        for record in consensus_record_list:
            size = record.annotations['size']
            new_record = SeqRecord(seq=Seq(record.seq), description = '',
                            id='centroid={0};size={1}'.format(record.id,size))
            SeqIO.write(new_record, out_handle, 'fasta')

def log_stats(log_file, inreads, groups, groups_more1, num_clus, start_time):
    with open(log_file, 'a') as out_handle:
        out_handle.write('\n'.join(['','--Log info for clustering--',
            'Total input reads: {0}'.format(inreads),
            'Derep groups: {0}'.format(groups),
            'Derep groups size >1: {0}'.format(groups_more1),
            'Number of clusters: {0}'.format(num_clus),
            'Time taken: {0:.1f} seconds'.format(time.time() - start_time)]))
