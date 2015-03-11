#!/usr/bin/env python
# -*- coding: utf-8 -*-

from clustering import distance

#
# Stage 2 functions
#

def write_matches(match_list, sample_obj, total_reads, handle):
    """Will print all the matches from a given sample
    """
    # Write header
    handle.write('\tMatches in {0} sample:\n'.format(sample_obj.title))
    # Write 2nd header
    handle.write('\t'.join(['\t',
                                'cluster_num',
                                'mis-matches',
                                'norm_size',
                                'actual_size',
                                'proportion',
                                'ID\n']))
    exact_sum = 0
    # Write all matches
    for match_entry in match_list:
        match_record = match_entry[0]
        mismatches = match_entry[1]
        norm_size = (match_record.size / (float(total_reads) / 1000000))
        if mismatches == 0:
            exact_sum += norm_size
        prop = float(match_record.size) / sample_obj.total_clusters_size
        handle.write('\t'.join(['\t',
                                str(match_record.clus_num),
                                str(mismatches),
                                str(norm_size),
                                str(match_record.size),
                                str(prop),
                                match_record.id]))
        handle.write('\n')
    handle.write('\t\tSum of clusters with 0 mismatches: {0}\n'.format(
                                                                exact_sum))
    handle.write('\n')
    
    return exact_sum

def find_matches(query_seq, sample_records, matches):
    """Given a list of records it will calculate the distance to each of the
       records and return a sorted list of the closest matches
    """
    distances = {}
    query_len = len(query_seq)
    # Calculate the distance to each target sequence in other sample
    for target in sample_records:
        distances[target] = distance(query_seq,
                                     target.seq[:query_len],
                                     range(query_len))
    # Sort distances into list
    distances_l = sorted(distances.iteritems(),
                         key=lambda x: (x[1], x[0].clus_num))
    # Get the top N matches. If more than one target has the same
    # distance then all will be reported.
    match_list = []
    while (len(match_list) < matches and
           len(distances_l) > 0):
        top_dist = distances_l[0][1]
        while distances_l[0][1] == top_dist:
            match_list.append(distances_l[0])
            del distances_l[0]
            
    return match_list
