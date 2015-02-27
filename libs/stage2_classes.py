#!/usr/bin/env python
# -*- coding: utf-8 -*-

import stage2_funcs as s2f

class Match:
    
    def __init__(self, query_record):
        self.query = query_record

class Results:
    # Class that will carry out a comparison two groups of records and
    # store the results in a useful way
    
    def __init__(self, sample1, sample2, top, matches):
        # Initiate variables
        self.sample1 = sample1
        self.sample2 = sample2
        self.top = top
        # Run comparison
        self.compare_results(self.sample1, self.sample2, top, matches)
    
    def compare_results(self, sample1, sample2, top, matches):
        """Will compare the each query sequence to the clusters in each of the
           samples to find matches. Also, calculates what cluster number each
           cluster is (for easier reference in the output).
        """
                
        # Need to work out what number cluster each record each sample is, and
        # store it as an attribute
        for record_list in [sample1.record_list, sample2.record_list]:
            cluster_num = 0
            for target_record in record_list:
                cluster_num += 1
                target_record.clus_num = cluster_num
        
        # Initiate variables
        self.results = []
        
        # For each of the top clusters (query)
        for query_record in sample1.record_list[:top]:
            
            # Find the starting position of the V segment
            try:
                v_pos = query_record.v[0].q_start - 1
            except IndexError:
                v_pos = None
            
            # Save query record to record class
            query_record.query_seq = query_record.seq[:v_pos]
            match_object = Match(query_record)
            
            # Look for matches in sample1
            match_object.s1_matches = s2f.find_matches(query_record.query_seq,
                                                       sample1.record_list,
                                                       matches)
            
            # Look for matches in sample2
            match_object.s2_matches = s2f.find_matches(query_record.query_seq,
                                                       sample2.record_list,
                                                       matches)
                                                   
            self.results.append(match_object)

    def write_shared_clusters(self, file_out, sample1_reads, sample2_reads):        
        with open(file_out,'a') as out_handle:
            
            # Write header
            line = ('#\n# Query sequences derived from the top {0} clusters' +
                    ' in the {1} sample\n#\n\n'.format(self.top,
                                                       self.sample1.title))
            out_handle.write(line)
            # For each of the queries
            query_num = 0
            for entry in self.results:
                query_num += 1
                
                # Print information about the query cluster
                query_record = entry.query
                line = '[{0} {1}] Query seq: {2} ({3} bp)\n'.format(
                        self.sample1.title,
                        query_num,
                        query_record.query_seq,
                        len(query_record.query_seq))
                out_handle.write(line)
                
                # Write matches in sample 1
                s1_size = s2f.write_matches(entry.s1_matches,
                                            self.sample1,
                                            sample1_reads,
                                            out_handle)
                
                # Write matches in sample 2
                s2_size = s2f.write_matches(entry.s2_matches,
                                            self.sample2,
                                            sample2_reads,
                                            out_handle)
                
                # Calculate and write fold-change
                if s1_size > 0 and s2_size > 0:
                    fc = s2_size/s1_size
                else:
                    fc = 'N/A'
                
                out_handle.write('\tFold-change between samples ' +
                                 '{0}/{1} = {2}\n'.format(s2_size,
                                                          s1_size,
                                                          fc))
                out_handle.write('\n')
