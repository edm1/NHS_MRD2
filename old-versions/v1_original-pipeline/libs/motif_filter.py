#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gzip
from contextlib import nested
import re
import time
import bio_file_parsers as parser

def return_handle(file_name):
    """If file_name has .gz then return uncompressed handle
    """
    if file_name.split('.')[-1] == 'gz':
        return gzip.open(file_name, 'r')
    else:
        return open(file_name,'r')

def filter_reads(j_reads,   # File containing J reads
                 out_id,    # Output identifier prefix
                 out_dir,   # Output directory
                 j_motif,   # J specific regex
                 germ_regex,    # Germline specific regex
                 j_after_bp,    # Number of bp after motif to keep
                 j_before_bp,   # Number of bp before motif to keep
                 max_dist,  # Max distance J motif can be from beginning
                 log_file = None):
    """Will do all the read filtering and return the path to the J and
       V end reads.
    """
    
    start_time = time.time()
    
    # Start counters
    match_pos, match_neg = 0,0
    match_far, match_close = 0,0
    germ_pos, germ_neg = 0,0
    correct_len, too_short = 0, 0
    total = 0
    
    # Calc how long J reads should be after clipping
    total_clip_len = j_after_bp + j_before_bp + len(j_motif)
    
    # Correct batch holder
    correct_batch = []
    
    # Open out/in handles
    with nested(open('{0}/{1}_Jreads.fa'.format(out_dir.rstrip('/'),
                                                out_id),'w'),
                open('{0}/{1}_germline.fa'.format(out_dir.rstrip('/'),
                                                  out_id), 'w'),
                open('{0}/{1}_no_motif.fa'.format(out_dir.rstrip('/'),
                                                  out_id), 'w'),
                open('{0}/{1}_motif_too_far.fa'.format(out_dir.rstrip('/'),
                                                  out_id), 'w'),
                open('{0}/{1}_motif_too_close.fa'.format(out_dir.rstrip('/'),
                                                  out_id), 'w'),
                return_handle(j_reads)
                ) as (out_J,
                      out_germ,
                      out_nomatch,
                      out_too_far,
                      out_too_close,
                      in_J):
        
        pattern = re.compile(j_motif)
        germ_pattern = re.compile(germ_regex)
        
        for j_record in parser.fastq_parser(in_J):
            record_header = j_record[0]
            record_seq = j_record[1]
            total += 1
            
            # Search to find required motif tag 
            match_obj = pattern.search(str(record_seq))
            if match_obj:
                match_pos += 1
                # Check the tag is within max distance
                if match_obj.start() < max_dist+1:
                    # Check tag is further than min distance
                    if match_obj.start() - j_before_bp > 0:
                        # Search sequence to see if it has germline tags
                        germ_match_obj = germ_pattern.search(str(record_seq))
                        if germ_match_obj:
                            germ_pos += 1
                            # Save germline seqs to separate file
                            parser.write_fasta(out_germ,
                                               record_header,
                                               record_seq)
                        else:
                            germ_neg += 1
                            sub_record_seq = record_seq[(match_obj.start() - 
                            j_before_bp):match_obj.end() + j_after_bp]
                            
                            # Check clip length is correct
                            if len(sub_record_seq) == total_clip_len:
                                correct_len += 1
                                
                                # Write the correct J end reads to file in
                                # batches (in order to improve IO time).
                                correct_batch.append((record_header,
                                                      sub_record_seq))
                                if len(correct_batch) == 10000:
                                    for corr_record in correct_batch:
                                        parser.write_fasta(out_J,
                                                           corr_record[0],
                                                           corr_record[1])
                                        correct_batch = []
                            else:
                                too_short += 1
                    else:
                        match_close += 1
                        parser.write_fasta(out_too_close,
                                           record_header,
                                           record_seq)
                else:
                    parser.write_fasta(out_too_far,
                                       record_header,
                                       record_seq)
                    match_far += 1
            else:
                parser.write_fasta(out_nomatch,
                                   record_header,
                                   record_seq)
                match_neg += 1
        
        # Clear the correct batch
        for corr_record in correct_batch:
            parser.write_fasta(out_J,
                               corr_record[0],
                               corr_record[1])
            
    if log_file:
        with open(log_file, 'a') as out_handle:
            out_handle.write('\n'.join(["\n--Log info for motif_filter--",
                "Sequences total: {0}".format(total),
                "  Sequences with tag: {0}".format(match_pos),
                "  Sequences with no tag: {0}".format(match_neg),
                "    Matches correct distance: {0}".format((match_pos - 
                    match_far - match_close)),
                "    Matches too far from 5': {0}".format(match_far),
                "    Matches too close to 5': {0}".format(match_close),
                "      Matches with no germline tag: {0}".format(germ_neg),
                "      Matches with germline tag: {0}".format(germ_pos),
                "        Correct length: {0}".format(correct_len),
                "        Too short: {0}".format(too_short),
                "Filtering took {0:.1f} seconds.".format((time.time() - 
                                                         start_time)),
                ""]))

    return '{0}/{1}_Jreads.fa'.format(out_dir.rstrip('/'), out_id), total
