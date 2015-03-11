#!/usr/bin/env pypy
# -*- coding: utf-8 -*-
#
#  Script to extract the J to V lengths from a <prefix.dat> file. Lenghts will
#  be sorted (descending) and printed to stout. This should be run from within
#  the extras folder.
#
#  Usage:
#    python JtoV_lengths.py <prefix.dat>
#  To save the lengths to a file:
#    python JtoV_lengths.py <prefix.dat> > <output_file_name>
#

# Import libraries
import os
import sys
import pickle

def main():
    
    # Import the required classes for reading .dat file
    cwd_folder = os.path.split(os.getcwd())[-1]
    if cwd_folder == 'extras':
        top_path = os.pardir
    elif cwd_folder == 'NHS_MRD':
        top_path = os.getcwd()
    else:
        print ('Current folder not recognised. Script should be ran from '
               'within the extras or NHS_MRD folder.')
    sys.path.append(top_path)
    
    from libs.stage1_classes import Records
    from libs.stage1_classes import Cluster
    from libs.stage1_classes import Segment
    
    # Get dat file name and load the pickle
    dat_file = sys.argv[1]
    dat_records = pickle.load(open(dat_file, 'r'))
    
    # For each record print the distance of V from beginnig of read
    lengths = []
    
    for record in dat_records.record_list:
        if len(record.v) > 0:
            lengths.append(record.v[0].q_start)
        elif len(record.j) > 0:
            lengths.append('NA')
    
    for length in sorted(lengths, reverse=True):
        print length
        
    return 0

if __name__ == '__main__':
    main()
