#!/usr/bin/env python
# -*- coding: utf-8 -*-

#~ # Import modules
#~ import os
#~ import time
#~ import pickle
#~ from shutil import rmtree
#~ 
#~ # Import functions
#~ from stage2_functions import write_matches
#~ from stage2_functions import find_matches
#~ 
#~ # Import classes
#~ from stage2_classes import Match
#~ from stage2_classes import Results

import time
import os
from shutil import rmtree
import pickle
import libs.stage2_classes as s2c

def main(args):
    
    if not args['--quiet']:
        print "\n{1} Job name {0} ---------------".format(
                args['<id>'],time.strftime('%H:%M:%S',time.localtime()))
    
    # Check if project folder already exists
    folder_name = './comparisons/{0}/'.format(args['<id>'])
    if os.path.exists(folder_name):
        resp = raw_input('This project folder already exists, would you ' + 
                         'like to overwrite it? (y/n): ')
        if resp == 'y':
            rmtree(folder_name)
        else:
            print 'Run again with different --id parameter'
            sys.exit()

    # Write new project folders
    os.mkdir(folder_name)

    # Load target and query results
    with open(args['<sample1.dat>'], 'r') as in_handle:
        diagnosis_records = pickle.load(in_handle)
        diagnosis_records.title = 'diagnostic'
    with open(args['<sample2.dat>'], 'r') as in_handle:
        followup_records = pickle.load(in_handle)
        followup_records.title = 'follow-up'

    # Run comparison
    # Diagnosis vs follow-up
    dia_size = diagnosis_records.total_input_seqs
    fol_size = followup_records.total_input_seqs
    
    results1 = s2c.Results(diagnosis_records,
                           followup_records,
                           args['--queries'],
                           args['--matches'])
    share_clus_file = './comparisons/{0}/{0}.shared'.format(args['<id>'])
    results1.write_shared_clusters(share_clus_file,
                                   dia_size,
                                   fol_size)
    del results1
    
    # Follow-up vs diagnosis
    results1 = s2c.Results(followup_records,
                           diagnosis_records,
                           args['--queries'],
                           args['--matches'])
    share_clus_file = './comparisons/{0}/{0}.shared'.format(args['<id>'])
    results1.write_shared_clusters(share_clus_file,
                                   fol_size,
                                   dia_size)
    del results1
    
    if not args['--quiet']:
        print "{0} Finished!\n".format(time.strftime('%H:%M:%S',
                                                     time.localtime()))
    
    return 0

if __name__ == '__main__':
	main()

