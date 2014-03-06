# python script <in_fasta> <out_fasta>

import os
import sys
import glob
from libs.bio_file_parsers import fasta_parser
from libs.bio_file_parsers import write_fasta

from sklearn.cluster import DBSCAN
from Levenshtein import ratio
import numpy as np

def main(args):
    
    cluster_groups(args['group_folder'], 'test/clusters')

def cluster_groups(groups_dir, cluster_dir):
    """ Main function.
    """
    i = 1 # DEBUG
    
    group_fasta_list = glob.glob(os.path.join(groups_dir, 'group_*.fasta'))
    for group_fasta in group_fasta_list:
        
        # DEBUG
        print 'Processed: {0} / {1}'.format(i, len(group_fasta_list))
        print os.path.split(group_fasta)[-1]
        i += 1
        
        # Get group number and make group folder
        group_name = os.path.split(group_fasta)[-1].split('.')[0]
        clus_group_dir = os.path.join(cluster_dir, group_name)
        os.mkdir(clus_group_dir)
        
        # Dereplicate the fasta file
        read_names, read_seqs = derep_fasta(group_fasta)
        
        # Calc dist matrix
        dist_matrix = calc_dist_matrix(read_seqs)
        
        # Do DBSCAN clustering
        db = DBSCAN(eps=0.05, min_samples=2, metric='precomputed')
        db_fit = db.fit(dist_matrix)
        
        # Write clusters
        write_clusters(read_names, read_seqs, db_fit.labels_, clus_group_dir)

def calc_dist_matrix(read_seqs):
    
    rl = len(read_seqs) # DEBUG
    
    # Populate dist matrix with Nones
    dist_matrix = []
    for i in range(len(read_seqs)):
        dist_matrix.append([None] * len(read_seqs))
    # Calculate distances
    for i in range(len(read_seqs)):
        
        print '\tClustering: {0} / {1}'.format(i, rl) # DEBUG
        
        for j in range(i, len(read_seqs)):
            dist = 1 - ratio(read_seqs[i], read_seqs[j])
            dist_matrix[i][j] = dist
            dist_matrix[j][i] = dist
    
    return np.array(dist_matrix)

def derep_fasta(fasta_file):
    """ Takes a fasta file, dereplicates it and returns a list of headers and
        a list of seqs.
    """
    seqs = {}
    with open(fasta_file) as fasta_in:
        # Parse the fasta
        for header, seq in fasta_parser(fasta_in):
            size = int(header.split(';size=')[-1])
            if not seq in seqs:
                seqs[seq] = {'title': header,
                             'size': size}
            else:
                seqs[seq]['size'] += size
    
    # Put headers and seqs into two lists
    read_names = []
    read_seqs = []
    for key in sorted(seqs.keys(), key=lambda x: seqs[x]['size'], reverse=True):
        title = '{0};size={1}'.format(seqs[key]['title'].split(';size=')[0],
                                      seqs[key]['size'])
        read_names.append(title)
        read_seqs.append(key)

    return read_names, read_seqs

def write_clusters(read_names, read_seqs, clusters, out_dir):
    """ Write the results of DBSCAN to fasta files.
    """
    # Convert clusters to ints
    clusters = [int(x) for x in clusters]
    
    for label in set(clusters):
        # Get indices for that label
        indicies = [i for i, x in enumerate(clusters) if x == label]
        # Open out fasta
        fasta_file = os.path.join(out_dir, 'cluster_{0}.fasta'.format(label))
        with open(fasta_file, 'w') as out_handle:
            for idx in indicies:
                write_fasta(out_handle, read_names[idx], read_seqs[idx])

if __name__ == '__main__':
    
    args = {}
    
    args['group_folder'] = sys.argv[1]
    
    main(args)
 
