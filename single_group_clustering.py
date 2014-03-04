# python script <in_fasta> <out_fasta>

import os
import sys
from libs.bio_file_parsers import fasta_parser
from libs.bio_file_parsers import write_fasta
from Levenshtein import ratio
from sklearn.cluster import DBSCAN
import numpy as np

import pylab as pl
from mpl_toolkits.mplot3d import Axes3D
from sklearn import decomposition

def main(args):
    
    # Parse the fasta into a list of headers and seqs
    read_names = []
    read_seqs = []
    with open(args['fasta_in']) as fasta_in:
        for header, seq in fasta_parser(fasta_in):
            read_names.append(header)
            read_seqs.append(seq)
    
    print 'Calculating distances...'
    # Populate dist matrix with Nones
    dist_matrix = []
    for i in range(len(read_seqs)):
        dist_matrix.append([None] * len(read_seqs))
    # Calculate distances
    for i in range(len(read_seqs)):
        for j in range(i, len(read_seqs)):
            dist = 1 - ratio(read_seqs[i], read_seqs[j])
            dist_matrix[i][j] = dist
            dist_matrix[j][i] = dist
    dist_matrix = np.array(dist_matrix)
    
    # Perform DBSCAN
    print 'Clusting by DBSCAN...'
    db = DBSCAN(eps=0.05, min_samples=2, metric='precomputed').fit(dist_matrix)
    print '- Number of clusters: {0}'.format(len(set(db.labels_)))
    
    write_clusters(read_names, read_seqs, db.labels_, 'results/clusters')
    sys.exit()
    
    # Do PCA
    print 'Doing PCA...'
    pca = decomposition.PCA(n_components=3).fit(dist_matrix)
    X = pca.transform(dist_matrix)
    
    # Plot PCA
    colours =  pl.cm.get_cmap("rainbow", len(set(db.labels_))) 
    fig = pl.figure()
    ax = Axes3D(fig)
    for row in range(X.shape[0]):
        # Get colour
        label = int(db.labels_[row])
        if label == -1:
            col = 'k'
        else:
            col = colours(label)
        # Plot
        ax.scatter(X[row, 0], X[row, 1], X[row, 2], color=col)
    pl.show()    

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
    
    args['fasta_in'] = sys.argv[1]
    
    main(args)
 
