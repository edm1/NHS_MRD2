# python script <in_fasta> <out_fasta>

import os
import sys
from libs.bio_file_parsers import fasta_parser
from Levenshtein import distance
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
    
    size = 300
    
    # Calculate distance matrix
    print 'Calculating distances...'
    dist_matrix = []
    for seq in read_seqs[:size]:
        seq_dists = []
        for seq2 in read_seqs[:size]:
            seq_dists.append(distance(seq, seq2))
        dist_matrix.append(seq_dists)
    dist_matrix = np.array(dist_matrix)
    
    # Perform DBSCAN
    print 'Clusting by DBSCAN...'
    db = DBSCAN(eps=2, min_samples=2, metric='precomputed').fit(dist_matrix)
    print '- Number of clusters: {0}'.format(len(set(db.labels_)))
    
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

if __name__ == '__main__':
    
    args = {}
    
    args['fasta_in'] = sys.argv[1]
    
    main(args)
