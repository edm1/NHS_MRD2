# pypy script <in_fasta> <out_fasta>

import sys
from libs.bio_file_parsers import fasta_parser
import cPickle
import numpy as np

def main(args):
    
    print 'Loading pickle...'
    
    #~ import cPickle
    #~ matrix = cPickle.load(open('results/kmer_obj.pickle', 'r'))
    
    kmer_array = np.load('results/kmer_obj.pickle')
    
    print 'Loading matrix...'
    #~ kmer_array = np.array(matrix)
    kmer_array = np.transpose(kmer_array)
    #~ del matrix
    
    #~ kmer_array = kmer_array[:1000]
    print kmer_array.shape
    
    print 'PCA...'
    import time
    #~ time.sleep(200)
    
    # Do PCA using scikit
    import pylab as pl
    from mpl_toolkits.mplot3d import Axes3D
    from sklearn import decomposition
    from sklearn.preprocessing import scale
    
    fig = pl.figure(1, figsize=(4, 3))
    pl.clf()
    ax = Axes3D(fig)
    
    kmer_array = scale(kmer_array, axis=0)
    
    pl.cla()
    pca = decomposition.PCA(n_components=400)
    pca.fit(kmer_array)
    
    print pca.explained_variance_ratio_
    cum = np.cumsum(pca.explained_variance_ratio_)
    print cum
    
    X = pca.transform(kmer_array)
    print X.shape
    
    print 'Plotting...'
    ax.scatter(X[:, 0], X[:, 1], X[:, 2], cmap=pl.cm.spectral)
    pl.show()
    
    #~ fig = pl.figure()
    #~ pl.plot(range(1, len(cum) + 1), cum)
    #~ pl.show()


if __name__ == '__main__':
    
    args = {}
    #~ args['in_fasta'] = sys.argv[1]
    
    main(args)
